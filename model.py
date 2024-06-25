import sys
sys.path.append('library/')
from transformers import BertForTokenClassification, BertModel
from typing import Optional, Union, Tuple
import torch
from utils import *
import torch.nn.functional as F
import pickle as pkl
import numpy as np


class Modified_BertModel(BertModel):
    def __init__(self, config, Class, add_pooling_layer=False):
        super().__init__(config)
        self.max_len = config.max_position_embeddings   # 56 by default
        self.hidden_size = config.hidden_size
        self.position_embeddings = torch.nn.Embedding(self.max_len, self.hidden_size)   # positional encodings for the 56 protein tokens
        self.LayerNorm = torch.nn.LayerNorm(self.hidden_size, eps=1e-12)
        self.dropout = torch.nn.Dropout(config.hidden_dropout_prob)
        self.raw_embed_dim = 1024
        self.fc = torch.nn.Linear(self.raw_embed_dim, config.hidden_size)
        del self.pooler, self.embeddings
        self.post_init()
        self.Class = Class

    def forward(
        self,
        proteins: Optional[torch.Tensor] = None,
        input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        encoder_hidden_states: Optional[torch.Tensor] = None,
        encoder_attention_mask: Optional[torch.Tensor] = None,
        past_key_values: Optional[List[torch.FloatTensor]] = None,
        use_cache: Optional[bool] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
        raw_embed: Optional[np.array] = None,
    ) -> Union[Tuple[torch.Tensor], BaseModelOutputWithPoolingAndCrossAttentions]:
        
        # customized part
        device = proteins.device
        batch_size, seq_length = proteins.size()
        input_shape = proteins.size()

        # global attention mask
        zero_tensor = torch.tensor(0, device=device)
        one_tensor = torch.tensor(1, device=device)
        attention_mask = torch.where(proteins == zero_tensor, zero_tensor, one_tensor)

        # assemble inputs_embeds
        proteins = proteins.view(-1, 1)
        raw_embedding_output = np.zeros((proteins.shape[0], self.raw_embed_dim))
        idx = -1
        for protein in proteins:
            idx+=1
            if(protein==0):
                continue 
            else:
                raw_embedding_output[idx] = raw_embed[protein-1]
        raw_embedding_output = torch.tensor(raw_embedding_output).to(device)
        raw_embedding_output = raw_embedding_output.to(torch.float32)
        raw_embedding_output = raw_embedding_output.view(batch_size, self.max_len, self.raw_embed_dim)
        embedding_output = self.fc(raw_embedding_output)    # token embedding

        # add positional embedding
        position_ids = torch.arange(self.max_len).expand((1, -1)).to(device)
        position_embeddings = self.position_embeddings(position_ids)
        position_embeddings = position_embeddings.repeat(batch_size, 1, 1)
        embedding_output = embedding_output + position_embeddings
        embedding_output = self.LayerNorm(embedding_output)
        embedding_output = self.dropout(embedding_output)

        # original part of huggingface 
        output_attentions = output_attentions if output_attentions is not None else self.config.output_attentions
        output_hidden_states = (
            output_hidden_states if output_hidden_states is not None else self.config.output_hidden_states
        )
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict
        use_cache = False
        past_key_values_length = past_key_values[0][0].shape[2] if past_key_values is not None else 0
        if attention_mask is None:
            attention_mask = torch.ones(((batch_size, seq_length + past_key_values_length)), device=device)
        extended_attention_mask: torch.Tensor = self.get_extended_attention_mask(attention_mask, input_shape)
        if self.config.is_decoder and encoder_hidden_states is not None:
            encoder_batch_size, encoder_sequence_length, _ = encoder_hidden_states.size()
            encoder_hidden_shape = (encoder_batch_size, encoder_sequence_length)
            if encoder_attention_mask is None:
                encoder_attention_mask = torch.ones(encoder_hidden_shape, device=device)
            encoder_extended_attention_mask = self.invert_attention_mask(encoder_attention_mask)
        else:
            encoder_extended_attention_mask = None
        head_mask = self.get_head_mask(head_mask, self.config.num_hidden_layers)

        # Transformer encoders
        encoder_outputs = self.encoder(
            embedding_output,
            attention_mask=extended_attention_mask,
            head_mask=head_mask,
            encoder_hidden_states=encoder_hidden_states,
            encoder_attention_mask=encoder_extended_attention_mask,
            past_key_values=past_key_values,
            use_cache=use_cache,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )
        sequence_output = encoder_outputs[0]
        sequence_output = sequence_output if self.Class in ['cc', 'mf'] else torch.cat((raw_embedding_output, sequence_output), dim=2)
        return sequence_output, attention_mask


def masked_BCEWithLogitsLoss(output, target, ignore_index): # mask padding tokens
    mask = target.ne(ignore_index)
    output = output[mask]
    target = target[mask]
    return F.binary_cross_entropy_with_logits(output, target)


class plasgo_model(BertForTokenClassification): # global BERT for multi-label token classification
    def __init__(self, config, Class, raw_embed):  # Class in ['mf', 'bp', 'cc']
        super().__init__(config)
        self.raw_embed_dim = 1024   # dimension of the original per-protein embeddings generated by ProtTrans-ProtT5
        self.importance_fc = torch.nn.Linear(self.raw_embed_dim+config.hidden_size, self.num_labels) if Class=='bp' else torch.nn.Linear(config.hidden_size, self.num_labels)        # learn confidence scores
        self.classifier = torch.nn.Linear(self.raw_embed_dim+config.hidden_size, config.num_labels) if Class=='bp' else torch.nn.Linear(config.hidden_size, config.num_labels)       # classifier
        self.sigmoid = torch.nn.Sigmoid()
        self.margin = 0.15  # for Loss_RR
        self.post_init()
        self.bert = Modified_BertModel(config, Class, add_pooling_layer=False)  # BERT model with modified token embedding layer
        self.raw_embed = raw_embed  # original per-protein embedding matrix generated by ProtT5

    def forward(
        self,
        proteins: Optional[torch.Tensor] = None,    # protein index from 1, 0 is padding token
        input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple[torch.Tensor], TokenClassifierOutput]:
        
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict
        
        # local BERT
        sequence_output, attention_mask = self.bert(
            proteins = proteins,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
            raw_embed=self.raw_embed,
        )

        importance_weight = self.importance_fc(sequence_output)
        importance_weight = self.sigmoid(importance_weight) # confidence scores
        sequence_output = self.dropout(sequence_output)
        logits = self.classifier(sequence_output)   # for multi-label prediction
        logits = torch.mul(logits, importance_weight)   # dot product

        loss = None
        if labels is not None:
            ignore_index = -100
            loss_bce = masked_BCEWithLogitsLoss(logits.view(-1, self.num_labels), labels.view(-1, self.num_labels), ignore_index)
            
            # loss_rr
            sorted_tensor = torch.flatten(importance_weight * attention_mask.unsqueeze(-1))
            nonzero_indices = torch.nonzero(sorted_tensor)
            sorted_tensor = sorted_tensor[nonzero_indices].squeeze()
            sorted_tensor, _ = torch.sort(sorted_tensor, descending=True)
            split_index = int(0.7 * len(sorted_tensor)) # 0.7 high and 0.3 low
            first_70_percent = sorted_tensor[:split_index]
            last_30_percent = sorted_tensor[split_index:]
            mean_first_70 = torch.mean(first_70_percent)
            mean_last_30 = torch.mean(last_30_percent)
            loss_rr = torch.max(self.margin-(mean_first_70-mean_last_30), torch.tensor(0))
            
            loss = loss_bce+loss_rr # final loss

        return TokenClassifierOutput(
            loss=loss,
            logits=logits,
            confidence_score = importance_weight
        )
