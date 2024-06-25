from transformers import BertTokenizer
from typing import List, Optional


class Modified_BertTokenizer(BertTokenizer):
    def build_inputs_with_special_tokens(
        self, token_ids_0: List[int], token_ids_1: Optional[List[int]] = None
    ) -> List[int]:
        if token_ids_1 is None:
            #return [self.cls_token_id] + token_ids_0 + [self.sep_token_id]
            return token_ids_0
        #cls = [self.cls_token_id]
        #sep = [self.sep_token_id]
        #return cls + token_ids_0 + sep + token_ids_1 + sep 
        return token_ids_0 + token_ids_1


    def create_token_type_ids_from_sequences(
        self, token_ids_0: List[int], token_ids_1: Optional[List[int]] = None
    ) -> List[int]:
        #sep = [self.sep_token_id]
        #cls = [self.cls_token_id]
        if token_ids_1 is None:
            #return len(cls + token_ids_0 + sep) * [0]
            return len(token_ids_0) * [0]
        #return len(cls + token_ids_0 + sep) * [0] + len(token_ids_1 + sep) * [1]
        return len(token_ids_0) * [0] + len(token_ids_1) * [1]
