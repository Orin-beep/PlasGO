import torch
from typing import Optional, Tuple, Any, List, Iterable
from collections import OrderedDict
from import_utils import is_torch_available
from dataclasses import is_dataclass, dataclass, fields


class ModelOutput(OrderedDict):
    def __init_subclass__(cls) -> None:
        """Register subclasses as pytree nodes.

        This is necessary to synchronize gradients when using `torch.nn.parallel.DistributedDataParallel` with
        `static_graph=True` with modules that output `ModelOutput` subclasses.
        """
        if is_torch_available():
            _torch_pytree._register_pytree_node(
                cls,
                _model_output_flatten,
                _model_output_unflatten,
            )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Subclasses of ModelOutput must use the @dataclass decorator
        # This check is done in __init__ because the @dataclass decorator operates after __init_subclass__
        # issubclass() would return True for issubclass(ModelOutput, ModelOutput) when False is needed
        # Just need to check that the current class is not ModelOutput
        is_modeloutput_subclass = self.__class__ != ModelOutput

        if is_modeloutput_subclass and not is_dataclass(self):
            raise TypeError(
                f"{self.__module__}.{self.__class__.__name__} is not a dataclasss."
                " This is a subclass of ModelOutput and so must use the @dataclass decorator."
            )

    def __post_init__(self):
        """Check the ModelOutput dataclass.

        Only occurs if @dataclass decorator has been used.
        """
        class_fields = fields(self)

        # Safety and consistency checks
        if not len(class_fields):
            raise ValueError(f"{self.__class__.__name__} has no fields.")
        if not all(field.default is None for field in class_fields[1:]):
            raise ValueError(f"{self.__class__.__name__} should not have more than one required field.")

        first_field = getattr(self, class_fields[0].name)
        other_fields_are_none = all(getattr(self, field.name) is None for field in class_fields[1:])

        if other_fields_are_none and not is_tensor(first_field):
            if isinstance(first_field, dict):
                iterator = first_field.items()
                first_field_iterator = True
            else:
                try:
                    iterator = iter(first_field)
                    first_field_iterator = True
                except TypeError:
                    first_field_iterator = False

            # if we provided an iterator as first field and the iterator is a (key, value) iterator
            # set the associated fields
            if first_field_iterator:
                for idx, element in enumerate(iterator):
                    if (
                        not isinstance(element, (list, tuple))
                        or not len(element) == 2
                        or not isinstance(element[0], str)
                    ):
                        if idx == 0:
                            # If we do not have an iterator of key/values, set it as attribute
                            self[class_fields[0].name] = first_field
                        else:
                            # If we have a mixed iterator, raise an error
                            raise ValueError(
                                f"Cannot set key/value for {element}. It needs to be a tuple (key, value)."
                            )
                        break
                    setattr(self, element[0], element[1])
                    if element[1] is not None:
                        self[element[0]] = element[1]
            elif first_field is not None:
                self[class_fields[0].name] = first_field
        else:
            for field in class_fields:
                v = getattr(self, field.name)
                if v is not None:
                    self[field.name] = v

    def __delitem__(self, *args, **kwargs):
        raise Exception(f"You cannot use ``__delitem__`` on a {self.__class__.__name__} instance.")

    def setdefault(self, *args, **kwargs):
        raise Exception(f"You cannot use ``setdefault`` on a {self.__class__.__name__} instance.")

    def pop(self, *args, **kwargs):
        raise Exception(f"You cannot use ``pop`` on a {self.__class__.__name__} instance.")

    def update(self, *args, **kwargs):
        raise Exception(f"You cannot use ``update`` on a {self.__class__.__name__} instance.")

    def __getitem__(self, k):
        if isinstance(k, str):
            inner_dict = dict(self.items())
            return inner_dict[k]
        else:
            return self.to_tuple()[k]

    def __setattr__(self, name, value):
        if name in self.keys() and value is not None:
            # Don't call self.__setitem__ to avoid recursion errors
            super().__setitem__(name, value)
        super().__setattr__(name, value)

    def __setitem__(self, key, value):
        # Will raise a KeyException if needed
        super().__setitem__(key, value)
        # Don't call self.__setattr__ to avoid recursion errors
        super().__setattr__(key, value)

    def __reduce__(self):
        if not is_dataclass(self):
            return super().__reduce__()
        callable, _args, *remaining = super().__reduce__()
        args = tuple(getattr(self, field.name) for field in fields(self))
        return callable, args, *remaining

    def to_tuple(self) -> Tuple[Any]:
        """
        Convert self to a tuple containing all the attributes/keys that are not `None`.
        """
        return tuple(self[k] for k in self.keys())


if is_torch_available():
    import torch.utils._pytree as _torch_pytree

    def _model_output_flatten(output: ModelOutput) -> Tuple[List[Any], "_torch_pytree.Context"]:
        return list(output.values()), (type(output), list(output.keys()))

    def _model_output_unflatten(values: Iterable[Any], context: "_torch_pytree.Context") -> ModelOutput:
        output_type, keys = context
        return output_type(**dict(zip(keys, values)))

    _torch_pytree._register_pytree_node(
        ModelOutput,
        _model_output_flatten,
        _model_output_unflatten,
    )


@dataclass
class TokenClassifierOutput(ModelOutput):
    loss: Optional[torch.FloatTensor] = None
    logits: torch.FloatTensor = None
    hidden_states: Optional[Tuple[torch.FloatTensor]] = None
    attentions: Optional[Tuple[torch.FloatTensor]] = None
    embedding_feat: Optional[torch.FloatTensor] = None
    confidence_score: Optional[torch.FloatTensor] = None
    learned_embed: Optional[torch.FloatTensor] = None


@dataclass
class BaseModelOutputWithPoolingAndCrossAttentions(ModelOutput):
    last_hidden_state: torch.FloatTensor = None
    pooler_output: torch.FloatTensor = None
    hidden_states: Optional[Tuple[torch.FloatTensor]] = None
    past_key_values: Optional[Tuple[Tuple[torch.FloatTensor]]] = None
    attentions: Optional[Tuple[torch.FloatTensor]] = None
    cross_attentions: Optional[Tuple[torch.FloatTensor]] = None
    embedding_feat: Optional[torch.FloatTensor] = None


@dataclass
class SequenceClassifierOutput(ModelOutput):
    loss: Optional[torch.FloatTensor] = None
    logits: torch.FloatTensor = None
    hidden_states: Optional[Tuple[torch.FloatTensor]] = None
    attentions: Optional[Tuple[torch.FloatTensor]] = None
    embedded_feat: Optional[torch.FloatTensor] = None


