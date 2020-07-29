from typing import Union


def params_except_keys(in_params: dict, except_keys: Union[str, list]) -> dict:
    """
    Remove specific item in dict by its key(s)
    :param in_params: input dict
    :param except_keys: key(s) value
    :return: dict
    """
    if isinstance(except_keys, str):
        except_keys = [except_keys]

    return {k: v for k, v in in_params.items() if k not in except_keys}
