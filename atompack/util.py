from typing import Any, Dict

import numpy as np


class AsDict(object):

    def as_dict(self) -> Dict[str, Any]:
        # collect all public attributes
        res = {k: getattr(self, k) for k in dir(self) if not k.startswith("_")}
        # remove callables
        res = {k: v for k, v in res.items() if not callable(v)}
        # convert all numpy arrays to list format
        res = {k: v.tolist() if type(v) is np.ndarray else v for k, v in res.items()}
        return res
