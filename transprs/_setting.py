from typing import Any, Union, Optional, Iterable, TextIO
from typing import Tuple, List, ContextManager


class TransPRSConfig:
    """\
    Config manager for trans_ethnic_prs.
    """

    def __init__(
        self,
        datasetdir: Union[str, Path] = "./data/",
    ):
        self.datasetdir = datasetdir

    @property
    def datasetdir(self) -> Path:
        """\
        Directory for example :mod:`~scanpy.datasets` (default `'./data/'`).
        """
        return self._datasetdir

    @datasetdir.setter
    def datasetdir(self, datasetdir: Union[str, Path]):
        _type_check(datasetdir, "datasetdir", (str, Path))
        self._datasetdir = Path(datasetdir).resolve()


def _type_check(var: Any, varname: str, types: Union[type, Tuple[type, ...]]):
    if isinstance(var, types):
        return
    if isinstance(types, type):
        possible_types_str = types.__name__
    else:
        type_names = [t.__name__ for t in types]
        possible_types_str = "{} or {}".format(
            ", ".join(type_names[:-1]), type_names[-1]
        )
    raise TypeError(f"{varname} must be of type {possible_types_str}")
