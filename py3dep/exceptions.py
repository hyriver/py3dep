"""Customized Py3DEP exceptions."""
from typing import List, Optional


class InvalidInputType(Exception):
    """Exception raised when a function argument type is invalid.

    Parameters
    ----------
    arg : str
        Name of the function argument
    valid_type : str
        The valid type of the argument
    example : str, optional
        An example of a valid form of the argument, defaults to None.
    """

    def __init__(self, arg: str, valid_type: str, example: Optional[str] = None) -> None:
        self.message = f"The {arg} argument should be of type {valid_type}"
        if example is not None:
            self.message += f":\n{example}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MissingColumns(Exception):
    """Exception raised when a required column is missing from a dataframe.

    Parameters
    ----------
    missing : list
        List of missing columns.
    """

    def __init__(self, missing: List[str]) -> None:
        self.message = "The following columns are missing:\n" + f"{', '.join(missing)}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MissingDependency(ImportError):
    """Exception raised when RichDEM is not installed."""

    def __init__(self) -> None:
        self.message = "Depression filling operation uses richdem package which is not installed."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MissingCRS(Exception):
    """Exception raised when input GeoDataFrame is missing CRS."""

    def __init__(self) -> None:
        self.message = "The input GeoDataFrame is missing CRS."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MissingAttribute(Exception):
    """Exception raised for missing attribute.

    Parameters
    ----------
    obj : object
        Object that is missing the attribute.
    attr : str
        Name of the missing attribute
    """

    def __init__(self, obj: str, attr: str) -> None:
        self.message = f"The ``{obj}`` object is missing {attr} attribute."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message
