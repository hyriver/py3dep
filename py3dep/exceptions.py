"""Customized Py3DEP exceptions."""
from __future__ import annotations

import async_retriever as ar
import pygeoogc as ogc


class InputTypeError(ar.InputTypeError):
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


class InputValueError(ar.InputValueError):
    """Exception raised for invalid input.

    Parameters
    ----------
    inp : str
        Name of the input parameter
    valid_inputs : tuple
        List of valid inputs
    """


class MissingColumnError(Exception):
    """Exception raised when a required column is missing from a dataframe.

    Parameters
    ----------
    missing : list
        List of missing columns.
    """

    def __init__(self, missing: list[str]) -> None:
        self.message = "The following columns are missing:\n" + f"{', '.join(missing)}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class DependencyError(ImportError):
    """Exception raised when RichDEM is not installed."""

    def __init__(self) -> None:
        self.message = "Depression filling operation uses richdem package which is not installed."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MissingCRSError(Exception):
    """Exception raised when input is missing CRS."""

    def __init__(self) -> None:
        self.message = "The input is missing CRS."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ServiceUnavailableError(ogc.ServiceUnavailableError):
    """Exception raised when the service is not available.

    Parameters
    ----------
    url : str
        The server url
    """
