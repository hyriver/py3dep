"""Customized Py3DEP exceptions."""

from __future__ import annotations

import async_retriever.exceptions as ar
import pygeoogc.exceptions as ogc


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
    given : str, optional
        The given input, defaults to None.
    """


class MissingColumnError(Exception):
    """Exception raised when a required column is missing from a dataframe.

    Parameters
    ----------
    missing : list
        List of missing columns.
    """

    def __init__(self, missing: list[str]) -> None:
        self.message = f"The following columns are missing:\n{', '.join(missing)}"
        super().__init__(self.message)

    def __str__(self) -> str:
        """Return the error message."""
        return self.message


class NoOutletError(Exception):
    """Exception raised when no outlet is found in a DEM."""

    def __init__(self) -> None:
        self.message = "No outlet was foundm, please provide an outlet point or adjust ``elv_max``."
        super().__init__(self.message)

    def __str__(self) -> str:
        """Return the error message."""
        return self.message


class MissingCRSError(Exception):
    """Exception raised when input is missing CRS."""

    def __init__(self) -> None:
        self.message = "The input is missing CRS."
        super().__init__(self.message)

    def __str__(self) -> str:
        """Return the error message."""
        return self.message


class ServiceUnavailableError(ogc.ServiceUnavailableError):
    """Exception raised when the service is not available.

    Parameters
    ----------
    url : str
        The server url
    """


class InputRangeError(Exception):
    """Exception raised when input is out of range.

    Parameters
    ----------
    inp : str
        Name of the input parameter
    valid_range : str
        The valid range of the input
    """

    def __init__(self, inp: str, valid_range: str) -> None:
        self.message = f"{inp} is out of range. Valid range is {valid_range}"
        super().__init__(self.message)

    def __str__(self) -> str:
        """Return the error message."""
        return self.message
