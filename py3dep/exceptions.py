"""Customized Py3DEP exceptions."""
from typing import List


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


class MissingOption(Exception):
    """Exception raised when layer is not provided."""

    def __init__(self) -> None:
        self.message = "layer option is required when target_type is geometry."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message
