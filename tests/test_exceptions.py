import pytest

from py3dep import InvalidInputType


def invalid_type():
    print(InvalidInputType("coords", "tuple", "(lon, lat)"))
    raise InvalidInputType("coords", "tuple", "(lon, lat)")


def test_invalid_type():
    with pytest.raises(InvalidInputType):
        invalid_type()
