import numpy as np
import pytest

from stix_parser.products.quicklook import get_energies_from_mask


def test_get_energies_from_mask():
    low, high = get_energies_from_mask()
    assert len(low) == 32
    assert len(high) == 32
    assert low[0] == 0.0
    assert high[0] == 4.0
    assert low[-1] == 150
    assert high[-1] == np.inf

    # Test mask
    mask = np.array([0]*32)
    mask[1:-1] = 1
    low, high = get_energies_from_mask(mask.tolist())
    assert low == 4.0
    assert high == 150.0

    # Test edge mask
    mask = np.array([0]*33)
    mask[[1, 7, 12, 17, 23, -2]] = 1
    low, high = get_energies_from_mask(mask)
    assert len(low) == len(high) == 5
    assert low[0] == 4.0
    assert low[-1] == 50.0
    assert high[0] == 10.0
    assert high[-1] == 150.0

    with pytest.raises(ValueError) as e_info:
        _ = get_energies_from_mask([1])
        assert e_info.message == 'Energy mask or edges must have a lenght of 32 or 33 not 1'
