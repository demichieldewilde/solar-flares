import pytest
import numpy as np
from use_nessi import monotopic_STDs

# Happy path tests with various realistic test values
@pytest.mark.parametrize("test_id, STD_Area, STD_Wavl, Deltas, expected_STD_Area, expected_STD_Wavl", [
    ("HP_1", np.array([1, 2, 3]), np.array([3, 2, 1]), np.array([0, 1, 2]), np.array([1, 2, 3]), np.array([3, 2, 1])),
    ("HP_2", np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 1, 2]), np.array([0, 0, 0]), np.array([0, 0, 0])),
    # Add more happy path cases with different realistic values
])
def test_monotopic_STDs_happy_path(test_id, STD_Area, STD_Wavl, Deltas, expected_STD_Area, expected_STD_Wavl):
    # Act
    result_STD_Area, result_STD_Wavl = monotopic_STDs(STD_Area, STD_Wavl, Deltas, show=False)

    # Assert
    np.testing.assert_array_equal(result_STD_Area, expected_STD_Area)
    np.testing.assert_array_equal(result_STD_Wavl, expected_STD_Wavl)

# Edge cases
@pytest.mark.parametrize("test_id, STD_Area, STD_Wavl, Deltas, expected_STD_Area, expected_STD_Wavl", [
    ("EC_1", np.array([np.nan, np.nan, np.nan]), np.array([np.nan, np.nan, np.nan]), np.array([0, 1, 2]), np.array([0, 0, 0]), np.array([0, 0, 0])),
    ("EC_2", np.array([np.inf, -np.inf, np.nan]), np.array([np.inf, -np.inf, np.nan]), np.array([0, 1, 2]), np.array([0, 0, 0]), np.array([0, 0, 0])),
    # Add more edge cases if necessary
])
def test_monotopic_STDs_edge_cases(test_id, STD_Area, STD_Wavl, Deltas, expected_STD_Area, expected_STD_Wavl):
    # Act
    result_STD_Area, result_STD_Wavl = monotopic_STDs(STD_Area, STD_Wavl, Deltas, show=False)

    # Assert
    np.testing.assert_array_equal(result_STD_Area, expected_STD_Area)
    np.testing.assert_array_equal(result_STD_Wavl, expected_STD_Wavl)

# Error cases
@pytest.mark.parametrize("test_id, STD_Area, STD_Wavl, Deltas, exception", [
    ("ERR_1", "not an array", np.array([1, 2, 3]), np.array([0, 1, 2]), TypeError),
    ("ERR_2", np.array([1, 2, 3]), "not an array", np.array([0, 1, 2]), TypeError),
    ("ERR_3", np.array([1, 2, 3]), np.array([1, 2]), np.array([0, 1, 2]), ValueError),
    # Add more error cases if necessary
])
def test_monotopic_STDs_error_cases(test_id, STD_Area, STD_Wavl, Deltas, exception):
    # Act & Assert
    with pytest.raises(exception):
        monotopic_STDs(STD_Area, STD_Wavl, Deltas, show=False)
