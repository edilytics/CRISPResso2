"""Tests for CRISPRessoMultiProcessing module."""

import os
import sys

import pandas as pd
import pytest

from CRISPResso2 import CRISPRessoMultiProcessing


# =============================================================================
# Tests for get_max_processes
# =============================================================================


def test_get_max_processes_returns_positive():
    """Test that get_max_processes returns a positive integer."""
    max_procs = CRISPRessoMultiProcessing.get_max_processes()
    assert max_procs >= 1


def test_get_max_processes_returns_int():
    """Test that get_max_processes returns an integer."""
    max_procs = CRISPRessoMultiProcessing.get_max_processes()
    assert isinstance(max_procs, int)


def test_get_max_processes_matches_cpu_count():
    """Test that get_max_processes matches os.cpu_count."""
    max_procs = CRISPRessoMultiProcessing.get_max_processes()
    assert max_procs == os.cpu_count()


# =============================================================================
# Tests for wrapper function
# =============================================================================


def test_wrapper_returns_index_and_result():
    """Test that wrapper returns tuple of (index, result)."""
    def add_one(x):
        return x + 1

    result = CRISPRessoMultiProcessing.wrapper(add_one, (0, 5))
    assert result == (0, 6)


def test_wrapper_preserves_index():
    """Test that wrapper preserves the index in the result."""
    def identity(x):
        return x

    result = CRISPRessoMultiProcessing.wrapper(identity, (42, "test"))
    assert result[0] == 42
    assert result[1] == "test"


def test_wrapper_with_different_indices():
    """Test wrapper with various indices."""
    def double(x):
        return x * 2

    for i in [0, 1, 10, 100]:
        result = CRISPRessoMultiProcessing.wrapper(double, (i, 5))
        assert result[0] == i
        assert result[1] == 10


# =============================================================================
# Tests for run_subprocess
# =============================================================================


def test_run_subprocess_echo():
    """Test run_subprocess with simple echo command."""
    result = CRISPRessoMultiProcessing.run_subprocess("echo hello")
    assert result == 0


def test_run_subprocess_true():
    """Test run_subprocess with true command (always succeeds)."""
    result = CRISPRessoMultiProcessing.run_subprocess("true")
    assert result == 0


def test_run_subprocess_false():
    """Test run_subprocess with false command (always fails)."""
    result = CRISPRessoMultiProcessing.run_subprocess("false")
    assert result != 0


def test_run_subprocess_invalid_command():
    """Test run_subprocess with invalid command returns non-zero."""
    result = CRISPRessoMultiProcessing.run_subprocess(
        "nonexistent_command_xyz_12345 2>/dev/null"
    )
    assert result != 0


def test_run_subprocess_exit_code():
    """Test run_subprocess captures specific exit codes."""
    # sh -c 'exit 5' will return exit code 5
    result = CRISPRessoMultiProcessing.run_subprocess("sh -c 'exit 5'")
    assert result == 5


# =============================================================================
# Tests for run_function_on_array_chunk_parallel
# =============================================================================


def test_run_function_on_array_chunk_parallel_single_process():
    """Test parallel function execution with single process."""
    def square_all(arr):
        return [x ** 2 for x in arr]

    input_array = [1, 2, 3, 4, 5]
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        input_array, square_all, n_processes=1
    )
    assert result == [1, 4, 9, 16, 25]


def test_run_function_on_array_chunk_parallel_empty_array():
    """Test parallel function with empty array."""
    def identity(arr):
        return arr

    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        [], identity, n_processes=1
    )
    assert result == []


def test_run_function_on_array_chunk_parallel_preserves_order_single():
    """Test that single process preserves order."""
    def add_index(arr):
        return [f"item_{x}" for x in arr]

    input_array = list(range(10))
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        input_array, add_index, n_processes=1
    )
    expected = [f"item_{i}" for i in range(10)]
    assert result == expected


def test_run_function_on_array_chunk_parallel_large_array():
    """Test parallel function with larger array."""
    def double_all(arr):
        return [x * 2 for x in arr]

    input_array = list(range(100))
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        input_array, double_all, n_processes=1
    )
    expected = [x * 2 for x in range(100)]
    assert result == expected


def test_run_function_on_array_chunk_parallel_with_strings():
    """Test parallel function with string data."""
    def uppercase_all(arr):
        return [s.upper() for s in arr]

    input_array = ["hello", "world", "test"]
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        input_array, uppercase_all, n_processes=1
    )
    assert result == ["HELLO", "WORLD", "TEST"]


# =============================================================================
# Tests for run_pandas_apply_parallel
# =============================================================================


def test_run_pandas_apply_parallel_single_process():
    """Test pandas parallel apply with single process."""
    def double_column(df_chunk):
        df_chunk["doubled"] = df_chunk["value"] * 2
        return df_chunk

    df = pd.DataFrame({"value": [1, 2, 3, 4, 5]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, double_column, n_processes=1
    )

    assert "doubled" in result.columns
    assert list(result["doubled"]) == [2, 4, 6, 8, 10]


def test_run_pandas_apply_parallel_preserves_columns():
    """Test that parallel apply preserves existing columns."""
    def add_column(df_chunk):
        df_chunk["new_col"] = "added"
        return df_chunk

    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, add_column, n_processes=1
    )

    assert "a" in result.columns
    assert "b" in result.columns
    assert "new_col" in result.columns


def test_run_pandas_apply_parallel_two_rows():
    """Test parallel apply with two-row dataframe."""
    def identity(df_chunk):
        return df_chunk

    df = pd.DataFrame({"a": [1, 2]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, identity, n_processes=1
    )

    assert len(result) == 2


def test_run_pandas_apply_parallel_single_row():
    """Test parallel apply with single row dataframe."""
    def add_computed(df_chunk):
        df_chunk["computed"] = df_chunk["x"] + df_chunk["y"]
        return df_chunk

    df = pd.DataFrame({"x": [5], "y": [3]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, add_computed, n_processes=1
    )

    assert len(result) == 1
    assert result["computed"].iloc[0] == 8


def test_run_pandas_apply_parallel_transform_values():
    """Test parallel apply transforming values."""
    def square_values(df_chunk):
        df_chunk["squared"] = df_chunk["num"] ** 2
        return df_chunk

    df = pd.DataFrame({"num": [1, 2, 3, 4, 5]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, square_values, n_processes=1
    )

    expected = [1, 4, 9, 16, 25]
    assert list(result["squared"]) == expected


# =============================================================================
# Tests for run_parallel_commands
# =============================================================================


def test_run_parallel_commands_single_process_success():
    """Test parallel commands with single process and successful commands."""
    commands = ["true", "true", "true"]
    # Should not raise
    CRISPRessoMultiProcessing.run_parallel_commands(
        commands, n_processes=1, descriptor="test"
    )


def test_run_parallel_commands_single_process_failure():
    """Test parallel commands with single process and failing command."""
    commands = ["true", "false", "true"]
    with pytest.raises(Exception, match="was failed"):
        CRISPRessoMultiProcessing.run_parallel_commands(
            commands, n_processes=1, descriptor="test"
        )


def test_run_parallel_commands_continue_on_fail():
    """Test parallel commands continues on failure when flag is set."""
    commands = ["true", "false", "true"]
    # Should not raise
    CRISPRessoMultiProcessing.run_parallel_commands(
        commands, n_processes=1, descriptor="test", continue_on_fail=True
    )


def test_run_parallel_commands_empty_list():
    """Test parallel commands with empty list."""
    commands = []
    # Should not raise and return immediately
    CRISPRessoMultiProcessing.run_parallel_commands(
        commands, n_processes=1, descriptor="test"
    )


def test_run_parallel_commands_single_command():
    """Test parallel commands with single command."""
    commands = ["echo hello"]
    CRISPRessoMultiProcessing.run_parallel_commands(
        commands, n_processes=1, descriptor="test"
    )


# =============================================================================
# Tests for run_crispresso_cmds with empty input
# =============================================================================


def test_run_crispresso_cmds_empty_list():
    """Test run_crispresso_cmds with empty list returns immediately."""
    # Should return without doing anything
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=[],
        n_processes="1",
        descriptor="test",
    )


# =============================================================================
# Edge case tests
# =============================================================================


def test_wrapper_with_none_result():
    """Test wrapper handles None result from function."""
    def return_none(x):
        return None

    result = CRISPRessoMultiProcessing.wrapper(return_none, (0, "anything"))
    assert result == (0, None)


def test_wrapper_with_complex_return():
    """Test wrapper handles complex return types."""
    def return_dict(x):
        return {"input": x, "processed": True}

    result = CRISPRessoMultiProcessing.wrapper(return_dict, (5, "data"))
    assert result[0] == 5
    assert result[1] == {"input": "data", "processed": True}


def test_run_function_on_array_chunk_parallel_exception_handling():
    """Test that exceptions in the function are propagated."""
    def raise_error(arr):
        raise ValueError("Test error")

    with pytest.raises(ValueError, match="Test error"):
        CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
            [1, 2, 3], raise_error, n_processes=1
        )
