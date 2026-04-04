import os
from multiprocessing import cpu_count
from unittest.mock import patch

import pytest

from prolif.parallel import MAX_JOBS, get_n_jobs


class TestGetNJobs:
    def test_explicit_n_jobs_returned_as_is(self) -> None:
        assert get_n_jobs(2) == 2
        assert get_n_jobs(2, capped=True) == 2
        assert get_n_jobs(2, capped=False) == 2

    def test_env_var_takes_priority_over_cpu_count(self) -> None:
        with patch.dict(os.environ, {"PROLIF_N_JOBS": "4"}):
            assert get_n_jobs(None, capped=True) == 4
            assert get_n_jobs(None, capped=False) == 4

    def test_env_var_zero_returns_none(self) -> None:
        with patch.dict(os.environ, {"PROLIF_N_JOBS": "0"}):
            assert get_n_jobs(None) is None

    def test_capped_limits_to_max_jobs(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            # Simulate many CPUs so the cap kicks in
            with patch("prolif.parallel.cpu_count", return_value=MAX_JOBS + 4):
                result = get_n_jobs(None, capped=True)
        assert result == MAX_JOBS

    def test_uncapped_returns_cpu_count(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            with patch("prolif.parallel.cpu_count", return_value=MAX_JOBS + 4):
                result = get_n_jobs(None, capped=False)
        assert result == MAX_JOBS + 4

    def test_capped_below_max_jobs_uses_cpu_count(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            with patch("prolif.parallel.cpu_count", return_value=2):
                result = get_n_jobs(None, capped=True)
        assert result == 2

    def test_no_env_var_no_explicit_uses_cpu_count(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            result_uncapped = get_n_jobs(None, capped=False)
        assert result_uncapped == cpu_count()
