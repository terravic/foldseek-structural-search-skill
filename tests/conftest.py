# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Shared pytest fixtures for foldseek-structural-search test suite."""

from __future__ import annotations

from pathlib import Path
import sys
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


@pytest.fixture()
def project_root() -> Path:
    """Return resolved path to the skill repository root."""
    return PROJECT_ROOT


@pytest.fixture()
def sample_pdb_path() -> Path:
    """Return path to the 100% synthetic sample_query.pdb fixture."""
    return PROJECT_ROOT / "sample_data" / "sample_query.pdb"


@pytest.fixture()
def sample_cif_path() -> Path:
    """Return path to the 100% synthetic sample_query.cif fixture."""
    return PROJECT_ROOT / "sample_data" / "sample_query.cif"


@pytest.fixture()
def sample_results_path() -> Path:
    """Return path to the synthetic sample_foldseek_results.json fixture."""
    return PROJECT_ROOT / "sample_data" / "sample_foldseek_results.json"


@pytest.fixture()
def isolated_workspace(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Provide an isolated workspace directory for testing license file creation."""
    ws = tmp_path / "workspace"
    ws.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("FOLDSEEK_WORKSPACE_ROOT", str(ws))
    return ws
