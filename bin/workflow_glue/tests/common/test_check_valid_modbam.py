"""Tests for modified-base BAM validation."""

import os

import pytest
from workflow_glue import check_valid_modbam


class FakeAlignment:
    """Minimal alignment stub exposing pysam's get_tags API."""

    def __init__(self, tags):
        """Store the synthetic SAM tags returned by ``get_tags``."""
        self._tags = tags

    def get_tags(self):
        """Return the synthetic tag list for this fake alignment."""
        return self._tags


def _args(*argv):
    return check_valid_modbam.argparser().parse_args(list(argv))


def test_main_accepts_bam_with_mm_and_ml_tags(monkeypatch):
    """A read carrying both MM and ML tags should pass validation."""
    monkeypatch.setattr(
        check_valid_modbam.pysam,
        "AlignmentFile",
        lambda _: [
            FakeAlignment([("MM", "A+a.,0;"), ("ML", [255])]),
        ],
    )

    check_valid_modbam.main(_args("input.bam"))


def test_main_rejects_bam_without_mod_tags(monkeypatch):
    """A BAM with no modified-base tags should exit with EX_DATAERR."""
    monkeypatch.setattr(
        check_valid_modbam.pysam,
        "AlignmentFile",
        lambda _: [
            FakeAlignment([]),
            FakeAlignment([("NM", 0)]),
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        check_valid_modbam.main(_args("input.bam"))

    assert exc_info.value.code == os.EX_DATAERR


def test_main_accepts_when_mod_tags_appear_later_in_scan(monkeypatch):
    """Validation should continue scanning until it finds a tagged read."""
    monkeypatch.setattr(
        check_valid_modbam.pysam,
        "AlignmentFile",
        lambda _: [
            FakeAlignment([("NM", 0)]),
            FakeAlignment([("mm", "C+m.,0;"), ("ml", [200])]),
        ],
    )

    check_valid_modbam.main(_args("input.bam"))
