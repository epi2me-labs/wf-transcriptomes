"""Check if a sample sheet is valid.

Loads validator classes from modules in the validators package to check various
aspects of a sample sheet. If any validator fails, the sample sheet is invalid.
"""
import codecs
import csv
from importlib import import_module
from inspect import getmembers, isclass
import json
import os
from pkgutil import iter_modules
import sys

from ..util import get_named_logger, wf_parser  # noqa: ABS101
from . import validators  # noqa: ABS101


# Some Excel users save their CSV as UTF-8 (and occasionally for a reason beyond my
# comprehension, UTF-16); Excel then adds a byte order mark (unnecessarily for UTF-8
# I should add). If we do not handle this with the correct encoding, the mark will
# appear in the parsed data, causing the header to be malformed.
# See CW-2310
def determine_codec(f):
    """Peek at a file and return an appropriate reading codec."""
    with open(f, 'rb') as f_bytes:
        # Could use chardet here if we need to expand codec support
        initial_bytes = f_bytes.read(8)

        for codec, encoding_name in [
            [codecs.BOM_UTF8, "utf-8-sig"],  # use the -sig codec to drop the mark
            [codecs.BOM_UTF16_BE, "utf-16"],  # don't specify LE or BE to drop mark
            [codecs.BOM_UTF16_LE, "utf-16"],
            [codecs.BOM_UTF32_BE, "utf-32"],  # handle 32 for completeness
            [codecs.BOM_UTF32_LE, "utf-32"],  # again skip LE or BE to drop mark
        ]:
            if initial_bytes.startswith(codec):
                return encoding_name
        return None  # will cause file to be opened with default encoding


def load_validators(modules, wf_params, options=None):
    """Load validator classes."""
    validator_classes = []
    for mod in modules:
        for _, validator_class in getmembers(mod, _is_sample_sheet_validator):
            validator_classes.append(validator_class(wf_params, options))
    return validator_classes


def load_validator_modules():
    """Load all modules from the validators package."""
    return [
        import_module(f"{validators.__name__}.{module_info.name}")
        for module_info in iter_modules(validators.__path__)
    ]


def _is_sample_sheet_validator(obj):
    """Return whether an object is a concrete sample sheet validator class.

    `issubclass(cls, Base)` is true when `cls is Base`, so the base validator
    class must be excluded explicitly.
    """
    return (
        isclass(obj) and issubclass(obj, validators.SampleSheetValidator)
        and obj is not validators.SampleSheetValidator
    )


def main(args):
    """Run the sample sheet checks."""
    logger = get_named_logger("checkSheet")
    with open(args.wf_params_path) as f:
        wf_params = json.load(f)
    if args.required_sample_types:
        wf_params['required_sample_types'] = args.required_sample_types

    validator_classes = load_validators(
        load_validator_modules(), wf_params, options={
            "no_barcode": args.no_barcode,
        }
    )

    rows = []

    if not os.path.exists(args.sample_sheet) or not os.path.isfile(args.sample_sheet):
        sys.stdout.write("Could not open sample sheet file.")
        sys.exit()

    try:
        encoding = determine_codec(args.sample_sheet)
        with open(args.sample_sheet, "r", encoding=encoding) as f:
            try:
                # Excel files don't throw any error until here
                csv.Sniffer().sniff(f.readline())
                f.seek(0)  # return to initial position again
            except Exception as e:
                # Excel fails with UniCode error
                sys.stdout.write(
                    "The sample sheet doesn't seem to be a CSV file.\n"
                    "The sample sheet has to be a CSV file.\n"
                    "Please verify that the sample sheet is a CSV file.\n"
                    f"Parsing error: {e}"
                 )

                sys.exit()
            csv_reader = csv.DictReader(f)
            columns = csv_reader.fieldnames
            rows = list(csv_reader)

    except Exception as e:
        sys.stdout.write(f"Parsing error: {e}")
        sys.exit()

    # Run all the validators.
    for v in validator_classes:
        v.on_header(columns)
    # Skip header row
    for lineno, row in enumerate(rows, start=1):
        for v in validator_classes:
            v.add_sheet_row(row, lineno)

    if not all((v.is_valid for v in validator_classes)):
        for v in validator_classes:
            for e in v.errors:
                sys.stdout.write(f"{e}\n")
        sys.exit()

    logger.info(f"Checked sample sheet {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_sample_sheet")
    parser.add_argument("sample_sheet", help="Sample sheet path to check")
    parser.add_argument("wf_params_path", help="Path to WF params JSON")
    parser.add_argument(
        "--required_sample_types",
        help="List of required sample types. Each sample type provided must "
             "appear at least once in the sample sheet",
        nargs="*"
    )
    parser.add_argument(
        "--no_barcode",
        action="store_true",
        help="Allow sample sheets without a barcode column "
        "and match rows by sample_name only",
    )
    return parser
