"""Scrape versions of conda packages."""

from collections import namedtuple
import subprocess


try:
    import pandas as pd
except ImportError:
    pass


PackageInfo = namedtuple(
    'PackageInfo', ('Name', 'Version', 'Build', 'Channel'))


def scrape_data(as_dataframe=False, include=None):
    """Return versions of conda packages in base environment."""
    cmd = """
. ~/conda/etc/profile.d/mamba.sh;
micromamba activate;
micromamba list;
    """
    proc = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    versions = dict()
    for line in proc.stdout.splitlines()[3:]:
        items = line.decode().strip().split()
        if len(items) == 3:
            # sometimes channel isn't listed :/
            items.append("")
        if include is None or items[0] in include:
            versions[items[0]] = PackageInfo(*items)
    if as_dataframe:
        versions = pd.DataFrame.from_records(
            list(versions.values()),
            columns=PackageInfo._fields)
    return versions
