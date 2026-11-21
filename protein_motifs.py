"""Sequence patterns verified against FeGenie and MHCscan; no activity calls.

Patterns are independently implemented, not a vendored copy of either program.
Counts deliberately use non-overlapping matches per pattern, as both sources do.
"""
import re

SOURCES = {
    'FeGenie': {'repository': 'https://github.com/Middle-Author-Bioinformatics/FeGenie',
                'commit': '034f554fb12528b43f7515b5c066de610a6ec6e8',
                'file': 'FeGenie.py', 'lines': '2474-2529', 'license': 'AGPL-3.0'},
    'MHCscan': {'repository': 'https://github.com/Middle-Author-Bioinformatics/MHCscan',
                'commit': '2fe42012d46f3c57281b856331dc51e3777b1872',
                'file': 'bin/cytoscan.py', 'lines': '376-380', 'license': 'GPL-3.0'},
}
PATTERNS = {
    'cxxch': r'C..CH', 'cxxxch': r'C...CH', 'cxxxxch': r'C....CH',
    'cx14ch': r'C.{14}CH', 'cx15ch': r'C.{15}CH',
    'heme_b': r'G.[HR].C[PLAV]G', 'hematite': r'[STC][AVILMFYWH][ST]P[ST]',
}
COUNT_FIELDS = [f'{name}_count' for name in PATTERNS] + ['mhcscan_hemes', 'fegenie_heme_c_total']


def characterize(sequence):
    """Return counts and 1-based inclusive positions for all pattern occurrences."""
    seq = str(sequence).upper()
    matches = [dict(motif=name, regex=regex, start=m.start()+1, end=m.end(), match=m.group())
               for name, regex in PATTERNS.items() for m in re.finditer(regex, seq)]
    counts = {f'{name}_count': sum(m['motif']==name for m in matches) for name in PATTERNS}
    counts['mhcscan_hemes'] = counts['cxxch_count'] + counts['cxxxch_count']
    counts['fegenie_heme_c_total'] = sum(counts[f'{name}_count'] for name in list(PATTERNS)[:5])
    return counts, matches
