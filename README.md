# MOSAEC
**M**etal **O**xidation **S**tate **A**utomated **E**rror **C**hecker: an algorithm to automatically detect chemically invalid crystal structures through analysis of their metal oxidation states.

This method operates on crystal structure data containing at least one metal atom defined in the CIF (metal-organic frameworks, coordination polymer, etc.) or MOL2 (metal complex, SBU, etc.) file formats. Validation of the calculated oxidation state and error status flagging outputs was mainly performed on metal-organic frameworks crystal structures.

## Installation
This code utilizes the [CSD Python API](https://www.ccdc.cam.ac.uk/solutions/csd-core/components/csd-python-api/), which must be installed  according to their [installation instructions](https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html).

Additional dependencies:
 - mendeleev
 - pandas

## Example — running MOSAEC
MOSAEC is designed to iterate over all structure files in the working directory. The following data files are also expected in the working directory for proper code operation:
 - Ionization_Energies.csv
 - KnownON.csv
 - Oxidation_Probabilities.csv

**A generic example of running MOSAEC on linux-based systems:**
```
cd <path to directory containing CIF/MOL2 files>
cp <path to MOSAEC directory>/*.csv .
python <path to MOSAEC directory>/mosaec.py
```

## Output
MOSAEC outputs a summary of the oxidation states and error flags in a .csv file. This file contains a separate row for each crystallographically unique metal atom site in the given structures. If **any** of the error flag columns for any metal in a given structure contains the "BAD" or "LIKELY_BAD" designation, the structure is likely to possess issues.

**OxStatesOutput.csv** — Description of Columns:
| Column Header | Description |
| -------------- | -------------- |
| `CIF` | structure file name |
| `Metal` | metal atom label |
| `ON_coordination_ONLY` | oxidation state computed by sharing across binding domains (localized distribution) |
| `EC_coordination_ONLY` | electron count computed by sharing across binding domains (localized distribution) |
| `ON_coordination+Outer_Sphere` | oxidation state computed by sharing across binding domains & distributing outer sphere charge (partially localized distribution) |
| `EC_coordination+Outer_Sphere` | electron count computed by sharing across binding domains & distributing outer sphere charge (partially localized distribution) |
| `ON_network_redistributed` | oxidation state computed by redistributing across metal networks according to IE (localized distribution) |
| `EC_network_redistributed` | electron count computed by redistributing across metal networks according to IE (localized distribution) |
| `ON_network+Outer_Sphere` | oxidation state computed by redistributing across metal networks according to IE & distributing outer sphere charge (partially localized distribution) |
| `EC_network+Outer_Sphere` | electron count computed by redistributing across metal networks according to IE & distributing outer sphere charge (partially localized distribution) |
| `ON_global` | oxidation state computed by global distribution of charges across all sites (delocalized distribution) |
| `EC_global` | electron count computed by global distribution of charges across all sites (delocalized distribution) |
| `Impossible` | **FLAG:** detected impossible oxidation state  i.e. exceeds available valence |
| `Unknown` | **FLAG:** detected unknown oxidation state i.e. never experimentally observed |
| `Zero_Valent` | **FLAG:** detected oxidation state of zero (0) — may be ignored in certain metal complexes & materials |
| `noint_flag` | **FLAG:** detected non-integer oxidation state  |
| `low_prob_1` | **FLAG:** detected oxidation state with a probability < 1% in CSD reporting |
| `low_prob_2` | **FLAG:** detected oxidation state with a probability < 0.1% in CSD reporting |
| `low_prob_3` | **FLAG:** detected oxidation state with a probability < 0.01% in CSD reporting |
| `low_prob_multi` | **FLAG:** detected multiple metal atoms with low probability (< 1%) |
| `high_count` | **FLAG:** detected electron counts of > 32 for f-block and > 20 for non-f-block elements |
| `low_count` | **FLAG:** detected electron counts of < 14  |

## Licensing
The [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) applies to MOSAEC. Follow the license guidelines regarding the use, sharing, adaptation, and attribution of this data.

## Citation
A [preprint version](https://doi.org/10.26434/chemrxiv-2024-ftsv3) of the manuscript outlining the application of the MOSAEC algorithm is available at the following citation. This section will be updated upon publication.

(1) White, A.; Gibaldi, M.; Burner, J.; Woo, T. K. Alarming Structural Error Rates in MOF Databases Used in Computational Screening Identified via a Novel Metal Oxidation State-Based Method. ChemRxiv 2024. https://doi.org/10.26434/chemrxiv-2024-ftsv3.

## Contact
Reach out to any of the following authors with any questions:

Andrew White: awhite2@uottawa.ca

Tom Woo: tom.woo@uottawa.ca