# SageMath code

## Introduction

This is a repository with the code that produces the
table in Chapter 12 of the paper

[NCCRs of cones over del Pezzo surfaces](https://arxiv.org/abs/2604.11319)

by Anya Nordskova and Michel Van den Bergh.

## Running the code

The code assumes that the [SageMath](https://www.sagemath.org) package is
installed.

If this is the case then two additional required Python packages should be installed
using the command:

```sh
sage_cmd --python -m pip install -r requirements.txt --upgrade
```

where here and below `sage_cmd` stands for the command that invokes `SageMath`.

## Internal consistency of the table

The table in the paper is automatically generated from the file
`minimal_collections.txt`.

Verifying the internal consistency of this file can be done with the
command:

```sh
sage_cmd validate_minimal_collections.py
```

The following checks are performed:

- `vtjson` is used to check that there is no missing data.
- The case identifiers are consecutive lists that have the correct length.
- The exceptional collections are really exceptional.
- The exceptional collections have the required Gram matrix.
- The exceptional collections are minimal.
- The exceptional collections are very strong.
- Paths between cases in the relation sections are valid.
- The relations graph for each del Pezzo is a tree.
- The quivers are correct.
- The ranks are correct.
- There is a certificate for each surface.
- Certificates are valid. Either:
   - The symmetry group is trivial or
   - all reflections give the same NCCR or
   - a valid certificate is supplied.
- All Gram matrices are distinct up to rotation.
- All wild card 3 block cases are valid.

## Generation scripts

## Minimal block-complete 4-block collections, case 1

The script

```sh
sage_cmd 4blocks_case_1.py
```

implements the algorithm 13.5.1 in the paper. It finishes
instantaneously and prints out two Gram matrices.

## Minimal block-complete 4-block collections, case 2

The script

```sh
sage_cmd 4blocks_case_2.py
```

implements the algorithm 13.5.2 in the paper. It finishes
instantaneously and prints out seven Gram matrices.

### Non-existence of minimal block-complete 5-block collections

The script

```sh
sage_cmd 5blocks_case_2.py
```

implements the algorithm 13.6.2 in the paper. It takes about
2000 seconds on a 4 year old laptop and ends with an empty file
`5blocks_case_2_data/5blocks_case_2_stage_3.txt`,
meaning that no solutions were found (as expected).

## Testing

The directory `test_data` contains data files with lists of very
strong exceptional collections. We use these for our own research.

The script

```
sage_cmd testing.py
```

applies a variant of the algorithm in the proof of Theorem 11.2 to samples
of exceptional collections extracted from these files and checks
that the resulting minimal block-complete collections are in the table.
