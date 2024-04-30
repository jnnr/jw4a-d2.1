# JW4A-D2.1

Workflow to run sector-coupled-eurocalliope with new wind regions and techs.

The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) or [conda](https://docs.conda.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml --no-default-packages

## Run the analysis

    snakemake

This will run all analysis steps to reproduce results.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake dag


## Repo structure

* `scripts`: contains the Python source code as scripts
* `rules`: contains Snakemake rule definitions
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `profiles`: Snakemake execution profiles
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
