from snakemake.utils import min_version
from snakemake.io import load_configfile
PANDOC = "pandoc --filter pantable --filter pandoc-crossref --citeproc"

configfile: "config/default.yaml"

module run_prebuilt:
    snakefile: "run-prebuilt-sector-coupled-euro-calliope/Snakefile"
    prefix: "run-prebuilt-sector-coupled-euro-calliope"
    config: load_configfile("config/default.yaml")

use rule * from run_prebuilt as run_prebuilt_*

include: "./rules/boundaries.smk"
include: "./rules/capacity_factors.smk"
include: "./rules/download_data.smk"
include: "./rules/model.smk"
include: "./rules/model_overrides.smk"
include: "./rules/potentials_offshore.smk"
include: "./rules/plot.smk"
include: "./rules/table.smk"
include: "./rules/model_assembly.smk"


min_version("7.8")


onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'jw4a-d2.1 succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'jw4a-d2.1 failed' {config[email]}")

rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.html",
        "build/test-report.html"


rule dag:
     message: "Plot dependency graph of the workflow."
     conda: "envs/dag.yaml"
     shell:
         """
         snakemake --rulegraph > build/dag.dot
         dot -Tpdf -o build/dag.pdf build/dag.dot
         """


rule clean: # removes all generated results
    message: "Remove all build results but keep downloaded data."
    run:
         import shutil

         shutil.rmtree("build")
         print("Data downloaded to data/ has not been cleaned.")


rule test:
    conda: "envs/test.yaml"
    output: "build/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
