from snakemake.utils import min_version
from snakemake.io import load_configfile
PANDOC = "pandoc --filter pantable --filter pandoc-crossref --citeproc"

configfile: "config/default.yaml"

wildcard_constraints:
    year = "2010|2011|2012|2013|2014|2015|2016|2017|2018",
    model_resolution = "res_1h|res_3h|res_6h|res_12h",
    scenario = "|".join(config["scenario"].keys())

module run_prebuilt:
    snakefile: "run-prebuilt-sector-coupled-euro-calliope/Snakefile"
    prefix: "run-prebuilt-sector-coupled-euro-calliope"
    config: load_configfile("config/default.yaml")

use rule * from run_prebuilt as run_prebuilt_*

include: "./rules/boundaries.smk"
include: "./rules/capacity_factors.smk"
include: "./rules/download_data.smk"
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
        rules.unzip_geodata_offshore.output,
        rules.map_eez_to_eurospores.output,
        rules.build_availabilitymatrix.output,
        expand(rules.prepare_capacity_factors.output, tech=config["prepare_capacity_factors"].keys()),
        # TODO Potentials onshore?
        rules.draw_map.output,
        expand(rules.draw_plots_capacityfactors.output[0], tech=config["prepare_capacity_factors"].keys()),
        rules.table_area_potential_offshore.output[0],
        rules.assemble_prebuild.output[0],
        expand(rules.plot_model_results.output, year="2016", resolution="eurospores", model_resolution="res_3h", scenario="noveltech"),
        expand("run-prebuilt-sector-coupled-euro-calliope/build/eurospores/outputs/{year}_{model_resolution}.nc", year="2016", model_resolution="res_3h")
        rules.draw_boxplot_capacityfactors.output,
        rules.plot_capacity_factor_distribution.output,

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
