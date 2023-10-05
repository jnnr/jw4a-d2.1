from snakemake.utils import min_version
PANDOC = "pandoc --filter pantable --filter pandoc-crossref --citeproc"

configfile: "config/default.yaml"

include: "./rules/boundaries.smk"
include: "./rules/capacity_factors.smk"
include: "./rules/download_data.smk"
include: "./rules/model.smk"
include: "./rules/model_overrides.smk"
include: "./rules/potentials_offshore.smk"
include: "./rules/plot.smk"
include: "./rules/table.smk"

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


rule run:
    message: "Runs the demo model."
    params:
        slope = config["slope"],
        x0 = config["x0"]
    output: "build/results.pickle"
    conda: "envs/default.yaml"
    script: "scripts/model.py"


rule plot:
    message: "Visualises the demo results."
    input:
        results = rules.run.output
    output: "build/plot.png"
    conda: "envs/default.yaml"
    script: "scripts/vis.py"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--embed-resources --standalone --to html5"
    elif suffix == "pdf":
        return "--pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "report/reset.css",
        "report/report.css",
        rules.plot.output
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


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
