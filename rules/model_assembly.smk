rule assemble_prebuild:
    input:
        prebuild_base=str(rules.run_prebuilt_download_pre_built.output),
        overrides="build/overrides"
    output:
        prebuild_assembled=directory(Path(str(rules.run_prebuilt_download_pre_built.output)).parent / "overrides")
    shell: "cp -r {input.overrides} {output.prebuild_assembled}"
