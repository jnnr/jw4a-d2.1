rule assemble_prebuild:
    input:
        prebuild_base=str(rules.run_prebuilt_download_pre_built.output),
        overrides="build/overrides"
    output:
        prebuild_assembled=directory(Path(str(rules.run_prebuilt_download_pre_built.output)).parent / "overrides")
    params:
        import_paths=[
                "'../../../overrides/locations-wind-offshore-deep.yaml'",
                "'../../../overrides/locations-wind-offshore-shallow.yaml'",
                "'../../../overrides/locations-wind-onshore.yaml'",
                "'../../../overrides/techs-novel-wind.yaml'",
            ],
        years=["2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018"]
    run:
        shell("cp -r {input.overrides} {output.prebuild_assembled}")
        for path in params.import_paths:
            for year in params.years:
                shell("""sed -i "2 i \ \ \ \ - {path}" {input.prebuild_base}/model/eurospores/model-{year}.yaml""")
