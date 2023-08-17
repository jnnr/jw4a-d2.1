rule extend_prebuild:
    input: 
        prebuild="data/{prebuild}",
        overrides="data/overrides"
    output: directory("build/{prebuild}_extended")
    shell: 
        """
        cp -r {input.prebuild} {output};
        cp -r {input.overrides} {output}/overrides
        """