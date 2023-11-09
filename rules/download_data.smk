from lib.download import download, unzip
rule download_geodata_offshore:
    params:
        water_depth="https://www.bodc.ac.uk/data/open_download/gebco/gebco_2023_sub_ice_topo/zip/",
        shipping_density_global="https://datacatalogfiles.worldbank.org/ddh-published/0037580/DR0045406/shipdensity_global.zip",
        natura2000_areas="https://sdi.eea.europa.eu/datashare/s/7WkpGo6K2jwFHGp/download"
    output:
        water_depth="data/potentials_offshore/gebco_2023_sub_ice_topo.zip",
        shipping_density_global="data/potentials_offshore/shipdensity_global.zip",
        natura2000_areas="data/potentials_offshore/natura2000_areas.zip"
    run: 
        download(params.water_depth, output.water_depth)
        download(params.shipping_density_global, output.shipping_density_global)
        download(params.natura2000_areas, output.natura2000_areas)

rule unzip_geodata_offshore:
    input: [item for item in rules.download_geodata_offshore.output if item.endswith(".zip")]
    output: [directory(os.path.splitext(item)[0]) for item in rules.download_geodata_offshore.output if item.endswith(".zip")] 
    run:
        for zipped, unzipped in zip(*{input}, *{output}):
            print(f"Unzipping {zipped} to {unzipped}")
            unzip(str(zipped), str(unzipped))

rule download_eez:
    # copied from euro-calliope
    message: "Download Exclusive Economic Zones as zip"
    output: protected("data/shapes/eez.zip")
    params: url = config["data-sources"]["eez"]
    # conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"

rule download_ERA5_cutout_modellevel:
    output: 
        target_dir="build/cutouts/cutout-era5-model-level.nc"
    script: "../scripts/download_cutout_era5_model_level.py"

rule download_ERA5_cutout_conventional:
    output:
        target_dir="build/cutouts/cutout-era5.nc"
    script: "../scripts/download_cutout_era5.py"
