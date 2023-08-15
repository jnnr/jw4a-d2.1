from lib.download import download, unzip
rule download_geodata_offshore:
    params:
        # eez="",
        # water_depth="",
        shipping_density_global="https://datacatalogfiles.worldbank.org/ddh-published/0037580/DR0045406/shipdensity_global.zip"
        #natura2000_areas="https://sdi.eea.europa.eu/datashare/s/7WkpGo6K2jwFHGp/download"
    output:
        # water_depth="data/potentials_offshore/water_depth.tif",
        shipping_density_global="data/potentials_offshore/shipping_density_global.zip"
        #natura2000_areas="data/potentials_offshore/natura2000_areas.zip",
    run: 
        # download(params.water_depth, output.water_depth)
        download(params.shipping_density_global, output.shipping_density_global)
        #download(params.natura2000_areas, output.natura2000_areas)

rule unzip_geodata_offshore:
    input: [item for item in rules.download_geodata_offshore.output if item.endswith(".zip")]
    output: [directory(os.path.splitext(item)[0]) for item in rules.download_geodata_offshore.output if item.endswith(".zip")] 
    run:
        for zipped, unzipped in zip({input}, {output}):
            print(f"Unzipping {zipped} to {unzipped}")
            unzip(str(zipped), str(unzipped))
            os.remove(str(zipped))
