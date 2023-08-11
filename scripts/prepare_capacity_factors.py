import pandas as pd
import netCDF4 as nc


def calculate_capacity_factors(windspeeds, powercurve):
    r"""
    Calculate capacity factors from wind speed data and a power curve.

    Parameters
    ----------
    windspeeds : :class:`netCDF4.Dataset`
        Wind speed data with a datetime index.
    powercurve : :class:`pandas.DataFrame`

    Returns
    -------
    capacity_factors : :class:`pandas.DataFrame`
        Capacity factors with a datetime index.
    """
    datetimeindex=pd.date_range("2013-01-01", periods= 52584, freq="H")

    df_list = []
    windspeed_avg_location={}

    for loc, cou, fil in zip(location, country, filename):
        data= nc.Dataset(output_directory + fil)

        u = data.variables['u'][:, :, :]
        v = data.variables['v'][:, :, :]

        # pythagoras
        u=u**2
        v=v**2
        windspeed=np.sqrt(u+v)

        # TODO: Ask Hidde about this.
        windspeed=windspeed.reshape(52584,4)
        windspeed=pd.DataFrame(windspeed)
        windspeed['windspeed']=windspeed.mean(axis=1)
        windspeed=windspeed['windspeed']
        windspeed=pd.DataFrame(windspeed)

        # round to 1 decimal to be able to merge with power curve
        windspeed=np.round(windspeed,decimals=1)
    
        # TODO: Have a look at this variable to see what is the use of this. 
        windspeedavg=windspeed.mean(axis=0)
        windspeed_avg_location[loc] = windspeedavg['windspeed']
    
        windspeed['location']=loc
        windspeed['country']=cou
        windspeed['time']=datetimeindex
        df_list.append(windspeed)
        df = pd.concat(df_list, ignore_index=False)

    # TODO: What happens here?
    capfactable=df.merge(power_curve, how='left', on='windspeed')
    capfactable=capfactable.drop(['windspeed','power'], axis=1)
    capfactable=capfactable.groupby(["country","time"])["capfac"].mean()
    capfactable=pd.DataFrame(capfactable)
    capfactable=capfactable.unstack(level=0)
    capfactable.columns=capfactable.columns.droplevel()

    # TODO: Generalize or drop this.
    #rename column to match North Sea Calliope 
    capfactable=capfactable.rename(columns={"Denmark":"DNK","France":"FRA","Ireland":"IRL","Norway":"NOR","Sweden":"SWE","UK":"GBR"}) #floating only
    capfactable=capfactable.rename(columns={"Netherlands":"NLD","Germany":"DEU","Belgium":"BEL","Luxembourg":"LUX"}) #fixed offshore

    return capfactable


if __name__ == "__main__":
    filepath_windspeeds = snakemake.inputs.windspeeds
    filepath_powercurve = snakemake.outputs.powercurve
    filepath_capacityfactors = snakemake.outputs.capacityfactors

    windspeeds = nc.Dataset(filepath_windspeeds)
    powercurve = pd.read_csv(filepath_powercurve)
    capacity_factors = calculate_capacity_factors(windspeeds, powercurve)

    capacity_factors.to_csv(filepath_capacityfactors)



    
