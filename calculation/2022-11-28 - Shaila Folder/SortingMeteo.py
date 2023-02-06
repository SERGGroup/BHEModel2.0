# %%-----------------------------------Path Setting-------------------------------------------------------

from pathlib import Path
import pandas as pd
import numpy as np
from pandas import Series,DataFrame
meteo_data=r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data\Munich ERA.csv'
df = pd.read_csv(meteo_data,usecols=['time','temperature'])
print(df.shape)
df.columns


# %%-----------------------------------Hourly Average-----------------------------------------------------
df['time']=pd.to_datetime(df['time'])
new_df=df[df["time"].between('2018-01-01T00:00','2018-12-31T00:00')]
print(new_df)
mean_df=new_df.groupby([new_df['time'].dt.month, new_df['time'].dt.hour]).temperature.mean()
print(mean_df)
hourly_mean_Munich=pd.DataFrame(mean_df)
filepath = Path(r'C:\Users\af_sh\PycharmProjects\BHEModel2.0\resources\meteo data.csv')
filepath.parent.mkdir(parents=True, exist_ok=True)
hourly_mean_Munich.to_csv(filepath)

