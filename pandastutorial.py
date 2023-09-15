import pandas as pd

print("hello")
df = pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_daily_reports_us/09-15-2021.csv')
df = df.set_index('Province_State')
print(df)