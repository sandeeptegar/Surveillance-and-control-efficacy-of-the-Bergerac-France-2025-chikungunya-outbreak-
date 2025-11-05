# Surveillance-and-control-efficacy-of-the-Bergerac-France-2025-chikungunya-outbreak-
This repository contains simulation code and analysis scripts used to evaluate the surveillance efficiency and control interventions during the 2025 chikungunya outbreak in Bergerac, France. The models simulate vector–host dynamics and disease transmission to assess the impact of different control strategies.

# Instructions to Run the Code

## 1. Chikungunya Outbreak Simulation

**Main code file:** chikv_bergerac.ipynb  
**Supporting function files:** ChikungunyaFunSR.jl and ErrorFun_bergerac.jl 

### Required Data:  
- Data files: `AMgam.csv`, `Car_Tol_1.csv`, `Fin_LSurv.csv`, `LDgam.csv`, `WL_re.csv`  
- Climate file: `latitude_longitude.csv`  
  (For Bergerac: latitude = 44.87, longitude = 0.49)
- Source: [ERA5-Land hourly data from 1950 to present](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview)

### Population Data:  
- Population data per sq. km around outbreak site in the city of Bergerac, France  
- Source: [Eurostat - Population density data]([https://ec.europa.eu/eurostat/statistics-explained/index.php?oldid=596753](https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/population-grids)
