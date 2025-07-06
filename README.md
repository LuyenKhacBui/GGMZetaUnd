## Getting Started

This Python tool computes height anomaly (ζ, ZETA) or geoid undulation (N, UND) from Global Gravitational Models (GGMs) using fully-normalized spherical harmonic coefficients. Input coordinates must be provided in geographic latitude and longitude.

By default, the tool computes ZETA. If Zeta-to-N conversion coefficients are supplied, it will also derive UND (geoid undulation).

Several high-degree GGM models are listed below as examples for reference and can be downloaded from 
🔗 https://icgem.gfz-potsdam.de/tom_longtime

- EGM2008 (The Earth Gravitational Model 2008) developed and released by the National Geospatial-Intelligence Agency (NGA) 🔗 https://doi.org/10.1029/2011JB008916
- EIGEN-6C4 (European Improved Gravity model of the Earth by New techniques-6C4) 🔗 https://doi.org/10.5880/icgem.2015.1
- SGG-UGM-2 🔗 http://doi.org/10.1016/j.eng.2020.05.008
- XGM2019e 🔗 https://doi.org/10.1007/s00190-020-01398-0

## Prerequisites
- pandas 🔗 https://pypi.org/project/pandas/
- ahrs 🔗 https://ahrs.readthedocs.io/en/latest/installation.html
- multiprocessing 🔗 https://docs.python.org/3/library/multiprocessing.html

## Usages
1. Edit the configuration file located at sample_data/config.cfg to match your own data (latitude/longitude file) and the desired GGM model.
2. Run the script using one of the available execution modes: sequential (a.k.a. serial) or parallel (i.e., multiprocessing).

Use the command below to display help information:

🪟 On Windows using PowerShell or Command Prompt:
```bash
python.exe GGMZetaUnd.py -h
```

🐧 On Linux using the terminal:
```bash
python GGMZetaUnd.py -h
```

To call the script in the sequential (serial) mode, use the following example command:

🪟 On Windows using PowerShell or Command Prompt.
```bash
python.exe GGMZetaUnd.py -g sample_data\eigen-6c4.gfc -i sample_data\LatLon.txt -o sample_data\ZetaUnd.xlsx -c sample_data\config.cfg
```

🐧 On Linux using the terminal:
```bash
python GGMZetaUnd.py -g sample_data/eigen-6c4.gfc -i sample_data/LatLon.txt -o sample_data/ZetaUnd.xlsx -c sample_data/config.cfg
```
To call the script in the parallel (multiprocessing) mode, use the following example command:

🪟 On Windows using PowerShell or Command Prompt.
```bash
python.exe GGMZetaUnd.py -g sample_data\eigen-6c4.gfc -i sample_data\LatLon.txt -o sample_data\ZetaUnd.xlsx -c sample_data\config.cfg -n 4
```

🐧 On Linux using the terminal:
```bash
python GGMZetaUnd.py -g sample_data/eigen-6c4.gfc -i sample_data/LatLon.txt -o sample_data/ZetaUnd.xlsx -c sample_data/config.cfg -n 4
```

## Author
* Luyen Bui: <bkluyen@gmail.com>