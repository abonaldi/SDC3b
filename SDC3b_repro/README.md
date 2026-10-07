# Dependencies
See SDC3b_repro.yml and SDC3b_repro_packages.txt for the environment used to run the simulations and the full list pf packages, respectively.

# Pipeline

1) run 21cmFAST to generate the density and emissivity fields, and the 21cm signal for 21cmFAST (Simulator 2)   
3) if using, run pyC2ray to generate the 21cm signal for pyc2ray (Simulator 1)
4) produce physical lightcones from the outputs of the simulators with run_lightcones.py
5) produce observed lightcones with run_obslc.py

# Full calling sequence

# PS1:
python run_py21cmfast_allmodels.py params/parameters_fullcubes.yml 68

python run_pyc2ray_allmodels.py params/parameters_fullcubes.yml 68

python run_lightcones.py params/parameters_fullcubes.yml_model88 with do_C2ray=true

python run_obslc.py params/parameters_fullcubes.yml_model88 with do_C2ray=true

# PS2: 
python run_py21cmfast_allmodels.py params/parameters_fullcubes.yml 56

python run_lightcones.py params/parameters_fullcubes.yml_model72 with do_C2ray=false

python run_obslc.py params/parameters_fullcubes.yml_model72 with do_C2ray=false

# PS3/IM1:
python run_py21cmfast_allmodels.py params/parameters_fullcubes.yml 27

python run_lightcones.py params/parameters_fullcubes.yml_model33 with do_C2ray=false

python run_obslc.py params/parameters_fullcubes.yml_model72 with do_C2ray=false
