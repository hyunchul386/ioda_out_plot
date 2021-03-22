# ioda_plot

This is an plotting suite for ioda format data, 
expanded from
/scratch1/NCEPDEV/da/Cory.R.Martin/JEDI/tutorial_dec2020/python/plot_jedi_obs.py

The expansions are mainly
1. A yaml file of plot_ioda_obs.yaml sets parameters of plot_ioda_obs.py.
2. It is able to handle for ioda format data JEDI/SOCA as well as GFS.
3. It is able to filter data by regional (horizonal/vertical) sets and Max/Min of data.

Remark: This plotting suit is currently only available in Orion 
        for the anaconda/PyYAML versions

To run the plot_ioda_obs.py in Orion,
- git clone
  git clone https://github.com/hyunchul386/ioda_plot.git
- set modules
  module use -a /home/cmartin/opt/modulefiles
  module load anaconda/anaconda3-2020.04.02
- run code
  cd ioda_plot
  edit plot_ioda_obs.yaml for the data
  python plot_ioda_obs.py
        
        
