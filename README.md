# ioda_plot
# module use /work/noaa/da/Hyun-Chul.Lee/modulefiles
# module load anaconda/1.7.2.1
<p>
This is a plotting suite for ioda format data, expanded from <br>
/scratch1/NCEPDEV/da/Cory.R.Martin/JEDI/tutorial_dec2020/python/plot_jedi_obs.py <br>

<br>

The expansions are mainly <br>
1. A yaml file of plot_ioda_data.yaml sets parameters of plot_ioda_data.py. <br>
2. It is able to handle ioda format data from JEDI/SOCA as well as GFS. <br>
3. It is able to filter data by regional (horizonal/vertical) sets and Max/Min of data. <br>

<b>Remark</b> This plotting suit is currently only available in Orion <br>
        for the anaconda/PyYAML versions <br>
<br>

To run the plot_ioda_data.py in Orion, <br>
- git clone <br>
  git clone https://github.com/hyunchul386/ioda_plot.git <br>
- set modules <br>
  module use -a /home/cmartin/opt/modulefiles <br>
  module load anaconda/anaconda3-2020.04.02 <br>
- run code <br>
  cd ioda_plot <br>
  edit plot_ioda_data.yaml for the data <br>
  python plot_ioda_data.py <br>
</p>
