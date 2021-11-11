#! /bin/csh

module load nco/4.9.3
#----------------
#set var_in = $1
#set grp_in = $2
#set dir_in = $3
#set file_in = $4
#set file_out = $5
#----------------
#set var_in = sea_water_temperature
#set grp_in = ombg
#set dir_in = /work/noaa/da/Hyun-Chul.Lee/s2s/3dvar_check/obs_out/2015/2015070500/ctrl
#set file_in = insitu_wod.20150705.nc
#set file_out = insitu_wod.20150705_ombg_T.nc
#set dim_data = 3
#----------------
#set var_in = sea_surface_temperature
#set grp_in = ombg
#set dir_in = /work/noaa/da/Hyun-Chul.Lee/s2s/3dvar_check/obs_out/2015/2015070500/ctrl
#set file_in = sst_noaa19.20150705.nc
#set file_out = sst_noaa19.20150705_ombg_T.nc
#set dim_data = 2
#----------------
set var_in = sea_surface_salinity
set grp_in = ombg
set dir_in = /work/noaa/da/Hyun-Chul.Lee/s2s/3dvar_check/obs_out/2015/2015070500/ctrl
set file_in = sss_smap.jpl.20150705.nc
set file_out = sss_smap.jpl.20150705_ombg_S.nc
set dim_data = 2
#----------------

echo "var_in = " $var_in
 echo "grp_in = " $grp_in
 echo "set dir_in = " $dir_in
 echo "set file_in = " $file_in 
 echo "set file_out = " $file_out

cd $dir_in
if (-e $file_out) rm $file_out

set inpfile = ${dir_in}/${file_in}
if ($dim_data == 3) then
  ncks -G :1 -g MetaData -v datetime,depth,latitude,longitude,record_number,time $file_in $file_out
else
  ncks -G :1 -g MetaData -v datetime,latitude,longitude,time $file_in $file_out
endif
ncrename -v /EffectiveQC0/$var_in,qc $file_in ${file_in}_Q
ncks -A -G :1 -g EffectiveQC0 -v qc ${file_in}_Q $file_out
ncks -A -G :1 -g $grp_in -v $var_in $file_in $file_out
rm ${file_in}_Q





