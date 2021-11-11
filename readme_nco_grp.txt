--Flattening Group
ncks -G :1 -g MetaData -v datetime,depth,latitude,longitude,record_number,time insitu_wod.20150705.nc insitu_wod.20150705_ombg_T.nc
ncks -A -G :1 -g ombg -v sea_water_temperature insitu_wod.20150705.nc insitu_wod.20150705_ombg_T.nc

-- Documents

http://nco.sourceforge.net/nco.html#Flattening-Groups

http://nco.sourceforge.net/nco.html#ncks

https://github.com/nco/nco
