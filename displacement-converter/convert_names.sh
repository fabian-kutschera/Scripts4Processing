ncrename -v elevation,z gebco_2022_n40.0_s36.0_w23.0_e29.0.nc
ncrename -d lat,y -d lon,x gebco_2022_n40.0_s36.0_w23.0_e29.0.nc
ncrename -v lat,y -v lon,x gebco_2022_n40.0_s36.0_w23.0_e29.0.nc
echo "Done. Names converted."
