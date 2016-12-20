f=../Data/abf_lanGrid0.colvars.traj
#awk '!/^#/ && NF>4 && $1%1000==0 {print $1,$2}' $f > traj.dat

out=thermostat0.colvars.traj
echo "#step         PosZ                 fs_PosZ                fa_PosZ" > $out
awk '!/^#/ && NF>4 && $1%1000==0 && $1<1e7 {print $1,$2,$3,$5}' $f >> $out
out=thermostat1.colvars.traj
echo "#step         PosZ                 fs_PosZ                fa_PosZ" > $out
awk '!/^#/ && NF>4 && $1%1000==0 && $1>1e7 {print $1,$2,$3,$5}' $f >> $out



