#DiffusionFusion/diffusionFusion smolCoarse_16ps.fusion > smolCoarse_16ps.log
bash costEvolution.sh smolCoarse_16ps.log
xmgrace smolCoarse_16ps.costUpdate.dat &

xmgrace smolCoarse_16ps.{0..9}.diffuse &

tclsh metroPosterior.tcl smolCoarse_16ps.traj coarse_oneHundred.dat 1000 1 post_smolCoarse_16ps
