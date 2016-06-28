for i in {1..5}
do
  DIRS=$(sed -n "${i}p" dir_list.txt)
  cd $DIRS
  cp full_dims_mhp1_occ2if.dcd /nfs/homes/tcolburn/Projects/Beckstein/Mhp1/simulations/dims/implicit_solvent/implicit_membrane/bench/1_core_psa/
  cd /nfs/homes/tcolburn/Projects/Beckstein/Mhp1/simulations/dims/implicit_solvent/implicit_membrane/bench/test_1
done
