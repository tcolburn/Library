for i in {1..16}
do
  DIR=$(sed -n "${i}p" analysis_dirs.txt)
  cd $DIR
  rm mhp1_analysis.py
  cp /nfs/homes/tcolburn/Projects/Beckstein/Mhp1/analysis/gen_tools/mhp1_analysis.py ./
  cd /nfs/homes/tcolburn/Projects/Beckstein/Mhp1/simulations/dims/implicit_solvent/implicit_membrane/bench/test_1
done
