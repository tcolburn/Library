find ./frywalker/003/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@frywalker {} \;
find ./chipbacca/003/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@chipbacca {} \;
find ./spudda/003/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@spudda {} \;
find ./yamsolo/003/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@yamsolo {} \;
