find ./mashteryoda/002/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@mashteryoda {} \;
find ./palpoutine/002/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@palpoutine {} \;
find ./spuddafett/002/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@spuddafett {} \;
find ./wedge/002/ -name "dims_seed_bench.ge" -type f -exec qsub_dependents.py -N 6 -- -q workstations.q@wedge {} \;
