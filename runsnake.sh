echo ""
echo "Running the snakemake workflow for the admixture process..."
echo ""
echo "Part 1 of 4 - Extract chromosomal DNA, remove duplicates, creates metrics file..."
/fh/fast/ha_g/app/bin/snakemake -s mixarray_p1.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} --requeue " --rerun-incomplete -j 40
echo ""
echo "Part 2 of 4 - Create the file samples2.yaml, look in logs/samples2yaml.txt for output ..."
python makeSamples2v2.py -i ./config/samples1.yaml -o ./config/samples2.yaml > ./logs/samples2yaml.txt
echo ""
echo "Part 3 of 4 - Downsample the source files as needed to prepare for mixing them..."
/fh/fast/ha_g/app/bin/snakemake -s mixarray_p2.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} --requeue " -j 40
echo ""
echo "Part 4 of 4 - Merge the files, create symlinks and metrics files..."
/fh/fast/ha_g/app/bin/snakemake -s mixarray_p3.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} --requeue " -j 40
