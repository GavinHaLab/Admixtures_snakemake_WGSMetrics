# Admixtures_snakemake_WGSMetrics

Adapted from the [original Admixture pipeline](https://github.com/GavinHaLab/Admixtures_snakemake), some new features were added:

- WGS-metrics based coverage assessments per input file, then downsample for a target unique fragment depth. Downsampling probabiltiies are adapted to the WGS-metrics based outputs.

- Modified pipeline to run on .cram files, seen in all three snakefiles. This removes barriers for when input data is downloaded from aws to be included in admixtures, should make it more user friendly long-term.

Full documentation can be found [here](https://fredhutch.atlassian.net/wiki/spaces/GHL/pages/3787260098/WGS+Metrics-Based+Approach). 

