#!/bin/env bash
snakemake  --rerun-incomplete --latency-wait 60 --jobs 100 --keep-going --cluster-config config/cluster_config.json \
  --drmaa " -cwd -j y -o logs/ -jc {cluster.jc} -pe smp {cluster.threads} -adds l_hard h_rt 43200 -mods l_hard mfree {cluster.mem_free} -adds l_hard local_free {cluster.local_free} -adds l_hard ramdisk {cluster.ramdisk}"
