# Falcon
## running
before running `fc_run.py`, switch to bash and then go to falcon directory and source `env.sh`

## change qsub
change `FALCON-integrate/pypeFlow/pwatcher/*.py`

## how to change already configured jobs?
find a json file, and modify it

# smrtportal, hgap
## how to run jobs manually?
go to `workflow` and look at `Workflow.details.svg` to figure out order of execution, and copy paste commands from `P_*/` folders.

## fix some robustness
### FilterM4
/home/yunfeiguo/Downloads/smrtanalysis2.3/install/smrtanalysis_2.3.0.140936/analysis/bin/filterm4.py
line 52, add
```
if len(rec.split()) != 13:
    continue
```
### 
/home/yunfeiguo/Downloads/smrtanalysis2.3/install/smrtanalysis_2.3.0.140936/analysis/bin/m4topre.py
line 195, add
```
if rec.qname in seqs:
    qseq = seqs[rec.qname][qst:qnd]
else:
    continue
```


# BUSCO
## installation
```
#install BUSCO
git clone https://gitlab.com/ezlab/busco.git
pushd busco
mkdir -p db
pushd db
wget http://busco.ezlab.org/v2/datasets/nematoda_odb9.tar.gz
popd
popd
#install HMMER
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
 tar zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
#install augustus
wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.3.tar.gz
tar zxvf augustus-3.2.3.tar.gz
pushd augustus-3.2.3
make
##add binaries, scripts/* to PATH
popd
##export AUGUSTUS_CONFIG_PATH=

#install NCBI blast
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
 tar zvxf ncbi-blast-2.6.0+-x64-linux.tar.gz

```
## analysis

# miniasm/minimap
## resource requirements
even use 24 threads, memory usage does not exceed 65GB, time is < 1hr, total CPU time is < 12hr.
## parameter tuning
* `-S` for minimap is required for de-novo assembly.
* combining raw reads with draft assembly will lead to very high N50 but very poor quality, e.g. very low GC content, at least with the default settings for minimap and miniasm.
* `-T1` will lead to very short assembly, but `-T0` and `-T10` are okay.
* current best parameter combination for minimap (w.r.t N50) is `minimap -k15 -w7 -f0.001 -r500 -m0 -c4 -L100 -g10000 -T0 -S`
## command
```
#qsub -q largemem -N mini3 -l walltime=48:0:0 -l nodes=1:ppn=24 -l mem=256GB -S /bin/bash -V -o stdout -e stderr -S /bin/bash cmd
INPUT=/home/rcf-02/yunfeigu/proj_dir/nematode/data/filtered_subreads.fasta
NCPU=24
# Overlap
minimap -Sw5 -L100 -m0 -t $NCPU $INPUT $INPUT | gzip -1 > out.paf.gz
# Layout
miniasm -f $INPUT out.paf.gz > out.gfa
gfa2fa.sh out.gfa > out.fa
```
## troubleshooting
### disk quota exceeded
it seems the intermediate file is very large, exceeding the remaining 1.x TB space.
