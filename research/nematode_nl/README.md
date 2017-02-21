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
during parameter tuning, max mem usage reached 88GB using 24 threads.
## parameter tuning
* `-S` for minimap is required for de-novo assembly.
* combining raw reads with draft assembly will lead to very high N50 but very poor quality, e.g. very low GC content, at least with the default settings for minimap and miniasm.
* `-T1` will lead to very short assembly, but `-T0` and `-T10` are okay.
* current best parameter combination for minimap (w.r.t N50) is `minimap -k15 -w7 -f0.001 -r600 -m0.05 -c4 -L100 -g10000 -T0 -S`
* increasing `-w` from 3 to 7 almost doubled assembly N50.
* increasing `-f` from 0.0005 to 0.002 increases N50
* increasing `-r` from 400 to 600 increases N50
* `-m0` < `-m0.1` > `-m0.2` w.r.t N50.
* increasing `-c` from 2 to 6 increases N50
* `-L70` `-L130` `-L100` do not seem to be significantly different
* `-g5000` doubles N50 at `-g10000`.
* current best parameter combination for miniasm (wrt N50): `-m100 -i0.05 -s2000 -c3 -h1000 -I0.8 -g1000 -d50000 -e4 -n3 -r0.7,0.5 -F0.8`
* `-R` decreases N50
* `-m50` `-m100` do not have difference, `-m150` reduces N50
* `-i0.05`->`-i0.02` or `-i0.08` halved N50.
* `-s2000`->`-s1500`or`-s2500` seems no difference
* `-c3`->`c2` or`-c4` reduced N50
* `-h1000`->`h600` significantly reduced N50, ->`-h1400` mildly reduced N50.
* `-I0.8`->`-I0.6` or `-I0.9` reduced N50
* `-g1000`->`-g500`,`-g1500` seems no effect
* `-d30000`->`-d50000`->`-d70000` increased N50
* `-e3,4,5` no difference
* `-n2,3,4` no difference
* `-r0.7,0.5`->`-r0.8,0.6` reduced N50
* `-F0.7,0.8,0.9` no difference
* `-b` will result in 0bp assembly.


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

# racon
## installation
I did fresh installation of gcc-6.3.0 to overcome issue of `cannot find -lstdc++`.
