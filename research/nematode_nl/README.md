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