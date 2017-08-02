#LINKS (scaffolder for long read)

```
#use the following command to locate perl CORE library
#perl -e 'use Config; print $Config{archlib};'
g++ -c BloomFilter_wrap.cxx -I /net/hippo/volumes/wadi/shared/prj/guoy28/Downloads/perl-5.24.0_install/lib/5.24.0/x86_64-linux/CORE -fPIC -Dbool=char -O3
g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so -O3
```
