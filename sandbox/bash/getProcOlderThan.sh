user=yunfeiguo
ps x -lf | grep $user | perl -ane '($h,$m,$s)=split /:/,$F[13]; !system("echo kill -SIGINT $F[14]") if ($h>12)'
