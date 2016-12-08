#include<stdio.h>
#include<stdlib.h>

#define MB 1024*1024
int main() {
    long total = 0;
    int i;
    for (i = 0; i < 100000000; i++) {
	char *i = malloc(sizeof(char)*100*MB);
	if(i == NULL) {
	    printf("malloc fails at %ld MB\n",total);
	    break;
	} else {
	    total += 100;
	}
    }
    return 0;
}
