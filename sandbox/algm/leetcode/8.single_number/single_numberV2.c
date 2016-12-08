#include <stdio.h>
int singleNumber(int* nums, int numsSize) {
    int ones = 0;
    int i = 0;
    for (i = 0; i < numsSize; i++) {
	ones = ones ^ nums[i];
    }
    return ones;
}

int main() {
    int x[3] = {1,2,1};
    printf("%d",singleNumber(x,3));
    return 0;
}
