#include <stdio.h>

main ()
{
    float celsius,fahr;
    int step,start,end;

    start=-10;
    end=100;
    step=5;

    celsius=start;
    printf("%-15s\t%-15s\n","Celsius","Fahrenheit");
    while(celsius<=end)
    {
	fahr=celsius*9/5.0+32;
	printf("%-15.0f\t%-15.1f\n",celsius,fahr);
	celsius=celsius+step;
    }
}
