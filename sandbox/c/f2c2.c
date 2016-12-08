#include <stdio.h>

/*print Fahrenheit-Celsius table
  for fahr = 0, 20, ..., 300 floating-point version */
main ()
{
    float fahr, celsius;
    int lower, upper, step;

    lower=0; /*lower limit of temp table*/
    upper=300; /*upper limit*/
    step=20; /*step size*/

    fahr=lower;
    printf("%-15s\t%-15s\n","Fahrenheit","Celsius");
    while (fahr <= upper)
    {
	celsius=5.0/9*(fahr-32);
	printf("%-15.0f\t%-15.2f\n",fahr,celsius);
	fahr=fahr+step;
    }
}
