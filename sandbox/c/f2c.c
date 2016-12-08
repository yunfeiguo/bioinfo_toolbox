#include <stdio.h>

/*print Fahrenheit-Celsius table
  for fahr = 0, 20, ..., 300 */
main ()
{
    int fahr, celsius;
    int lower, upper, step;

    lower=0; /*lower limit of temp table*/
    upper=300; /*upper limit*/
    step=20; /*step size*/

    fahr=lower;
    while (fahr <= upper)
    {
	celsius=5*(fahr-32)/9;
	printf("%5d\t%5d\n",fahr,celsius);
	fahr=fahr+step;
    }
}
