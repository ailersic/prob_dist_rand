#include "config.h"
#include "probdist.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


real randomUniform() /* uniform distribution, (0..1] */
{
	return (rand() + 1.0) / (RAND_MAX + 1.0);
}

real randomDist() /* normal distribution, centered on 0, std dev 1 */
{
	weibullcd test (5, 10);
	
	real i = test.randval();
	printf("%f\n", i);
	
	return i;
}
