#include "config.h"
#include "func.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>


int main(int argc, char* argv[])
{
	srand(time(NULL));
	
	real * randnums = (real *) malloc(500 * sizeof(real));

	for (whole i = 0; i < 500; i++)
	{
		randnums[i] = randomDist();
	}
	
	whole * randbuckets = (whole *) malloc(40 * sizeof(whole));
	
	for (whole i = 0; i < 40; i++)
	{
		randbuckets[i] = 0;
	
		for (whole j = 0; j < 500; j++)
		{
			if (randnums[j] >= (real) i && randnums[j] < (real) i + 1) randbuckets[i]++;
		}
		
		printf("%d", i);
		
		for (whole j = 0; j < randbuckets[i]; j++)
		{
			printf("#");
		}
		
		printf("\n");
	}
	
	return 0;
}