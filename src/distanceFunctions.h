#include <R.h>
#include <Rinternals.h>
#include <math.h>

void distanceOnTorus_(int *d, double *x, double *y, double *distance)
{
	int k;
	double squaredDistance, addedTerm, altTerm;
	squaredDistance = 0;
	for (k = 0; k < *d; ++k)
	{
		addedTerm = fabs(x[k] - y[k]);
		altTerm = fabs(x[k] + (1 - y[k]));
		if (altTerm < addedTerm)
		{
			addedTerm = altTerm;
		}
		altTerm = fabs(y[k] + (1 - x[k]));
		if (altTerm < addedTerm)
		{
			addedTerm = altTerm;
		}
		squaredDistance += addedTerm*addedTerm;
	}
	*distance = sqrt(squaredDistance);
}

void euclideanDistance_(int *d, double *x, double *y, double *distance)
{
	int k;
	double squaredDistance, addedTerm;
	squaredDistance = 0;
	for (k = 0; k < *d; ++k)
	{
		addedTerm = x[k] - y[k];
		squaredDistance += addedTerm*addedTerm;
	}
	*distance = sqrt(squaredDistance);
}