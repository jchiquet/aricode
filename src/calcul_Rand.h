#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#define MIN(X,Y) ( ( (X) < (Y) ) ? (X) : (Y) )
#define MAX(X,Y) ( ( (X) > (Y) ) ? (X) : (Y) )

// suppose the first class are from [0 to N1] and [0 to N2]
void c_SortPairs(int * c1, int * c2, int * new_c1, int * new_c2,
	       int *pair_c1, int *pair_c2, int *pair_count,
	       int *count1, int *count2, int *n_, int *nzero) ;



