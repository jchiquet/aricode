#include "calcul_Rand.h"
#include <R_ext/Rdynload.h>

// suppose the first class are from [0 to N1] and [0 to N2]
void c_SortPairs(int * c1, int * c2, int * new_c1, int * new_c2,
	       int *pair_c1, int *pair_c2, int *pair_count,
	       int *count1, int *count2, int *n_, int *nzero){

	int n  = n_[0];
	int N1 = 0;
	int N2 = 0;

	// Number of class in c1 and c2
	for(int i=0; i <n; i++) N1 = MAX(N1,c1[i]);
	N1++;
	for(int i=0; i <n; i++) N2 = MAX(N2,c2[i]);
	N2++;

	// Number of individuals for each c1 and c2 class
	for(int i=0; i < n; i++) count1[i]=0;
	for(int i=0; i < n; i++) count2[i]=0;


	for(int i=0; i < n; i++) count1[c1[i]]++;
	for(int i=0; i < n; i++) count2[c2[i]]++;

	int * tmp_c2 = malloc(n*sizeof(int));
	int * tmp_c1 = malloc(n*sizeof(int));

	// Sort along c2
	int * shift2 = malloc((N2+1)*sizeof(int));
	shift2[0] = 0;
	for(int i=1; i < N2+1; i++) shift2[i]=shift2[i-1]+count2[i-1];


	for(int i=0; i <n; i++)
	{
		tmp_c1[shift2[c2[i]]] = c1[i];
		tmp_c2[shift2[c2[i]]] = c2[i];
		shift2[c2[i]]++;
	}

	// Sort along c1
	int * shift1 = malloc((N1+1)*sizeof(int));
	shift1[0] = 0;
	for(int i=1; i < N1+1; i++) shift1[i]=shift1[i-1]+count1[i-1];

	for(int i=0; i <n; i++)
	{
		new_c1[shift1[tmp_c1[i]]] = tmp_c1[i];
		new_c2[shift1[tmp_c1[i]]] = tmp_c2[i];
		shift1[tmp_c1[i]]++;
	}

	// Calcul of pairs
	for(int i=0; i < n; i++) pair_count[i]=0;
	int pair_cur_c1, pair_cur_c2, i_index;
	i_index=0;

	pair_cur_c1 = new_c1[0];
	pair_cur_c2 = new_c2[0];
	pair_c1[i_index] = pair_cur_c1;
	pair_c2[i_index] = pair_cur_c2;

	for(int i=0;i <n;i++)
	{
		if((new_c1[i] == pair_cur_c1) & (new_c2[i] == pair_cur_c2))
		{
			pair_count[i_index]++;
		} else
		{
			pair_cur_c1 = new_c1[i];
			pair_cur_c2 = new_c2[i];
			i_index++;
			pair_count[i_index]++;
			pair_c1[i_index] = pair_cur_c1;
			pair_c2[i_index] = pair_cur_c2;
		}
	}
	nzero[0] = i_index ;

	//pair_c1, pair_c2, pair_count,
	free(shift1);
	free(tmp_c1);
	free(tmp_c2);
}


