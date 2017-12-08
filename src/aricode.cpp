#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List countPairs(IntegerVector classi1, IntegerVector classi2, IntegerVector order) {
  // first path to count pairs
  int n = classi1.size();

  // count per classification
  IntegerVector count1(n, 0);
  for(int i = 0; i < n; i++) count1[classi1[i]-1]++;

  IntegerVector count2(n, 0);
  for(int i = 0; i < n; i++) count2[classi2[i]-1]++;

  // count per pairs
  int count = 1;
  int class1_cur = classi1[order[0]];
  int class2_cur = classi2[order[0]];

  for(int i = 1; i < n; i++){
    if( (class1_cur != classi1[order[i]]) | (class2_cur != classi2[order[i]]) ){
      count += 1;
      class1_cur = classi1[order[i]];
      class2_cur = classi2[order[i]];
    }
  }

  // create output Integer Vector for pairs and initialize
  IntegerVector nameClassi1(count, 0);
  IntegerVector nameClassi2(count, 0);
  IntegerVector numberPair(count, 0);

  int current_position = 0;
  nameClassi1[0] = classi1[order[0]];
  nameClassi2[0] = classi2[order[0]];
  numberPair[0] = 1;

  // count pairs
  for(int i = 1; i < n; i++){
    if( ( nameClassi1[current_position] == classi1[order[i]]) & (nameClassi2[current_position] == classi2[order[i]]) ){
      numberPair[current_position]++;
    } else {
      current_position += 1;
      nameClassi1[current_position] = classi1[order[i]];
      nameClassi2[current_position] = classi2[order[i]];
      numberPair[current_position]  = 1;
    }
  }

  // output as a list
  List ListOut;
  ListOut["pair_nb"] = numberPair;
  ListOut["pair_c1"] = nameClassi1;
  ListOut["pair_c2"] = nameClassi2;
  ListOut["c1_nb"]   = count1[count1 > 0];
  ListOut["c2_nb"]   = count2[count2 > 0];
  return(ListOut);
}
