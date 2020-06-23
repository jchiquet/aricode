#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List countPairs(IntegerVector classi1, IntegerVector classi2, IntegerVector order) {
  // first path to count pairs
  int n = classi1.size();

  // count per classification
  IntegerVector count1(n, 0);
  for(int i = 0; i < n; i++) count1[classi1[i]]++;

  IntegerVector count2(n, 0);
  for(int i = 0; i < n; i++) count2[classi2[i]]++;

  // count per pairs
  int count = 1;
  int class1_cur = classi1[order[0]];
  int class2_cur = classi2[order[0]];

  for(int i = 1; i < n; i++){
    if( (class1_cur != classi1[order[i]]) | (class2_cur != classi2[order[i]]) ){
      count++;
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


// [[Rcpp::export]]
double expected_MI(IntegerVector ni_, IntegerVector n_j) {

  int N = sum(ni_) ;

  double emi = 0.0 ;

  NumericVector ni_f = lfactorial(ni_) ;
  NumericVector nj_f = lfactorial(n_j) ;
  NumericVector Nmni_f = lfactorial(N - ni_) ;
  NumericVector Nmnj_f = lfactorial(N - n_j) ;
  double N_f = lgamma(N + 1) ;

  for (int i=0; i< ni_.size(); i++) {
    for (int j=0; j< n_j.size(); j++) {

      int start_nij = std::max(1, ni_[i] + n_j[j] - N) ;
      int end_nij = std::min(ni_[i], n_j[j]) ;

      for (int nij = start_nij; nij <= end_nij; nij++ ) {

          double t1 = ((float) nij / (float) N) * std::log((float)(nij * N) / (float)(ni_[i]*n_j[j])) ;

          double t2 = std::exp((ni_f[i] + nj_f[j] + Nmni_f[i] + Nmnj_f[j] - N_f - lgamma(1 + nij) - lgamma(1 + ni_[i] - nij) - lgamma(1 + n_j[j] - nij) - lgamma(1 + N - ni_[i] - n_j[j] + nij))) ;

          emi += t1*t2;
      }
    }
  }
  return emi;

}

// [[Rcpp::export]]
List getRank(IntegerVector classi){
   int maxi = max(classi);
   int mini = min(classi);

   // Present
   LogicalVector present(maxi - mini + 1);
   for(int i=0; i< classi.size(); i++) present[classi[i]-mini] = TRUE;
   
   // Count
   IntegerVector translator(maxi - mini + 1);
   int nbIndex = 0;
   for(int i=0; i< present.size(); i++) {
     if(present[i]) nbIndex++;
   }

   // Translator and Index Vector
   IntegerVector index(nbIndex);
   int indexCur = 0;
   for(int i=0; i< present.size(); i++) {
     if(present[i]) { 
        translator[i] = indexCur;
	index[indexCur] = i+mini;
	indexCur++;
     } else {
        translator[i] = R_NaN;
     }
   }
   // Converted Vector
   IntegerVector translated(classi.size());
   for(int i=0; i< classi.size(); i++) translated[i] = translator[classi[i] - mini];

   // output as a list
   List ListOut;
   ListOut["index"] =  index;
   ListOut["translator"] = translator;
   ListOut["translated"]   = translated;
   return ListOut;
} 
  

   
