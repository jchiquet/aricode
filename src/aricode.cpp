#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List cpp_SortPairs(IntegerVector c1, IntegerVector c2) {

  int n = c1.size();
  if (n == 0) return List::create();

  int N1 = 0;
  int N2 = 0;

  // Nombre de classes (N1 et N2)
  N1 = Rcpp::max(c1);
  N2 = Rcpp::max(c2);
  N1++; // Car les classes commencent à 0
  N2++;

  // Comptage des individus par classe
  IntegerVector count1(N1);
  IntegerVector count2(N2);

  for(int i = 0; i < n; i++) {
    count1[c1[i]]++;
    count2[c2[i]]++;
  }

  // Allocation temporaire (automatiquement libérée en fin de scope)
  IntegerVector tmp_c1(n);
  IntegerVector tmp_c2(n);

  // Tri selon c2 (Counting Sort logic)
  IntegerVector shift2(N2 + 1);
  for(int i = 1; i <= N2; i++) {
    shift2[i] = shift2[i-1] + count2[i-1];
  }

  for(int i = 0; i < n; i++) {
    int idx = shift2[c2[i]];
    tmp_c1[idx] = c1[i];
    tmp_c2[idx] = c2[i];
    shift2[c2[i]]++;
  }

  // Tri final selon c1
  IntegerVector new_c1(n);
  IntegerVector new_c2(n);
  IntegerVector shift1(N1 + 1);
  for(int i = 1; i <= N1; i++) {
    shift1[i] = shift1[i-1] + count1[i-1];
  }

  for(int i = 0; i < n; i++) {
    int idx = shift1[tmp_c1[i]];
    new_c1[idx] = tmp_c1[i];
    new_c2[idx] = tmp_c2[i];
    shift1[tmp_c1[i]]++;
  }

  // Calcul des paires uniques et comptage
  // On initialise avec la taille max possible (n)
  IntegerVector pair_c1(n);
  IntegerVector pair_c2(n);
  IntegerVector pair_count(n);

  int i_index = 0;
  int pair_cur_c1 = new_c1[0];
  int pair_cur_c2 = new_c2[0];

  pair_c1[0] = pair_cur_c1;
  pair_c2[0] = pair_cur_c2;

  for(int i = 0; i < n; i++) {
    if((new_c1[i] == pair_cur_c1) && (new_c2[i] == pair_cur_c2)) {
      pair_count[i_index]++;
    } else {
      i_index++;
      pair_cur_c1 = new_c1[i];
      pair_cur_c2 = new_c2[i];
      pair_c1[i_index] = pair_cur_c1;
      pair_c2[i_index] = pair_cur_c2;
      pair_count[i_index]++;
    }
  }

  // Comme i_index est un index basé sur 0, le nombre de paires est i_index + 1
  int num_pairs = i_index + 1;

  // Retourne une liste de résultats (équivalent à nzero et aux pointeurs modifiés)
  return List::create(
    Named("pair_nb") = pair_count[seq(0, i_index)],
    Named("pair_c1") = pair_c1[seq(0, i_index)],
    Named("pair_c2") = pair_c2[seq(0, i_index)],
    Named("c1_nb")   = count1[count1 > 0],
    Named("c2_nb")   = count2[count2 > 0]
  );
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
        translator[i] = NA_INTEGER;
     }
   }
   // Converted Vector
   IntegerVector translated(classi.size());
   for(int i=0; i< classi.size(); i++) translated[i] = translator[classi[i] - mini];

   // output as a list
   List ListOut;
   ListOut["index"] = index;
   ListOut["translator"] = translator;
   ListOut["translated"] = translated;
   return ListOut;
}




//
// // [[Rcpp::export]]
// List countPairs(IntegerVector classi1, IntegerVector classi2, IntegerVector order) {
//   // first path to count pairs
//   unsigned int n = classi1.size();
//
//   // count per classification
//   IntegerVector count1(n, 0);
//   IntegerVector count2(n, 0);
//   for(unsigned int i = 0; i < n; i++) {
//     count1[classi1[i]]++;
//     count2[classi2[i]]++;
//   }
//
//   // count per pairs
//   unsigned int count = 1;
//   unsigned int class1_cur = classi1[order[0]];
//   unsigned int class2_cur = classi2[order[0]];
//
//   for(unsigned int i = 1; i < n; i++){
//     if( (class1_cur != classi1[order[i]]) || (class2_cur != classi2[order[i]]) ){
//       count++;
//       class1_cur = classi1[order[i]];
//       class2_cur = classi2[order[i]];
//     }
//   }
//
//   // create output Integer Vector for pairs and initialize
//   IntegerVector nameClassi1(count, 0);
//   IntegerVector nameClassi2(count, 0);
//   IntegerVector numberPair(count, 0);
//
//   unsigned int current_position = 0;
//   nameClassi1[0] = classi1[order[0]];
//   nameClassi2[0] = classi2[order[0]];
//   numberPair[0] = 1;
//
//   // count pairs
//   for(unsigned int i = 1; i < n; i++){
//     if( ( nameClassi1[current_position] == classi1[order[i]]) && (nameClassi2[current_position] == classi2[order[i]]) ){
//       numberPair[current_position]++;
//     } else {
//       current_position += 1;
//       nameClassi1[current_position] = classi1[order[i]];
//       nameClassi2[current_position] = classi2[order[i]];
//       numberPair[current_position]  = 1;
//     }
//   }
//
//   // output as a list
//   return List::create(
//     Named("pair_nb") = numberPair,
//     Named("pair_c1") = nameClassi1,
//     Named("pair_c2") = nameClassi2,
//     Named("c1_nb")   = count1[count1 > 0],
//                              Named("c2_nb")   = count2[count2 > 0]
//   ) ;
// }
