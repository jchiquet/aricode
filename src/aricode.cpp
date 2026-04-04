/*
 *   Copyright (c) 2026 INRAE
 *   guillem.rigaill@inrae.fr, julien.chiquet@inrae.fr
 *   INRAE
 */

#include <Rcpp.h>

#include <algorithm>
#include <vector>

// [[Rcpp::export]]
Rcpp::List cpp_SortPairs(Rcpp::IntegerVector& c1, Rcpp::IntegerVector& c2,
                         unsigned int N1, unsigned int N2) {
  unsigned int n = c1.size();
  if (n == 0) return Rcpp::List::create();

  // Comptage des individus par classe
  Rcpp::IntegerVector count1(N1);
  Rcpp::IntegerVector count2(N2);

  for (unsigned int i = 0; i < n; i++) {
    count1[c1[i]]++;
    count2[c2[i]]++;
  }

  // Allocation temporaire (automatiquement libérée en fin de scope)
  Rcpp::IntegerVector tmp_c1(n), tmp_c2(n);

  // Tri selon c2 (Counting Sort logic)
  Rcpp::IntegerVector shift2(N2 + 1);
  for (unsigned int i = 1; i <= N2; i++) {
    shift2[i] = shift2[i - 1] + count2[i - 1];
  }

  for (unsigned int i = 0; i < n; i++) {
    int idx = shift2[c2[i]];
    tmp_c1[idx] = c1[i];
    tmp_c2[idx] = c2[i];
    shift2[c2[i]]++;
  }

  // Tri final selon c1
  Rcpp::IntegerVector new_c1(n), new_c2(n), shift1(N1 + 1);
  for (unsigned int i = 1; i <= N1; i++) {
    shift1[i] = shift1[i - 1] + count1[i - 1];
  }

  for (unsigned int i = 0; i < n; i++) {
    int idx = shift1[tmp_c1[i]];
    new_c1[idx] = tmp_c1[i];
    new_c2[idx] = tmp_c2[i];
    shift1[tmp_c1[i]]++;
  }

  // Calcul des paires uniques et comptage
  Rcpp::IntegerVector pair_c1(n), pair_c2(n), pair_count(n);

  unsigned int i_index = 0;
  int pair_cur_c1 = new_c1[0];
  int pair_cur_c2 = new_c2[0];

  pair_c1[0] = pair_cur_c1;
  pair_c2[0] = pair_cur_c2;

  for (unsigned int i = 0; i < n; i++) {
    if ((new_c1[i] == pair_cur_c1) && (new_c2[i] == pair_cur_c2)) {
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

  Rcpp::IntegerVector ind = Rcpp::seq(0, i_index);
  return Rcpp::List::create(Rcpp::Named("pair_nb") = pair_count[ind],
                            Rcpp::Named("pair_c1") = pair_c1[ind],
                            Rcpp::Named("pair_c2") = pair_c2[ind],
                            Rcpp::Named("c1_nb") = count1[count1 > 0],
                            Rcpp::Named("c2_nb") = count2[count2 > 0]);
}

// [[Rcpp::export]]
Rcpp::List std_SortPairs(Rcpp::IntegerVector c1_in, Rcpp::IntegerVector c2_in,
                         unsigned int N1, unsigned int N2) {
  unsigned int n = c1_in.size();
  if (n == 0) return Rcpp::List::create();

  // Conversion des entrées en std::vector
  std::vector<int> c1 = Rcpp::as<std::vector<int>>(c1_in);
  std::vector<int> c2 = Rcpp::as<std::vector<int>>(c2_in);

  // Comptage des individus par classe
  std::vector<int> count1(N1, 0);
  std::vector<int> count2(N2, 0);

  for (unsigned int i = 0; i < n; i++) {
    count1[c1[i]]++;
    count2[c2[i]]++;
  }

  // Allocation temporaire
  std::vector<int> tmp_c1(n), tmp_c2(n);

  // Tri selon c2 (Counting Sort logic)
  std::vector<int> shift2(N2 + 1, 0);
  for (unsigned int i = 1; i <= N2; i++) {
    shift2[i] = shift2[i - 1] + count2[i - 1];
  }

  for (unsigned int i = 0; i < n; i++) {
    int idx = shift2[c2[i]];
    tmp_c1[idx] = c1[i];
    tmp_c2[idx] = c2[i];
    shift2[c2[i]]++;
  }

  // Tri final selon c1
  std::vector<int> new_c1(n), new_c2(n);
  std::vector<int> shift1(N1 + 1, 0);
  for (unsigned int i = 1; i <= N1; i++) {
    shift1[i] = shift1[i - 1] + count1[i - 1];
  }

  for (unsigned int i = 0; i < n; i++) {
    int idx = shift1[tmp_c1[i]];
    new_c1[idx] = tmp_c1[i];
    new_c2[idx] = tmp_c2[i];
    shift1[tmp_c1[i]]++;
  }

  // Calcul des paires uniques et comptage
  std::vector<int> pair_c1(n), pair_c2(n), pair_count(n, 0);

  unsigned int i_index = 0;
  int pair_cur_c1 = new_c1[0];
  int pair_cur_c2 = new_c2[0];

  pair_c1[0] = pair_cur_c1;
  pair_c2[0] = pair_cur_c2;

  for (unsigned int i = 0; i < n; i++) {
    if ((new_c1[i] == pair_cur_c1) && (new_c2[i] == pair_cur_c2)) {
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

  // Préparation du retour : on retaille les vecteurs à la taille réelle i_index
  // + 1
  pair_c1.resize(i_index + 1);
  pair_c2.resize(i_index + 1);
  pair_count.resize(i_index + 1);

  // Filtrage des counts > 0
  std::vector<int> res_count1, res_count2;
  for (int c : count1)
    if (c > 0) res_count1.push_back(c);
  for (int c : count2)
    if (c > 0) res_count2.push_back(c);

  return Rcpp::List::create(
      Rcpp::Named("pair_nb") = pair_count, Rcpp::Named("pair_c1") = pair_c1,
      Rcpp::Named("pair_c2") = pair_c2, Rcpp::Named("c1_nb") = res_count1,
      Rcpp::Named("c2_nb") = res_count2);
}

// [[Rcpp::export]]
double expected_MI(Rcpp::IntegerVector ni_, Rcpp::IntegerVector n_j) {
  int n_rows = ni_.size();
  int n_cols = n_j.size();
  int N = sum(ni_);
  double N_double = static_cast<double>(N);

  double emi = 0.0;
  double log_N = std::log(N_double);

  // Pré-calcul des log-factorielles pour gagner en performance
  Rcpp::NumericVector log_fact = lfactorial(Rcpp::seq(0, N));
  double log_fact_N = log_fact[N];

  for (int i = 0; i < n_rows; i++) {
    int ni = ni_[i];
    if (ni == 0) continue;

    double term_ni = log_fact[ni] + log_fact[N - ni];

    for (int j = 0; j < n_cols; j++) {
      int nj = n_j[j];
      if (nj == 0) continue;

      int start_nij = std::max(1, ni + nj - N);
      int end_nij = std::min(ni, nj);

      double common_term =
          term_ni + log_fact[nj] + log_fact[N - nj] - log_fact_N;
      double log_ni_nj = std::log(static_cast<double>(ni) * nj);

      for (int nij = start_nij; nij <= end_nij; nij++) {
        // Probabilité hypergéométrique P(nij)
        double log_p_nij =
            common_term - (log_fact[nij] + log_fact[ni - nij] +
                           log_fact[nj - nij] + log_fact[N - ni - nj + nij]);

        double p_nij = std::exp(log_p_nij);

        // Calcul du terme MI : (nij / N) * log( (nij * N) / (ni * nj) )
        // Équivalent à : (nij / N) * (log(nij) + log(N) - log(ni * nj))
        double mi_term =
            (static_cast<double>(nij) / N_double) *
            (std::log(static_cast<double>(nij)) + log_N - log_ni_nj);

        emi += p_nij * mi_term;
      }
    }
  }

  return emi;
}
// [[Rcpp::export]]
Rcpp::List getRank(Rcpp::IntegerVector classi) {
  int maxi = Rcpp::max(classi);
  int mini = Rcpp::min(classi);

  // Present
  Rcpp::LogicalVector present(maxi - mini + 1);
  for (int i = 0; i < classi.size(); i++) present[classi[i] - mini] = TRUE;

  // Count
  Rcpp::IntegerVector translator(maxi - mini + 1);
  int nbIndex = 0;
  for (int i = 0; i < present.size(); i++) {
    if (present[i]) nbIndex++;
  }

  // Translator and Index Vector
  Rcpp::IntegerVector index(nbIndex);
  int indexCur = 0;
  for (int i = 0; i < present.size(); i++) {
    if (present[i]) {
      translator[i] = indexCur;
      index[indexCur] = i + mini;
      indexCur++;
    } else {
      translator[i] = NA_INTEGER;
    }
  }

  // Converted Vector
  Rcpp::IntegerVector translated(classi.size());
  for (int i = 0; i < classi.size(); i++)
    translated[i] = translator[classi[i] - mini];

  // output as a list
  return Rcpp::List::create(Rcpp::Named("index") = index,
                            Rcpp::Named("translator") = translator,
                            Rcpp::Named("translated") = translated);
}

//
// // [[Rcpp::export]]
// List countPairs(IntegerVector classi1, IntegerVector classi2, IntegerVector
// order) {
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
//     if( (class1_cur != classi1[order[i]]) || (class2_cur !=
//     classi2[order[i]]) ){
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
//     if( ( nameClassi1[current_position] == classi1[order[i]]) &&
//     (nameClassi2[current_position] == classi2[order[i]]) ){
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
