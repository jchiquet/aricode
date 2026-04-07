/*
 *   Copyright (c) 2026 INRAE
 *   guillem.rigaill@inrae.fr, julien.chiquet@inrae.fr
 *   INRAE
 */

#include <Rcpp.h>

#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>

// ====================================================================
//
// Function get_rank
//
// This function performs a rank transformation (dense encoding)
// on an integer vector.
//
// [[Rcpp::export]]
Rcpp::List get_rank(Rcpp::IntegerVector& classif) {
  unsigned int n = classif.size();

  // 1. Find min and max in a single pass using std::minmax_element
  auto mm = std::minmax_element(classif.begin(), classif.end());
  int mini = *mm.first, maxi = *mm.second;
  int range = maxi - mini + 1;

  // Mark presence (Use std::vector<char> instead of bool/LogicalVector for
  // speed/space)
  std::vector<char> present(range, 0);
  for (int val : classif) {
    present[val - mini] = 1;
  }

  // Translator and Index Vector
  Rcpp::IntegerVector translator(range);
  std::vector<int> index_vec;
  index_vec.reserve(range);  // avoid frequent reallocs

  // We can use the presence array to calculate ranks
  int current_rank = 0;
  for (int i = 0; i < range; ++i) {
    if (present[i]) {
      translator[i] = current_rank;
      index_vec.push_back(i + mini);
      current_rank++;
    } else {
      translator[i] = NA_INTEGER;
    }
  }

  // Translate the original vector
  Rcpp::IntegerVector translated(n);
  for (unsigned int i = 0; i < n; ++i) {
    translated[i] = translator[classif[i] - mini];
  }

  return Rcpp::List::create(Rcpp::Named("index") = index_vec,
                            Rcpp::Named("translator") = translator,
                            Rcpp::Named("translated") = translated);
}

// ====================================================================
//
// Sort pairs by bucket sorting : O(n) implementation
//
//

// Use a simple struct for better readability than std::pair
struct Pair {
  int c1, c2;
};

// [[Rcpp::export]]
Rcpp::List std_sort_pairs(Rcpp::IntegerVector& c1_in,
                          Rcpp::IntegerVector& c2_in, unsigned int N1,
                          unsigned int N2) {
  // ----------------------------------------------------------
  // Initialization
  //
  unsigned int n = c1_in.size();

  // Counting (Single pass over input)
  std::vector<int> count_c1(N1, 0), count_c2(N2, 0);
  for (unsigned int i = 0; i < n; ++i) {
    count_c1[c1_in[i]]++;
    count_c2[c2_in[i]]++;
  }

  // ----------------------------------------------------------
  // Sorting c1 and c2 with Radix
  //

  // 2. Offsets (Prefix Sum)
  std::vector<int> shift_c1(N1), shift_c2(N2);
  std::exclusive_scan(count_c1.begin(), count_c1.end(), shift_c1.begin(), 0);
  std::exclusive_scan(count_c2.begin(), count_c2.end(), shift_c2.begin(), 0);

  // Radix Sort using a single vector of Structs for intermediate
  // This improves cache locality significantly during the sort moves
  std::vector<Pair> buffer1(n), buffer2(n);

  // Pass 1: Sort by c2 into buffer1
  for (unsigned int i = 0; i < n; ++i) {
    buffer1[shift_c2[c2_in[i]]++] = {c1_in[i], c2_in[i]};
  }

  // Pass 2: Sort by c1 into buffer2
  for (unsigned int i = 0; i < n; ++i) {
    buffer2[shift_c1[buffer1[i].c1]++] = buffer1[i];
  }

  // Counting Unique Pairs (Single pass over sorted buffer2)
  std::vector<int> pair_c1, pair_c2, count_pair;
  pair_c1.reserve(n);  // Reserve to avoid reallocations
  pair_c2.reserve(n);
  count_pair.reserve(n);

  for (unsigned int i = 0; i < n; ++i) {
    if (i == 0 || buffer2[i].c1 != buffer2[i - 1].c1 ||
        buffer2[i].c2 != buffer2[i - 1].c2) {
      pair_c1.push_back(buffer2[i].c1);
      pair_c2.push_back(buffer2[i].c2);
      count_pair.push_back(1);
    } else {
      count_pair.back()++;
    }
  }

  // ----------------------------------------------------------
  // Preparing Output
  //

  std::vector<int> count_c1_f, count_c2_f;
  for (int c : count_c1)
    if (c > 0) count_c1_f.push_back(c);
  for (int c : count_c2)
    if (c > 0) count_c2_f.push_back(c);

  return Rcpp::List::create(Rcpp::Named("count_pair") = count_pair,
                            Rcpp::Named("count_c1") = count_c1_f,
                            Rcpp::Named("count_c2") = count_c2_f,
                            Rcpp::Named("pair_c1") = pair_c1,
                            Rcpp::Named("pair_c2") = pair_c2);
}

// ====================================================================
//
// Compute the expected Mutual Information between two classification
//

// [[Rcpp::export]]
double expected_MI(const Rcpp::IntegerVector& ni_r,
                   const Rcpp::IntegerVector& nj_r) {
  int n_rows = ni_r.size();
  int n_cols = nj_r.size();

  //  Direct pointer access (bypasses Rcpp operator[] bounds checking)
  const int* ni_ptr = INTEGER(ni_r);
  const int* nj_ptr = INTEGER(nj_r);

  int64_t N = 0;
  for (int i = 0; i < n_rows; ++i) N += ni_ptr[i];

  if (N == 0) return 0.0;
  double N_double = static_cast<double>(N);
  double log_N = std::log(N_double);

  // Pre-calculate logs AND log-factorials in one pass using std::vector
  std::vector<double> log_fact(N + 1);
  std::vector<double> log_ints(N + 1);
  log_fact[0] = 0.0;
  log_ints[0] = -INFINITY;  // log(0) is undefined, but helps catch errors

  for (int i = 1; i <= N; ++i) {
    double val = static_cast<double>(i);
    log_ints[i] = std::log(val);
    log_fact[i] = log_fact[i - 1] + log_ints[i];
  }

  double emi = 0.0;
  double log_fact_N = log_fact[N];

  for (int i = 0; i < n_rows; ++i) {
    int ni = ni_ptr[i];
    if (ni <= 0) continue;

    double term_ni = log_fact[ni] + log_fact[N - ni];

    for (int j = 0; j < n_cols; ++j) {
      int nj = nj_ptr[j];
      if (nj <= 0) continue;
      int start_nij = std::max(1, ni + nj - static_cast<int>(N));
      int end_nij = std::min(ni, nj);

      // Move as much as possible out of the nij loop
      double common_term =
          term_ni + log_fact[nj] + log_fact[N - nj] - log_fact_N;
      double log_ni_nj = log_ints[ni] + log_ints[nj];

      for (int nij = start_nij; nij <= end_nij; ++nij) {
        // Hypergeometric probability using pre-calculated STL vector
        double log_p_nij =
            common_term -
            (log_fact[nij] + log_fact[ni - nij] + log_fact[nj - nj] +
             log_fact[nj - nij] + log_fact[N - ni - nj + nij]);

        double p_nij = std::exp(log_p_nij);

        // Optimized MI term: No calls to std::log here
        double mi_term = (static_cast<double>(nij) / N_double) *
                         (log_ints[nij] + log_N - log_ni_nj);

        emi += p_nij * mi_term;
      }
    }
  }

  return emi;
}
