/*
 * This file is based on the fgsea package (https://github.com/ctlab/fgsea)
 * Copyright (c) 2016-2024 Alexey Sergushichev
 *
 * It has been adapted for the enrichit package.
 */

#ifndef ESCALCULATION_H
#define ESCALCULATION_H

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

using namespace std;

// Assuming 64-bit environment with __int128 support
// If compilation fails on 32-bit, we might need a fallback or Boost dependency
__extension__ typedef __int128 int128;
__extension__ typedef unsigned __int128 uint128;

struct score_t {
  //  score = coef_ns / NS - coef_const / diff
  int64_t NS;
  int64_t coef_NS;
  int64_t diff;   //  n - k
  int64_t coef_const;

  double getDouble() const {
      return 1.0 * coef_NS / NS - 1.0 * coef_const / diff;
  }

  int64_t getNumerator() const {
    return coef_NS * diff - coef_const * NS;
  }

  int128 Compare(score_t const& other) const {
    int64_t P1 = coef_NS * diff + NS * (other.coef_const - coef_const);
    int64_t Q1 = NS * diff;
    int64_t P2 = other.coef_NS;
    int64_t Q2 = other.NS;
    return int128(P1) * Q2 - int128(P2) * Q1;
  }

  bool operator<(score_t const& other) const {
    return Compare(other) < 0;
  }

  bool operator<=(score_t const& other) const {
    return Compare(other) <= 0;
  }

  bool operator>(score_t const& other) const {
    return Compare(other) > 0;
  }

  bool operator>=(score_t const& other) const {
    return Compare(other) >= 0;
  }

  bool operator==(score_t const& other) const {
    return Compare(other) == 0;
  }

  bool operator!=(score_t const& other) const {
    return Compare(other) != 0;
  }

  static int64_t getMaxNS() {
      return int64_t(1 << 30);
  }

  score_t operator-() const {
    return score_t{NS, -coef_NS, diff, -coef_const};
  }

  score_t abs() const {
    return std::max(*this, -(*this));
  }
};

score_t calcES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS);
score_t calcES(const vector<int64_t> &ranks, const vector<int> &p);
score_t calcPositiveES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS);
score_t calcPositiveES(const vector<int64_t> &ranks, const vector<int> &p);
pair<score_t, score_t> calcSignedES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS);
pair<score_t, score_t> calcSignedES(const vector<int64_t> &ranks, const vector<int> &p);
bool compareStat(const vector<double> &ranks, const vector<int> &p, double NS, double bound);

#endif // ESCALCULATION_H
