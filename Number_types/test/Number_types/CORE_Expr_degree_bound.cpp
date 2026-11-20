#include <CGAL/config.h>

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#include <mpfr.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

// Adapted from Phillip Keldenich's 66-term reproducer in CGAL issue #9708.
// https://github.com/CGAL/cgal/issues/9708
using Expr = CORE::Expr;
using Point = std::pair<double, double>;

struct Term
{
  int weight; // +1 or -1
  Point s, t;
};

// BEGIN DATA
static const std::vector<Term> terms = {
    {+1, {0x1.c6d04eafa8499p+25, 0x1.e4e444feca712p+28}, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}},
    {+1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {+1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.2e6d642c33e46p+26, 0x1.006dfac05dffbp+29}},
    {+1, {0x1.2e6d642c33e46p+26, 0x1.006dfac05dffbp+29}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {+1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca766p+28}},
    {+1, {0x1.2e6d642c33e46p+26, 0x1.006dfac05dffbp+29}, {0x1.3aaf598c2635fp+26, 0x1.004dde8aeb0f0p+29}},
    {+1, {0x1.3aaf598c2635fp+26, 0x1.004dde8aeb0f0p+29}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {+1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.c6d04eafa84a9p+25, 0x1.e4e444feca748p+28}},
    {+1, {0x1.c6d04eafa84a9p+25, 0x1.e4e444feca748p+28}, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca766p+28}},
    {+1, {0x1.3aaf598c2637bp+26, 0x1.004dde8aeb0f2p+29}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {+1, {0x1.3aaf598c2637bp+26, 0x1.004dde8aeb0f2p+29}, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}},
    {+1, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {+1, {0x1.bdbf2307987f1p+26, 0x1.d328de30f9011p+28}, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}},
    {+1, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}, {0x1.1fde988161026p+25, 0x1.8ccadae164bfep+28}},
    {+1, {0x1.bdbf2307987f1p+26, 0x1.d328de30f9011p+28}, {0x1.bdbf2307987ecp+26, 0x1.d328de30f8fdep+28}},
    {+1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161026p+25, 0x1.8ccadae164bfep+28}},
    {+1, {0x1.bdbf2307987f1p+26, 0x1.d328de30f9011p+28}, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}},
    {+1, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}, {0x1.bdbf2307987ecp+26, 0x1.d328de30f8fdep+28}},
    {+1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161035p+25, 0x1.8ccadae164bfdp+28}},
    {+1, {0x1.bdbf2307987f6p+26, 0x1.d328de30f900fp+28}, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}},
    {+1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf7p+28}},
    {+1, {0x1.bdbf2307987f6p+26, 0x1.d328de30f900fp+28}, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}},
    {+1, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}},
    {+1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf3p+28}},
    {+1, {0x1.bdbf230798825p+26, 0x1.d328de30f900fp+28}, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}},
    {+1, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}, {0x1.bdbf23079880fp+26, 0x1.d328de30f8fffp+28}},
    {+1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}},
    {+1, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf3p+28}},
    {+1, {0x1.84179d78b8a40p+26, 0x1.6ffadedc148bap+28}, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}},
    {+1, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}},
    {+1, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf3p+28}},
    {+1, {0x1.1fde98816101ep+25, 0x1.8ccadae164bc6p+28}, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}},
    {+1, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}, {0x1.1fde98816103bp+25, 0x1.8ccadae164be9p+28}},
    {-1, {0x1.c6d04eafa84a2p+25, 0x1.e4e444feca725p+28}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {-1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {-1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.2e6d642c33e46p+26, 0x1.006dfac05dffbp+29}},
    {-1, {0x1.2e6d642c33e46p+26, 0x1.006dfac05dffbp+29}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {-1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca766p+28}},
    {-1, {0x1.2e6d642c33e46p+26, 0x1.006dfac05dffbp+29}, {0x1.3aaf598c2635fp+26, 0x1.004dde8aeb0f0p+29}},
    {-1, {0x1.3aaf598c2635fp+26, 0x1.004dde8aeb0f0p+29}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {-1, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca737p+28}, {0x1.c6d04eafa84a9p+25, 0x1.e4e444feca748p+28}},
    {-1, {0x1.c6d04eafa84a9p+25, 0x1.e4e444feca748p+28}, {0x1.c6d04eafa84abp+25, 0x1.e4e444feca766p+28}},
    {-1, {0x1.3aaf598c2637bp+26, 0x1.004dde8aeb0f2p+29}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {-1, {0x1.3aaf598c2637bp+26, 0x1.004dde8aeb0f2p+29}, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}},
    {-1, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}, {0x1.c6d04eafa8497p+25, 0x1.e4e444feca705p+28}},
    {-1, {0x1.bdbf2307987f1p+26, 0x1.d328de30f9011p+28}, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}},
    {-1, {0x1.bdbf2307987e9p+26, 0x1.d328de30f8fccp+28}, {0x1.1fde988161026p+25, 0x1.8ccadae164bfep+28}},
    {-1, {0x1.bdbf2307987f1p+26, 0x1.d328de30f9011p+28}, {0x1.bdbf2307987ecp+26, 0x1.d328de30f8fdep+28}},
    {-1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161026p+25, 0x1.8ccadae164bfep+28}},
    {-1, {0x1.bdbf2307987f1p+26, 0x1.d328de30f9011p+28}, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}},
    {-1, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}, {0x1.bdbf2307987ecp+26, 0x1.d328de30f8fdep+28}},
    {-1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161035p+25, 0x1.8ccadae164bfdp+28}},
    {-1, {0x1.bdbf2307987f6p+26, 0x1.d328de30f900fp+28}, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}},
    {-1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf7p+28}},
    {-1, {0x1.bdbf2307987f6p+26, 0x1.d328de30f900fp+28}, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}},
    {-1, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}, {0x1.bdbf2307987f4p+26, 0x1.d328de30f8ffbp+28}},
    {-1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf3p+28}},
    {-1, {0x1.bdbf230798825p+26, 0x1.d328de30f900fp+28}, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}},
    {-1, {0x1.bdbf230798820p+26, 0x1.d328de30f9003p+28}, {0x1.bdbf23079880fp+26, 0x1.d328de30f8fffp+28}},
    {-1, {0x1.84179d78b8a66p+26, 0x1.6ffadedc148bdp+28}, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}},
    {-1, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf3p+28}},
    {-1, {0x1.84179d78b8a40p+26, 0x1.6ffadedc148bap+28}, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}},
    {-1, {0x1.1fde988161045p+25, 0x1.8ccadae164bbap+28}, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}},
    {-1, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}, {0x1.1fde988161046p+25, 0x1.8ccadae164bf3p+28}},
    {-1, {0x1.1fde98816101ep+25, 0x1.8ccadae164bc6p+28}, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}},
    {-1, {0x1.1fde98816102fp+25, 0x1.8ccadae164bd4p+28}, {0x1.1fde98816103bp+25, 0x1.8ccadae164be9p+28}}};
// END DATA

static Expr length(const Term &term)
{
  const Expr dx = Expr(term.s.first) - Expr(term.t.first);
  const Expr dy = Expr(term.s.second) - Expr(term.t.second);
  return CGAL::sqrt(dx * dx + dy * dy);
}

static Expr signed_length(const Term &term) { return term.weight > 0 ? length(term) : -length(term); }

static Expr pairwise_sum(std::vector<Expr> values)
{
  while (values.size() > 1)
  {
    std::vector<Expr> next;
    for (std::size_t i = 0; i + 1 < values.size(); i += 2)
    {
      next.push_back(values[i] + values[i + 1]);
    }
    if (values.size() % 2 == 1)
    {
      next.push_back(values.back());
    }
    values.swap(next);
  }
  return values.empty() ? Expr(0) : values.front();
}

static int proven_sign(mpfr_prec_t precision)
{
  mpfr_t dx, dy, squared, root_down, root_up, lower, upper;
  mpfr_inits2(precision, dx, dy, squared, root_down, root_up, lower, upper, (mpfr_ptr) nullptr);
  mpfr_set_zero(lower, 1);
  mpfr_set_zero(upper, 1);
  auto exact = [](int ternary)
  {
    if (ternary != 0)
    {
      std::cerr << "inexact intermediate result; increase the precision" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  };
  for (const Term &term : terms)
  {
    exact(mpfr_set_d(dx, term.s.first, MPFR_RNDN));
    exact(mpfr_sub_d(dx, dx, term.t.first, MPFR_RNDN));
    exact(mpfr_set_d(dy, term.s.second, MPFR_RNDN));
    exact(mpfr_sub_d(dy, dy, term.t.second, MPFR_RNDN));
    exact(mpfr_sqr(dx, dx, MPFR_RNDN));
    exact(mpfr_sqr(dy, dy, MPFR_RNDN));
    exact(mpfr_add(squared, dx, dy, MPFR_RNDN));
    mpfr_sqrt(root_down, squared, MPFR_RNDD); // root_down <= length <= root_up
    mpfr_sqrt(root_up, squared, MPFR_RNDU);
    if (term.weight > 0)
    {
      mpfr_add(lower, lower, root_down, MPFR_RNDD);
      mpfr_add(upper, upper, root_up, MPFR_RNDU);
    }
    else
    {
      mpfr_sub(lower, lower, root_up, MPFR_RNDD);
      mpfr_sub(upper, upper, root_down, MPFR_RNDU);
    }
  }
  std::cout << "MPFR enclosure (" << precision << " bits, directed rounding): [";
  mpfr_out_str(stdout, 10, 12, lower, MPFR_RNDD);
  std::cout << std::flush << ", ";
  mpfr_out_str(stdout, 10, 12, upper, MPFR_RNDU);
  std::cout << "]" << std::endl;
  const int sign = mpfr_sgn(lower) > 0 ? 1 : (mpfr_sgn(upper) < 0 ? -1 : 0);
  mpfr_clears(dx, dy, squared, root_down, root_up, lower, upper, (mpfr_ptr) nullptr);
  return sign;
}

static void check_expression(const Expr &sum, int expected)
{
  const Expr expression = (Expr(1) + sum) - Expr(1);
  const int wrapped_sign = expression.Rep()->getExactSign();
  assert(wrapped_sign == expected);
  CORE::ExprRep *rep = expression.Rep();
  std::cout << "wrapped degree: " << rep->d_e() << ", measure: " << rep->measure() << ", bound: " << rep->computeBound()
            << std::endl;
  assert(!rep->measure().isNaN());
  assert(!rep->computeBound().isNaN());
  const int sum_sign = sum.Rep()->getExactSign();
  assert(sum_sign == expected);
}
#endif

int main()
{
#ifdef CGAL_USE_CORE
  const int expected = proven_sign(1024);
  assert(expected == 1);
  Expr left_deep(0);
  std::vector<Expr> signed_lengths, positive, negative;
  for (const Term &term : terms)
  {
    left_deep += signed_length(term);
    signed_lengths.push_back(signed_length(term));
    (term.weight > 0 ? positive : negative).push_back(length(term));
  }
  check_expression(left_deep, expected);
  check_expression(pairwise_sum(signed_lengths), expected);
  check_expression(pairwise_sum(positive) - pairwise_sum(negative), expected);
#endif
}
