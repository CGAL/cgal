// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

// Note: whenever a comment refers to "Khachiyan's paper" then the
// paper "Rounding of polytopes in the real number model of
// computation" is ment (Mathematics of Operations Research, Vol. 21,
// No. 2, May 1996).  Nontheless, most comments refer to the
// accompanying documentation sheet (and not to the above paper), see
// the file(s) in documentation/.

#ifndef CGAL_KHACHIYAN_APPROXIMATION_H
#define CGAL_KHACHIYAN_APPROXIMATION_H

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Approximate_min_ellipsoid_d/Approximate_min_ellipsoid_d_configure.h>

// If CGAL is not being used, we need to define certain things:
#ifdef CGAL_JUST_MINIBALL
  namespace CGAL {
    struct Tag_true {};
    struct Tag_false {};
  }
#endif

namespace CGAL {

  // In order not to pollute the CGAL namespace with implementation
  // fuzz, we put such things into the following attic:
  namespace Appel_impl {

    template<typename FT>
    struct Selector {
      typedef Tag_true Is_exact;
    };

    template<>
    struct Selector<float> {
      typedef Tag_false Is_exact;
    };

    template<>
    struct Selector<double> {
      typedef Tag_false Is_exact;
    };

  } // namespace Appel_impl

  template<bool Embed,class Traits>
  class Khachiyan_approximation {
  public: // types:
    typedef typename Traits::FT FT;
    typedef typename Traits::ET ET;
    typedef typename Traits::Point Point;
    typedef typename Traits::Cartesian_const_iterator C_it;
    typedef typename Appel_impl::Selector<FT>::Is_exact Exact_flag;

  private: // member variables and invariants for them:
    Traits& tco;                    // traits class object
    std::vector<const Point *> P;   // input points
    int n;                          // number of input points, i.e., P.size()

    // This class comes in two flavours:
    //
    //  (i) When Embed is false, the input points are taken to be
    //    ordinary points in R^{d_P}, where d_P is the dimension of the
    //    input points (as obtained by Traits::dimension()).  In this
    //    case the "dimension of the ambient space" (see variable d
    //    below) simply equals d_P.
    //
    //  (ii) When Embed is true, the input points are embedded by
    //    means of p -> (p,1) into R^{d_P+1}.  In this case, the
    //    dimension of the amient space is d_P+1.
    //
    // Notice that (at almost all places in this code) whenever a
    // point "p" is mentioned in a comment then the embedded point is
    // ment in case (ii).
    const int d_P;                  // dimension of the input points
    const int d;                    // dimension of the ambient space

    // An instance of this class will compute (provided the points are
    // not degenerate, see is_degenerate() below) a solution x to
    // Khachiyan's program (D) such that the "relaxed optimality
    // conditions"
    //
    //      p_i^T M(x)^{-1} p_i <= (1+desired_eps) d             (*)
    //
    // hold for all input points p_i.  (Notice that the p_i in this
    // case are the embedded points when Embed is true.)  Khachiyan's
    // lemma then guarantees that the ellipsoid defined by the matrix
    // M(x)^{-1} is a "good" approximation in some sense (for further
    // details on "good" refer to is_valid() below).
    //
    // The current epsilon for which (*) holds is stored in the
    // following variable; the variable is only defined if is_deg is
    // false.
    FT eps;

    // In addition to the current epsilon (see variable eps) above, we
    // also maintain the exact epsilon for which (*) holds; this exact
    // epsilon is eps_exact. (Notice that (*) for the variable eps
    // only holds modulo rounding-errors if FT is an inexact type.)
    // In order not to unnecessarily recompute eps_exact, we keep the
    // flag is_exact_eps_uptodate which is true if and only if the
    // value in eps_exact is correct.
    double eps_exact;
    bool is_exact_eps_uptodate;

    // Khachian's algorithm is only guaranteed to work when the input
    // points linearly span the whole space (i.e., if dim(span(P)) =
    // d, or, equivalently, when the smallest enclosing ellipsoid has
    // a volume > 0).  The user of this class must call
    // is_degenerate() which returns false iff dim(span(P)) = d;
    // is_degenerate() simply looks up the following variable:
    bool is_deg;

    #ifdef CGAL_APPEL_ASSERTION_MODE
    // We store the current solution x to program (D) as a n-vector.
    // If is_deg is true, the variable x has no meaning.
    std::vector<FT> x;
    #endif // CGAL_APPEL_ASSERTION_MODE

    // The inverse of the (d x d)-matrix M(x) is stored in the
    // variable mi.  The element (M(x)^{-1})_{ij} is stored in
    // mi[i+j*d].  (It would actually be enough to only store the
    // upper half of the matrix since it's symmetric -- but we don't
    // care.)  If is_deg is true, the variable mi is undefined.
    //
    // Todo: optimize by only storing one half of mi.
    std::vector<FT> mi;

    // The following variable sum is used to build, during the
    // initialization phase, the (d x d)-matrix M(x).  The variable is
    // only defined if is_deg is true; once is_deg becomes false, the
    // variable sum isn't used any more (except in routines for
    // debugging and statistics).  If is_deg is true then sum represents
    // the martix
    //
    //    sum(p_i p_i^T,i=0..(n-1));
    //
    // the (i,j)th element of which is stored as sum[i+j*d].  (It
    // would actually be enough to only store the upper half of the
    // matrix since it's symmetric -- but we don't care.)
    //
    // Todo: optimize by only storing one half of sum.
    std::vector<FT> sum;

    // Khachiyan's algorithm heavily relies on the excess ex[i] of
    // input point P[i] w.r.t. the current solution x (see routine
    // excess() below), i in {0,...,n-1}.  For the first iteration, we
    // compute these numbers explicitly, which is quite costly, namely
    // O(d^2) per excess, for a total of O(n d^2).  However, in all
    // remaining iterations we only update the excesses which can be
    // done in O(n d).  For this we need the excesses of the previous
    // round, and that's why we store them.  Thus, if is_deg is false
    // then ex[i] is the excess of point P[i] w.r.t. the current
    // solution x:
    std::vector<FT> ex;

    // We remember the index (into ex) of a largest element in ex: If
    // is_deg is false then
    //
    //    ex_max = argmax_{0 <= i < n} ex[i].
    //
    int ex_max;

    // We sometimes need temporary storage.  We allocate this once at
    // instance construction time:
    mutable std::vector<FT> tmp;
    mutable std::vector<FT> t;

  public: // member variables for statistics

    #ifdef CGAL_APPEL_STATS_MODE
    // The value max_error_m is the maximal absolute value in an
    // entry of (M(x) M(x)^{-1} - I), where I is the identity matrix.
    // Of course, in theory, M(x) M(x)^{-1} - I should be the zero
    // matrix, but if FT is an inexact type (i.e., double), then this
    // need not be the case anymore, due to rounding errors.
    //
    // The variable max_error_m_all is the maximum over all
    // max_error values ever encountered during the computation.
    FT max_error_m, max_error_m_all;
    #endif // CGAL_APPEL_STATS_MODE

  private: // internal helper routines:

    template<typename NumberType, typename InputIterator>
    NumberType excess(InputIterator p)
      // Computes p^T M(x)^{-1} p, where x is the current solution.
      // (When Embed is true, then the p in the previous sentence's
      // formula is an embedded point...)
      //
      // Complexity: O(d^2)
      //
      // Note: this routine is parametrized by a number-type argument
      // because we will use it in exact_epsilon() and is_valid()
      // below with the exact number-type ET (which may be different
      // from the possibly inexact number-type FT).
    {
      typedef NumberType NT;
      NT result(0);
      InputIterator qi(p);

      for (int i=0; i<d; ++i) {

        // compute i-th entry of the vector M(x)^{-1} p into tmp:
        NT tmp(0);
        InputIterator q(p);
        for (int j=0; j<d_P; ++j, ++q)
          tmp += NT(*q) * NT(mi[i+j*d]);
        if (Embed)
          tmp += NT(mi[i+d_P*d]);

        // add tmp*p_i to result:
        if (!Embed || i < d_P) {
          result += NT(*qi) * NT(tmp);
          ++qi;
        } else
          result += NT(tmp);
      }

      return result;
    }

    void update_sum(int start)
      // Adds to the variable sum the matrices P[i] P[i]^T where i
      // ranges from start to n-1.
      //
      // Complexity: O((n-start) d^2)
      //
      // Todo: maybe use something like Kahan Summation here? See
      // <https://en.wikipedia.org/wiki/Kahan_Summation_Algorithm>.
    {
      for (int k=start; k<n; ++k) {
        C_it pi = tco.cartesian_begin(*P[k]);
        for (int i=0; i<d_P; ++i, ++pi) {
          C_it pj = tco.cartesian_begin(*P[k]);
          for (int j=0; j<d_P; ++j, ++pj)
            sum[i+j*d] += (*pi) * (*pj);
          if (Embed)
            sum[i+d_P*d] += (*pi);
        }

        if (Embed) {
          C_it pj = tco.cartesian_begin(*P[k]);
          for (int j=0; j<d_P; ++j, ++pj)
            sum[d_P+j*d] += (*pj);
          sum[d_P+d_P*d] += 1;
        }
      }
    }

  public: // construction & destruction:

    Khachiyan_approximation(int dim,int n_est,Traits& tco) :
      // An instance of this class represents, for some given point
      // set P, a solution x to program (D) which satisfies the
      // relaxed optimality conditions (*) for some desired_epsilon
      // (to be specified by the user at a later point).  After
      // calling this constructor, the instance will represent the
      // (degenerate) solution to program (D) for P={}.  By adding
      // points to P (see routine add()), you can compute an
      // approximate solution to (D) for any point set P.
      //
      // (A solution does not always exist; the set P must linearly
      // span the whole space in order for a solution to exist.  See
      // routine is_degenerate() and the result of routine add() for
      // more information.)
      //
      // The number n_est is an estimate on the number of points you
      // are going to add; n_est need not coincide with the actual
      // number of points you will eventually have added.  It merely
      // gives a hint on how much storage is needed.
      tco(tco), n(0),
      d_P(dim), d(Embed? d_P+1 : d_P), is_deg(true),
      #ifdef CGAL_APPEL_ASSERTION_MODE
      x(n_est),
      #endif // CGAL_APPEL_ASSERTION_MODE
      mi(d*d), sum(d*d), ex(n_est), tmp(d), t(d*d)
    {
      CGAL_APPEL_LOG("appel","Entering Khachiyan with d=" << d << " (" <<
                (Embed? "" : "not ") << "embedded)." << std::endl);
      CGAL_APPEL_TIMER_START("khachiyan");

      // In order to satisfy the invariant on m, we have to initalize
      // m with the zero matrix:
      for (int i=0; i<d; ++i)
        for (int j=0; j<d; ++j)
          sum[i+j*d] = FT(0);
    }

    ~Khachiyan_approximation();

    template<typename InputIterator>
    bool add(InputIterator first,InputIterator last,double desired_eps)
      // Adds the points from the range [first,last) to the instance's
      // set P and computes a solution x to program (D) satisfying the
      // relaxed optimality conditions (*) for the given value
      // desired_eps.  (Khachiyan's lemma then guarantees that the
      // ellipsoid defined by the matrix M(x)^{-1} is a "good"
      // approximation in some sense ...)
      //
      // Returned is false if and only if the smallest enclosing
      // ellipsoid of all points added so far is degenerate (see
      // routine is_degenerate()), in which case a solution to program
      // (D) has not been computed (because it doesn't exist).
      //
      // Observe that the points you add must have dimension d in any
      // case.  If Embed is true, then an input point p is internally
      // interpreted as (p,1) (a (d+1)-dimensional point); if Embed is
      // false, the point is simply left as it is.
    {
      // In the following we denote the points which have been added
      // so far (not including the points we are about to add) by
      // P_old.
      CGAL_APPEL_LOG("appel","  add() called with desired_eps=" <<
                desired_eps << std::endl);

      // We first add the new points:
      const int n_old = n;
      for (; first != last; first++)
        P.push_back(&*first);
      n = static_cast<int>(P.size());

      // debugging output:
      CGAL_APPEL_LOG("appel","  Add()'ing " << n-n_old << " points to a " <<
                (is_deg? "" : "not ") << "degenerate set of " <<
                n_old << " points." << std::endl);
      #ifdef CGAL_APPEL_ASSERTION_MODE
      x.resize(n);
      #endif // CGAL_APPEL_ASSERTION_MODE
      ex.resize(n);

      if (is_deg) {
        // Here, we are still collecting points until is_deg gets true.

        // update variable sum:
        update_sum(n_old);

        // Next, we need to check whether the matrix M(x)^{-1} for
        // x=(1/n,...,1/n) exists (in which case we could start
        // Khachiyan's algorithm): compute_initial_inverse() computes
        // from the matrix sum the inverse of M(x), if possible:
        is_deg = !compute_initial_inverse_from_sum();
      }

      else {
        // Here, is_deg is false, i.e., the points P_old already span
        // the whole space (and by our invariants the variable x
        // therefore already represents a solution to the program (D)
        // for the points P_old).

        // The only thing we need to update is the array ex and ex_max:
        for (int i=n_old; i<n; ++i)
          if ((ex[i] = excess<FT>(tco.cartesian_begin(*P[i]))) > ex[ex_max])
              ex_max = i;
        CGAL_APPEL_LOG("appel","  Maximal excess is " << to_double(ex[ex_max])
                  << "." << std::endl);

        // (Here, we could restart with x = (1/n,...,1/n).  But this
        // is probably a bad idea because (i) computing the initial
        // inverse of M(x) is expensive and (ii) because our current
        // ellipsoid is, hopefully, quite good already.)

        // Todo: If the needed eps is larger than n-1, then
        // restarting with x = (1/n,...,1/n) is better, probably (at
        // least, we have a better bound then).
      }

      // Finally, if all points are non-degenerate now, we run
      // Khachiyan's algorithm:
      CGAL_APPEL_LOG("appel","  The points are " << (is_deg? "" : "not ") <<
                "degenerate." << std::endl);
      if (!is_deg)
        run(desired_eps);

      return !is_deg;
    }

  private:
    void run(const double desired_eps)
      // Runs Khachiyans algorithm, provided the points added so far
      // are non-degnerate.  On return, the variable x will be a
      // solution to program (D) for the points P such that the
      // relaxed optimality conditions (*) hold.
      //
      // Precondition: !is_deg
    {
      CGAL_APPEL_ASSERT(!is_deg);

      // compute upper bound on needed iterations; the bound only
      // holds from the moment on when eps is <= 1.
      //
      // Todo: the bound does not really say anything, because all it
      // says is that under exact computation and this many iterations
      // we have reached the desired approximation ratio. But as we
      // are (possibly) using inexact arithmetic, the bound does not
      // hold. However, we do not promise anything in the
      // specification, so we can stop after this many iterations (and
      // the user will query exact_epsilon() to see whether or not we
      // reached the desired approximation ratio). So this
      // implementation DOES meet the specification, but its sort of
      // clumsy.
      const int max_iterations = static_cast<int>(2*d/desired_eps/
                                         (std::log(1.5)-1.0/3.0));

      // run Khachiyan's algorithm until we find a good enough eps:
      //
      // Todo: add a better termination guarantee
      int iterations = 1;
      bool counting = eps <= 1.0;
      while (!improve(desired_eps) &&        // not good enough and ...
             (!counting ||                   // ... (eps not yet <= 1 or ...
              iterations <= max_iterations)) // ... not yet enough iterations)
        if (counting)
          ++iterations;
        else if (eps <= 1.0)
          counting = true;

      // output stats:
      CGAL_APPEL_LOG("appel",
                "  Took " << iterations << " iterations (upper " <<
                "bound was " << max_iterations << "." << std::endl);
      CGAL_APPEL_TIMER_PRINT("appel","khachiyan","  Time at end: ");
      CGAL_APPEL_IF_STATS(CGAL_APPEL_LOG("appel",
                  "  The overall represenation error in m is:   " <<
                     max_error_m_all << "." << std::endl);)
    }

  public: // access:

    bool is_degenerate() const
    // Returns true if and only if dim(span(P)) < d, with P being the
    // input points (i.e., the points added so far) and d the
    // dimension of the ambient space.  In other words, the routine
    // returns false if and only if the smallest enclosing ellipsoid
    // of all points added so far has a volume > 0 in R^d.
    {
      return is_deg;
    }

    const FT& matrix(int i,int j) const
    // Returns the (i,j)-element of the matrix M(x)^{-1}.
    //
    // Precondition: 0<=i<d && 0<=j<d
    {
      CGAL_APPEL_ASSERT(0<=i && i<d && 0<=j && j<d);
      return mi[i+d*j];
    }

  public: // convenience members:

    // Computes the inverse of the (d-1xd-1)-matrix defined by the above
    // member matrix(i,j),
    //
    //   [ matrix(i,j) ]_{0<=i,j<d-1},
    //
    // into the matrix given by the Iterator inverse. On exit,
    // inverse[i+(d-1)*j] contains the (i,j)-element of the
    // inverse. If the return value is false, stability problems
    // prevented the computation from being successfully completed.
    //
    // Precondition: !is_deg
    template<typename Iterator>
    bool compute_inverse_of_submatrix(Iterator inverse);

  public: // validity:

    // Returns a double-precision upper bound on the exact
    // approximation ratio. I.e., if FT is an exact number type, this
    // routine returns eps; if FT is inexact, it computes (see (*))
    //
    //    tmp:= max p_i^T M(x)^{-1} p_i
    //
    // using exact arithmetic and returns
    //
    //    epsilon = rnd_down(tmp/d - 1)
    //
    // where rnd_down(x) is the largest double smaller or equal to x.
    //
    // In other words, Eq. (*) will hold in exact arithmetic for the
    // epsilon returned by this routine (something which is not true
    // for the current epsilon eps, where (*) only holds modulo
    // rounding errors).
    //
    // Note: if recompute is false, the routine is allowed to use the
    // cached result from a previous invocation (if still valid).  In
    // is_valid() (see below) we will call the routine with recompute
    // equal to true (simply because a validation routine should not
    // rely on cached results...).
    //
    // Note: it turns out that under floating-point arithmetic (i.e.,
    // when FT is double), the algorithm does not always correctly
    // decide whether or not the instance is degenerate. For example,
    // it may happen that the instance is "seen" to be nondegenerate
    // (is_degenerate() returns false) but still the points are
    // degenerate; I have seen this in 2D (embedded in 3D) with two
    // points, each having multiplicity two (meaning n=4). In such a
    // situation, exact_epsilon() may return a negative number (and we
    // should then conclude that something went wrong, of course).
    //
    // Complexity: O(n d^2)
    //
    // Precondition: !is_deg
    double exact_epsilon(bool recompute = false);

    // Checks that the computed solution indeed satisfies (*) for the
    // epsilon obtained by exact_epsilon() above.  In particular, this
    // shows that the computed ellipsoid is indeed a
    // (1+exact_epsilon())-approximation to the smallest centrally
    // symmetric ellipsoid of P.
    //
    // The test is carried out using exact arithmetic.
    //
    // Complexity: O(n d^2)
    bool is_valid(bool verbose);

  private: // internal assertion routines (only avalable in debug mode):

    #ifdef CGAL_APPEL_ASSERTION_MODE
    void compute_M_of_x();
    // Given the current solution x, this routine computes the matrix
    // M(x) into the matrix represented by the variable t.

    FT representation_error_of_mi();
    // Computes a number indicating whether the matrix stored in the
    // member variable mi is the inverse of the matrix M(x): If FT is
    // an exact number type, then this is the case if and only if
    // representation_error_of_mi() returns 0.  For an inexact number
    // type FT, the returned number gives only a "feeling" for the
    // error.
    //
    // Note: This is a _very_ slow routine and should only be used
    // for debugging purposes.
    #endif // CGAL_APPEL_ASSERTION_MODE

  private: // internal assertion checking routines:

    void print(std::ostream& o);
    // Outputs the internally used matrices to stream o.

  public: // internal debugging routines:

    int write_eps() const;
    // Writes the current approximation and the point set P to an EPS
    // file.  Returned is the id under which the file was stored
    // (filename 'id.eps').
    //
    // Precondition: d == 2.

  private: // internal routines:

    bool compute_inverse_of_t_into_mi(const Tag_true exact);
    // Computes the inverse of the matrix t into variable mi and
    // returns true iff the inverse exists.
    //
    // Note that the algorithm doing the inversion (its Gauss, more or
    // less) is not well-suited for inexact number types FT, e.g., it
    // should not be used with double.
    //
    // Complexity: O(d^3)

    bool compute_inverse_of_t_into_mi(const Tag_false exact);
    // Computes the inverse of the matrix t into variable mi and
    // returns true iff the inverse exists.
    //
    // Note (todo): This routine works well for double arithmetic, I
    // haven't tested it for anything else (and some epsilons should
    // probably be adjusted for other floating-point types).  It uses
    // the (very) stable Cholesky decomposition internally.
    //
    // Complexity: O(d^3)

    bool compute_initial_inverse_from_sum();
    // Initializes the current solution x to (1/n,...,1/n), computes
    // into mi the inverse of M(x) and sets eps to some (trivial)
    // value such that the relaxed optimiality conditions (*) hold.
    // In addition, this routine intializes the array ex to hold the
    // excesses (w.r.t. this trivial solution) of the input points and
    // initializes ex_max to the index into ex of a largest element.
    // The routine returns true if and only if M(x)^{-1} exists.
    //
    // Complexity: O(d^3 + n d^2)

    void rank_1_update(int k,const FT& tau);
    // Efficiently computes into mi the inverse of M(x'), where
    //
    //   x' = (1 - tau) x + tau e_k
    //
    // with tau = eps/((1+eps)d-1) and where eps satisfies
    //
    //   p_k^T M(x)^{-1} p_k = (1+eps) d.                       (**)
    //
    // In addition, the routine efficiently updates the array ex to
    // reflect the excesses w.r.t. the new solution x.
    //
    // Complexity: O(d^2 + n d)
    //
    // Precondition: tau = eps/((1+eps)d-1) and (**) holds.

    bool improve(const double desired_eps);
    // Given the current solution x to program (D), improve() finds
    // better solution to (D) by applying one iteration of the
    // barycentric coordinate descent method.  The routine returns
    // true if the solution satisfies the relavex optimality
    // conditions (*) for some eps <= desired_eps.
    //
    // Complexity: O(d^2 + n d)
  };

}

#include <CGAL/Approximate_min_ellipsoid_d/Khachiyan_approximation_impl.h>

#endif // CGAL_KHACHIYAN_APPROXIMATION_H
