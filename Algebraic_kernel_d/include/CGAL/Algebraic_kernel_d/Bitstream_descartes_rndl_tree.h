// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Bitstream_descartes_rndl_tree.h
    \brief Definition of \c NiX::Bitstream_descartes_rndl_tree.
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_H
#define CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_H

#include <vector>
#include <list>
#include <utility>
#include <iterator>
#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/tss.h>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
/*#include <CGAL/Handle.h>
#include <CGAL/Random.h>*/

/*
 *  AUXILIARY CLASSES AND FUNCTIONS
 */

namespace CGAL {

namespace internal {

// TODO: Copied from CGAL/enums.h:
enum Three_valued_estimate {
    CLEARLY_NEGATIVE         = -1, //!< = -1. Sign of a value is clearly
                                   //!<       negative.
    UNCLEAR_SIGN             =  0, //!< =  0. It is unclear whether value
                                   //!<       has negative, zero, or positive
                                   //!<       sign.
    CLEARLY_POSITIVE         =  1, //!< = +1. Sign of a value is clearly
                                   //!<       positive.

    CLEARLY_LESS             = -1, //!< = -1. First value is clearly less
                                   //!<       than second value.
    UNCLEAR_COMPARISON       =  0, //!< =  0. It is unclear whether first value
                                   //!<       is less than, equal to, or
                                   //!<       greater than second value.
    CLEARLY_GREATER          =  1, //!< = +1. First value is clearly greater
                                   //!<       than second value.

    CLEARLY_CLOCKWISE        = -1, //!< = -1. Points in the plane or space are
                                   //!<       clearly oriented clockwise.
    CLEAR_RIGHT_TURN         = -1, //!< = -1. Points in the plane or space
                                   //!<       clearly form a right turn.
    UNCLEAR_ORIENTATION      =  0, //!< =  0. It is not clear whether points in
                                   //!<       the plane or space are oriented
                                   //!<       clockwise, counterclockwise, or
                                   //!<       are in degenerate position.

    UNCLEAR_TURN             =  0, //!< =  0. It is not clear whether points in
                                   //!<       the plane or space form a left
                                   //!<       turn, a right turn, or are in
                                   //!<       degenerate position.
    CLEAR_LEFT_TURN          =  1, //!< = +1. Points in the plane or space
                                   //!<       clearly form a left turn.
    CLEARLY_COUNTERCLOCKWISE =  1, //!< = +1. Points in the plane or space are
                                   //!<       clearly oriented
                                   //!<       counterclockwise.

    CLEARLY_ON_NEGATIVE_SIDE  = -1, //!< = -1. Point is clearly on the negative
                                    //!<       side of an oriented surface.
    ON_UNCLEAR_ORIENTED_SIDE  =  0, //!< =  0. It is unclear whether point is
                                    //!<       outside, on, or inside an
                                    //!<       oriented surface.
    CLEARLY_ON_POSITIVE_SIDE  =  1, //!< = +1. Point is clearly on the positive
                                    //!<       side of an oriented surface.

    CLEARLY_ON_UNBOUNDED_SIDE = -1, //!< = -1. Point is clearly in the
                                    //!<       exterior of a closed surface.
    ON_UNCLEAR_SIDE           =  0, //!< =  0. It is unclear whether point is
                                    //!<       outside, on, or inside a closed
                                    //!<       surface.
    CLEARLY_ON_BOUNDED_SIDE   =  1  //!< =  1. Point is clearly in the
                                    //!<       interior of a closed surface.
};

//! reverses the \a sign (from plus to minus and minus to plus ;-)
inline Three_valued_estimate operator- ( Three_valued_estimate sign ) {
    return Three_valued_estimate( - int( sign));
}


// END: Copied from CGAL/enums.h

/*
 * Helper functions
 */

// compute g(x) = 2^p f((ax+b)/c) for c = 2^log_c with absolute error <= 1
template <class Integer, class RandomAccessIterator, class OutputIterator,
class Approximator, class CeilLog2AbsInteger, class CeilLog2AbsLong>
OutputIterator
polynomial_affine_transform_approx_log_denom(
        RandomAccessIterator first, RandomAccessIterator beyond,
        OutputIterator out,
        Integer a, Integer b, long log_c,
        long p,
        Approximator approx,
        CeilLog2AbsInteger log, CeilLog2AbsLong logl
) {
    // degree of input polynomial
    const int n = int((beyond-first)-1);
    CGAL_precondition(n >= 0);

    // u[i][j] = \binom{i}{j} a^j b^{i-j}  (0 <= i <= n; 0 <= j <= i)
    long max_log_u; // = \max_{i,j} \log\abs{u[i][j]}
    std::vector< std::vector< Integer > > u(n+1);
    u[0].push_back(Integer(1));
    max_log_u = 0; // = log(1)
    for (int i = 1; i <= n; ++i) {
        u[i].reserve(i+1);
        u[i].push_back(b*u[i-1][0]);
        if (u[i][0] != 0) max_log_u = (std::max)(max_log_u, log(u[i][0]));
        for (int j = 1; j <= i-1; ++j) {
            u[i].push_back(a*u[i-1][j-1] + b*u[i-1][j]);
            if (u[i][j] != 0) max_log_u = (std::max)(max_log_u, log(u[i][j]));
        }
        u[i].push_back(a*u[i-1][i-1]);
        if (u[i][i] != 0) max_log_u = (std::max)(max_log_u, log(u[i][i]));
    }

    long q = logl(n+1) + max_log_u + 1;
    Integer half = Integer(1) << q-1;
    std::vector< Integer > h(n+1);
    RandomAccessIterator it = first;
    for (int i = 0; i <= n; ++i, ++it) {
        h[i] = approx(*it, -i*log_c + p + q);
    }
    Integer sum;
    for (int j = 0; j <= n; ++j) {
        sum = 0;
        for (int i = j; i <= n; ++i) {
            sum += u[i][j]*h[i];
        }
        sum += half; sum >>= q;  // round to nearest
        *out++ = sum;
    }
    return out;
} // polynomial_affine_transform_approx_log_denom()


template <class Integer>
Integer caching_factorial(int n) {
    CGAL_precondition(n >= 0);

    // table of factorials; augment if necessary
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(std::vector< Integer >, factorial);
    factorial.reserve(n+1);
    if (factorial.empty()) {
        factorial.push_back(Integer(1)); // 0! = 1
        factorial.push_back(Integer(1)); // 1! = 1
    }
    for (int i = int(factorial.size()); i <= n; ++i) {
        factorial.push_back(i*factorial[i-1]);
    }

    return factorial[n];
}


// Bernst coeff of 2^p n! f(x) wrt [lwr_num,upr_num]/2^log_den w/ abs err <= 1
template <class Integer,
    class RandomAccessIterator, class OutputIterator,
    class Approximator, class CeilLog2AbsInteger, class CeilLog2AbsLong
>
OutputIterator
polynomial_power_to_bernstein_approx(
        RandomAccessIterator first, RandomAccessIterator beyond,
        OutputIterator out,
        Integer lower_num, Integer upper_num, long log_denom,
        long p,
        Approximator approx, CeilLog2AbsInteger log, CeilLog2AbsLong logl
) {
    // degree of input polynomial
    const int n = int((beyond-first)-1);
    CGAL_precondition(n >= 0);

    long q = log(caching_factorial<Integer>(n)) + logl(n+1) + 1;
    Integer half = Integer(1) << q-1;
    std::vector<Integer> f(n+1);
    polynomial_affine_transform_approx_log_denom(
            first, beyond, f.begin(),
            upper_num - lower_num, lower_num, log_denom,
            p+q,
            approx, log, logl
    );
    Integer sum, lprod;
    for (int l = 0; l <= n; ++l) {
        sum = 0;
        lprod = 1; // = l*(l-1)*(l-2)*...*(l-(k-1))
        for (int k = 0; k <= l; ++k) {
            sum += lprod * caching_factorial<Integer>(n-k) * f[k];
            lprod *= l-k;
        }
        sum += half; sum >>= q;  // round to nearest
        *out++ = sum;
    }
    return out;
} // polynomial_power_to_bernstein_approx()


// min/max number of variations in epsilon-sign
template <class InputIterator, class UnaryFunction>
void var_eps( 
        InputIterator first, InputIterator beyond,
        int& min_var, int& max_var,
        const UnaryFunction& sign_eps
) {
    min_var = max_var = 0;
    InputIterator it = first;
    CGAL_precondition(it != beyond);

    internal::Three_valued_estimate last_sign_min, last_sign_max; // always non-zero
    last_sign_min = last_sign_max = sign_eps(*it);
    CGAL_assertion(last_sign_min != internal::UNCLEAR_SIGN);

    while (++it != beyond) {
        internal::Three_valued_estimate cur_sign = sign_eps(*it);
        if (cur_sign == internal::UNCLEAR_SIGN) {
            last_sign_max = -last_sign_max;
            ++max_var;
        } else {
            if (last_sign_max != cur_sign) ++max_var;
            if (last_sign_min != cur_sign) ++min_var;
            last_sign_min = last_sign_max = cur_sign;
        }
    }
} // var_eps()

/*
 * The generic de Casteljau method
 */

template <class NT_>
class Convex_combinator_generic {
public:
    typedef NT_ NT;
private:
    NT alpha_, beta_;
public:
    Convex_combinator_generic(NT alpha) : alpha_(alpha), beta_(NT(1)-alpha) { }
    void into_first (NT& a, const NT& b) const { a *= alpha_; a +=  beta_*b; }
    void into_second(const NT& a, NT& b) const { b *=  beta_; b += alpha_*a; }
    void into_third (const NT& a, const NT& b, NT& c) const {
        c = a; c *= alpha_; c += beta_*b; // c might alias a but not b
    }
};

template <class NT_>
class Convex_combinator_approx_long_log {
public:
    typedef NT_ NT;

private:
    long alpha_num_, beta_num_, half_;
    int log_denom_;

public:
    Convex_combinator_approx_long_log(
            long alpha_num = 1, int log_denom = 1
    ) : alpha_num_(alpha_num),
        beta_num_((1L<<log_denom) - alpha_num),
        half_((log_denom > 0) ? (1L << (log_denom-1)) : 0),
        log_denom_(log_denom)
    {
        CGAL_precondition(log_denom_ >= 0);
    }
    void into_first(NT& a, const NT& b) const {
        a *= alpha_num_; a += beta_num_*b;
        a += half_; a >>= log_denom_;  // round to nearest
    }
    void into_second(const NT& a, NT& b) const {
        b *= beta_num_; b += alpha_num_*a;
        b += half_; b >>= log_denom_;  // round to nearest
    }
    void into_third(const NT& a, const NT& b, NT& c) const {
        c = a; c *= alpha_num_; c += beta_num_*b; // c might alias a but not b
        c += half_; c >>= log_denom_;  // round to nearest
    }
};

template <class NT_>
class Convex_combinator_approx_midpoint {
public:
    typedef NT_ NT;

public:
    void into_first(NT& a, const NT& b) const {
        a += b;
        a += 1; a >>= 1;  // round to nearest
    }
    void into_second(const NT& a, NT& b) const {
        b += a;
        b += 1; b >>= 1;  // round to nearest
    }
    void into_third(const NT& a, const NT& b, NT& c) const {
        c = a; c += b;  // c might alias a but not b
        c += 1; c >>= 1;  // round to nearest
    }
};

template <class ForwardIterator1, class ForwardIterator2,
    class ForwardIterator3, class Combinator
>
std::pair<ForwardIterator2, ForwardIterator3>
de_casteljau_generic(
        ForwardIterator1 first, ForwardIterator1 beyond,
        ForwardIterator2 left, ForwardIterator3 right,
        const Combinator& combine
) {
    // left may not alias any of the other two
    CGAL_assertion((void*)(&(*left)) != (void*)(&(*first)));
    CGAL_assertion((void*)(&(*left)) != (void*)(&(*right)));

    /* In the sequel, we think of de Casteljau's algorithm as
     * filling out a triangular array of numbers:
     *
     *     * * * .. *   row 0
     *      * * .. *    row 1
     *       .   .       .
     *        . .        .
     *         *
     *
     * The inputs [first, beyond) form row 0 (at the top).
     * Row 1 (below row 0) consists of combinations of any two
     * adjacent elements from row 0.
     * Inductively, row i+1 consists of combinations of any two
     * adjacent elements from row i, until we arrive at a row
     * of length 1.
     *
     * We output through iterator left the values appearing
     * on the triangle's left side.
     *
     * We use the container pointed to by iterator right as
     * storage for the rows (one after the other).
     * Since each row is one element shorter than the preceding
     * one, towards the end of this container the values appearing
     * on the triangle's right accumulate. When we terminate,
     * the container pointed to by iterator right contains
     * the triangle's right side, from bottom to top.
     */

    *left = *first; ++left; // output leftmost element of row 0

    // to compute row 1, we iterate over pairs (*iit1, *iit2)
    // of adjacent elements of row 0 and combine them into *rit1
    ForwardIterator1 iit1 = first, iit2 = first; ++iit2; // source
    ForwardIterator3 rit1 = right;                       // target
    for (;;) {
        combine.into_third(*iit1, *iit2, *rit1);
        ++iit1; ++iit2;
        if (iit2 == beyond) break;
        ++rit1;
    }
    ForwardIterator3 right_end = rit1;   // point to rightmost element of row 1
    *++rit1 = *iit1;                     // output rightmost element of row 0
    ForwardIterator3 right_beyond = ++rit1; // past-the-end
    *left = *right; ++left;              // output leftmost element of row 1

    // compute rows 2 and later in the same style
    // invariant: right_end is the rightmost element of the previous row
    ForwardIterator3 rit2;
    while (right != right_end) { // prev row was longer than 1 element
        rit1 = rit2 = right; ++rit2;
        combine.into_first(*rit1, *rit2);
        *left = *rit1; ++left;           // output leftmost element
        while (rit2 != right_end) {
            ++rit1; ++rit2;
            combine.into_first(*rit1, *rit2);
        } 
        right_end = rit1;
    }

    return std::make_pair(left, right_beyond);
} // de_casteljau_generic()



/*
 * Helper functors
 */

template <class CeilLog2Abs>
class Abs_le_pow2 {
public:
    typedef CeilLog2Abs Ceil_log2_abs;
    typedef bool result_type;
    typedef typename Ceil_log2_abs::argument_type first_argument_type;
    typedef typename Ceil_log2_abs::result_type   second_argument_type;
    result_type operator() (first_argument_type x, second_argument_type p) {
        return x == 0 || Ceil_log2_abs()(x) <= p;
    }
};

template <class Integer_, class AbsLePow2, class Sign_>
class Sign_eps_log2 {
private:
    long log_eps_;

public:
    typedef Integer_ Integer;
    typedef AbsLePow2 Abs_le_pow2;
    typedef Sign_ Sign;
    typedef internal::Three_valued_estimate result_type;
    typedef Integer argument_type;

    Sign_eps_log2(long log_eps = 0) : log_eps_(log_eps) { }
    long log_eps() const { return log_eps_; }
    void set_log_eps(long log_eps) { log_eps_ = log_eps; }

    result_type operator() (argument_type x) const {
        if (Abs_le_pow2()(x, log_eps_)) {
            return internal::UNCLEAR_SIGN;
        } else {
            return internal::Three_valued_estimate(Sign()(x));
        }
    }
}; // class Sign_eps_log2

} // namespace internal


namespace internal {
    template <class BitstreamDescartesRndlTreeTraits>
    class Bitstream_descartes_rndl_tree;

    template <class BitstreamDescartesRndlTreeTraits>
    struct Bitstream_descartes_rndl_node;

    template <class BitstreamDescartesRndlTreeTraits>
    class Bitstream_descartes_rndl_tree_rep;
} // namespace internal

/* The template argument supplied as BitstreamDescartesRndlTreeTraits
 * shall be a class containing the following types in its scope:
 *   Coefficient:           caller-supplied coefficient type
 *   Bound:              type for interval bound output (exact)
 *   Integer:               integer type for actual calculations, needs >>, <<
 *   Approximator:          functor to get Integer approx to x*2^p from coeff x
 *   Lower_bound_log2_abs:  functor for lower bound to log|x| for coeff x
 *   Bound_creator:      functor to create bound x*2^p from x and p
 *   Sign:                  functor to get sign of Integer x
 *   Ceil_log2_abs_Integer: functor to get smallest long >= log|x| for Integer
 *   Ceil_log2_abs_long:    functor to get smallest long >= log|x| for long
 */

/*
 * macros for common typedefs
 */

// bring types from traits into local scope
#define CGAL_SNAP_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_TYPEDEFS(TRAITS) \
  typedef typename TRAITS::Coefficient           Coefficient;           \
  typedef typename TRAITS::Bound              Bound;                    \
  typedef typename TRAITS::Integer               Integer;               \
  typedef typename TRAITS::Approximator          Approximator;          \
  typedef typename TRAITS::Lower_bound_log2_abs  Lower_bound_log2_abs;  \
  typedef typename TRAITS::Bound_creator      Bound_creator;            \
  typedef typename TRAITS::Sign                  Sign;                  \
  typedef typename TRAITS::Ceil_log2_abs_Integer Ceil_log2_abs_Integer; \
  typedef typename TRAITS::Ceil_log2_abs_long    Ceil_log2_abs_long     \

// end #define

// common typedefs for all Bitstream_descartes_rndl_* classes
#define CGAL_BITSTREAM_DESCARTES_RNDL_TREE_COMMON_TYPEDEFS              \
  typedef BitstreamDescartesRndlTreeTraits TRAITS;                      \
  typedef TRAITS Bitstream_descartes_rndl_tree_traits;                  \
  CGAL_SNAP_BITSTREAM_DESCARTES_RNDL_TREE_TRAITS_TYPEDEFS(TRAITS);      \
  typedef std::vector<Coefficient> Coefficient_vector;                  \
  typedef std::vector<Integer> Integer_vector;                          \
  typedef internal::Abs_le_pow2<Ceil_log2_abs_Integer> Abs_le_pow2;     \
  typedef internal::Sign_eps_log2<Integer, Abs_le_pow2, Sign>           \
  Sign_eps_log2                                                         \
  
// end #define

// typedefs for Bitstream_descartes_rndl_tree{,_rep}
#define CGAL_BITSTREAM_DESCARTES_RNDL_TREE_TYPEDEFS                   \
  CGAL_BITSTREAM_DESCARTES_RNDL_TREE_COMMON_TYPEDEFS;                 \
  typedef internal::Bitstream_descartes_rndl_node<TRAITS> Node;       \
  typedef std::list<Node> Node_list                                   \

// end #define


namespace internal {

/*
 * class Bitstream_descartes_rndl_node
 */

template <class BitstreamDescartesRndlTreeTraits>
struct Bitstream_descartes_rndl_node {
public:
    typedef Bitstream_descartes_rndl_node Self;
    CGAL_BITSTREAM_DESCARTES_RNDL_TREE_COMMON_TYPEDEFS;

    friend class internal::Bitstream_descartes_rndl_tree<TRAITS>;
    friend class internal::Bitstream_descartes_rndl_tree_rep<TRAITS>;

private:
    // "node data" (set individually in subdivision)
    Integer lower_num_, upper_num_; // TODO use lower_num_, width_num_ instead
    long log_bdry_den_;
    Integer_vector coeff_; // wrt [lower_, upper_], approximate
    int min_var_, max_var_;
    // "state data" (copied en bloc by .copy_state_from())
    long subdiv_tries_, subdiv_fails_;
    long recdepth_;
    long log_sep_, delta_log_sep_, log_eps_, log_C_eps_;

    Bitstream_descartes_rndl_node(int degree = -1,
            Integer lower_num = Integer(0), Integer upper_num = Integer(0),
            long log_bdry_den = 0, int min_var = -1, int max_var = -1
    ) : lower_num_(lower_num), upper_num_(upper_num),
          log_bdry_den_(log_bdry_den),
          coeff_(degree+1),
          min_var_(min_var), max_var_(max_var),
          subdiv_tries_(0), subdiv_fails_(0),
          recdepth_(-1),
          log_sep_(0), delta_log_sep_(0), log_eps_(0), log_C_eps_(0)
    { }

    void copy_state_from(const Self& n) {
        subdiv_tries_  = n.subdiv_tries_;
        subdiv_fails_  = n.subdiv_fails_;
        recdepth_      = n.recdepth_;
        log_sep_       = n.log_sep_;
        delta_log_sep_ = n.delta_log_sep_;
        log_eps_       = n.log_eps_;
        log_C_eps_     = n.log_C_eps_;
    }

    // const Self& operator= (const Self&); // assignment is forbidden
}; // struct Bitstream_descartes_rndl_node


/*
 * class Bitstream_descartes_rndl_tree_rep
 */

template <class BitstreamDescartesRndlTreeTraits>
class Bitstream_descartes_rndl_tree_rep {
public:
    typedef Bitstream_descartes_rndl_tree_rep Self;
    CGAL_BITSTREAM_DESCARTES_RNDL_TREE_TYPEDEFS;

    class Monomial_basis_tag { };

    friend class internal::Bitstream_descartes_rndl_tree<TRAITS>;

private:
    Coefficient_vector input_monomial_coeff_;
    int degree_;
    long ceil_log_degree_;
    long lbd_log_lcoeff_;
    Node_list node_list_;

    // temporary data fields for subdivision
    Integer_vector tmp1_coeff_, tmp2_coeff_;
    Integer splitpoint_num_;
    long log_splitpoint_den_;

    // function objects
    Approximator         approximator_;
    Lower_bound_log2_abs lower_bound_log2_abs_;

public:
    Bitstream_descartes_rndl_tree_rep() : degree_(-1) { }

    template <class InputIterator>
    Bitstream_descartes_rndl_tree_rep(
            Integer lower_num,  Integer upper_num, long log_bdry_den,
            InputIterator first, InputIterator beyond, Monomial_basis_tag,
            const TRAITS& traits
    ) : input_monomial_coeff_(first, beyond),
        splitpoint_num_(0), log_splitpoint_den_(0),
        approximator_(traits.approximator_object()),
        lower_bound_log2_abs_(traits.lower_bound_log2_abs_object())
    {
        degree_ = int(input_monomial_coeff_.size() - 1);
        CGAL_precondition(degree_ >= 0);
        ceil_log_degree_ = (degree_ > 0) ? Ceil_log2_abs_long()(degree_) : -1;
        lbd_log_lcoeff_
            = lower_bound_log2_abs_(input_monomial_coeff_[degree_]);
        node_list_.push_front(
                Node(degree_, lower_num, upper_num, log_bdry_den)
        );
        tmp1_coeff_.resize(degree_ + 1);
        tmp2_coeff_.resize(degree_ + 1);
    }
}; // class Bitstream_descartes_rndl_tree_rep

/*
 * class Bitstream_descartes_rndl_tree
 */

/*! \ingroup NiX_Bitstream_descartes_tree
    \brief Subdivision tree of the BitstreamDescartes method (rndl variant)

    Before you try to understand this class fully, you might want
    to have a look at the paper on the BitstreamDescartes method
    mentioned \link NiX_Bitstream_descartes here \endlink.
    The next paragraph gives a brief summary;
    the description of this class follows after it.

    <b>The BitstreamDescartes method</b>

    The BitstreamDescartes method searches the real roots of
    a polynomial in some initial interval by subdividing this
    interval recursively into open subintervals.
    Each subinterval is subjected to the Descartes Test,
    which gives an integer that is an upper bound on the
    number of real roots in the interval.
    For efficiency, the BitstreamDescartes method does not compute
    with the coefficients given by the caller, only approximations of them.
    Therefore, the result of the Descartes Test may only be known
    in the form of lower and upper bounds on the exact test,
    called min_var and max_var.
    However, the approximation quality of the input coefficients
    and the choice of the subdivision points are automatically
    controlled in a way that allows the following conclusions:
     - If min_var 0 for some interval, this interval
       does not contain any real root.
     - If max_var is 1 for some interval, this interval
       contains exactly one simple real root.
     - If all real roots of the input polynomial are simple,
       repeated subdivision will eventually produce intervals
       that all have min_var equal to max_var equal to 0 or 1.

    Hence we have an algorithm for isolating the real roots
    of square-free polynomials.  We think of it as constructing
    a binary tree:
    Each subinterval considered by the algorithm is a node of the tree.
    The children of a node are the two subintervals created by subdivision.
    The root of the tree is the initial interval.
    At each stage of the algorithm, the interesting nodes are
    the leaves of the tree.  Subdivision of a leaf turns it
    into an internal node whose children are leaves.

    <b>Description of class</b>

    This class lets you interactively explore the subdivision
    tree of the BitstreamDescartes method (or more precisely,
    its variant called "rndL" in the paper).
    An object \c T of this class is constructed from an
    initial interval and the polynomial.
    The polynomial is read from an iterator range,
    whose first element is the constant coefficient and whose
    last element is the leading coefficient (which has to be
    non-zero).

    After construction, \c T represents the list of leaves
    in the current subdivision tree.  (Initially, there is
    only one leaf: the node representing the initial interval.)
    You can iterate through this list in the style of an
    \c std::list; that is, a \c Node_iterator is a
    \c BidirectionalIterator.
    At any time, the nodes are sorted in the order of the intervals
    they stand for.

    You cannot dereference an iterator into anything meaningful;
    however, you can pass it to various member functions
    that tell you the relevant data about the interval
    represented by the node, e.g., \c lower() and \c upper()
    bound and \c min_var() and \c max_var().

    Most importantly, you can \c subdivide() a node.
    Conceptually, this replaces the interval by two subintervals.
    However, subintervals with \c min_var() equal to 0
    are immediately discarded; so \c subdivide() may
    actually replace one node by zero, one or two new nodes.
    In fact, if \c min_var() is 0 already for the initial interval,
    a newly constructed object \c T has an empty list of leaves.
    (In particular, this happens if the polynomial supplied in
    construction is constant.)

    Unlike STL containers, this class is implemented using
    \c CGAL::Handle.  That means, an object is just a ref-counted
    pointer to the actual representation.  Thus, copying an
    object is cheap.  However, all copies of an object alias
    each other.  If you modify one, this changes the state
    of all copies; including invalidation of iterators pointing
    to destroyed nodes.

    A \c Node_iterator remains valid until the node it points to
    is destroyed by \c subdivide() or \c erase().  Destruction
    of one node does not affect validity of iterators pointing
    to other nodes.

    <b>Example</b> (root isolation (square-free case))

    Once you have constructed \c T, implementing the BitstreamDescartes
    method by exploring \c T is a matter of a few lines:
    \code
    Node_iterator it = T.begin();
    Node_iterator chld_first, chld_beyond;
    while (it != T.end()) {
        if (T.max_var(it) == 1) {
            cout << "found [" << T.lower(it) << ", " << T.upper(it) << "]\n"; 
            ++it;
        } else {
            T.subdivide(it, chld_first, chld_beyond);
            it = chld_first;
        }
    }
    \endcode

    <b>Supplying a traits class</b>

    This class is actually a class template.
    To use it, you need to instanciate it with a traits class
    that defines the following three types and the various
    functors on them listed below.
     - \c Coefficient: The type of coefficients supplied
       during construction. Must be \c Assignable .
     - \c Integer: A type for infinite-precision integer arithmetic
       (such as \c leda::integer or \c CORE::BigInt ).
       All internal computations are done using this type.
       Must be a model of \c Ring and additionally provide
       operators \c >> and \c << with the usual semantics.
     - \c Bound:  \c lower() and \c upper() return
       interval boundaries in this type.  Must be \c Assignable.
       The canonical choice is \c NiX::Exact_float_number<Integer>.
       If you never instanciate \c lower() and \c upper()
       (maybe use \c boundaries() instead), you might be lucky
       and get away with typedef'ing this to \c void.

    The traits class must also contain the following functors
    and member functions for their construction:
     - \c Approximator: A \c BinaryFunction with signature
       <tt>Integer y = Approximator()(Coefficient x, long p)</tt>
       that computes an \c Integer approximation
       to 2<sup><i>p</i></sup>&nbsp;<tt>*</tt>&nbsp;<i>x</i> satisfying
       |<i>y</i> - 2<sup><i>p</i></sup>&nbsp;<tt>*</tt>&nbsp;<i>x</i>| <= 1.
     - \c approximator_object(): A \c const member function
       taking no arguments and returning a function object
       of class \c Approximator.  This function is called once
       at construction of \c T to initialize one \c Approximator
       that is used for all subsequent coefficient approximations.
       It is only applied to arguments \c x that have had one
       of the coefficients assigned to them that were supplied
       during construction of \c T. Hence it can
       keep state and maybe cache some knowledge about coefficients.
     - \c Lower_bound_log2_abs: A \c UnaryFunction with signature
       <tt>long l = Lower_bound_log2_abs()(Coefficient x)</tt>.
       The result \c l must be a lower bound to log<sub>2</sub>(|<i>x</i>|).
       If \c Coefficient posesses \c NiX::NT_traits::Floor_log2_abs,
       you can simply use that.
     - \c lower_bound_log2_abs_object(): A \c const member function
       taking no arguments and returning a function object
       of class \c Lower_bound_log2_abs.  This function is called once
       at construction of \c T to get a \c Lower_bound_log2_abs
       on the polynomial's leading coefficient.
     - \c Bound_creator: A functor with signature
       <tt>Bound b = Bound_creator()(Integer x, long p)</tt>
       to construct \c b with value
       <i>x</i>&nbsp;<tt>*</tt>&nbsp;2<sup><i>p</i></sup>.
       If \c Bound has a matching constructor
       (as \c NiX::Exact_float_number<Integer> does), you can simply
       <tt>typedef CGAL::Creator_2 <Integer, long, Bound>
       Bound_creator;</tt>.
     - \c Sign: A functor working identically to
       \c NiX::NT_traits::Sign for \c NT equal to \c Integer.
       (You can just typedef to that.)
     - \c Ceil_log2_abs_Integer: A functor working identically to
       \c NiX::NT_traits::Ceil_log2_abs for \c NT equal to \c Integer.
       (You can just typedef to that.)
     - \c Ceil_log2_abs_long: A functor working identically to
       \c NiX::NT_traits::Ceil_log2_abs for \c NT equal to \c long.
       (You can just typedef to that.)

    In brief, the core requirement is that you can approximate
    \c Coefficient to any arbitrarily small absolute error
    2<sup><i>-p</i></sup> (for integral <i>p</i>)
    and deliver that approximation scaled with 2<sup><i>p</i></sup>
    as an \c Integer.
    For the leading coefficient, you also need to be able to locate
    the leading 1-bit in its binary expansion.
    The functors dealing with \c Coefficient are accessed
    through \c _object() member functions so that the user of
    this class can supply them with an internal state, because
    \c Coefficient might hide some non-trivial
    approximation or evaluation process.
    For the functors dealing with \c Integer, this does not
    seem necessary.

    <b>Example</b> (traits class)

    You can use the BitstreamDescartes method for polynomials with
    integer coefficients.  If the coefficients are very long, this
    saves time over the exact Descartes method, because they are
    only needed in truncated form. A suitable traits class looks
    like this:
    \code
    template <class Integer_>
    class Bitstream_descartes_rndl_tree_traits_from_Integer_coeff {
    public:
        typedef Integer_ Coefficient;
        typedef Integer_ Integer;
        typedef NiX::Exact_float_number<Integer> Bound;

        class Approximator {
        public:
            Integer operator() (Coefficient x, long p) {
                if (p >= 0) return x << p; else return x >> -p;
            }
        };
        Approximator approximator_object() const { return Approximator(); }

        typedef typename NiX::NT_traits<Coefficient>::Floor_log2_abs Lower_bound_log2_abs;
        Lower_bound_log2_abs lower_bound_log2_abs_object() const { return Lower_bound_log2_abs(); }

        typedef CGAL::Creator_2<Integer, long, Bound> Bound_creator;
        typedef typename NiX::NT_traits<Integer>::Sign Sign;
        typedef typename NiX::NT_traits<Integer>::Ceil_log2_abs Ceil_log2_abs_Integer;
        typedef typename NiX::NT_traits<long>::Ceil_log2_abs Ceil_log2_abs_long;
    };
    \endcode

    <b>Technical remarks</b>

    This class implements essentially the "rndL" variant
    of the BitstreamDescartes method. The following remarks apply:

    It is an invariant that for all nodes, the first and
    last Bernstein coefficient are larger in magnitude than C eps.
    Consequently, Lemma 5 of the paper allow us to conclude
    right away from \c min_var() being 0 what would happen
    after one further subdivision (namely: there are no roots).
    We don't have to do this extra subdivision, we know right away.
    Unfortunately, the analogous argument for Lemma 6 doesn't work,
    because we still would have to verify that the value at the
    bisection point is large.

    Subdivision points with value larger than C eps are found
    by trying randomly and checking.  This randomization means
    the same polynomial and same initial interval may give rise
    to different intervals each time this class is used.
    As indicated in the paper, we favour subdivision ratios
    with a small denominator. Hence we first try denominator
    2 (subdivision at midpoint), then denominator 16, and
    only then the "proper" denominator prescribed by theory.
    Failures are only counted for the "proper" tries.

    Unlike the algorithm in the paper, we do not have one global
    estimate for the separation of roots, and we do not restart
    globally if that estimate turns out wrong.  Instead, each node
    maintains an estimate.  Upon subdivision, its children
    inherit it, including the counts of tried and failed
    subdivisions. If fails/tries >= 1/2 and tries >= 2,
    the estimate of separation (and all other parameters
    coming out of it) are updated only for this one node.
    The node's interval does not change; subdivision
    resumes from this interval with an improved estimate of
    separation.  All other nodes are unaffected.
    Global restart is not an option for this class, because the user has
    already observed the subintervals found up to this point,
    so we cannot simply switch to other intervals.

    This implementation relies on the assumption that the
    <i>logarithms</i> of certain relevant quantities, in particular
    the degree times the logarithm of the estimated separation,
    are small enough to be representable in a <tt>long int</tt>.
    Also, the degree of the input polynomial and related quantities
    are assumed to fit into an \c int.
 */
template <class BitstreamDescartesRndlTreeTraits>
class Bitstream_descartes_rndl_tree
// TODO: Replaced CGAL::Handle by following CGAL::Handle_with_policy, is this correct?
    : public ::CGAL::Handle_with_policy< Bitstream_descartes_rndl_tree_rep<
        BitstreamDescartesRndlTreeTraits
    >, ::CGAL::Handle_policy_no_union >
{
public:
    typedef Bitstream_descartes_rndl_tree Self;
    CGAL_BITSTREAM_DESCARTES_RNDL_TREE_TYPEDEFS;
    typedef Bitstream_descartes_rndl_tree_rep<TRAITS> Rep;
    typedef ::CGAL::Handle_with_policy< Rep, ::CGAL::Handle_policy_no_union > Base;

    //! node iterator.
    typedef typename Node_list::iterator       Node_iterator;
    //! node iterator (for STL compatibility only).
    typedef typename Node_list::iterator       iterator;
    //! const node iterator.
    typedef typename Node_list::const_iterator Node_const_iterator;
    //! const node iterator (for STL compatibility only).
    typedef typename Node_list::const_iterator const_iterator;

    //! tag type to distinguish a certain constructor.
    typedef typename Rep::Monomial_basis_tag Monomial_basis_tag;

public:
    //! default constructor (makes <tt>degree() == -1</tt>)
    Bitstream_descartes_rndl_tree() : Base(Rep()) { }

    //! copy constructor
    Bitstream_descartes_rndl_tree(const Self& p)
        : Base(static_cast<const Base&>(p))
    { }

    //! Internal function called by constructor. Avoids code duplication
    void init_tree() {
        Node_iterator n = this->ptr()->node_list_.begin();
        if (this->ptr()->degree_ > 0) {
            initial_guess_sep(n);
            while (!reinit_from_sep(n)) next_guess_sep(n);
            if (n->min_var_ == 0) this->ptr()->node_list_.erase(n);
        } else {
            this->ptr()->node_list_.erase(n);
        }
    }
        

    /*! \brief construct from initial interval and coefficients
     *
     *  The initial interval is
     *  [\c lower_num, \c upper_num] / 2^(\c log_bdry_den ).
     *
     *  The iterator range [\c first, \c beyond ) gives the
     *  coefficients of 1, <i>x</i>, <i>x</i><sup>2</sup>, ...
     *  The leading coefficient (last in sequence) must be non-zero.
     *
     *  The \c Monomial_basis_tag is required for the benefit of
     *  future extensions to coefficients w.r.t. other bases.
     */
    template <class InputIterator>
    Bitstream_descartes_rndl_tree(
            Integer lower_num,  Integer upper_num, long log_bdry_den,
            InputIterator first, InputIterator beyond, Monomial_basis_tag tag,
            const BitstreamDescartesRndlTreeTraits& traits
                                        = BitstreamDescartesRndlTreeTraits()
    ) : Base(Rep(lower_num, upper_num, log_bdry_den,
                    first, beyond, tag, traits))
    {
        CGAL_precondition(lower_num < upper_num);
        init_tree();
        
    }

    /*! 
     * This is needed for compatibility with other tree implementations
     * The initial interval is
     *  [-1, 1] / 2^(\c -log_bdry_den ).
     * Be aware that log_bdry_den is negated here!
     */
    template <class InputIterator>
    Bitstream_descartes_rndl_tree(
            long log_bdry_den,
            InputIterator first, InputIterator beyond, Monomial_basis_tag tag,
            const BitstreamDescartesRndlTreeTraits& traits
                                        = BitstreamDescartesRndlTreeTraits()
    )
        : Base(Rep(Integer(-1), Integer(1), -log_bdry_den, 
                   first, beyond, tag, traits))
    {
        init_tree();
    }

    //! return degree of polynomial
    int degree() const { return this->ptr()->degree_; }

    //! iterator to first node
    Node_iterator begin() {
        return this->ptr()->node_list_.begin();
    }
    //! iterator beyond last node
    Node_iterator end() {
        return this->ptr()->node_list_.end();
    }
    //! const iterator to first node
    Node_const_iterator begin() const {
        return this->ptr()->node_list_.begin();
    }
    //! const iterator beyond last node
    Node_const_iterator end() const {
        return this->ptr()->node_list_.end();
    }

    //! get lower bound of interval at node \c n.
    Bound lower(Node_iterator n) const {
        return Bound_creator()(n->lower_num_, -n->log_bdry_den_);
    }
    //! get lower bound of interval at node \c n.
    Bound lower(Node_const_iterator n) const {
        return Bound_creator()(n->lower_num_, -n->log_bdry_den_);
    }
    //! get upper bound of interval at node \c n.
    Bound upper(Node_iterator n) const {
        return Bound_creator()(n->upper_num_, -n->log_bdry_den_);
    }
    //! get upper bound of interval at node \c n.
    Bound upper(Node_const_iterator n) const {
        return Bound_creator()(n->upper_num_, -n->log_bdry_den_);
    }

    //! get boundaries: interval at node \c n is
    //! [\c lower_num, \c upper_num] / 2^(\c log_bdry_den ).
    void boundaries(Node_iterator n,
            Integer& lower_num, Integer& upper_num, long& log_bdry_den
    ) {
        lower_num    = n->lower_num_;
        upper_num    = n->upper_num_;
        log_bdry_den = n->log_bdry_den_;
    }

    //! get minimum number of sign variations in Descartes Test
    //! for approximate polynomial at node \c n
    int min_var(Node_const_iterator n) const { return n->min_var_; }
    //! get maximum number of sign variations in Descartes Test
    //! for approximate polynomial at node \c n
    int max_var(Node_const_iterator n) const { return n->max_var_; }

    /*! \brief subdivide interval at node \c n.
     *
     *  The node representing interval \c n is replaced in the list
     *  of nodes by 0, 1, or 2 nodes that represent those among the
     *  two subintervals that have \c min_var() greater than 0.
     *  The number of new nodes is returned as result.
     *  The subrange of nodes consisting of the new nodes is
     *  returned in the arguments \c first and \c beyond
     *
     *  Subdividing a node invalidates all iterators to it;
     *  both for the object where \c subdivide() was called
     *  and all copies of it (since they share the same
     *  representation and state).
     *
     *  The parameter \c n is passed by value.  Hence you can
     *  implement depth-first search for isolating intervals like this:
     *  \code
     *    Node_iterator dummy, curr = tree.begin();
     *    while (curr != tree.end()) {
     *        if (tree.max_var(curr) == 1) ++curr;
     *        else tree.subdivide(curr, curr, dummy);
     *    }
     *  \endcode
     */
    int subdivide(
            Node_iterator n, Node_iterator& first, Node_iterator& beyond
    );

    /*! \brief erase node \c n.
     *
     *  Erasing a node invalidates all iterators to it;
     *  both for the object where \c subdivide() was called
     *  and all copies of it (since they share the same
     *  representation and state).
     */
    void erase(Node_iterator n) {
        this->ptr()->node_list_.erase(n);
    }

    /*! \brief Replace traits class
     */
    void set_traits(TRAITS& traits) {

      this->ptr()->approximator_ 
        = traits.approximator_object();
      this->ptr()->lower_bound_log2_abs_ 
        = traits.lower_bound_log2_abs_object();

    }

    /*! \brief Returns a copy of this with its own representation
     */
    Self make_unique() const {
      Self tmp = *this;
      tmp.copy_on_write();
      return tmp;
    }


protected:
    int subdivide_at_midpoint(
            Node_iterator n, Node_iterator& first, Node_iterator& beyond
    );

    int subdivide_at(
            Node_iterator n, Node_iterator& first, Node_iterator& beyond,
            long alpha_num, int log_alpha_den
    );

private:
    int replace_by_tmp(
            Node_iterator n, Node_iterator& first, Node_iterator& beyond
    );

    void initial_guess_sep(Node_iterator n) {
        Ceil_log2_abs_Integer log;
        long log_I = log(n->upper_num_-n->lower_num_) - n->log_bdry_den_;
        n->delta_log_sep_ = -5;
        n->log_sep_ = log_I + n->delta_log_sep_;
    }

    void next_guess_sep(Node_iterator n) {
        n->delta_log_sep_ *= 2;
        CGAL_warning_msg(-n->delta_log_sep_ < 1L<<24, "delta_log_sep >= 1L<<24");
        n->log_sep_ += n->delta_log_sep_;
    }

    bool reinit_from_sep(Node_iterator n);

}; // class Bitstream_descartes_rndl_tree


/*
 * Non-inline member functions of class Bitstream_descartes_rndl_tree
 */

template <class BitstreamDescartesRndlTreeTraits>
bool
Bitstream_descartes_rndl_tree<BitstreamDescartesRndlTreeTraits>
::reinit_from_sep(Node_iterator n) {
    n->subdiv_tries_ = n->subdiv_fails_ = 0;

    Ceil_log2_abs_Integer log;
    /* We want to set recdepth to
     *   floor( (log(|I|/sep) / log(4/3)) + 5/2 )
     * or something slightly larger.
     * Using the continued fractions expansion [2,2,2,3],
     * we find an upper bound of 41/17 = 82/34 = 2.41176...
     * for the exact multiplier     1/log(4/3) = 2.40942...
     * which is off by less than 0.1%, namely    0.00234...
     * Noting 5/2 = 85/34, we hence set recdepth to
     *   floor( (log|I| - log(sep))*82 + 85) / 34 )
     */
    n->recdepth_ = (
            (log(n->upper_num_ - n->lower_num_) - n->log_bdry_den_ // log|I|
            - n->log_sep_
        )*82 + 85) / 34;
    if (n->recdepth_ < 6) n->recdepth_ = 6; // TODO find rationale

    n->log_eps_ = this->ptr()->ceil_log_degree_
        + Ceil_log2_abs_long()(n->recdepth_);
    n->log_C_eps_ = n->log_eps_ + 4*degree(); // C = 16^n
    long target_log_lcf = (4 - n->log_sep_)*this->ptr()->degree_
        + n->log_eps_ + 2*this->ptr()->ceil_log_degree_ + 10;
    long log_lcf_scale = target_log_lcf - this->ptr()->lbd_log_lcoeff_
        - (log(caching_factorial<Integer>(this->ptr()->degree_)) - 1);

    polynomial_power_to_bernstein_approx(
            this->ptr()->input_monomial_coeff_.begin(),
            this->ptr()->input_monomial_coeff_.end(),
            n->coeff_.begin(),
            n->lower_num_, n->upper_num_, n->log_bdry_den_,
            log_lcf_scale - (n->log_eps_ - 1),
            this->ptr()->approximator_,
            Ceil_log2_abs_Integer(), Ceil_log2_abs_long()
    );
    for (int i = 0; i <= degree(); ++i) {
        n->coeff_[i] <<= (n->log_eps_ - 1); // TODO avoid preceding rshift
    }

    if (Abs_le_pow2()(n->coeff_[ degree()], n->log_C_eps_)
            ||  Abs_le_pow2()(n->coeff_[0], n->log_C_eps_)
    ) {
        return false;
    } else {
        var_eps(n->coeff_.begin(), n->coeff_.end(),
                n->min_var_, n->max_var_, Sign_eps_log2(n->log_eps_)
        );
        return true;
    }
} // Bitstream_descartes_rndl_tree::reinit_from_sep()

template <class BitstreamDescartesRndlTreeTraits>
int
Bitstream_descartes_rndl_tree<BitstreamDescartesRndlTreeTraits>
::subdivide_at_midpoint(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond
) {
    de_casteljau_generic(n->coeff_.begin(), n->coeff_.end(),
            this->ptr()->tmp1_coeff_.begin(), this->ptr()->tmp2_coeff_.begin(),
            Convex_combinator_approx_midpoint<Integer>()
    );
    this->ptr()->splitpoint_num_     = n->lower_num_ + n->upper_num_;
    this->ptr()->log_splitpoint_den_ = n->log_bdry_den_ + 1;

    if (Abs_le_pow2()(this->ptr()->tmp2_coeff_[0], n->log_C_eps_)) {
        return -1;
    } else {
        return replace_by_tmp(n, first, beyond);
    }
} // Bitstream_descartes_rndl_tree::subdivide_at()

template <class BitstreamDescartesRndlTreeTraits>
int
Bitstream_descartes_rndl_tree<BitstreamDescartesRndlTreeTraits>
::subdivide_at(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond,
        long alpha_num, int log_alpha_den
) {
    de_casteljau_generic(n->coeff_.begin(), n->coeff_.end(),
        this->ptr()->tmp1_coeff_.begin(), this->ptr()->tmp2_coeff_.begin(),
        Convex_combinator_approx_long_log<Integer>(alpha_num, log_alpha_den)
    );
    this->ptr()->splitpoint_num_ =
        alpha_num * n->lower_num_
            + ((1L << log_alpha_den) - alpha_num) * n->upper_num_;
    this->ptr()->log_splitpoint_den_ = n->log_bdry_den_ + log_alpha_den;

    if (Abs_le_pow2()(this->ptr()->tmp2_coeff_[0], n->log_C_eps_)) {
        return -1;
    } else {
        return replace_by_tmp(n, first, beyond);
    }
} // Bitstream_descartes_rndl_tree::subdivide_at()


template <class BitstreamDescartesRndlTreeTraits>
int
Bitstream_descartes_rndl_tree<BitstreamDescartesRndlTreeTraits>
::replace_by_tmp(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond
) {
    --(n->recdepth_);

    long delta_log_bdry_den =
        this->ptr()->log_splitpoint_den_ - n->log_bdry_den_;
    CGAL_assertion(delta_log_bdry_den >= 0);

    int l_min_var, l_max_var, r_min_var, r_max_var;
    var_eps(this->ptr()->tmp1_coeff_.begin(),
            this->ptr()->tmp1_coeff_.end(),
            l_min_var, l_max_var, Sign_eps_log2(n->log_eps_)
    );
    var_eps(this->ptr()->tmp2_coeff_.begin(),
            this->ptr()->tmp2_coeff_.end(),
            r_min_var, r_max_var, Sign_eps_log2(n->log_eps_)
    );
    CGAL_assertion(l_min_var >= 0 && l_max_var >= 0);
    CGAL_assertion(r_min_var >= 0 && r_max_var >= 0);

    beyond = first = n;
    ++beyond;

    if (l_min_var > 0) {
        int children = 1;
        if (r_min_var > 0) {
            // create new node for right child
            Node_iterator r = 
                this->ptr()->node_list_.insert(beyond, Node(degree(),
                            this->ptr()->splitpoint_num_,        // lower
                            n->upper_num_ << delta_log_bdry_den, // upper
                            this->ptr()->log_splitpoint_den_,
                            r_min_var, r_max_var
                ));
            r->coeff_.swap(this->ptr()->tmp2_coeff_);
            r->copy_state_from(*n);
            ++children;
        }
        // put left child into n
        n->lower_num_  <<= delta_log_bdry_den;
        n->upper_num_    = this->ptr()->splitpoint_num_;
        n->log_bdry_den_ = this->ptr()->log_splitpoint_den_;
        n->min_var_      = l_min_var;
        n->max_var_      = l_max_var;
        n->coeff_.swap(this->ptr()->tmp1_coeff_);
        return children;
    } else if (r_min_var > 0) {
        // put right child into n
        n->lower_num_    = this->ptr()->splitpoint_num_;
        n->upper_num_  <<= delta_log_bdry_den;
        n->log_bdry_den_ = this->ptr()->log_splitpoint_den_;
        n->min_var_      = r_min_var;
        n->max_var_      = r_max_var;
        n->coeff_.swap(this->ptr()->tmp2_coeff_);
        return 1;
    } else /* l_min_var == 0 && r_min_var == 0 */ {
        // delete n
        first = beyond;
        this->ptr()->node_list_.erase(n);
        return 0;
    }
} // Bitstream_descartes_rndl_tree::replace_by_tmp()

template <class BitstreamDescartesRndlTreeTraits>
int
Bitstream_descartes_rndl_tree<BitstreamDescartesRndlTreeTraits>
::subdivide(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond
) {
    long alpha_num;
    int  log_alpha_den;
    long alpha_den_4;
    int  ret;

    for (;;) {
        if (n->recdepth_ > 0) { // TODO decouple recdepth from guess_sep
            // first try heuristic alpha = 1/2 (failures don't count)
            ++(n->subdiv_tries_);
            ret = subdivide_at_midpoint(n, first, beyond);
            if (ret >= 0) { return ret; } else { --(n->subdiv_tries_); }

            // next try heuristic alpha with small denom (failures don't count)
            log_alpha_den = 4;
            alpha_den_4 = 1L << (log_alpha_den - 2);
            alpha_num = CGAL::get_default_random().get_int(  // TODO .get_long
                    alpha_den_4, 3*alpha_den_4 + 1
            );
            ++(n->subdiv_tries_);
            ret = subdivide_at(n, first, beyond, alpha_num, log_alpha_den);
            if (ret >= 0) { return ret; } else { --(n->subdiv_tries_); }

            // now try alpha properly randomized, counting failure rate
            log_alpha_den = 5 + this->ptr()->ceil_log_degree_;
            alpha_den_4 = 1L << (log_alpha_den - 2);
            do {
                alpha_num = CGAL::get_default_random().get_int(  // TODO .get_long
                        alpha_den_4, 3*alpha_den_4 + 1
                );
                ++(n->subdiv_tries_);
                ret = subdivide_at(n, first, beyond, alpha_num, log_alpha_den);
                if (ret >= 0) {
                    return ret;
                } else {
                    ++(n->subdiv_fails_);
                }
            } while (n->subdiv_fails_ < 2  // Arno says: 2 (not 6) is enough
                    || 2 * n->subdiv_fails_ < n->subdiv_tries_);
        } // if (n->recdepth_ > 0)

        // if failure rate too high or recdepth exceeded,
        // decrease guess of sep and restart
        next_guess_sep(n);
        bool reinit_success = reinit_from_sep(n);
        CGAL_assertion(reinit_success); (void)reinit_success;
    } // for (;;)
} // Bitstream_descartes_rndl_tree::subdivide()

} // namespace internal


/*
 * FUJIWARA ROOT BOUND
 */

namespace internal {

struct Fujiwara_root_bound_queue_entry {
    typedef Fujiwara_root_bound_queue_entry Self;

    int n_minus_i;
    long ub_log2_qi;
    bool is_certainly_zero;
    bool is_tight;

    bool operator < (Self& rhs) {
        if (is_certainly_zero) return !rhs.is_certainly_zero;
        if (rhs.is_certainly_zero) return false;
        return ub_log2_qi * rhs.n_minus_i < rhs.ub_log2_qi * n_minus_i;
    }
};

class Fujiwara_root_bound_queue_entry_ptr_less {
public:
    typedef bool result_type;
    typedef Fujiwara_root_bound_queue_entry* first_argument_type;
    typedef Fujiwara_root_bound_queue_entry* second_argument_type;
    result_type operator() (first_argument_type a, second_argument_type b) {
        return *a < *b;
    }
};

template <class CeilLog2Abs>
class Upper_bound_log2_abs_approximator_from_ceil_log2_abs {
public:
    typedef Upper_bound_log2_abs_approximator_from_ceil_log2_abs Self;
    typedef CeilLog2Abs Ceil_log2_abs;
    typedef typename Ceil_log2_abs::argument_type NT;
    bool initial_upper_bound(NT x, long& ub_log2, bool& is_certainly_zero) {
        is_certainly_zero = (x == NT(0));
        if (!is_certainly_zero) ub_log2 = Ceil_log2_abs()(x);
        return true; // reported bound is tight
    }
    bool improve_upper_bound(NT, long&, bool&) { return true; }
};

/*! \ingroup NiX_Bitstream_descartes_tree
    \brief Fujiwara root bound (as logarithm wrt base 2)

    All complex roots of a polynomial
    <i>A</i>(<i>X</i>)&nbsp;=&nbsp;<i>a<sub>n</sub>X<sup>n</sup></i>&nbsp;+&nbsp;...&nbsp;+&nbsp;<i>a</i><sub>0</sub>
    are bounded in magnitude by
    F(<i>A</i>)&nbsp;=&nbsp;2&nbsp;max<sub><i>i</i>&lt;<i>n</i></sub>&nbsp;<i>q<sub>i</sub></i><sup>1/(<i>n</i>-<i>i</i>)</sup>
    where <i>q<sub>i</sub></i>&nbsp;=&nbsp;|<i>a<sub>i</sub></i>/<i>a<sub>n</sub></i>|
    for 0 < <i>i</i> < <i>n</i>
    and <i>q</i><sub>0</sub>
    =&nbsp;|<i>a</i><sub>0</sub>/(2<i>a<sub>n</sub></i>)|.
    This bound goes back to M. Fujiwara
    [<i>Tohoku Math. J.</i> <b>10</b> (1916) 167-171],
    cited here after P. Batra's PhD thesis [TU Hamburg-Harburg, Germany, 1999].

    This function computes an integer upper bound for
    log<sub>2</sub>(F(<i>A</i>)) from an iterator range [first, beyond)
    over values of some type \c Coefficient
    where <tt>*(first+i)</tt> is <i>a<sub>i</sub></i>.
    This function does not operate on type \c Coefficient
    except through the two function objects passed to it:

    <tt>LowerBoundLog2Abs lblog2</tt> has to be a function object
    with a function call operator taking an argument
    <tt>*(beyond-1)</tt> of type \c Coefficient
    and returning a \c long which is a lower bound for
    log<sub>2</sub>(|<i>a<sub>n</sub></i>|).
    This functor is called once and should return a bound as
    good as possible.

    <tt>UpperBoundLog2AbsApproximator ublog2apx</tt> has to be an
    object with two member functions
    \c initial_upper_bound() and \c improve_upper_bound().
    Both of them take three arguments:
    <tt>(Coefficient ai, long& ub_log2, bool& is_certainly_zero)</tt>
    and return <tt>bool</tt>. \c ai is one of the coefficients.
    First, the member function \c initial_upper_bound() is invoked
    on \c ai and uninitialized arguments \c ub_log2 and \c is_certainly_zero.
    Then, \c improve_upper_bound() is invoked repeatedly on \c ai and
    arguments \c ub_log2 and \c is_certainly_zero as set by the
    previous call.  This sequence of calls has to create a
    sequence of estimates of upper bounds for
    log<sub>2</sub>(|<tt>ai</tt>|).  As long as
    a later call will return a better approximation,
    the function returns \c false. If the current bound is to be
    regarded as best possible, the function returns \c true.
    If, while improving the bound, it is discovered that
    <tt>ai</tt> is zero, in which case the correct estimate would be
    "minus infinity", \c is_certainly_zero is set to \c true;
    otherwise it is always \c false.
    If \c is_certainly_zero is set to \c true,
    the value of \c ub_log2 is ignored.

    Internally, this function works as follows:
    It improves upper bounds of log<sub>2</sub>(|<i>a<sub>i</sub></i>|)
    for all <i>i</i> &lt; <i>n</i> until it has found one that
    is designated as best possible (i.e., has returned \c true)
    and that realizes the maximum in the definition
    of the Fujiwara bound. Unlike \c LowerBoundLog2Abs, which is
    called only once, namely for the leading coefficient, and should
    provide the best possible bound at once,
    \c UpperBoundLog2AbsApproximator is called repeatedly,
    and each call should perform only a limited amount of work
    (e.g., one round of interval refinement for its argument)
    so that this function can distribute the approximation work
    evenly over all coefficients until the maximum is found.

    The leading coefficient <tt>*(beyond-1)</tt> has to be non-zero.
    In the special case that all other coefficients are zero,
    the result 0 (standing for 2<sup>0</sup>=1) is returned.
    Obviously, this cannot happen for a square-free polynomial
    of degree larger than 1.

    <b>Warning:</b> This bound is tight for certain polynomials
    (e.g., (<i>X</i>-2)(<i>X<sup>n</sup></i>-1)/(<i>X</i>-1)).
    The BitstreamDescartes method, however, needs an initial interval
    whose boundaries are far away from the zeroes. Hence you should
    add 1 to the result of this function before using it there.
    Also bear in mind that this function returns a logarithm of
    the bound, whereas certain other functions may expect a
    logarithm of a denominator, necessitating a negation.

    For coefficients that are known explicitly (e.g., big integers)
    and possess <tt>NiX::NT_traits::{Floor,Ceil}_log2_abs</tt> functors,
    the lower and upper bound function objects for this function
    should be constructed from them.  This is done automatically
    by an overloaded version of this functions that takes only the
    first two arguments.
 */
template<class RandomAccessIterator,
    class LowerBoundLog2Abs,
    class UpperBoundLog2AbsApproximator
>
long Fujiwara_root_bound_log(
    RandomAccessIterator first, RandomAccessIterator beyond,
    LowerBoundLog2Abs lblog2, UpperBoundLog2AbsApproximator ublog2apx
) {
    int n = beyond - first - 1; // degree
    if (n < 1) return 0;
    long lblog2_lcoeff = lblog2(*(beyond - 1));

    Fujiwara_root_bound_queue_entry_ptr_less less;
    typedef Fujiwara_root_bound_queue_entry QE;
    std::vector<QE>  entries(n); // entries are never copied
    std::vector<QE*> heap(n);    // heap is built from pointers to them
    for (int i = 0; i < n; ++i) {
        QE& entry = entries[i];
        entry.n_minus_i = n - i;
        entry.is_tight = ublog2apx.initial_upper_bound(
                *(first + i), entry.ub_log2_qi, entry.is_certainly_zero
        );
        CGAL_assertion(entry.is_tight || !entry.is_certainly_zero);
        if (!entry.is_certainly_zero) entry.ub_log2_qi -= lblog2_lcoeff;
        heap[i] = &(entry);
    }
    entries[0].ub_log2_qi -= 1;

    std::make_heap(heap.begin(), heap.end(), less);
    while (!heap[0]->is_tight) {
        std::pop_heap(heap.begin(), heap.end(), less);
        QE& popped = **(heap.end() - 1);
        int i = n - popped.n_minus_i;
        CGAL_assertion(i >= 0 && i < n);
        CGAL_assertion(&popped == &(entries[i]));
        popped.is_tight = ublog2apx.improve_upper_bound(
                *(first + i),
                popped.ub_log2_qi, popped.is_certainly_zero
        );
        if (!popped.is_certainly_zero) {
            popped.ub_log2_qi -= lblog2_lcoeff;
            if (i == 0) popped.ub_log2_qi -= 1;
        }
        std::push_heap(heap.begin(), heap.end(), less);
    }
    QE& maxi = *(heap[0]);
    if (maxi.is_certainly_zero) return 0;
    long max_log2_qi_div_n_minus_i = maxi.ub_log2_qi / maxi.n_minus_i;
    if (maxi.ub_log2_qi % maxi.n_minus_i > 0) ++max_log2_qi_div_n_minus_i;
    return 1 + max_log2_qi_div_n_minus_i; // = log ( 2 * max_i q_i^{1/(n-i)} )
}

template<class RandomAccessIterator>
inline
long Fujiwara_root_bound_log(
    RandomAccessIterator first, RandomAccessIterator beyond
) {
    typedef typename RandomAccessIterator::value_type NT;
    typedef typename internal::Real_embeddable_extension<NT>::Floor_log2_abs Lbd;
    typedef Upper_bound_log2_abs_approximator_from_ceil_log2_abs<
                typename internal::Real_embeddable_extension<NT>::Ceil_log2_abs
            > Ubd;
    return Fujiwara_root_bound_log(first, beyond, Lbd(), Ubd());
}


} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_RNDL_TREE_H

// EOF
