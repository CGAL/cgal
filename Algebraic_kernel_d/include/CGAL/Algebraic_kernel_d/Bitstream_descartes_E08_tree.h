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


#ifndef CGAL_BITSTREAM_DESCARTES_E08_TREE_H
#define CGAL_BITSTREAM_DESCARTES_E08_TREE_H

#include <vector>
#include <list>
#include <utility>
#include <iterator>
#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree.h> // TODO remove
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Random.h>
#include <CGAL/tss.h>

#include <boost/optional.hpp>

/*
 *  AUXILIARY CLASSES AND FUNCTIONS
 */

namespace CGAL {

namespace internal {

template <class Integer>
Integer caching_binomial(int n, int k) {
    CGAL_precondition(n >= 0);
    if (k < 0 || k > n) return Integer(0);

    // Pascal's triangle; augment if necessary
    // TODO flat array with manual index computation should be slightly faster
    typedef std::vector< Integer > Row;
    typedef std::vector< Row > Triangle;
    // MSVC uses "pascal" as a keyword or defines it as a macro!
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Triangle, my_pascal);

    int old_size = int(my_pascal.size());
    if (n >= old_size) {
        my_pascal.resize(n+1);
        if (old_size == 0) {
            my_pascal[0].push_back(Integer(1));
            old_size = 1;
        }
        for (int i = old_size; i <= n; ++i) {
            Row& prev = my_pascal[i-1];
            Row& curr = my_pascal[i];
            curr.reserve(i+1);
            curr.push_back(Integer(1));
            for (int j = 1; j < i; ++j) {
                curr.push_back(prev[j-1] + prev[j]);
            }
            curr.push_back(Integer(1));
        }
    }
    return my_pascal[n][k];
}


template <class Integer_>
class Power_to_Bernstein_pm1_nofrac_matrix {
public:
    typedef Integer_ Integer;

private:
    int dim_;
    std::vector<Integer> m_;
    int ij_to_idx(int i, int j) const { return i*dim_ + j; }
    void init_m();

public:
    Power_to_Bernstein_pm1_nofrac_matrix(int degree = -1)
        : dim_(degree + 1), m_(dim_ * dim_)
    {
        CGAL_assertion(degree >= -1);
        init_m();
    }

    void set_degree(int degree) {
        CGAL_assertion(degree >= -1);
        dim_ = degree + 1;
        m_.resize(dim_ * dim_);
        init_m();
    }

    Integer operator()(int i, int j) const {
        CGAL_assertion(0 <= i && i < dim_);
        CGAL_assertion(0 <= j && j < dim_);
        return m_[ij_to_idx(i,j)];
    }

    int degree() const { return dim_ - 1; }
}; // class Power_to_Bernstein_pm1_nofrac_matrix

template <class Integer_>
void Power_to_Bernstein_pm1_nofrac_matrix<Integer_>::init_m() {
    // TODO: this implements the definition, but [E08] describes a faster way
    int degree = dim_ - 1;
    for (int i = 0; i < dim_; ++i) {
        for (int j = 0; j < dim_; ++j) {
            Integer sum(0);
            int nu_lo = (std::max)(0, i+j-degree);
            int nu_hi = (std::min)(i, j);
            for (int nu = nu_lo; nu <= nu_hi; ++nu) {
                Integer term = caching_binomial<Integer>(j, nu)
                             * caching_binomial<Integer>(degree-j, i-nu);
                if ((j - nu) & 1) { sum -= term; }
                else              { sum += term; }
            }
            m_[ij_to_idx(i,j)] =
                sum 
              * CGAL::internal::caching_factorial<Integer>(i)
              * CGAL::internal::caching_factorial<Integer>(degree - i);
        }
    }
}

template <class Traits_>
class Bitstream_bernstein_from_power {
public:
    typedef Traits_  Traits;
    typedef typename Traits::Coefficient           Coefficient;
    typedef typename Traits::Integer               Integer;
    typedef typename Traits::Approximator          Approximator;
    typedef typename Traits::Lower_bound_log2_abs  Lower_bound_log2_abs;
    typedef typename Traits::Ceil_log2_abs_Integer Ceil_log2_abs_Integer;
    typedef typename Traits::Ceil_log2_abs_long    Ceil_log2_abs_long;

    typedef std::vector<Coefficient> Coefficient_vector;
    typedef std::vector<Integer>     Integer_vector;
    typedef Power_to_Bernstein_pm1_nofrac_matrix<Integer> Mn;
    typedef std::map<int, Mn> Mn_cache;

private:
    Coefficient_vector input_power_coeff_;
    int degree_;        // $ =   n $ in [E08]
    long log_radius_;   // $ = r+1 $ in [E08]
    Approximator approximator_;
    Lower_bound_log2_abs lower_bound_log2_abs_;
    long ceil_log_degfac_, floor_log_degfac_, ceil_log_degp1_;
    long lbd_log_lcoeff_;
    long log_lcfscale_; // $ =   l $ in [E08]

    Integer_vector bernstein_coeff_; // empty if uninitialized
    long bernstein_coeff_prec_; // $ = p+1 $ in [E08]


    static const Mn& get_mn(int n) {
      CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Mn_cache, mn_cache_);
        Mn& mn = mn_cache_[n];
        if (mn.degree() == -1) mn.set_degree(n);
        return mn;
    }

public:
    Bitstream_bernstein_from_power() : degree_(-1) { }

    template <class InputIterator>
    Bitstream_bernstein_from_power(
            InputIterator first, InputIterator beyond,
            long log_radius,
            const Traits_& traits = Traits_()
    );

    template <class OutputIterator>
    OutputIterator operator()(OutputIterator oi, long prec);

    int degree() const { return degree_; }

    void set_traits(Traits& traits) {
        approximator_ = traits.approximator_object();
        lower_bound_log2_abs_ = traits.lower_bound_log2_abs_object();
    }
}; // class Bitstream_bernstein_from_power


// non-member functions
template <class Traits_>
template <class InputIterator>
Bitstream_bernstein_from_power<Traits_>::Bitstream_bernstein_from_power(
        InputIterator first, InputIterator beyond,
        long log_radius,
        const Traits_& traits // = Traits_()
) : input_power_coeff_(first, beyond),
    degree_(int(input_power_coeff_.size() - 1)),
    log_radius_(log_radius),
    approximator_(traits.approximator_object()),
    lower_bound_log2_abs_(traits.lower_bound_log2_abs_object())
{ 
    CGAL_assertion(degree_ >= 0);
    lbd_log_lcoeff_ = lower_bound_log2_abs_(input_power_coeff_[degree_]);
    if (degree_ == 2) {
        ceil_log_degfac_ = floor_log_degfac_ = 1;
        ceil_log_degp1_ = 2;
    } else {
        ceil_log_degfac_
            = Ceil_log2_abs_Integer()
                (CGAL::internal::
                 caching_factorial<Integer>(degree_));
        floor_log_degfac_
            = ceil_log_degfac_ - 1;
        ceil_log_degp1_
            = Ceil_log2_abs_long()(degree_ + 1);
    }
    log_lcfscale_
        = lbd_log_lcoeff_ + floor_log_degfac_ + degree_ * (log_radius+1);
}

template <class Traits_>
template <class OutputIterator>
OutputIterator
Bitstream_bernstein_from_power<Traits_>::operator()(
        OutputIterator oi,
        long prec // $ = p+1$ in [E08]
) {
    CGAL_precondition(degree_ >= 0);

    // try to produce output from existing coefficient approximations
    if (!bernstein_coeff_.empty()) {
        long delta_prec = bernstein_coeff_prec_ - prec;
        if (delta_prec == 0) {
            for (int i = 0; i <= degree_; ++i) {
                *oi = bernstein_coeff_[i];
                ++oi;
            }
            return oi;
        } else if (delta_prec > 0) {
            Integer half = Integer(1) << (delta_prec-1);
            for (int i = 0; i <= degree_; ++i) {
                *oi = (bernstein_coeff_[i] + half) >> delta_prec;
                ++oi;
            }
            return oi;
        }
        // else delta_prec < 0: fall through to precision increase
    } else {
        bernstein_coeff_.resize(degree_ + 1);
    }

    // compute coefficients at new precision
    Integer_vector c(degree_ + 1);
    long q = prec + ceil_log_degfac_ + ceil_log_degp1_ + 1;
    for (int j = 0; j <= degree_; ++j) {
        c[j] = approximator_(
                input_power_coeff_[j],
                j*log_radius_ - log_lcfscale_ + q
        );
    }
    const Mn& mn = get_mn(degree_);
    long shift = q - prec;
    Integer half = Integer(1) << (shift-1);
    for (int i = 0; i <= degree_; ++i) {
        Integer bi(0);
        for (int j = 0; j <= degree_; ++j) bi += mn(i,j) * c[j];
        bernstein_coeff_[i] = (bi + half) >> shift;
    }
    bernstein_coeff_prec_ = prec;

    // output new approximations
    for (int i = 0; i <= degree_; ++i) {
        *oi = bernstein_coeff_[i];
        ++oi;
    }
    return oi;
}

} // namespace internal


/*
 * The generic de Casteljau method
 */

template <class Integer_>
class Convex_combinator_approx_Integer_log {
public:
    typedef Integer_ Integer;

private:
    Integer alpha_num_, beta_num_, half_;
    int log_denom_;

public:
    Convex_combinator_approx_Integer_log(
            Integer alpha_num = Integer(1), int log_denom = 1
    ) : alpha_num_(alpha_num),
        beta_num_((Integer(1) << log_denom) - alpha_num),
        half_((log_denom > 0) ? (Integer(1) << log_denom-1) : 0),
        log_denom_(log_denom)
    {
        CGAL_precondition(log_denom_ >= 0);
    }
    void into_first(Integer& a, const Integer& b) const {
        a *= alpha_num_; a += beta_num_*b;
        a += half_; a >>= log_denom_;  // round to nearest
    }
    void into_second(const Integer& a, Integer& b) const {
        b *= beta_num_; b += alpha_num_*a;
        b += half_; b >>= log_denom_;  // round to nearest
    }
    void into_third(const Integer& a, const Integer& b, Integer& c) const {
        c = a; c *= alpha_num_; c += beta_num_*b; // c might alias a but not b
        c += half_; c >>= log_denom_;  // round to nearest
    }
};

template <class NT_>
class Convex_combinator_approx_fraction {
public:
    typedef NT_ NT;
private:
    NT alpha_num_, beta_num_, denom_, half_;
public:
    Convex_combinator_approx_fraction(NT alpha_num, NT denom)
        : alpha_num_(alpha_num), beta_num_(denom - alpha_num),
          denom_(denom), half_(denom >> 1)
    { }
    void into_first(NT& a, const NT& b) const {
        a *= alpha_num_; a += beta_num_*b;
        a += half_; a /= denom_;  // round to nearest
    }
    void into_second(const NT& a, NT& b) const {
        b *= beta_num_; b += alpha_num_*a;
        b += half_; b /= denom_;  // round to nearest
    }
    void into_third(const NT& a, const NT& b, NT& c) const {
        c = a; c *= alpha_num_; c += beta_num_*b; // c might alias a but not b
        c += half_; c /= denom_;  // round to nearest
    }
}; // class Convex_combinator_approx_fraction

namespace internal {

/*
 * THE ACTUAL TREE CLASSES
 */
template <class BitstreamDescartesE08TreeTraits>
class Bitstream_descartes_E08_tree;

template <class BitstreamDescartesE08TreeTraits>
        struct Bitstream_descartes_E08_node;

        template <class BitstreamDescartesE08TreeTraits>
        class Bitstream_descartes_E08_tree_rep;
} // namespace internal

} //namespace CGAL

/* The template argument supplied as BitstreamDescartesE08TreeTraits
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
#define CGAL_SNAP_BITSTREAM_DESCARTES_E08_TREE_TRAITS_TYPEDEFS(TRAITS)  \
  typedef typename TRAITS::Coefficient           Coefficient;           \
  typedef typename TRAITS::Bound                 Bound;                 \
  typedef typename TRAITS::Integer               Integer;               \
  typedef typename TRAITS::Approximator          Approximator;          \
  typedef typename TRAITS::Lower_bound_log2_abs  Lower_bound_log2_abs;  \
  typedef typename TRAITS::Bound_creator      Bound_creator;            \
  typedef typename TRAITS::Sign                  Sign;                  \
  typedef typename TRAITS::Ceil_log2_abs_Integer Ceil_log2_abs_Integer; \
  typedef typename TRAITS::Ceil_log2_abs_long    Ceil_log2_abs_long     \

// end #define

// common typedefs for all Bitstream_descartes_E08_* classes
#define CGAL_BITSTREAM_DESCARTES_E08_TREE_COMMON_TYPEDEFS               \
  typedef BitstreamDescartesE08TreeTraits TRAITS;                       \
  typedef TRAITS Bitstream_descartes_E08_tree_traits;                   \
  CGAL_SNAP_BITSTREAM_DESCARTES_E08_TREE_TRAITS_TYPEDEFS(TRAITS);       \
  typedef CGAL::internal::Bitstream_bernstein_from_power<TRAITS> B_from_p; \
  typedef std::vector<Integer> Integer_vector;                          \
  typedef CGAL::internal::Abs_le_pow2<Ceil_log2_abs_Integer>            \
  Abs_le_pow2;                                                          \
  typedef CGAL::internal::Sign_eps_log2                                 \
  <Integer, Abs_le_pow2, Sign>                                          \
  Sign_eps_log2                                                         \
  
// end #define

// typedefs for Bitstream_descartes_E08_tree{,_rep}
#define CGAL_BITSTREAM_DESCARTES_E08_TREE_TYPEDEFS                     \
  CGAL_BITSTREAM_DESCARTES_E08_TREE_COMMON_TYPEDEFS;                   \
  typedef CGAL::internal::Bitstream_descartes_E08_node<TRAITS> Node;   \
  typedef std::list<Node> Node_list                                    \
  
// end #define


namespace CGAL {

namespace internal {

/*
 * class Bitstream_descartes_E08_node
 */

template <class BitstreamDescartesE08TreeTraits>
struct Bitstream_descartes_E08_node {
public:
  typedef Bitstream_descartes_E08_node Self;
  CGAL_BITSTREAM_DESCARTES_E08_TREE_COMMON_TYPEDEFS;
  
  friend class CGAL::internal::Bitstream_descartes_E08_tree<TRAITS>;
  friend class CGAL::internal::Bitstream_descartes_E08_tree_rep<TRAITS>;

private:
    // "node data" (set individually in subdivision)
    Integer lower_num_, upper_num_; // TODO use lower_num_, width_num_ instead
    long log_bdry_den_;
    Integer_vector coeff_; // wrt [lower_, upper_], approximate
    int min_var_, max_var_;
    bool coeff_update_delayed_;
    // "state data" (copied en bloc by .copy_state_from())
    long subdepth_bound_, subdepth_current_;
    long log_eps_;   // $q - p$
    long log_C_eps_; // $q - p + 4n$

    Bitstream_descartes_E08_node(int degree = -1,
            Integer lower_num = Integer(0), Integer upper_num = Integer(0),
            long log_bdry_den = 0, int min_var = -1, int max_var = -1
    ) : lower_num_(lower_num), upper_num_(upper_num),
          log_bdry_den_(log_bdry_den),
          coeff_(degree+1),
          min_var_(min_var), max_var_(max_var),
          coeff_update_delayed_(false),
          subdepth_bound_(0), subdepth_current_(0),
           log_eps_(0), log_C_eps_(0)
    { }

    void copy_state_from(const Self& n) {
        subdepth_bound_   = n.subdepth_bound_;
        subdepth_current_ = n.subdepth_current_;
        log_eps_          = n.log_eps_;
        log_C_eps_        = n.log_C_eps_;
    }

    // const Self& operator= (const Self&); // assignment is forbidden
}; // struct Bitstream_descartes_E08_node


/*
 * class Bitstream_descartes_E08_tree_rep
 */

template <class BitstreamDescartesE08TreeTraits>
class Bitstream_descartes_E08_tree_rep {
public:
    typedef Bitstream_descartes_E08_tree_rep Self;
    CGAL_BITSTREAM_DESCARTES_E08_TREE_TYPEDEFS;

    class Monomial_basis_tag { };

    friend class CGAL::internal::Bitstream_descartes_E08_tree<TRAITS>;

private:
    B_from_p b_from_p_;
    long log_radius_;
    int degree_;
    long ceil_log_degree_;
    Node_list node_list_;

    long payload_prec_;
    int subdiv_tries_, subdiv_fails_;
    int bisect_tries_, bisect_fails_;

    // temporary data fields for subdivision
    Integer_vector tmp1_coeff_, tmp2_coeff_;
    Integer splitpoint_num_;
    long log_splitpoint_den_;

public:
    Bitstream_descartes_E08_tree_rep() : degree_(-1) { }

    template <class InputIterator>
    Bitstream_descartes_E08_tree_rep(
            long log_radius,
            InputIterator first, InputIterator beyond, Monomial_basis_tag,
            const TRAITS& traits
    ) : b_from_p_(first, beyond, log_radius, traits),
        log_radius_(log_radius),
        subdiv_tries_(0), subdiv_fails_(0),
        bisect_tries_(0), bisect_fails_(0),
        splitpoint_num_(0), log_splitpoint_den_(0)
    {
        degree_ = b_from_p_.degree();
        CGAL_precondition(degree_ >= 0);
        ceil_log_degree_ = (degree_ > 0) ? Ceil_log2_abs_long()(degree_) : -1;
        node_list_.push_front(
                Node(degree_, Integer(-1), Integer(1), -log_radius)
        );
        payload_prec_ = 6 * degree_ + 20;
        tmp1_coeff_.resize(degree_ + 1);
        tmp2_coeff_.resize(degree_ + 1);
    }
}; // class Bitstream_descartes_E08_tree_rep

/*
 * class Bitstream_descartes_E08_tree
 */

/*! \ingroup CGAL_Bitstream_descartes_tree
 *  \brief Subdivision tree of the BitstreamDescartes method (E08 variant)
 */
template <class BitstreamDescartesE08TreeTraits>
class Bitstream_descartes_E08_tree
    : public
        ::CGAL::Handle_with_policy<
                 internal::Bitstream_descartes_E08_tree_rep<
                         BitstreamDescartesE08TreeTraits
                 >
         >
{
public:
    typedef Bitstream_descartes_E08_tree Self;
    CGAL_BITSTREAM_DESCARTES_E08_TREE_TYPEDEFS;
    typedef internal::Bitstream_descartes_E08_tree_rep<TRAITS> Rep;
    typedef ::CGAL::Handle_with_policy<Rep> Base;

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
    Bitstream_descartes_E08_tree() : Base(Rep()) { }

    //! copy constructor
    Bitstream_descartes_E08_tree(const Self& p)
        : Base(static_cast<const Base&>(p))
    { }

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
    Bitstream_descartes_E08_tree(
            long log_radius,
            InputIterator first, InputIterator beyond, Monomial_basis_tag tag,
            const BitstreamDescartesE08TreeTraits& traits
                                        = BitstreamDescartesE08TreeTraits()
    ) : Base(Rep(log_radius, first, beyond, tag, traits))
    {
        Node_iterator n = this->ptr()->node_list_.begin();
        if (this->ptr()->degree_ > 0) {
            bool init_ok = reinit_from_prec(n);
            CGAL_assertion(init_ok); (void)init_ok;
            if (n->min_var_ == 0) this->ptr()->node_list_.erase(n);
        } else {
            this->ptr()->node_list_.erase(n);
        }
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
        CGAL_assertion(is_iterator_valid(n));
        return Bound_creator()(n->lower_num_, -n->log_bdry_den_);
    }
    //! get lower bound of interval at node \c n.
    Bound lower(Node_const_iterator n) const {
        CGAL_assertion(is_iterator_valid(n));
        return Bound_creator()(n->lower_num_, -n->log_bdry_den_);
    }
    //! get upper bound of interval at node \c n.
    Bound upper(Node_iterator n) const {
        CGAL_assertion(is_iterator_valid(n));
        return Bound_creator()(n->upper_num_, -n->log_bdry_den_);
    }
    //! get upper bound of interval at node \c n.
    Bound upper(Node_const_iterator n) const {
        CGAL_assertion(is_iterator_valid(n));
        return Bound_creator()(n->upper_num_, -n->log_bdry_den_);
    }

    //! get boundaries: interval at node \c n is
    //! [\c lower_num, \c upper_num] / 2^(\c log_bdry_den ).
    void boundaries(Node_iterator n,
            Integer& lower_num, Integer& upper_num, long& log_bdry_den
    ) {
        CGAL_assertion(is_iterator_valid(n));
        lower_num    = n->lower_num_;
        upper_num    = n->upper_num_;
        log_bdry_den = n->log_bdry_den_;
    }

    //! get minimum number of sign variations in Descartes Test
    //! for approximate polynomial at node \c n
    int min_var(Node_const_iterator n) const {
        CGAL_assertion(is_iterator_valid(n));
        return n->min_var_;
    }
    //! get maximum number of sign variations in Descartes Test
    //! for approximate polynomial at node \c n
    int max_var(Node_const_iterator n) const {
        CGAL_assertion(is_iterator_valid(n));
        return n->max_var_;
    }

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
        CGAL_assertion(is_iterator_valid(n));
        this->ptr()->node_list_.erase(n);
    }

    /*! \brief Replace traits class
     */
    void set_traits(TRAITS& traits) {
        this->ptr()->b_from_p_.set_traits(traits);
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

    bool is_iterator_valid(Node_const_iterator n) const {
        Node_const_iterator it  = this->ptr()->node_list_.begin();
        Node_const_iterator end = this->ptr()->node_list_.end();
        while (it != end) {
            if (it == n) return true;
            ++it;
        }
        return false;
    }

private:
    int replace_by_tmp(
            Node_iterator n, Node_iterator& first, Node_iterator& beyond
    );

    static const long log_subdepth_bound_init_ = 6;

    bool reinit_from_prec(Node_iterator n);
    void global_prec_increase(Node_iterator n);

}; // class Bitstream_descartes_E08_tree

template <class BitstreamDescartesE08TreeTraits>
const long
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::log_subdepth_bound_init_;



/*
 * Non-inline member functions of class Bitstream_descartes_E08_tree
 */

template <class BitstreamDescartesE08TreeTraits>
bool
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::reinit_from_prec(Node_iterator n) {

    n->subdepth_bound_   = 1L << log_subdepth_bound_init_;
    n->subdepth_current_ = 0;
    n->log_eps_ = this->ptr()->ceil_log_degree_
                + log_subdepth_bound_init_
                + 1;
    n->log_C_eps_ = n->log_eps_ + 4*this->degree();

    this->ptr()->b_from_p_(
            this->ptr()->tmp1_coeff_.begin(),
            this->ptr()->payload_prec_ + 1
    );
    for (int i = 0; i <= degree(); ++i) { // TODO avoid preceding rshift
        this->ptr()->tmp1_coeff_[i] <<= (n->log_eps_ - 1);
    }

    Integer alpha_num =
        (Integer(1) << (n->log_bdry_den_ + this->ptr()->log_radius_))
        -  n->lower_num_;
    long log_alpha_den = n->log_bdry_den_ + this->ptr()->log_radius_ + 1;
    if (alpha_num != Integer(1) << log_alpha_den) {
        de_casteljau_generic(
            this->ptr()->tmp1_coeff_.begin(), this->ptr()->tmp1_coeff_.end(),
            this->ptr()->tmp2_coeff_.begin(), this->ptr()->tmp1_coeff_.begin(),
            Convex_combinator_approx_Integer_log<Integer>(
                    alpha_num, log_alpha_den
            )
        );
        ++(n->subdepth_current_);
    }
    alpha_num =
        (Integer(1) << (n->log_bdry_den_ + this->ptr()->log_radius_))
        -  n->upper_num_;
    if (alpha_num != Integer(0)) {
        Integer alpha_den =
            (Integer(1) << (n->log_bdry_den_ + this->ptr()->log_radius_))
            - n->lower_num_;
        de_casteljau_generic(
            this->ptr()->tmp1_coeff_.begin(), this->ptr()->tmp1_coeff_.end(),
            n->coeff_.begin(), this->ptr()->tmp1_coeff_.begin(),
            Convex_combinator_approx_fraction<Integer>(alpha_num, alpha_den)
        );
        ++(n->subdepth_current_);
    } else {
        n->coeff_.swap(this->ptr()->tmp1_coeff_);
    }

    if (Abs_le_pow2()(n->coeff_[degree()], n->log_C_eps_)
            || Abs_le_pow2()(n->coeff_[0], n->log_C_eps_)
    ) {
        return false;
    } else {
        internal::var_eps(n->coeff_.begin(), n->coeff_.end(),
                n->min_var_, n->max_var_, Sign_eps_log2(n->log_eps_)
        );
        return true;
    }
} // Bitstream_descartes_E08_tree::reinit_from_sep()

template <class BitstreamDescartesE08TreeTraits>
int
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::subdivide_at_midpoint(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond
) {
    de_casteljau_generic(n->coeff_.begin(), n->coeff_.end(),
            this->ptr()->tmp1_coeff_.begin(), this->ptr()->tmp2_coeff_.begin(),
            CGAL::internal::Convex_combinator_approx_midpoint<Integer>()
    );
    this->ptr()->splitpoint_num_     = n->lower_num_ + n->upper_num_;
    this->ptr()->log_splitpoint_den_ = n->log_bdry_den_ + 1;

    if (Abs_le_pow2()(this->ptr()->tmp2_coeff_[0], n->log_C_eps_)) {
        return -1;
    } else {
        return replace_by_tmp(n, first, beyond);
    }
} // Bitstream_descartes_E08_tree::subdivide_at()

template <class BitstreamDescartesE08TreeTraits>
int
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::subdivide_at(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond,
        long alpha_num, int log_alpha_den
) {
    de_casteljau_generic(n->coeff_.begin(), n->coeff_.end(),
        this->ptr()->tmp1_coeff_.begin(), this->ptr()->tmp2_coeff_.begin(),
        CGAL::internal::Convex_combinator_approx_long_log<Integer>
            (alpha_num, log_alpha_den)
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
} // Bitstream_descartes_E08_tree::subdivide_at()


template <class BitstreamDescartesE08TreeTraits>
int
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::replace_by_tmp(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond
) {

    ++(n->subdepth_current_);

    long delta_log_bdry_den =
        this->ptr()->log_splitpoint_den_ - n->log_bdry_den_;
    CGAL_assertion(delta_log_bdry_den >= 0);

    int l_min_var, l_max_var, r_min_var, r_max_var;
    internal::var_eps(this->ptr()->tmp1_coeff_.begin(),
                                     this->ptr()->tmp1_coeff_.end(),
                                     l_min_var, l_max_var, 
                                     Sign_eps_log2(n->log_eps_)
    );
    internal::var_eps(this->ptr()->tmp2_coeff_.begin(),
            this->ptr()->tmp2_coeff_.end(),
            r_min_var, r_max_var, Sign_eps_log2(n->log_eps_)
    );
    CGAL_assertion(0 <= l_min_var && l_min_var <= l_max_var);
    CGAL_assertion(0 <= r_min_var && r_min_var <= r_max_var);

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
} // Bitstream_descartes_E08_tree::replace_by_tmp()

template <class BitstreamDescartesE08TreeTraits>
int
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::subdivide(
        Node_iterator n, Node_iterator& first, Node_iterator& beyond
) {
    long alpha_num;
    int  log_alpha_den;
    long alpha_den_4;
    int  ret;

    CGAL_assertion(is_iterator_valid(n));

    // check for delayed update
    if (n->coeff_update_delayed_) {
        reinit_from_prec(n);
        n->coeff_update_delayed_ = false;
    }

    // apply Zeno trap
    if (n->subdepth_current_ == n->subdepth_bound_) {
        for (int i = 0; i <= degree(); ++i)  n->coeff_[i] <<= 2;
        // n->working_prec_ += 2;
        n->log_eps_ += 2;
        n->log_C_eps_ += 2;
        n->subdepth_bound_ <<= 1;
        n->subdepth_current_ = 0;
    }

    for (;;) {
        if (true) { // used to be recdepth > 0
            // first try heuristic alpha = 1/2 (failures don't count)
            ++(this->ptr()->bisect_tries_);
            ret = subdivide_at_midpoint(n, first, beyond);
            if (ret >= 0) { return ret; } else { ++(this->ptr()->bisect_fails_); }

            // next try heuristic alpha with small denom (failures don't count)
            log_alpha_den = 4;
            alpha_den_4 = 1L << (log_alpha_den - 2);
            alpha_num = CGAL::get_default_random().get_int(  // TODO .get_long
                    alpha_den_4, 3*alpha_den_4 + 1
            );
            ++(this->ptr()->subdiv_tries_);
            ret = subdivide_at(n, first, beyond, alpha_num, log_alpha_den);
            if (ret >= 0) { return ret; } else { --(this->ptr()->subdiv_tries_); }

            // now try alpha properly randomized, counting failure rate
            log_alpha_den = 4 + this->ptr()->ceil_log_degree_;
            alpha_den_4 = 1L << (log_alpha_den - 2);
            do {
                alpha_num = CGAL::get_default_random().get_int(  // TODO .get_long
                        alpha_den_4, 3*alpha_den_4 + 1
                );
                ++(this->ptr()->subdiv_tries_);
                ret = subdivide_at(n, first, beyond, alpha_num, log_alpha_den);
                if (ret >= 0) {
                    return ret;
                } else {
                    ++(this->ptr()->subdiv_fails_);
                }
            } while (
(this->ptr()->subdiv_fails_ < 2 || 2 * this->ptr()->subdiv_fails_ < this->ptr()->subdiv_tries_)
&& (this->ptr()->bisect_fails_ < degree())
            );
        } // if (true) // used to be recdepth > 0

        // if failure rate too high, decrease guess of prec and restart
        global_prec_increase(n);
        // TODO what if now  n->max_var_ == 1  ??
    } // for (;;)
} // Bitstream_descartes_E08_tree::subdivide()

template <class BitstreamDescartesE08TreeTraits>
void
Bitstream_descartes_E08_tree<BitstreamDescartesE08TreeTraits>
::global_prec_increase(Node_iterator n)
{
    this->ptr()->payload_prec_ *= 2;

    this->ptr()->subdiv_tries_ = this->ptr()->subdiv_fails_ = 0;
    this->ptr()->bisect_tries_ = this->ptr()->bisect_fails_ = 0;

    Node_iterator it  = this->ptr()->node_list_.begin();
    Node_iterator end = this->ptr()->node_list_.end();
    while (it != end) {
        if (it != n && it->min_var_ == it->max_var_) {
            it->coeff_update_delayed_ = true;
        } else {
            bool reinit_ok = reinit_from_prec(it);
            CGAL_assertion(reinit_ok); (void)reinit_ok;
            it->coeff_update_delayed_ = false;
        }
        ++it;
    }
} // Bitstream_descartes_E08_tree::global_prec_increase()

} // namespace internal

} //namespace CGAL

#endif // CGAL_BITSTREAM_DESCARTES_E08_TREE_H

// EOF
