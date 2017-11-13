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
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_ACK_CURVE_PAIR_ANALYSIS_H
#define CGAL_ACK_CURVE_PAIR_ANALYSIS_H 1

#include <vector>
#include <algorithm>

#include <boost/optional.hpp>

#include <CGAL/Handle_with_policy.h>
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_d/Shear_controller.h>
#include <CGAL/Algebraic_kernel_d/Shear_transformation.h>
#include <CGAL/Algebraic_kernel_d/enums.h>
#include <CGAL/Algebraic_kernel_d/exceptions.h>
#include <CGAL/Algebraic_kernel_d/Status_line_CPA_1.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4290)
#endif

namespace CGAL {

namespace internal {
    
template<class AlgebraicReal_1>
class Distinct_compare {

public:

    typedef AlgebraicReal_1 Algebraic_real_1;

    typedef ::CGAL::Comparison_result result_type;
    typedef Algebraic_real_1                  first_argument_type;
    typedef Algebraic_real_1                  second_argument_type;

    ::CGAL::Comparison_result operator() 
        (Algebraic_real_1 a,Algebraic_real_1 b) {
        return a.compare_distinct(b);
    }

};

}// namespace internal

//////////////////////////////////////////////////////////////////////////////
// Curve_pair_2

// Forwards
template < typename AlgebraicKernelWithAnalysis_2 >
class Curve_pair_analysis_2;

template<typename AlgebraicKernelWithAnalysis_2>
std::ostream& operator<< 
    (std::ostream&,const Curve_pair_analysis_2
                           <AlgebraicKernelWithAnalysis_2>&);

namespace internal {

// Internally used enums and structs

enum Slice_type {
    FIRST_CURVE = 0,
    SECOND_CURVE = 1,
    INTERSECTION = 2,
    CANDIDATE = 3
};


/*!
 * An x-event of the curve pair is either a root of a dicriminant of a single
 * curve, or a root of the resultant of both curves, or both.
 * The \c Event_indices vector stores a triple <tt>(fg,ffy,ggy)</tt> denoting
 * that some event is the <tt>fg</tt> root of <tt>res(f,g,y)</tt>, 
 * the <tt>ffy</tt>th root of <tt>disc(f,y)</tt> and
 * the <tt>ggy</tt>th root of <tt>disc(g,y)</tt>.
 */
template<typename size_type>
struct Event_indices {

    size_type fg;
    size_type ffy;
    size_type ggy;
    Event_indices(size_type fg,size_type ffy, size_type ggy) 
    : fg(fg), ffy(ffy), ggy(ggy) {}
};

// Representation class for curve pairs
template < class AlgebraicKernelWithAnalysis_2 >
class Curve_pair_analysis_2_rep {
public:

    //! \name public typedefs
    //! @{
    typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

    typedef Curve_pair_analysis_2_rep<Algebraic_kernel_with_analysis_2> Self;

    typedef Curve_pair_analysis_2<Algebraic_kernel_with_analysis_2> Handle;

    typedef typename Algebraic_kernel_with_analysis_2::Curve_analysis_2 
        Curve_analysis_2;

    typedef typename Curve_analysis_2::size_type size_type;

    typedef typename Curve_analysis_2::Polynomial_2 Polynomial_2;
        
    typedef typename Curve_analysis_2::Algebraic_real_1 Algebraic_real_1;
        
    typedef typename Polynomial_2::NT Polynomial_1;

    typedef typename Curve_analysis_2::Bound Bound;

    typedef CGAL::internal::Status_line_CPA_1<Handle> Status_line_CPA_1;

    typedef std::pair<Slice_type,size_type> Slice_element;

    typedef std::vector<Slice_element> Slice_info;

    typedef boost::optional<Slice_info> Lazy_slice_info;

    typedef boost::optional<Bound> Lazy_bound;

    typedef CGAL::internal::Event_indices<size_type> Event_indices;

    struct Intersection_info {
        typename Curve_analysis_2::Status_line_1 ev;
        size_type index;
        size_type mult;
    };

    typedef std::vector<std::vector<Intersection_info> > 
        Intersection_info_container;
        
    typedef boost::optional<Intersection_info_container> 
        Lazy_intersection_info_container;

    // For lazy evaluation of Status_line_CPA_1s.
    typedef boost::optional<Status_line_CPA_1> Lazy_status_line_CPA_1;

    //! @}

    //! \name Constructors
    //! @{
    
    // DefaultConstructible
    Curve_pair_analysis_2_rep() :
        c1_(), c2_() {
    }

    Curve_pair_analysis_2_rep(Algebraic_kernel_with_analysis_2 *kernel,
                              Curve_analysis_2 c1, Curve_analysis_2 c2,
                              CGAL::Degeneracy_strategy strategy) :
        _m_kernel(kernel),
        c1_(c1), c2_(c2), f(c1.polynomial_2()), g(c2.polynomial_2()),
        degeneracy_strategy(strategy) {
    }

    //! @}

private:

    //! \name members
    //! @{

    Algebraic_kernel_with_analysis_2* _m_kernel;
    
    Curve_analysis_2 c1_;
    Curve_analysis_2 c2_;

    Polynomial_2 f;
    Polynomial_2 g;

    
    mutable boost::optional<std::vector<Polynomial_2> > subresultants;

    mutable boost::optional<std::vector<Polynomial_1> > 
        principal_subresultants;
    mutable boost::optional<std::vector<Polynomial_1> > 
        coprincipal_subresultants;
        
    mutable boost::optional<Polynomial_1> resultant;

    mutable boost::optional<std::vector<Algebraic_real_1> > resultant_roots;
    mutable boost::optional<std::vector<Algebraic_real_1> > 
        event_x_coordinates;
    mutable boost::optional<std::vector<size_type> > 
        multiplicities_of_resultant_roots;

    mutable boost::optional<std::vector<Bound> > stripe_values;
        
    mutable std::vector< Lazy_status_line_CPA_1 > event_slices;

    mutable boost::optional<std::vector< Lazy_bound > > intermediate_values;

    mutable boost::optional< std::vector< Lazy_status_line_CPA_1 > >
        intermediate_slices;

    mutable boost::optional<std::vector<Event_indices> > event_indices;

    mutable Lazy_intersection_info_container intersection_info_container;

    typedef typename Curve_analysis_2::Integer Integer;

    CGAL::Degeneracy_strategy degeneracy_strategy;

    mutable CGAL::internal::Shear_controller<Integer> shear_controller;

    //! @}

    //! \name friends
    //! @{
    
    friend class Curve_pair_analysis_2<Algebraic_kernel_with_analysis_2>;

    //!@}

};

} // namespace internal

/*!
 * A model for <tt>AlgebraicKernelWithAnalysis_2::CurvePairAnalysis_2</tt>
 * It provides topological-geometric information about the intersection
 * points, and the vertical order of arcs of two algebraic plane curves.
 *
 * The curve pair is passed by two \c Curve_analysis_2 instances.
 * It is required that they do not share a component, i.e., the number
 * of common points must be finite. Note that overlapping curves are handled
 * by \c Algebraic_curve_kernel_2::Construct_curve_pair_2.
 * Also for caching reasons, it is recommended to construct curve pairs
 * always with this method.
 *
 * As for the single-curve analysis, the curve pair analysis is implemented
 * in a "lazy" fashion. That means, any computation is triggered when
 * the result is actually queried by the user. This prevents
 * expensive symbolic computations in some cases.
 *
 * For all algorithmic details of the curve pair analysis, we refer to
 * Arno Eigenwillig, Michael Kerber: Exact and Efficient 2D-Arrangements 
 * of Arbitrary Algebraic Curves. Proceedings of the Nineteenth Annual 
 * ACM-SIAM Symposium on Discrete Algorithms (SODA 2008), pp. 122-131
 */
template < typename AlgebraicKernelWithAnalysis_2 >
class Curve_pair_analysis_2 : 
    public ::CGAL::Handle_with_policy
        < CGAL::internal::Curve_pair_analysis_2_rep
              < AlgebraicKernelWithAnalysis_2 > > {
    

public:

    //! \name typedefs
    //! @{

    //! The algebraic kernel that uses the curve pair analysis
    typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

private:
    
    //! Representation class
    typedef CGAL::internal::Curve_pair_analysis_2_rep
      < Algebraic_kernel_with_analysis_2 > Rep;
    
    //! Base class
    typedef ::CGAL::Handle_with_policy< Rep >        Base;

public:
    //! The Curve_pair_analysis_2 itself
    typedef Curve_pair_analysis_2<Algebraic_kernel_with_analysis_2> Self;
    
    //! The corresponding Curve_analysis_2 class
    typedef typename Rep::Curve_analysis_2 Curve_analysis_2;

    //! Index type
    typedef typename Rep::size_type size_type;

    //! Univariate polynomials
    typedef typename Rep::Polynomial_1 Polynomial_1;

    //! Bivariate polynomials
    typedef typename Rep::Polynomial_2 Polynomial_2;

    //! Type for algebraic numbers (one dimension)
    typedef typename Rep::Algebraic_real_1 Algebraic_real_1;

    //! Type for points with algebraic coordinates
    typedef typename Algebraic_kernel_with_analysis_2::Algebraic_real_2 
        Algebraic_real_2;

    //! Bound type (for rational numbers)
    typedef typename Rep::Bound Bound;

private:
    // Optional for boundaries
    typedef typename Rep::Lazy_bound Lazy_bound;

    // Object to store information about intersection points 
    typedef typename Rep::Intersection_info_container 
        Intersection_info_container;

    // Its lazy version
    typedef typename Rep::Lazy_intersection_info_container 
        Lazy_intersection_info_container;

    // Type for indices of events.
    typedef typename Rep::Event_indices Event_indices;

    // Integer type
    typedef typename Curve_analysis_2::Integer Integer;

    // Status line of single curve analysis
    typedef typename Curve_analysis_2::Status_line_1 Status_line_CA_1;

    // Coefficient type
    typedef typename Curve_analysis_2::Coefficient Coefficient;
    
    // Polynomial traits class
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    // Polynomial traits class
    typedef CGAL::Polynomial_traits_d<Polynomial_1> Polynomial_traits_1;
  
public:
    
    //! The event slice object type
    typedef typename Rep::Status_line_CPA_1 Status_line_CPA_1;
    
    /*!
     * Required by the concept. The name is not used internally 
     * to distinguish from one curve status_lines syntactically
     */
    typedef Status_line_CPA_1 Status_line_1;

private:

    // Lazy version of status lines
    typedef typename Rep::Lazy_status_line_CPA_1 Lazy_status_line_CPA_1;

    // Coercion between Bound and Coefficient type
    typedef CGAL::Coercion_traits<Bound, Coefficient> Coercion;
    
    // The common supertype
    typedef typename Coercion::Type Coercion_type;

    // Polynomials over that supertype
    typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
        ::template Rebind<Coercion_type,1>::Other::Type Poly_coer_1;

    // Functor to isolate real roots of univariate polynomials
    typedef typename Algebraic_kernel_with_analysis_2::Solve_1 Solve_1;

    // Slice info objects
    typedef typename Rep::Slice_info Slice_info;

    // Lazy version
    typedef typename Rep::Lazy_slice_info Lazy_slice_info;

    //! @}

private:

    //! \name Internal structs
    //! @{

    struct Curves_at_event_functor {
        
        typedef size_type argument_type;
        typedef CGAL::internal::Slice_type result_type;

        Curves_at_event_functor(const Status_line_CPA_1& status_line) 
            : status_line(status_line)
        {}

        CGAL::internal::Slice_type operator() (size_type i) const {
            typedef typename Status_line_CPA_1::size_type 
                Status_line_size_type;
            std::pair<Status_line_size_type,Status_line_size_type> pair =
                status_line.curves_at_event(i);
            CGAL_assertion(pair.first>=0 || pair.second >=0);
            if(pair.first==-1) {
                return CGAL::internal::SECOND_CURVE;
            }
            if(pair.second==-1) {
                return CGAL::internal::FIRST_CURVE;
            }
            return CGAL::internal::INTERSECTION;
        }
    private:
        
        const Status_line_CPA_1& status_line;
        
    };

    typedef boost::transform_iterator<Curves_at_event_functor, 
                              boost::counting_iterator<size_type> > 
        Status_line_CPA_iterator;

    struct Xval_of_status_line_CA_1 {
        typedef Status_line_CA_1 argument_type;
        typedef Algebraic_real_1 result_type;
        Algebraic_real_1 operator() (const Status_line_CA_1& status_line) 
            const {
            return status_line.x();
        }
    };

    // @}

    //! \name Constructors
    //! @{
    
public:

    //! DefaultConstructible
    Curve_pair_analysis_2() :
        Base(Rep()) {
    };

    //! \brief Copy constructor
    Curve_pair_analysis_2(const Self& alg_curve_pair)
        : Base(static_cast<const Base&>(alg_curve_pair)) 
    {
    }

    // Assignable
    
    /*!
     * \brief Constructable from two curves
     *
     * Create a curve pair object for the two curves \c c1 and \c c2, 
     * given by their curve analysis object. The two curves are checked
     * to have no common vertical line component (if they have, an
     * exception of type \c CGAL::internal::Non_generic_position_exception
     * is thrown), no further computation is performed.
     *
     * \param strategy If a degenerate situation (e.g., two covertical 
     * intersection at the same x-coordinate) occurs during the analysis, 
     * this value controls the strategy to handle it. If set to
     * CGAL::EXCEPTION_STRATEGY, an exception of type
     * \c CGAL::internal::Non_generic_position_exception is thrown whenever
     * such a degeneracy occurs. If set to \c CGAL::SHEAR_STRATEGY, a shear
     * transformation is performed, and the sheared curve pair is used
     * to handle degenerate situations. Finally, if set to 
     * CGAL::SHEAR_ONLY_AT_IRRATIONAL_STRATEGY, degeneracies at rational
     * x-ccordinates are handled directly, and a shear is only applied
     * in other situations. The default argument for \c strategy is
     * \c CGAL::SHEAR_ONLY_AT_IRRATIONAL_STRATEGY.
     */
    Curve_pair_analysis_2(Algebraic_kernel_with_analysis_2* kernel,
                          Curve_analysis_2 c1, 
                          Curve_analysis_2 c2,
                          CGAL::Degeneracy_strategy strategy
                              = CGAL_ACK_DEFAULT_DEGENERACY_STRATEGY) 
        : Base(Rep(kernel,c1, c2, strategy)) 
    {
        
#if CGAL_ACK_DEBUG_FLAG
        CGAL::set_pretty_mode(CGAL_ACK_DEBUG_PRINT);
#endif
            
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Check content for squarefreeness.." 
                             << std::flush;
#endif
        if(CGAL::degree(this->ptr()->c1_.content())>0 &&
           CGAL::degree(this->ptr()->c2_.content())>0) {
            typename Polynomial_traits_1::Gcd_up_to_constant_factor gcd_utcf;
            if(CGAL::degree(gcd_utcf
                            (this->ptr()->c1_.content(), 
                             this->ptr()->c2_.content())) >= 1) {
                
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Common vertical line discovered" 
                                     << std::endl;
#endif
                throw CGAL::internal::Non_generic_position_exception();
            } else {
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
            }
        }
        
    }
    
    //! @}


private:

    // Computes the resultant of the defining polynomials wrt \c y
    void compute_resultant() const;

    // Computes the subresultant coefficients of the defining polynomials
    void compute_subresultants() const;

    /* 
     * Computes the roots of the resultants (via isolation) and their
     * multiplicities
     */
    void compute_resultant_roots_with_multiplicities() const;
    
    /*
     * Computes all x-events of the curve pair, 
     * together with their event indices 
     */
    void compute_event_x_coordinates_with_event_indices() const;

    /*
     * \brief Computes the intermediate x-coordinates and their status lines
     *
     * In fact, it only fills the data fields with boost::none instances,
     * according to the lazy philosophy of the whole class.
     */
    void compute_intermediate_values_and_slices() const;

public:

    Algebraic_kernel_with_analysis_2* kernel() const {
        return this->ptr()->_m_kernel;
    }

    //! Returns the resultant of the defing polynomials wrt \c y
    Polynomial_1 resultant() const {
        if(! this->ptr()->resultant) {
            compute_resultant();
        }
        CGAL_assertion(bool(this->ptr()->resultant));
        return this->ptr()->resultant.get();
    }

    std::vector<Algebraic_real_1>& resultant_roots() const {
        if(! this->ptr()->resultant_roots) {
            compute_resultant_roots_with_multiplicities();
        }
        CGAL_assertion(bool(this->ptr()->resultant_roots));
        return this->ptr()->resultant_roots.get();
    }

    Algebraic_real_1& resultant_roots(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>(resultant_roots().size()));
        return resultant_roots()[i];
    }

    std::vector<size_type>& multiplicities_of_resultant_roots() const {
        if(! this->ptr()->multiplicities_of_resultant_roots) {
            compute_resultant_roots_with_multiplicities();
        }
        CGAL_assertion(bool(this->ptr()->multiplicities_of_resultant_roots));
        return this->ptr()->multiplicities_of_resultant_roots.get();
    }    
    
    size_type multiplicities_of_resultant_roots(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>
                           (multiplicities_of_resultant_roots().size()));
        return multiplicities_of_resultant_roots()[i];
    }

    std::vector<Bound>& stripe_values() const {
        if(! this->ptr()->stripe_values) {
            this->ptr()->stripe_values = std::vector<Bound>();
            find_intermediate_values
	      (kernel(),
	       resultant_roots().begin(),
	       resultant_roots().end(),
	       std::back_inserter(this->ptr()->stripe_values.get()));
        }
        CGAL_assertion(bool(this->ptr()->stripe_values));
        return this->ptr()->stripe_values.get();
    }

    std::vector<Algebraic_real_1>& event_x_coordinates() const {
        if(! this->ptr()->event_x_coordinates) {
            compute_event_x_coordinates_with_event_indices();
        }
        CGAL_assertion(bool(this->ptr()->event_x_coordinates));
        return this->ptr()->event_x_coordinates.get();
    }

    std::vector<Event_indices>& event_indices() const {
        if(! this->ptr()->event_indices) {
            compute_event_x_coordinates_with_event_indices();
        }
        CGAL_assertion(bool(this->ptr()->event_indices));
        return this->ptr()->event_indices.get();
    } 

public:

    /*
     * \brief Returns the indices of the <tt>i</tt>th event value
     *
     * Returns a Event_indices <tt>(fg,ffy,ggy)</tt> such that
     * the <tt>i</tt>th event root is the <tt>fg</tt>th root of the 
     * resultant of \c f and \c g, the <tt>ffy</tt>th root of the 
     * discriminant of \c f, and  the <tt>ggy</tt>th root of the 
     * discriminant of \c g.
     */
    Event_indices event_indices(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>
                           (event_indices().size()));
        return event_indices()[i];
    } 

private:

    std::vector<Lazy_bound>& intermediate_values() const {
        if(! this->ptr()->intermediate_values) {
            compute_intermediate_values_and_slices();
        }
        CGAL_assertion(bool(this->ptr()->intermediate_values));
        return this->ptr()->intermediate_values.get();
    }

    std::vector<Lazy_status_line_CPA_1>& intermediate_slices() const {
        if(! this->ptr()->intermediate_slices) {
            compute_intermediate_values_and_slices();
        }
        CGAL_assertion(bool(this->ptr()->intermediate_slices));
        return this->ptr()->intermediate_slices.get();
    }


private:

    std::vector<Polynomial_2>& subresultants() const {
        if(! this->ptr()->subresultants) {
            compute_subresultants();
        }
        CGAL_assertion(bool(this->ptr()->subresultants));
        return this->ptr()->subresultants.get();
    }

    Polynomial_2& subresultants(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>(subresultants().size()));
        return subresultants()[i];
    }

    std::vector<Polynomial_1>& principal_subresultants() const {
        if(! this->ptr()->principal_subresultants) {
            compute_subresultants();
        }
        CGAL_assertion(bool(this->ptr()->principal_subresultants));
        return this->ptr()->principal_subresultants.get();
    }

    Polynomial_1& principal_subresultants(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>
                           (principal_subresultants().size()));
        return principal_subresultants()[i];
    }

    std::vector<Polynomial_1>& coprincipal_subresultants() const {
        if(! this->ptr()->coprincipal_subresultants) {
            compute_subresultants();
        }
        CGAL_assertion(bool(this->ptr()->coprincipal_subresultants));
        return this->ptr()->coprincipal_subresultants.get();
    }

    Polynomial_1& coprincipal_subresultants(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i < static_cast<size_type>
                           (coprincipal_subresultants().size()));
        return coprincipal_subresultants()[i];
    }
    


private:

    /* 
     * Refines the isolating intervals until they are disjoint
     * Returns CGAL::SMALLER, if the y-coordinate defined by <tt>(e1,i1)</tt>
     * is smaller than the y-coordinate <tt>(e2,i2)</tt>,
     * and CGAL::GREATER otherwise
     *
     * If both y-coordinates are equal, this method does not terminate
     */
    CGAL::Sign split_compare(Status_line_CA_1& e1, size_type i1,
                       Status_line_CA_1& e2, size_type i2) const {
        while(overlap(e1,i1,e2,i2)) {
            if(e1.interval_length(i1)<e2.interval_length(i2)) {
                e2.refine(i2);
            }
            else {
                e1.refine(i1);
            } 
        }
        return (e1.lower_bound(i1) < e2.lower_bound(i2)) 
            ? CGAL::SMALLER
            : CGAL::LARGER;
    }

private:

    /*!
     * TODO doc
     */
    Status_line_CPA_1 create_event_slice(size_type i) 
        const {
#if !CGAL_ACK_NO_ARC_FLIP
        size_type index_in_fg = event_indices(i).fg;
        if(index_in_fg == -1 ) {
            return create_slice_with_multiplicity_zero_or_one(i);
        } else {
            size_type mult_of_alpha 
                = multiplicities_of_resultant_roots(index_in_fg);
            if(mult_of_alpha == 1) {
                return create_slice_with_multiplicity_zero_or_one(i);
            } else {
#endif
                return create_slice_of_higher_multiplicity(i);
#if !CGAL_ACK_NO_ARC_FLIP
            }
        }
#endif
    }

    Status_line_CPA_1 create_slice_of_higher_multiplicity(size_type i) 
        const {
        bool is_resultant_root = event_indices(i).fg >=0;
        if(is_resultant_root &&
           this->ptr()->intersection_info_container) {
            return create_event_slice_with_shear(i);
        }
        try {
            Status_line_CPA_1 slice = construct_generic_case(i);

            return slice;
        } catch(CGAL::internal::Non_generic_position_exception ex) {
            switch(this->ptr()->degeneracy_strategy) {
            case(CGAL::EXCEPTION_STRATEGY): {
                throw ex;
                break;
            }
            case(CGAL::SHEAR_ONLY_AT_IRRATIONAL_STRATEGY): {
                if(event_x(i).is_rational()) {
                    return create_event_slice_at_rational(i);
                }

                CGAL_FALLTHROUGH;
            }
            case(CGAL::SHEAR_STRATEGY): {
                return create_event_slice_with_shear(i);
            }
            }
            
            
            // NEVER HAPPENS
            return Status_line_CPA_1();
            
        }
    }

private:
    Status_line_CPA_1 create_event_slice_at_rational(size_type i) const {
        
        Algebraic_real_1& x = event_x(i);

        CGAL_precondition(x.is_rational());
        Bound r = x.rational();

        int k = degree_of_local_gcd(event_indices(i).fg,x);
        Polynomial_2 sres = subresultants(k);
        
        Polynomial_1 gcd = kernel()->evaluate_utcf_2_object() 
	  (typename Polynomial_traits_2::Swap() (sres,0,1),r);
        std::vector<Algebraic_real_1> gcd_roots;
        kernel()->solve_1_object()(gcd,std::back_inserter(gcd_roots),false);
        int m = gcd_roots.size();

        Slice_info slice_info = construct_slice_info(x);
        reduce_number_of_candidates_and_intersections_to
            (m,
             this->ptr()->c1_.status_line_at_exact_x(x),
             this->ptr()->c2_.status_line_at_exact_x(x),
             slice_info);
        for(typename Slice_info::iterator it=slice_info.begin();
            it!=slice_info.end();
            it++) {

            if(it->first==CGAL::internal::CANDIDATE) {
                it->first=CGAL::internal::INTERSECTION;
            }
        }
        
        return create_slice_from_slice_info(i,slice_info,true);
    }

private:
    
    /*!
     * TODO doc
     */
    Status_line_CPA_1 create_slice_with_multiplicity_zero_or_one(size_type i) 
        const;
       
private:

    // Creates an intermediate slice at a rational value
    Status_line_CPA_1 create_intermediate_slice_at(int i) const;
    
private:

    // Create a slice with id \c id from the Slice_info object
    Status_line_CPA_1 create_slice_from_slice_info(size_type id,
                                                   const Slice_info& slice,
                                                   bool event_flag) const;

private:

    // Computes a slice_info object at Algebraic_real_1 \c alpha
    Slice_info construct_slice_info(Algebraic_real_1 alpha) const;
      
private:
      
    Status_line_CPA_1 construct_generic_case(size_type i) const;        

private:
      
    bool check_candidate_by_arc_pattern(size_type index,
                                        Status_line_CA_1& e1,
                                        size_type i1,
                                        Status_line_CA_1& e2,
                                        size_type i2) const;
private:

    /*
     * TODO update doc
     * Checks the point on e1 with index i1, and
     * the point on e2 with index i2 really intersect. The \c slice_info
     * is updated accordingly: If not intersecting, the corresponding 
     * points are refined until they can be arranged in the correct order.
     * If intersecting, the corresponding Slice_info element is set to 
     * INTERSECTION.
     */
    template<typename InputIterator>
    void check_candidate(Status_line_CA_1& e1,size_type i1,
                         Status_line_CA_1& e2,size_type i2,
                         size_type k,
                         Slice_info& slice_info,
                         InputIterator slice_it,
                         size_type root_index) const;

private:

    /*
     * Checks intersection with symbolic methods
     */
    bool check_candidate_symbolically(Status_line_CA_1& e1,size_type ,
                                      Status_line_CA_1& CGAL_precondition_code(e2),size_type ,
                                      size_type k) const {
        Polynomial_1 p = -coprincipal_subresultants(k-1);
        Polynomial_1 q = principal_subresultants(k)*Coefficient(k);
        Algebraic_real_1 alpha = e1.x();
        CGAL_assertion(alpha==e2.x());
        if(CGAL::internal::zero_test_bivariate
	   <Algebraic_kernel_with_analysis_2>
	     (kernel(),alpha,this->ptr()->f,p,q) && 
           CGAL::internal::zero_test_bivariate
	   <Algebraic_kernel_with_analysis_2>
	     (kernel(),alpha,this->ptr()->g,p,q)) {
            return true;
        }
        else {
            throw CGAL::internal::Non_generic_position_exception();
        }
        return false; // never happens
    }

private:

    /*
     * Checks whether the isolting intervals for the point on \c e1 with
     * index \c index1, and for the point on \c e2 with index \c index2
     * overlap
     */
    bool overlap(Status_line_CA_1& e1, 
                 size_type index1,
                 Status_line_CA_1& e2, 
                 size_type index2) const {
        if(e1.lower_bound(index1) > e2.upper_bound(index2)) {
            return false;
        }
        else if(e1.upper_bound(index1) < e2.lower_bound(index2)) {
            return false;
        }
        else {
            return true;
        }
    }

    /*
     * For the point \c p on \c e1 with index \c index1, find the 
     * unique point on \c e2 which might be equal to \c p. If no point
     * can be equal, -1 is returned.
     */
    size_type find_possible_matching(Status_line_CA_1& e1, 
                                     size_type index1,
                                     Status_line_CA_1& e2) const;


    size_type degree_of_local_gcd(size_type index_of_fg,
                            Algebraic_real_1 alpha) const {
        
        if(multiplicities_of_resultant_roots(index_of_fg) == 1) {
            return 1;
        } else {
            size_type k=1;
            while(kernel()->is_zero_at_1_object()
                  (principal_subresultants(k),alpha)) {
                k++;
            }
            return k;
        }
    }

public:

    //! Returns curve analysis for the cth curve
    Curve_analysis_2 curve_analysis(bool c) const {
        return c ? this->ptr()->c2_ : this->ptr()->c1_;
    }

    size_type event_of_curve_analysis(size_type i, bool c) const {
        Event_indices& ev_ind = event_indices(i);
        return c ? ev_ind.ggy : ev_ind.ffy;
    }

    size_type event_of_curve_analysis(size_type i, 
                                      const Curve_analysis_2& c) const {
        CGAL_assertion(c.id()==curve_analysis(false).id() ||
                       c.id()==curve_analysis(true).id());
        Event_indices& ev_ind = event_indices(i);
        return (c.id()==curve_analysis(false).id()) ? ev_ind.ffy : ev_ind.ggy;
    }

    /*! 
     * \brief Returns the number of event slices
     *
     * Precisely, this is the number of points which are either root of
     * the resultant of the two curves, or root of discriminant of one
     * of the curves
     */
    size_type number_of_status_lines_with_event() const {
        return static_cast<size_type>(event_x_coordinates().size());
    }

    //! Returns the x-coordinate of the <tt>i</tt>th event
    Algebraic_real_1& event_x(size_type i) const {
        CGAL_assertion(i>=0 && 
                       i<static_cast<size_type>(event_x_coordinates().size()));
        return event_x_coordinates()[i];
    }

    /*!
     * \brief The index of the x-coordinate
     *
     * For x-value \c x, the index of the suitable slice is computed. For
     * event value, the \c event flag is set to true, otherwise to false
     * and the slice of the interval to which \c x belongs is returned
     */
    void x_to_index(Algebraic_real_1 x, 
                    size_type& idx, bool& event) const {
        const std::vector<Algebraic_real_1>& sl = event_x_coordinates();
        idx = std::lower_bound(sl.begin(),
                               sl.end(),
                               x) - sl.begin();
        event = (idx < static_cast<size_type>(sl.size()) && (sl[idx] == x));

    }


    Status_line_CPA_1 status_line_for_x(Algebraic_real_1 x, 
                                        CGAL::Sign perturb = CGAL::ZERO) 
    const {
        size_type index;
        bool evt;
        x_to_index(x,index,evt);
        if(evt) {
            switch(perturb) {
            case(CGAL::ZERO): return status_line_at_event(index);
            case(CGAL::NEGATIVE): return status_line_of_interval(index);
            case(CGAL::POSITIVE): return status_line_of_interval(index+1);
            }
        } // else: 
        return status_line_of_interval(index);
        
        
    }
        

    Status_line_CPA_1 status_line_at_exact_x(Algebraic_real_1 x) {
        return status_line_for_x(x);
    }


public:
      
    //! Returns the Status_line_CPA_1 at the <tt>i</tt>th event
    const Status_line_CPA_1& status_line_at_event(size_type i) const {
        if(! this->ptr()->event_slices[i]) {
            this->ptr()->event_slices[i] = create_event_slice(i);
        }
        CGAL_assertion(bool(this->ptr()->event_slices[i]));
        return this->ptr()->event_slices[i].get();
    }
      


    //! Returns the Status_line_CPA_1 at the <tt>i</tt>th interval
    const Status_line_CPA_1& status_line_of_interval(size_type i) const {

        if(! intermediate_slices()[i]) {

            intermediate_slices()[i] 
                = create_intermediate_slice_at(i);
          
        }

        return intermediate_slices()[i].get();
    }
        
    //!  Returns bound representative value at the <tt>i</tt>th interval
    const Bound bound_value_in_interval(size_type i) const {

        const std::vector<Algebraic_real_1>& events = event_x_coordinates(); 

        if(! intermediate_values()[i]) {
            // Create the intermediate x-coordinate first
            if(events.size()==0) {
                CGAL_assertion(i==0);
                intermediate_values()[0]=Bound(0);
            } else {
                if(i==0) {
                    intermediate_values()[i] 
		      = bound_left_of(kernel(),events[i]);
                } else if(i == static_cast<size_type>(events.size())) {
                    intermediate_values()[i]
		      = bound_right_of(kernel(),events[i-1]);

                } else {
                    intermediate_values()[i]
		      = kernel()->bound_between_1_object()
		          (events[i-1],events[i]);
                }
            }
        }
        CGAL_assertion(bool(intermediate_values()[i]));
        return intermediate_values()[i].get();

    }
        
private:
      
    struct Bound_to_coercion_functor {
        
        typedef Bound argument_type;
        typedef Coercion_type result_type;

        result_type operator() (argument_type x) const {
            typename CGAL::Coercion_traits<Bound,Coefficient>::Cast cast;
            return cast(x);
        }
    };

    struct Coefficient_to_coercion_functor {
        
        typedef Coefficient argument_type;
        typedef Coercion_type result_type;

        result_type operator() (argument_type x) const {
            typename CGAL::Coercion_traits<Bound,Coefficient>::Cast cast;
            return cast(x);
        }
    };


    // If a new shear was used, update intersection multiplicities
    void merge_new_intersection_info
    (const Intersection_info_container& new_info_container) const {
        if(! this->ptr()->intersection_info_container) {
            // ok, nothing existed, so take the new intersection info
            this->ptr()->intersection_info_container
                = new_info_container;
            return;
        }
        Intersection_info_container& old_info_container
            = *(this->ptr()->intersection_info_container);
        size_type n = old_info_container.size();
        CGAL_assertion(n == static_cast<size_type>
                              ( new_info_container.size()));
        //iterate through the vector and update 
        // (-1 stands for "multiplicity unknown")
        for(size_type i=0;i<n;i++) {
            size_type m = old_info_container[i].size();
            CGAL_assertion(m == static_cast<size_type>\
                                  (new_info_container[i].size()));
            for(size_type j=0;j<m;j++) {
                old_info_container[i][j].mult
                  = (std::max)(new_info_container[i][j].mult,
                               old_info_container[i][j].mult);
            }
          
        }
    }

    void new_shear_for_intersection_info
        (Intersection_info_container& info_container) const;

    Status_line_CPA_1 create_event_slice_with_shear(size_type i) const {
        while(true) { // we know that it works at some point
            try {
                if(! this->ptr()->intersection_info_container) {
                    Intersection_info_container info_container;
                    new_shear_for_intersection_info(info_container);
                    merge_new_intersection_info(info_container);
                }		  
                Status_line_CPA_1 slice 
                    = create_event_slice_from_current_intersection_info(i);

                return slice;
            } catch(CGAL::internal::Non_generic_position_exception /* ex */) {
                // just try the next one
                Intersection_info_container info_container;
                new_shear_for_intersection_info(info_container);
                merge_new_intersection_info(info_container);
            }
        }
    }

      
    Status_line_CPA_1 
        create_event_slice_from_current_intersection_info (size_type i) const;

    Bound x_sheared(Bound x, Bound y,Integer sh) const {
        return x-sh*y;
    }

    void update_intersection_info(Intersection_info_container& 
                                  info_container,
                                  Self& sh_pair,
                                  Status_line_CPA_1 slice,
                                  size_type i,
                                  size_type j,
                                  Integer s) const;
      
    /*
     * \brief Reduces the number of possible intersections
     *
     * At the position given by the event lins \c e1 and \c e2 and the slice 
     * info object \c slice, the points on the event lines are further refined
     * until there are only \c n possible intersection points. The method can
     * be interrupted if all possible intersection points are known to have
     * a maximal intersection mulipicity smaller \c k, and a 
     * Non_generic_position_exception is thrown then.
     */
    size_type reduce_number_of_candidates_and_intersections_to
        (size_type n,
         Status_line_CA_1& e1,
         Status_line_CA_1& e2,
         Slice_info& slice,
         size_type k=-1) const;
        
    // Handle provides
    // .id()
    // .is_identical

    friend std::ostream& operator<< <>
               (std::ostream& out, 
                const Self& curve_pair);

}; // end of Curve_pair_analysis_2

//! \brief Prints the objects.
template<typename AlgebraicKernelWithAnalysis_2>
std::ostream& operator<< 
    (std::ostream& out, 
     const Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>& curve_pair) {
    typedef Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2> 
        Curve_pair_analysis_2;
    typedef typename Curve_pair_analysis_2::size_type size_type;
    typedef typename Curve_pair_analysis_2::Event_indices Event_indices;
    typedef typename Curve_pair_analysis_2::Status_line_CPA_1 Slice;
    out << "--------------- Analysis results ---------------" << std::endl;
    out << "Number of constructed event lines: " 
        << curve_pair.number_of_status_lines_with_event() 
        << std::endl;
      
    out << "Intermediate line: "  << std::flush;
    Slice slice=curve_pair.status_line_of_interval(0);
    
    out << slice.number_of_events() << " passing arcs" << std::endl ;
    out << "in order: " << std::flush;
    for(size_type i=0;i<slice.number_of_events();i++) {
        CGAL_assertion(slice.curves_at_event(i).first==-1 || 
                       slice.curves_at_event(i).second==-1 );
        if(slice.curves_at_event(i).second==-1) {
            out << "First " <<std::flush;
        } else {
            out << "Second " <<std::flush;
        }
    }          
    out << std::endl << std::endl;
    for(size_type j = 0; 
        j < curve_pair.number_of_status_lines_with_event();
        j++) {

        out << "Event line at " << CGAL::to_double(curve_pair.event_x(j))
            << ": " << std::endl;
        out << "Indices: ";
        Event_indices ev_ind = curve_pair.event_indices(j);
        out << "fg: " << ev_ind.fg << ", ffy: " 
            << ev_ind.ffy <<", ggy: " << ev_ind.ggy 
            << std::endl; 
        slice = curve_pair.status_line_at_event(j);
        out << slice.number_of_events() << " passing arcs" << std::endl ;
        out << "in order: " << std::flush;
        for(size_type i=0;i<slice.number_of_events();i++) {
            if(slice.curves_at_event(i).second==-1) {
                out << "First " <<std::flush;
            } else if(slice.curves_at_event(i).first==-1) {
                out << "Second " <<std::flush;
            } else {
                out << "Inter," << slice.multiplicity_of_intersection(i)
                    << " " << std::flush;
            }
          
        }   
        out << std::endl << std::endl;

        out << "Intermediate line:"  << std::flush;
        Slice slice=curve_pair.status_line_of_interval(j+1);
        out << slice.number_of_events() << " passing arcs" << std::endl ;
        out << "in order: " << std::flush;
        for(size_type i=0;i<slice.number_of_events();i++) {
            CGAL_assertion(slice.curves_at_event(i).first==-1 || 
                       slice.curves_at_event(i).second==-1 );
            if(slice.curves_at_event(i).second==-1) {
                out << "First " <<std::flush;
            } else {
                out << "Second " <<std::flush;
            }
            
        }          
        out << std::endl << std::endl;
    }
    out << "------------------------------------------------" << std::endl;
    
    return out;
}



// Implementation of functions from Curve_pair_analysis class

//////////////////// compute_resultant()

template <typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::compute_resultant()
    const {
    
#if CGAL_ACK_RESULTANT_FIRST_STRATEGY
#ifndef CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD
    bool speed_up = true;
#else
    bool speed_up = (std::min)
        (CGAL::degree(curve_analysis(false).polynomial_2(),1),
         CGAL::degree(curve_analysis(true).polynomial_2(),1)) >= 
         CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD;
#endif
#else
    bool speed_up=false;
#endif
    if(speed_up) {
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute the resultant of f and g..." 
                             << std::flush;
#endif
        this->ptr()->resultant 
            = CGAL::resultant(this->ptr()->f,this->ptr()->g);
    } else {
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute the subres-seq of f and g..." 
                             << std::flush;
#endif
        compute_subresultants();
        
        this->ptr()->resultant 
            = this->ptr()->principal_subresultants.get()[0];
    }
    
    
    if(this->ptr()->resultant.get().is_zero()) {
        throw CGAL::internal::Zero_resultant_exception<Polynomial_2>
            (this->ptr()->f,
             this->ptr()->g);
        }
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    
}

//////////////////// compute_resultant_roots_with_multiplicities()

template<typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
compute_resultant_roots_with_multiplicities() const {
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Isolate the real roots of resultant..." 
                         << std::flush;
#endif
    Solve_1 solve_1;
    this->ptr()->resultant_roots = std::vector<Algebraic_real_1>();
    this->ptr()->multiplicities_of_resultant_roots
        = std::vector<size_type>();
    std::vector<std::pair<Algebraic_real_1, size_type> > res_pairs;
    solve_1(resultant(), std::back_inserter(res_pairs));

    for(int i=0; i < static_cast<int>(res_pairs.size()); i++ ) {
        this->ptr()->resultant_roots.get().push_back(res_pairs[i].first);
        this->ptr()->multiplicities_of_resultant_roots.get()
            .push_back(res_pairs[i].second);
    }
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    
#if CGAL_ACK_DEBUG_FLAG
    for(size_type i = 0;
        i<static_cast<size_type>
            (this->ptr()->resultant_roots.get().size());
        i++) {
        CGAL_ACK_DEBUG_PRINT 
            << "Root at " 
            << CGAL::to_double(this->ptr()->resultant_roots.get()[i])
            << " with multiplicity "
            << this->ptr()->multiplicities_of_resultant_roots.get()[i]
            << std::endl;
    }
#endif
    
}

//////////////////// compute_event_x_coordinates_with_event_indices

template<typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
compute_event_x_coordinates_with_event_indices() const {
    
    Xval_of_status_line_CA_1 xval;
    const Curve_analysis_2& c1=this->ptr()->c1_, c2=this->ptr()->c2_;
    
    std::vector<Algebraic_real_1> one_curve_events;
    
    std::vector<CGAL::internal::Three_valued> one_curve_events_type;
    
    typename CGAL::Real_embeddable_traits<Algebraic_real_1>::Compare compare;
    
    CGAL::internal::set_union_with_source
        (::boost::make_transform_iterator(c1.event_begin(),xval),
         ::boost::make_transform_iterator(c1.event_end(),xval),
         ::boost::make_transform_iterator(c2.event_begin(),xval),
         ::boost::make_transform_iterator(c2.event_end(),xval),
         std::back_inserter(one_curve_events),
         std::back_inserter(one_curve_events_type),
         compare);
    
    this->ptr()->event_x_coordinates = std::vector<Algebraic_real_1>();
    std::vector<CGAL::internal::Three_valued> events_type;
    CGAL::internal::set_union_with_source
        (one_curve_events.begin(),
         one_curve_events.end(),
         resultant_roots().begin(),
         resultant_roots().end(),
         std::back_inserter(this->ptr()->event_x_coordinates.get()),
         std::back_inserter(events_type),
         compare);
    std::vector<Algebraic_real_1>& events 
        = this->ptr()->event_x_coordinates.get();
    
    typename std::vector<CGAL::internal::Three_valued>::iterator one_curve_it
        =one_curve_events_type.begin();
    size_type inter_count=0, f_count=0,g_count=0;
    this->ptr()->event_indices = std::vector<Event_indices>();
    std::vector<Event_indices>& event_indices
        = this->ptr()->event_indices.get();
    for(size_type i=0;i<static_cast<size_type>(events.size());i++) {
/*
        #if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << CGAL::to_double(events[i]) << std::flush;
        #endif
*/
        switch(events_type[i]) {
        case(CGAL::internal::ROOT_OF_FIRST_SET): {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << " one curve event" << std::endl;
#endif
*/
            this->ptr()->event_slices.push_back(Lazy_status_line_CPA_1());
            switch(*(one_curve_it++)) {
            case(CGAL::internal::ROOT_OF_FIRST_SET): {
                event_indices.push_back(Event_indices(-1,f_count,-1));
                f_count++;
                break;
            }
            case(CGAL::internal::ROOT_OF_SECOND_SET): {
                event_indices.push_back(Event_indices(-1,-1,g_count));
                g_count++;
                break;
            }
            case(CGAL::internal::ROOT_OF_BOTH_SETS): {
                event_indices.push_back(Event_indices(-1,f_count,g_count));
                f_count++;
                g_count++;
                break;
            }
            }
            break;
        }
        case(CGAL::internal::ROOT_OF_SECOND_SET): {
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << " two curve event" << std::endl;
#endif
*/
            this->ptr()->
                event_slices.push_back(Lazy_status_line_CPA_1());
            
            event_indices.push_back
                (Event_indices(inter_count,-1,-1));
            inter_count++;
            break;
        }
        case(CGAL::internal::ROOT_OF_BOTH_SETS): {
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << " one and two curve event" 
                                     << std::endl;
#endif
*/
            this->ptr()->event_slices.push_back(Lazy_status_line_CPA_1());
            
            
            switch(*(one_curve_it++)) {
            case(CGAL::internal::ROOT_OF_FIRST_SET): {
                event_indices.push_back
                    (Event_indices(inter_count,f_count,-1));
                f_count++;
                break;
            }
            case(CGAL::internal::ROOT_OF_SECOND_SET): {
                event_indices.push_back
                    (Event_indices(inter_count,-1,g_count));
                g_count++;
                break;
            }
            case(CGAL::internal::ROOT_OF_BOTH_SETS): {
                event_indices.push_back
                    (Event_indices(inter_count,f_count,g_count));
                f_count++;
                g_count++;
                break;
            }
            }
            inter_count++;
            break;
        }
        }
    }
    CGAL_assertion(inter_count 
                   == static_cast<size_type>
                   (resultant_roots().size()));
    CGAL_assertion(one_curve_it==one_curve_events_type.end());
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "done" << std::endl; 
#endif
    
}

//////////////////// compute_intermediate_values_and_slices()

template<typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
compute_intermediate_values_and_slices() const {
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Prepare intermediate slices.." << std::flush;
#endif
    this->ptr()->intermediate_values=std::vector<Lazy_bound>();
    this->ptr()->intermediate_slices=std::vector<Lazy_status_line_CPA_1>();
    
    for(size_type i=0;
        i<=static_cast<size_type>(event_x_coordinates().size());
        i++) {
        this->ptr()->intermediate_values.get().push_back(Lazy_bound());
        this->ptr()->intermediate_slices.get().push_back
            (Lazy_status_line_CPA_1());
    }
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
}

//////////////////// compute_subresultants

template<typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
compute_subresultants() const {
    typedef std::vector<Polynomial_1> Polynomial_container;
    this->ptr()->principal_subresultants = Polynomial_container();
    this->ptr()->coprincipal_subresultants = Polynomial_container();
    const Polynomial_2& f = this->ptr()->f, g=this->ptr()->g;
    this->ptr()->subresultants = std::vector<Polynomial_2>();
    if(CGAL::degree(f,1)<CGAL::degree(g,1)) {
#if CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS 
        CGAL::internal::bezout_polynomial_subresultants
            (g,f,std::back_inserter(this->ptr()->subresultants.get()));
#else
        typename CGAL::Polynomial_traits_d<Polynomial_2>
            ::Polynomial_subresultants()
            (g,f,std::back_inserter(this->ptr()->subresultants.get()));
#endif
    } else {
#if CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS 
        CGAL::internal::bezout_polynomial_subresultants
            (f,g,std::back_inserter(this->ptr()->subresultants.get()));
#else
        typename CGAL::Polynomial_traits_d<Polynomial_2>
            ::Polynomial_subresultants()
            (f,g,std::back_inserter(this->ptr()->subresultants.get()));
#endif
    }
    
    std::vector<Polynomial_2>& subresultants 
        = this->ptr()->subresultants.get();
    
    size_type n = static_cast<size_type>(subresultants.size());
    
    for(size_type i=0;i<n;i++) {
        if(CGAL::degree(subresultants[i]) < i) {
            this->ptr()->principal_subresultants->
                push_back(Polynomial_1(0));
        }
        else {
            this->ptr()->principal_subresultants->
                push_back(subresultants[i][i]);
        }
    }
    for(size_type i=1;i<n;i++) {
        if(CGAL::degree(subresultants[i]) < i-1) {
            this->ptr()->coprincipal_subresultants->
                push_back(Polynomial_1(0));
        }
        else {
            this->ptr()->coprincipal_subresultants->
                push_back(subresultants[i][i-1]);
        }
    }
    // This must be corrected, if f and g have same degree:
    if(CGAL::degree(f,1) == CGAL::degree(g,1)) {
        if(n>=1) {
            this->ptr()->principal_subresultants.get()[n-1] 
                = Polynomial_1(CGAL::leading_coefficient(g));
        }
        if(n>=2) {
            this->ptr()->coprincipal_subresultants.get()[n-2] 
                = Polynomial_1(g[CGAL::degree(g,1)-1]);
        }
    }
    
}

//////////////////// create_slice_with_multiplicity_zero_or_one

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>
    ::Status_line_CPA_1 
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
create_slice_with_multiplicity_zero_or_one(size_type i) const {
    
    const std::vector<Algebraic_real_1>& events 
        = event_x_coordinates();
    Algebraic_real_1 alpha = events[i];
    const Curve_analysis_2& c1=curve_analysis(false), c2=curve_analysis(true);
    size_type i1,i2;
    bool flag1,flag2;
    c1.x_to_index(alpha,i1,flag1);
    c2.x_to_index(alpha,i2,flag2);
    
    bool exactly_at_alpha_1 = flag1, exactly_at_alpha_2 = flag2;
    Status_line_CA_1 e1=flag1 ? c1.status_line_at_event(i1)
        : c1.status_line_of_interval(i1);
    Status_line_CA_1 e2=flag2 ? c2.status_line_at_event(i2)
        : c2.status_line_of_interval(i2);
    
    Status_line_CPA_1 left_slice = this->status_line_of_interval(i),
        right_slice = this->status_line_of_interval(i+1);
    
    Status_line_CPA_iterator left_it 
        = ::boost::make_transform_iterator
        (::boost::counting_iterator<size_type>(0),
         Curves_at_event_functor(left_slice));
    Status_line_CPA_iterator right_it 
        = ::boost::make_transform_iterator
        (::boost::counting_iterator<size_type>(0),
         Curves_at_event_functor(right_slice));

    Status_line_CPA_iterator left_end 
        = ::boost::make_transform_iterator
        (::boost::counting_iterator<size_type>
         (left_slice.number_of_events()),
         Curves_at_event_functor(left_slice));
    Status_line_CPA_iterator right_end 
        = ::boost::make_transform_iterator
        (::boost::counting_iterator<size_type>
         (right_slice.number_of_events()),
         Curves_at_event_functor(right_slice));

    // Take out asymptotes
    size_type asym_lm_1 
        = e1.number_of_branches_approaching_minus_infinity().first;
    size_type asym_rm_1 
        = e1.number_of_branches_approaching_minus_infinity().second;
    size_type asym_lp_1 
        = e1.number_of_branches_approaching_plus_infinity().first;
    size_type asym_rp_1 
        = e1.number_of_branches_approaching_plus_infinity().second;
    size_type asym_lm_2 
        = e2.number_of_branches_approaching_minus_infinity().first;
    size_type asym_rm_2 
        = e2.number_of_branches_approaching_minus_infinity().second;
    size_type asym_lp_2 
        = e2.number_of_branches_approaching_plus_infinity().first;
    size_type asym_rp_2 
        = e2.number_of_branches_approaching_plus_infinity().second;
            
    while(asym_lm_1 != 0 || asym_lm_2 != 0) {
        CGAL_assertion(*left_it != CGAL::internal::INTERSECTION);
        if(*left_it == CGAL::internal::FIRST_CURVE) {
            CGAL_assertion(asym_lm_1!=0);
            asym_lm_1--;
        }
        if(*left_it == CGAL::internal::SECOND_CURVE) {
            CGAL_assertion(asym_lm_2!=0);
            asym_lm_2--;
        }
        left_it++;
    }
    while(asym_rm_1 != 0 || asym_rm_2 != 0) {
        if(*right_it == CGAL::internal::FIRST_CURVE) {
            CGAL_assertion(asym_rm_1!=0);
            asym_rm_1--;
        }
        if(*right_it == CGAL::internal::SECOND_CURVE) {
            CGAL_assertion(asym_rm_2!=0);
            asym_rm_2--;
        }
        right_it++;
    }
    while(asym_lp_1 != 0 || asym_lp_2 != 0) {
        left_end--;
        if(*left_end == CGAL::internal::FIRST_CURVE) {
            CGAL_assertion(asym_lp_1!=0);
            asym_lp_1--;
        }
        if(*left_end == CGAL::internal::SECOND_CURVE) {
            CGAL_assertion(asym_lp_2!=0);
            asym_lp_2--;
        }
    }
    while(asym_rp_1 != 0 || asym_rp_2 != 0) {
        right_end--;
        if(*right_end == CGAL::internal::FIRST_CURVE) {
            CGAL_assertion(asym_rp_1!=0);
            asym_rp_1--;
        }
        if(*right_end == CGAL::internal::SECOND_CURVE) {
            CGAL_assertion(asym_rp_2!=0);
            asym_rp_2--;
        }
    }
    // Now, the iterator ranges [left_it,left_end)
    // and [right_it,right_end) give the arcs really
    // going into the event line

    Slice_info slice_info;
    CGAL::internal::Slice_type curr_lowest_arc;
    size_type curr_multiplicity;

    size_type event_index_1=0, event_index_2=0;

    while(event_index_1 != e1.number_of_events() || 
          event_index_2 != e2.number_of_events()) {

        CGAL_assertion(event_index_1 != e1.number_of_events() || 
                       event_index_2 != e2.number_of_events());
        if(event_index_1==e1.number_of_events()) {
            curr_lowest_arc=CGAL::internal::SECOND_CURVE;
        } 
        else if(event_index_2==e2.number_of_events()) {
            curr_lowest_arc=CGAL::internal::FIRST_CURVE;
        } 
        else if((e1.number_of_incident_branches(event_index_1).first>0 && 
                 e2.number_of_incident_branches(event_index_2).first>0)) {
            // The next arc on the left must come as next:
            curr_lowest_arc=*left_it;
        } 
        else if((e1.number_of_incident_branches(event_index_1).second>0 && 
                 e2.number_of_incident_branches(event_index_2).second>0)) {
            // The next arc on the right must come as next:
            curr_lowest_arc=*right_it;
        } 
        else {
            // We cannot decide it from the arcs, so we have to compare
            // isolating intervals
            if(! exactly_at_alpha_1) {
                e1 = c1.status_line_at_exact_x(alpha);
                CGAL_assertion(e1.number_of_events()>event_index_1);
            }
            if(! exactly_at_alpha_2) {
                e2 = c2.status_line_at_exact_x(alpha);
                CGAL_assertion(e2.number_of_events()>event_index_2);
            }
            CGAL::Sign e1_smaller 
                = split_compare(e1,event_index_1,e2,event_index_2);
            curr_lowest_arc 
                = (e1_smaller==CGAL::SMALLER) 
                ? CGAL::internal::FIRST_CURVE : CGAL::internal::SECOND_CURVE;
        }

        curr_multiplicity = -1;

        // Move the iterators
        size_type arcs_of_other_curve_left=0, arcs_of_other_curve_right=0;
        if(curr_lowest_arc==CGAL::internal::FIRST_CURVE) {
            size_type j=0;
            while(j<e1.number_of_incident_branches(event_index_1).first) {
                if(*left_it==CGAL::internal::FIRST_CURVE) {
                    j++;
                } else {
                    CGAL_assertion(event_index_2 < e2.number_of_events());
                    arcs_of_other_curve_left++;
                }
                left_it++;
            }

            j=0;
            while(j<e1.number_of_incident_branches(event_index_1).second) {
                if(*right_it==CGAL::internal::FIRST_CURVE) {
                    j++;
                } else {
                    CGAL_assertion(event_index_2 < e2.number_of_events());
                    arcs_of_other_curve_right++;
                }
                right_it++;
            }
            event_index_1++;
            if(arcs_of_other_curve_left+arcs_of_other_curve_right>0) {
                // Intersection! Iterate over the remaining arcs
                // on both sides belonging to this intersection
                for(size_type j=arcs_of_other_curve_left;
                    j<e2.number_of_incident_branches(event_index_2).first;
                    j++) {
                    CGAL_assertion(*left_it==CGAL::internal::SECOND_CURVE);
                    left_it++;
                }
                for(size_type j=arcs_of_other_curve_right;
                    j<e2.number_of_incident_branches(event_index_2).second;
                    j++) {
                    CGAL_assertion(*right_it==CGAL::internal::SECOND_CURVE);
                    right_it++;
                }
                event_index_2++;
                curr_lowest_arc=CGAL::internal::INTERSECTION;
                curr_multiplicity=1;
            }
        } else { // curr_lowest_arc=CGAL::internal::SECOND_CURVE
            size_type j=0;
            while(j<e2.number_of_incident_branches(event_index_2).first) {
                if(*left_it==CGAL::internal::SECOND_CURVE) {
                    j++;
                } else {
                    CGAL_assertion(event_index_1 < e1.number_of_events());
                    arcs_of_other_curve_left++;
                }
                left_it++;
            }
              
            j=0;
            while(j<e2.number_of_incident_branches(event_index_2).second) {
                if(*right_it==CGAL::internal::SECOND_CURVE) {
                    j++;
                } else {
                    CGAL_assertion(event_index_1 < e1.number_of_events());
                    arcs_of_other_curve_right++;
                }
                right_it++;
            }
            event_index_2++;
            if(arcs_of_other_curve_left+arcs_of_other_curve_right>0) {
                // Intersection! Iterate over the remaining arcs
                // on both sides belonging to this intersection
                for(size_type j=arcs_of_other_curve_left;
                    j<e1.number_of_incident_branches(event_index_1).first;
                    j++) {
                    CGAL_assertion(*left_it==CGAL::internal::FIRST_CURVE);
                    left_it++;
                }
                for(size_type j=arcs_of_other_curve_right;
                    j<e1.number_of_incident_branches(event_index_1).second;
                    j++) {
                    CGAL_assertion(*right_it==CGAL::internal::FIRST_CURVE);
                    right_it++;
                }
                event_index_1++;
                curr_lowest_arc=CGAL::internal::INTERSECTION;
                curr_multiplicity=1;
            }
            
        }
        slice_info.push_back(std::make_pair(curr_lowest_arc,
                                            curr_multiplicity));
    }
    CGAL_assertion(left_it == left_end && 
                   right_it == right_end);

    return create_slice_from_slice_info(i,slice_info,true);
}

//////////////////// create_intermediate_slice_at

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>
    ::Status_line_CPA_1 
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
create_intermediate_slice_at(int i) const {
    
    Bound r = bound_value_in_interval(i);

    std::vector<Algebraic_real_1> p1_roots,p2_roots;

    this->ptr()->c1_.get_roots_at_rational(r,std::back_inserter(p1_roots));
    this->ptr()->c2_.get_roots_at_rational(r,std::back_inserter(p2_roots));

    size_type number_of_roots 
        = static_cast<size_type>(p1_roots.size() + p2_roots.size());
    std::vector<Algebraic_real_1> p12_roots;
    p12_roots.reserve(number_of_roots);
    std::vector<CGAL::internal::Three_valued> p12_order;
    p12_order.reserve(number_of_roots);
    
    CGAL::internal::Distinct_compare<Algebraic_real_1> distinct_compare;
    set_union_with_source(p1_roots.begin(),
                          p1_roots.end(),
                          p2_roots.begin(),
                          p2_roots.end(),
                          std::back_inserter(p12_roots),
                          std::back_inserter(p12_order),
                          distinct_compare);
    
    Slice_info slice_info;
    
    for(typename std::vector<CGAL::internal::Three_valued>::const_iterator 
            it = p12_order.begin();
        it!=p12_order.end();
        it++) {
        switch(*it){
        case(CGAL::internal::ROOT_OF_FIRST_SET): {
            slice_info.push_back
                (std::make_pair(CGAL::internal::FIRST_CURVE,-1));
            break;
        }
        case(CGAL::internal::ROOT_OF_SECOND_SET): {
            slice_info.push_back
                (std::make_pair(CGAL::internal::SECOND_CURVE,-1));
            break;
        }
        case(CGAL::internal::ROOT_OF_BOTH_SETS): {  
            CGAL_assertion(false);
            break;
        }
        }
    }
    
    Status_line_CPA_1 new_slice 
        = create_slice_from_slice_info(i,slice_info,false);
    
    return new_slice;          
}

//////////////////// create_slice_from_slice_info

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>
    ::Status_line_CPA_1 
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
create_slice_from_slice_info(size_type id,
                             const Slice_info& slice,
                             bool event_flag) const {
    typedef typename Status_line_CPA_1::Arc_pair Arc_pair;
    typedef typename Status_line_CPA_1::Arc_container Arc_container;
    typedef typename Status_line_CPA_1::Int_container Int_container;
    Arc_container arc_container;
    Int_container int_container;
    
    for(typename Slice_info::const_iterator it = slice.begin();
        it!=slice.end();
        it++) {
        CGAL_assertion(it->first != CGAL::internal::CANDIDATE);
        switch(it->first) {
        case(CGAL::internal::FIRST_CURVE): {
            if(event_flag) {
                arc_container.push_back(Arc_pair(0,it->second));
            } else {
                int_container.push_back(0);
            }
            break;
        }
        case(CGAL::internal::SECOND_CURVE): {
            if(event_flag) {
                arc_container.push_back(Arc_pair(1,it->second));
            } else {
                int_container.push_back(1);
            }
            break;
        }
        case(CGAL::internal::INTERSECTION): {
            CGAL_assertion(event_flag);
            arc_container.push_back(Arc_pair(2,it->second));
            break;
        }
        case(CGAL::internal::CANDIDATE): {
            CGAL_assertion(false);
            break;
        }
        }
    }
    
    return event_flag 
        ? Status_line_CPA_1(id,arc_container,*this)
        : Status_line_CPA_1(id,int_container,*this);  
}

//////////////////// construct_slice_info

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::Slice_info 
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
construct_slice_info(Algebraic_real_1 alpha) const {
    
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "Consider alpha=" << CGAL::to_double(alpha) 
  << std::endl;
  #endif
*/      
    
    Status_line_CA_1 e1 = this->ptr()->c1_.status_line_at_exact_x(alpha);
    
    Status_line_CA_1 e2 = this->ptr()->c2_.status_line_at_exact_x(alpha);
    
    std::vector<std::pair<size_type,size_type> > matchings;
    for(size_type i=0;i<e1.number_of_events();i++) {
        size_type match=find_possible_matching(e1,i,e2);
        if(match==-1) {
            continue;
        }
        if(find_possible_matching(e2,match,e1) != i) {
            continue;
        }
/*
  #if CGAL_ACK_DEBUG_FLAG
  CGAL_ACK_DEBUG_PRINT << "New matching: (" << i 
  << "," << match << ")" << std::endl;
  #endif
*/
        matchings.push_back(std::make_pair(i,match));
    }
    size_type i1=0, i2=0,
        n1=e1.number_of_events(), n2=e2.number_of_events();
    
    typename std::vector<std::pair<size_type,size_type> >::const_iterator 
        match = matchings.begin();
    Slice_info slice_info;
    while(i1<n1 || i2<n2) {
        if(i1==n1) {
            slice_info.push_back
                (std::make_pair(CGAL::internal::SECOND_CURVE,-1));
            i2++;
            continue;
        }
        if(i2==n2) {
            slice_info.push_back
                (std::make_pair(CGAL::internal::FIRST_CURVE,-1));
            i1++;
            continue;
        }
        if(match!=matchings.end() && 
           i1==match->first && 
           i2==match->second) {
            slice_info.push_back(std::make_pair(CGAL::internal::CANDIDATE,1));
            i1++;
            i2++;
            match++;
            continue;
        }
        CGAL_assertion(!overlap(e1,i1,e2,i2));
        if(e1.lower_bound(i1) < e2.lower_bound(i2)) {
            slice_info.push_back
                (std::make_pair(CGAL::internal::FIRST_CURVE,-1));
            i1++;
            continue;
        } else {
            slice_info.push_back
                (std::make_pair(CGAL::internal::SECOND_CURVE,-1));
            i2++;
            continue;
        }
    }
    CGAL_assertion(match==matchings.end());
    return slice_info;
}

//////////////////// construct_generic_case

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>
    ::Status_line_CPA_1 
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
construct_generic_case(size_type i) const {
    
    Algebraic_real_1 alpha = event_x(i);
    
    Slice_info slice_info;
    
    size_type index_of_fg = event_indices(i).fg;
    size_type index_of_ffy =event_indices(i).ffy;
    size_type index_of_ggy =event_indices(i).ggy;
    if(index_of_fg>=0) {
        if(kernel()->is_zero_at_1_object() 
             (CGAL::leading_coefficient
              (this->ptr()->c1_.polynomial_2()),alpha)
           ||
           kernel()->is_zero_at_1_object() 
             (CGAL::leading_coefficient
              (this->ptr()->c2_.polynomial_2()),alpha)) {
            throw CGAL::internal::Non_generic_position_exception();
        }
        size_type k = -1; // not yet computed
        if(index_of_ffy==-1 && index_of_ggy==-1) {
            // this means, we need the multiplicity of the intersections
            if(kernel()->is_zero_at_1_object() 
               (principal_subresultants(1),alpha)) {
                // multiplicity cannot be determined, throw exception
                throw CGAL::internal::Non_generic_position_exception();
            } else {
                k=1;
            }
        } else {
            k = degree_of_local_gcd(index_of_fg,alpha);
        }
        Status_line_CA_1 e1 
            = this->ptr()->c1_.status_line_at_exact_x(alpha);
        Status_line_CA_1 e2 
            = this->ptr()->c2_.status_line_at_exact_x(alpha);
        slice_info = construct_slice_info(alpha);
        size_type no_candidates=
            reduce_number_of_candidates_and_intersections_to
            (1,e1,e2,slice_info,k);
        CGAL_assertion(no_candidates==0 || no_candidates==1);
        if(no_candidates==1) {
            typename Slice_info::iterator slice_it
                = slice_info.begin();
            size_type i1=0,i2=0;
            while(slice_it->first!=CGAL::internal::CANDIDATE) {
                if(slice_it->first==CGAL::internal::FIRST_CURVE) {
                    i1++;
                }
                if(slice_it->first==CGAL::internal::SECOND_CURVE) {
                    i2++;
                }
                if(slice_it->first==CGAL::internal::INTERSECTION) {
                    i1++;
                    i2++;
                }
                slice_it++;
            }
            check_candidate(e1,i1,e2,i2,k,
                            slice_info,
                            slice_it,i);
        }
    } else {
        Status_line_CA_1 e1 
            = this->ptr()->c1_.status_line_at_exact_x(alpha);
        Status_line_CA_1 e2 
            = this->ptr()->c2_.status_line_at_exact_x(alpha);
        slice_info = construct_slice_info(alpha);
        reduce_number_of_candidates_and_intersections_to
            (0,e1,e2,slice_info,0);
    }
    return create_slice_from_slice_info(i,slice_info,true);
}

//////////////////// check_candidate_by_arc_pattern

template<typename AlgebraicKernelWithAnalysis_2>

bool Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
check_candidate_by_arc_pattern(size_type index,
                               Status_line_CA_1& e1,
                               size_type i1,
                               Status_line_CA_1& e2,
                               size_type i2) const {
    
    Status_line_CPA_1 left_slice = status_line_of_interval(index),
        right_slice = status_line_of_interval(index+1);
    size_type left_index=0,right_index=0;
    for(size_type i=0;i<i1;i++) {
        left_index+=e1.number_of_incident_branches(i).first;
        right_index+=e1.number_of_incident_branches(i).second;
    }
    for(size_type i=0;i<i2;i++) {
        left_index+=e2.number_of_incident_branches(i).first;
        right_index+=e2.number_of_incident_branches(i).second;
    }
    // left_index and right_index now point at the position
    // of the first arc into the candidate
    size_type num_of_arcs_to_candidate_left
        = e1.number_of_incident_branches(i1).first 
        + e2.number_of_incident_branches(i2).first,
        num_of_arcs_to_candidate_right
        = e1.number_of_incident_branches(i1).second 
        + e2.number_of_incident_branches(i2).second;
    CGAL_assertion(left_index + num_of_arcs_to_candidate_left <= 
                   left_slice.number_of_events());
    CGAL_assertion(right_index + num_of_arcs_to_candidate_right <=
                   right_slice.number_of_events());
    
    CGAL::internal::Slice_type curr;
    Curves_at_event_functor left_functor(left_slice);
    size_type number_of_changes;
    if(left_index < left_slice.number_of_events()) {
        curr = left_functor(left_index);
        number_of_changes=0;
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT <<  num_of_arcs_to_candidate_left 
                             << num_of_arcs_to_candidate_right 
                             << left_index << right_index << std::endl;
#endif
*/          
        for(size_type i=1;i<num_of_arcs_to_candidate_left;i++) {
            if(curr != left_functor(left_index+i)) {
                curr = left_functor(left_index+i);
                number_of_changes++;
                if(number_of_changes>=2) {
                    return true;
                }
            }
        }
    }
    
    Curves_at_event_functor right_functor(right_slice);
    
    if(right_index < right_slice.number_of_events()) {
        curr = right_functor(right_index);
        number_of_changes=0;
        for(size_type i=1;i<num_of_arcs_to_candidate_right;i++) {
            if(curr != right_functor(right_index+i)) {
                curr = right_functor(right_index+i);
                number_of_changes++;
                if(number_of_changes>=2) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

//////////////////// check_candidate

template<typename AlgebraicKernelWithAnalysis_2>
template<typename InputIterator>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
check_candidate(Status_line_CA_1& e1,size_type i1,
                Status_line_CA_1& e2,size_type i2,
                size_type k,
                Slice_info& slice_info,
                InputIterator slice_it,
                size_type root_index) const {
    
    Algebraic_real_1 xval = e1.x();
    CGAL_assertion(xval==e2.x());
    bool is_intersection=false;
    size_type index_in_fg = event_indices(root_index).fg;
    size_type mult_of_resultant
        = multiplicities_of_resultant_roots(index_in_fg);
    
    if(k%2==1 || mult_of_resultant%2==1) {
        is_intersection=true;
    } else {
        if(check_candidate_by_arc_pattern(root_index,e1,i1,e2,i2)) {
            is_intersection=true;
        } else {
            is_intersection=check_candidate_symbolically(e1,i1,e2,i2,k);
        }
    }
    if(is_intersection) {
        slice_it=slice_info.erase(slice_it);
        size_type mult_of_intersection;
        if(k==1) {
            mult_of_intersection = mult_of_resultant;
        } else {
            mult_of_intersection = -1; 
        }
        slice_info.insert(slice_it,
                          std::make_pair(CGAL::internal::INTERSECTION,
                                         mult_of_intersection));
    }
    else {
        CGAL::Sign e1_smaller=split_compare(e1,i1,e2,i2);
        
        slice_it=slice_info.erase(slice_it);
        if(e1_smaller==CGAL::SMALLER) {
            slice_it = slice_info.insert
                (slice_it,std::make_pair(CGAL::internal::FIRST_CURVE,-1));
            slice_it++;
            slice_it = slice_info.insert
                (slice_it,std::make_pair(CGAL::internal::SECOND_CURVE,-1));
        } else {
            slice_it = slice_info.insert
                (slice_it,std::make_pair(CGAL::internal::SECOND_CURVE,-1));
            slice_it++;
            slice_it = slice_info.insert
                (slice_it,std::make_pair(CGAL::internal::FIRST_CURVE,-1));
        }
    }
}

//////////////////// find_possible_matching

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::size_type
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
find_possible_matching(Status_line_CA_1& e1, 
                       size_type index1,
                       Status_line_CA_1& e2) const {
    
    std::vector<size_type> possible_overlaps;
    for(size_type i=0;i<e2.number_of_events();i++) {
        if(overlap(e1,index1,e2,i)) {
            possible_overlaps.push_back(i);
        }
    }
    while(possible_overlaps.size()>1) {
        if(possible_overlaps.size()==2) {
            // Prevent that both intervals touch in a bound
            while(overlap(e2,possible_overlaps[0],
                          e2,possible_overlaps[1])) {
                e2.refine(possible_overlaps[0]);
                e2.refine(possible_overlaps[1]);
            }
        }
        e1.refine(index1);
        
        typename std::vector<size_type>::iterator it
            = possible_overlaps.begin();
        while(it!=possible_overlaps.end()) {
            if(!overlap(e1,index1,e2,*it)) {
                it=possible_overlaps.erase(it);
            }
            else {
                it++;
            }
        }
    }
    if(possible_overlaps.size()==0) {
        return -1;
    }
    else {
        return possible_overlaps[0];
    }
}


//////////////////// new_shear_for_intersection_info

template<typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
new_shear_for_intersection_info(Intersection_info_container& info_container) 
    const {
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Use shear for intersections.." << std::endl;
#endif
    bool good_direction_found=false;
    Integer s;
    
    while(! good_direction_found) {
        try {
            info_container.clear();
            info_container.resize(resultant_roots().size());               
            s = this->ptr()->shear_controller.get_shear_factor();
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Try shear factor " << s << std::endl;
            CGAL_ACK_DEBUG_PRINT 
                << ">>>>>>>>>>> Transform first curve"  << std::endl;
#endif
            
            Curve_analysis_2 sh1 
                = this->ptr()->c1_.shear_primitive_part(s);
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT 
                << "<<<<<<<<<<< End of transform first curve" << std::endl;
            CGAL_ACK_DEBUG_PRINT << ">>>>>>>>>>> Transform second curve" 
                                 << std::endl;
#endif
            Curve_analysis_2 sh2 = this->ptr()->c2_.shear_primitive_part(s);
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT 
                << "<<<<<<<<<<< End of transform second curve" 
                << std::endl;
#endif
            Self sh_pair(kernel(),sh1,sh2,CGAL::EXCEPTION_STRATEGY);
            
#if CGAL_ACK_DEBUG_FLAG 
            CGAL_ACK_DEBUG_PRINT << "Shear back intersection points..." 
                                 << std::flush;
#endif
            for(size_type i=0;
                i<static_cast<size_type>
                    (sh_pair.event_x_coordinates().size());
                i++) {
                if(sh_pair.event_indices(i).fg==-1) {
                    continue;
                }
                Status_line_CPA_1 slice 
                    = sh_pair.status_line_at_event(i);
                Curves_at_event_functor functor(slice);
                for(size_type j=0;j<slice.number_of_events();j++) {
                    if(functor(j) == CGAL::internal::INTERSECTION) {
                        this->update_intersection_info(info_container,
                                                       sh_pair,
                                                       slice,
                                                       i,j,s);
                    }
                }
            }
            good_direction_found=true;
        }
        catch(CGAL::internal::Non_generic_position_exception /* ex */) {
            this->ptr()->shear_controller.report_failure(s);
        }
    }
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    return;
}

//////////////////// create_event_slice_from_current_intersection_info

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>
    ::Status_line_CPA_1 
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
create_event_slice_from_current_intersection_info (size_type i) const{
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Reduce the candidates.." << std::flush;
#endif
    Event_indices ev_ind = event_indices(i);
    size_type index_of_fg = ev_ind.fg;
    Intersection_info_container& intersection_info_container
        = *(this->ptr()->intersection_info_container);
    Algebraic_real_1 alpha = event_x(i);
    CGAL_assertion(index_of_fg <
                   static_cast<size_type>
                   (intersection_info_container.size()));
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << i << "th slice has " 
                         << intersection_info_container[index_of_fg].size()
                         << " intersections" << std::endl;
#endif
    Status_line_CA_1 e1=this->ptr()->c1_.
        status_line_at_exact_x(resultant_roots(index_of_fg)),
        e2=this->ptr()->c2_.
        status_line_at_exact_x(resultant_roots(index_of_fg));
    Slice_info slice=construct_slice_info(alpha);
    CGAL_assertion_code(size_type no_intersections=)
        reduce_number_of_candidates_and_intersections_to
        (static_cast<size_type>
         (intersection_info_container[index_of_fg].size()),
         e1,
         e2,
         slice,
         -1);
    CGAL_assertion(no_intersections==static_cast<size_type>
                   (intersection_info_container[index_of_fg].size()));
    typename std::vector<typename Rep::Intersection_info>::iterator 
        inter_info_it 
        = intersection_info_container[index_of_fg].begin();
    for(size_type j=0;j<static_cast<size_type>(slice.size());j++) {
        if(slice[j].first==CGAL::internal::INTERSECTION) {
            inter_info_it++;
        }
        if(slice[j].first==CGAL::internal::CANDIDATE) {
            slice[j].first=CGAL::internal::INTERSECTION;
            if(ev_ind.ffy==-1 && ev_ind.ggy==-1 && inter_info_it->mult==-1) {
                // Multiplicity unknown for case where we need it
                throw CGAL::internal::Non_generic_position_exception();
            }
            slice[j].second=inter_info_it->mult;              
            inter_info_it++;
        }
    }
    
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    return create_slice_from_slice_info(i,slice,true);
    
}

//////////////////// update_intersection_info

template<typename AlgebraicKernelWithAnalysis_2>
void Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
update_intersection_info(Intersection_info_container& 
                         info_container,
                         Self& sh_pair,
                         Status_line_CPA_1 slice,
                         size_type i,
                         size_type j,
                         Integer s) const {
    typedef typename Rep::Intersection_info Intersection_info;
    const Algebraic_real_1& xval = sh_pair.event_x(i);
    CGAL_assertion(Curves_at_event_functor(slice)(j)
                   ==CGAL::internal::INTERSECTION);
    Status_line_CA_1 ev = sh_pair.ptr()->c1_.status_line_at_exact_x(xval);
    // x_coordinate is given by xval
    // y_coordinate by ev[index]
    Intersection_info intersection_info;
    intersection_info.ev=ev;
    int index = slice.curves_at_event(j).first;
    intersection_info.index=index;
    intersection_info.mult=slice.multiplicity_of_intersection(j);
    // Find the right position to insert the object
    // first the "x-coordiante"
    size_type left_index = -1, 
        right_index = static_cast<size_type>(stripe_values().size()-1);
    Algebraic_real_1 xv = ev.x();
    Bound lx = xv.low(), rx=xv.high(),
        x_iv_size = rx-lx;
    Bound ly = ev.lower_bound(index),
        ry = ev.upper_bound(index);;
    while(left_index < right_index) {
        if(x_iv_size > ry-ly) {
            xv.refine();
            lx = xv.low();
            rx=xv.high();
            x_iv_size=rx-lx;
            continue;
        }
        ev.refine(index);
        ly = ev.lower_bound(index);
        ry = ev.upper_bound(index);
        Bound left=(s<0) ? x_sheared(lx,ry,-s): x_sheared(lx,ly,-s);
        Bound right = (s<0) ? x_sheared(rx,ly,-s) : x_sheared(rx,ry,-s);
        CGAL_assertion(left<right);
        while(left_index<right_index && 
              stripe_values()[left_index+1]<left) {
            ++left_index;
        }
        while(left_index<right_index && 
              right<stripe_values()[right_index]) {
            --right_index;
        }
    }
    CGAL_assertion(left_index==right_index);
        
    // Now, the "y-coordinate"
    typedef std::vector<Intersection_info> Intersection_info_vector;
    Intersection_info_vector& info_vec 
        = info_container[left_index];
    typename Intersection_info_vector::iterator info_it=info_vec.begin();
    while(info_it!=info_vec.end()) {
        Status_line_CA_1& comp_ev=info_it->ev;
        size_type comp_index = info_it->index;
        CGAL::Sign index_smaller 
            = split_compare(ev,index,comp_ev,comp_index);
        if(index_smaller==CGAL::LARGER) {
            info_it++;
        } else {
            break;
        }
    }
    info_vec.insert(info_it,intersection_info);        
}

//////////////////// reduce_number_of_candidates_and_intersections_to

template<typename AlgebraicKernelWithAnalysis_2>
typename Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::size_type
Curve_pair_analysis_2<AlgebraicKernelWithAnalysis_2>::
reduce_number_of_candidates_and_intersections_to(size_type n,
                                                 Status_line_CA_1& e1,
                                                 Status_line_CA_1& e2,
                                                 Slice_info& slice,
                                                 size_type k) const {
/*
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Reduce: " << n << " " 
                         << CGAL::to_double(e1.x()) << " " << k 
                         << std::endl;
#endif
*/
    size_type number_of_intersections=0;
    size_type number_of_candidates=0;
    for(size_type i=0;i<static_cast<size_type>(slice.size());i++) {
        if(slice[i].first==CGAL::internal::CANDIDATE) {
            number_of_candidates++;
        }
        if(slice[i].first==CGAL::internal::INTERSECTION) {
            number_of_intersections++;
        }
    }
    CGAL_assertion(number_of_intersections<=n);
    
    typename Slice_info::iterator slice_it=slice.begin();
    size_type i1=0,i2=0;
    size_type max_candidate_mult=0;
    while(n<number_of_candidates+number_of_intersections) {
        if(slice_it==slice.end()) {
            CGAL_assertion(e1.number_of_events()==i1 && 
                           e2.number_of_events()==i2);
            if(max_candidate_mult<k) {
                throw CGAL::internal::Non_generic_position_exception();
            } else {
                slice_it=slice.begin();
                max_candidate_mult=0;
                i1=i2=0;
            }
        }
        switch(slice_it->first) {
        case(CGAL::internal::FIRST_CURVE): {
            i1++;
            break;
        }
        case(CGAL::internal::SECOND_CURVE): {
            i2++;
            break;
        }
        case(CGAL::internal::CANDIDATE): {
            if(e1.interval_length(i1)<e2.interval_length(i2)) {
                e2.refine(i2);
            }
            else {
                e1.refine(i1);
            }
            if(! overlap(e1,i1,e2,i2)) {
                number_of_candidates--;
                slice_it=slice.erase(slice_it);
                if(e1.lower_bound(i1)<e2.lower_bound(i2)) {
                    slice_it=slice.insert
                        (slice_it,std::make_pair(CGAL::internal::FIRST_CURVE,-1));
                    slice_it++;
                    slice_it=slice.insert
                        (slice_it,std::make_pair
                             (CGAL::internal::SECOND_CURVE,-1));
                } else {
                    slice_it=slice.
                        insert(slice_it,std::make_pair
                                   (CGAL::internal::SECOND_CURVE,-1));
                    slice_it++;
                    slice_it=slice.
                        insert(slice_it,std::make_pair
                                   (CGAL::internal::FIRST_CURVE,-1));
                }
            } else {
                size_type m1 = e1.get_upper_bound_for_multiplicity(i1),
                    m2 = e2.get_upper_bound_for_multiplicity(i2),
                    min_m = m1<m2 ? m1 : m2;
                max_candidate_mult = max_candidate_mult>min_m 
                    ? max_candidate_mult : min_m;
            }
            i1++;
            i2++;
            break;
        }
        case(CGAL::internal::INTERSECTION): {
            i1++;
            i2++;
            break;
        }
        }
        slice_it++;
    }
    return number_of_intersections+number_of_candidates;
}      

} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif


#endif
