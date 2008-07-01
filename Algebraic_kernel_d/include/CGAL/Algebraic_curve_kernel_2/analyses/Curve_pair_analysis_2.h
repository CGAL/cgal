// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
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

#include <CGAL/Algebraic_curve_kernel_2/analyses/Zero_resultant_exception.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/Shear_controller.h>
#include <CGAL/Algebraic_curve_kernel_2/analyses/Shear_transformation.h>
#include <CGAL/Algebraic_curve_kernel_2/Bitstream_descartes_at_x/Non_generic_position_exception.h>

#include <CGAL/Algebraic_curve_kernel_2/Status_line_CPA_1.h>

#if CGAL_USE_LEDA
#if CGAL_LEDA_VERSION < 500
  #include <LEDA/list.h>
#else
  #include <LEDA/core/list.h>
#endif
#endif

CGAL_BEGIN_NAMESPACE

// Needed size_typeernally
  
namespace CGALi {
    
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

}// namespace CGALi

//////////////////////////////////////////////////////////////////////////////
// Curve_pair_2

template < typename AlgebraicKernel_2 >
struct Curve_pair_analysis_2;

namespace CGALi {

enum Slice_type {
    FIRST_CURVE = 0,
    SECOND_CURVE = 1,
    INTERSECTION = 2,
    CANDIDATE = 3
};

template<typename size_type>
struct Index_triple {

    size_type fg;
    size_type ffy;
    size_type ggy;
    Index_triple(size_type fg,size_type ffy, size_type ggy) : fg(fg), ffy(ffy), ggy(ggy) {}
};


template < class AlgebraicKernel_2 >
struct Curve_pair_analysis_2_rep {

    typedef AlgebraicKernel_2 Algebraic_kernel_2;

    typedef Curve_pair_analysis_2_rep<Algebraic_kernel_2> Self;

    typedef Curve_pair_analysis_2<Algebraic_kernel_2> Handle;

    typedef typename Algebraic_kernel_2::Curve_analysis_2 Curve_analysis_2;

    typedef typename Curve_analysis_2::size_type size_type;

    typedef typename Curve_analysis_2::Polynomial_2 Polynomial_2;
        
    typedef typename Curve_analysis_2::X_coordinate_1 X_coordinate_1;
        
    typedef typename Polynomial_2::NT Polynomial_1;

    typedef typename Curve_analysis_2::Boundary Boundary;

    typedef CGAL::CGALi::Status_line_CPA_1<Handle> Status_line_CPA_1;

    Curve_analysis_2 c1_;
    Curve_analysis_2 c2_;

    Polynomial_2 f;
    Polynomial_2 g;
    
    // DefaultConstructible
    Curve_pair_analysis_2_rep() :
        c1_(), c2_() {
    }

    Curve_pair_analysis_2_rep(Curve_analysis_2 c1, Curve_analysis_2 c2) :
        c1_(c1), c2_(c2) {
    }
    
    mutable boost::optional<std::vector<Polynomial_1> > 
    principal_subresultants;
    mutable boost::optional<std::vector<Polynomial_1> > 
    coprincipal_subresultants;
        
    Polynomial_1 resultant;

    std::vector<X_coordinate_1> resultant_roots, event_x_coordinates;
    std::vector<size_type> multiplicities_of_resultant_roots;

    mutable std::vector<Boundary> stripe_values;
        
    typedef std::pair<Slice_type,size_type> Slice_element;

    typedef std::vector<Slice_element> Slice_info;

    typedef boost::optional<Slice_info> Lazy_slice_info;

    std::vector< Lazy_slice_info > slice_info_container;

    // For lazy evaluation of Status_line_CPA_1s.
    typedef boost::optional<Status_line_CPA_1> Lazy_status_line_1;

    mutable std::vector< Lazy_status_line_1 > event_slices;

    typedef boost::optional<Boundary> Lazy_boundary;

    mutable std::vector< Lazy_boundary > intermediate_values;

    mutable std::vector< Lazy_status_line_1 > intermediate_slices;

    typedef CGAL::CGALi::Index_triple<size_type> Index_triple;

    std::vector<Index_triple> index_triples;

    struct Intersection_info {
        typename Curve_analysis_2::Status_line_1 ev;
        size_type index;
        size_type mult;
    };
        
    typedef std::vector<std::vector<Intersection_info> > 
    Intersection_info_container;
        
    typedef boost::optional<Intersection_info_container> 
    Lazy_intersection_info_container;

    mutable Lazy_intersection_info_container intersection_info_container;

    typedef typename Curve_analysis_2::Integer Integer;

    mutable CGAL::CGALi::Shear_controller<Integer> shear_controller;

    friend class Curve_pair_analysis_2<Self>;

};

} // namespace CGALi

  /*!
   * \brief A non-trivial model of an CurvePairAnalysis_2. Analysis is 
   * performed
   */
template < typename AlgebraicKernel_2 >
class Curve_pair_analysis_2 : 
    public ::CGAL::Handle_with_policy< CGAL::CGALi::Curve_pair_analysis_2_rep< AlgebraicKernel_2 > > {

public:

    typedef AlgebraicKernel_2 Algebraic_kernel_2;

    typedef CGAL::CGALi::Curve_pair_analysis_2_rep< Algebraic_kernel_2 > Rep;
    typedef ::CGAL::Handle_with_policy< Rep >        Base;

    typedef Curve_pair_analysis_2<Algebraic_kernel_2> Self;

    typedef typename Rep::Curve_analysis_2 Curve_analysis_2;
    typedef typename Rep::size_type size_type;
    typedef typename Rep::Polynomial_1 Polynomial_1;
    typedef typename Rep::Polynomial_2 Polynomial_2;
    typedef typename Rep::X_coordinate_1 X_coordinate_1;

    typedef X_coordinate_1 Algebraic_real_1;
    typedef typename Algebraic_kernel_2::Algebraic_real_2 Algebraic_real_2;
    typedef typename Algebraic_kernel_2::Xy_coordinate_2 Xy_coordinate_2;

    typedef typename Rep::Boundary Boundary;
    typedef typename Rep::Lazy_boundary Lazy_boundary;
    typedef typename Rep::Intersection_info_container 
        Intersection_info_container;
    typedef typename Rep::Lazy_intersection_info_container 
        Lazy_intersection_info_container;
    typedef typename Rep::Index_triple Index_triple;
    typedef typename Curve_analysis_2::Integer Integer;
    typedef typename Curve_analysis_2::Status_line_1 Status_line_CA_1;

    typedef typename Curve_analysis_2::Coefficient Coefficient;
    
    typedef typename Curve_analysis_2::Polynomial_traits_2 Polynomial_traits_2;
  
    typedef typename Rep::Lazy_status_line_1 Lazy_status_line_1;

    typedef typename Curve_analysis_2::Coercion Coercion;
    
    typedef typename Curve_analysis_2::Coercion_type Coercion_type;

    typedef typename Curve_analysis_2::Poly_coer_1 Poly_coer_1;

    typedef typename Curve_analysis_2::Solve_1 Solve_1;

    typedef typename Rep::Slice_info Slice_info;
    typedef typename Rep::Lazy_slice_info Lazy_slice_info;

    //! The event slice object type
    typedef CGAL::CGALi::Status_line_CPA_1< Self > 
        Status_line_CPA_1;

    // Required by the concept. It is not used internally to distinguish
    // from one curve status_lines
    typedef Status_line_CPA_1 Status_line_1;

    struct Curves_at_event_functor {
        
        typedef size_type argument_type;
        typedef CGAL::CGALi::Slice_type result_type;

        Curves_at_event_functor(const Status_line_CPA_1& status_line) 
            : status_line(status_line)
        {}

        CGAL::CGALi::Slice_type operator() (size_type i) const {
            typedef typename Status_line_CPA_1::size_type 
                Status_line_size_type;
            std::pair<Status_line_size_type,Status_line_size_type> pair =
                status_line.curves_at_event(i);
            CGAL_assertion(pair.first>=0 || pair.second >=0);
            if(pair.first==-1) {
                return CGAL::CGALi::SECOND_CURVE;
            }
            if(pair.second==-1) {
                return CGAL::CGALi::FIRST_CURVE;
            }
            return CGAL::CGALi::INTERSECTION;
        }

    private:

        const Status_line_CPA_1& status_line;

    };

    typedef boost::transform_iterator<Curves_at_event_functor, 
                              boost::counting_iterator<size_type> > 
        Status_line_CPA_iterator;


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
     * If \c full_analysis is set, the curves are possibly sheared, if not
     * in generic position. Also, the complete Status_line_CPA_1 objects are
     * constructed. If not set, only internal pre-information is computed
     * which cannot be used from outside the object, or a 
     * Non_generic_position_exception is thrown.
     */
    Curve_pair_analysis_2(Curve_analysis_2 c1, Curve_analysis_2 c2) 
        throw(CGAL::CGALi::Zero_resultant_exception<Polynomial_2>,
              CGAL::CGALi::Non_generic_position_exception)
        : Base(Rep(c1, c2)) {

#if CGAL_ACK_DEBUG_FLAG
            CGAL::set_pretty_mode(CGAL_ACK_DEBUG_PRINT);
#endif
            
            
            this->ptr()->f = this->ptr()->c1_.f();
            this->ptr()->g = this->ptr()->c2_.f();
            
#if CGAL_ACK_RESULTANT_FIRST_STRATEGY
#ifndef CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD
        bool speed_up = true;
#else
        bool speed_up = std::min(curve_analysis(false).degree(),
                                 curve_analysis(true).degree()) >= 
            CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD;
#endif
#else
        bool speed_up=false;
#endif

#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Check content for squarefreeness.." 
                                 << std::flush;
#endif
            if(this->ptr()->c1_.content().degree()>0 &&
               this->ptr()->c2_.content().degree()>0) {
                if(CGAL::CGALi::gcd_utcf
                   (this->ptr()->c1_.content(), 
                    this->ptr()->c2_.content()).degree() >= 1) {

#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "Common vertical line discovered" 
                                         << std::endl;
#endif
                    throw CGAL::CGALi::Non_generic_position_exception();
                } else {
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
                }
            }

            if(speed_up) {
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Compute the resultant of f and g..." 
                                     << std::flush;
#endif
                this->ptr()->resultant 
                    = CGAL::CGALi::resultant(this->ptr()->f,this->ptr()->g);
            } else {
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Compute the subres-seq of f and g..." 
                                     << std::flush;
#endif
                compute_subresultant_coefficients();
                     
                this->ptr()->resultant 
                    = this->ptr()->principal_subresultants.get()[0];
            }


            if(this->ptr()->resultant.is_zero()) {
                throw CGAL::CGALi::Zero_resultant_exception<Polynomial_2>
                    (this->ptr()->f,
                     this->ptr()->g);
            }
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
            //CGAL_ACK_DEBUG_PRINT << "R=" << resultant << std::endl;
            CGAL_ACK_DEBUG_PRINT << "Isolate the real roots of resultant..." 
                                 << std::flush;
#endif
            Solve_1 solve_1;
            solve_1(this->ptr()->resultant, 
                    std::back_inserter(this->ptr()->resultant_roots),
                    std::back_inserter(this->ptr()->
                                       multiplicities_of_resultant_roots) );
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

            for(size_type i = 0;
                i<static_cast<size_type>(this->ptr()->resultant_roots.size());
                i++) {
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT 
                    << "Root at " 
                    << CGAL::to_double(this->ptr()->resultant_roots[i])
                    << " with multiplicity "
                    << this->ptr()->multiplicities_of_resultant_roots[i]
                    << std::endl;
#endif
            }
        
            compute_event_and_intermediate_values();

        }

private:

    // Refines the isolating intervals until they are disjoint
    // Returns true, if the first root is lower
    bool split_compare(Status_line_CA_1& e1, size_type i1,
                       Status_line_CA_1& e2, size_type i2) const {
        while(overlap(e1,i1,e2,i2)) {
            if(e1.interval_length(i1)<e2.interval_length(i2)) {
                e2.refine(i2);
            }
            else {
                e1.refine(i1);
            } 
        }
        return e1.lower_boundary(i1) < e2.lower_boundary(i2);
    }


    Status_line_CPA_1 create_event_slice(size_type i,bool use_shear=true) const {
#if !CGAL_ACK_NO_ARC_FLIP
        size_type index_in_fg = event_indices(i).fg;
        if(index_in_fg == -1 ) {
            return create_slice_with_multiplicity_zero_or_one(i);
        } else {
            size_type mult_of_alpha 
                = this->ptr()->multiplicities_of_resultant_roots[index_in_fg];
            if(mult_of_alpha == 1) {
                return create_slice_with_multiplicity_zero_or_one(i);
            } else {
#endif
                return create_slice_of_higher_multiplicity(i,use_shear);
#if !CGAL_ACK_NO_ARC_FLIP
            }
        }
#endif
    }

    Status_line_CPA_1 create_slice_of_higher_multiplicity(size_type i,
                                                     bool use_shear=true)  const {
        bool is_resultant_root = this->ptr()->index_triples[i].fg >=0;
        if(is_resultant_root &&
           this->ptr()->intersection_info_container) {
            return create_event_slice_with_shear(i);
        }
        try {
            Status_line_CPA_1 slice = construct_generic_case(i);

            return slice;
        } catch(CGAL::CGALi::Non_generic_position_exception ex) {

            if(! use_shear) {
                throw;
            } else {
                return create_event_slice_with_shear(i);
            }
        }

        // NEVER HAPPENS
        return Status_line_CPA_1();

    }


    Status_line_CPA_1 create_slice_with_multiplicity_zero_or_one(size_type i) const {

        const std::vector<X_coordinate_1>& events 
            = this->ptr()->event_x_coordinates;
        X_coordinate_1 alpha = events[i];
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
            CGAL_assertion(*left_it != CGAL::CGALi::INTERSECTION);
            if(*left_it == CGAL::CGALi::FIRST_CURVE) {
                CGAL_assertion(asym_lm_1!=0);
                asym_lm_1--;
            }
            if(*left_it == CGAL::CGALi::SECOND_CURVE) {
                CGAL_assertion(asym_lm_2!=0);
                asym_lm_2--;
            }
            left_it++;
        }
        while(asym_rm_1 != 0 || asym_rm_2 != 0) {
            if(*right_it == CGAL::CGALi::FIRST_CURVE) {
                CGAL_assertion(asym_rm_1!=0);
                asym_rm_1--;
            }
            if(*right_it == CGAL::CGALi::SECOND_CURVE) {
                CGAL_assertion(asym_rm_2!=0);
                asym_rm_2--;
            }
            right_it++;
        }
        while(asym_lp_1 != 0 || asym_lp_2 != 0) {
            left_end--;
            if(*left_end == CGAL::CGALi::FIRST_CURVE) {
                CGAL_assertion(asym_lp_1!=0);
                asym_lp_1--;
            }
            if(*left_end == CGAL::CGALi::SECOND_CURVE) {
                CGAL_assertion(asym_lp_2!=0);
                asym_lp_2--;
            }
        }
        while(asym_rp_1 != 0 || asym_rp_2 != 0) {
            right_end--;
            if(*right_end == CGAL::CGALi::FIRST_CURVE) {
                CGAL_assertion(asym_rp_1!=0);
                asym_rp_1--;
            }
            if(*right_end == CGAL::CGALi::SECOND_CURVE) {
                CGAL_assertion(asym_rp_2!=0);
                asym_rp_2--;
            }
        }
        // Now, the iterator ranges [left_it,left_end)
        // and [right_it,right_end) give the arcs really
        // going into the event line

        Slice_info slice_info;
        CGAL::CGALi::Slice_type curr_lowest_arc;
        size_type curr_multiplicity;

        size_type event_index_1=0, event_index_2=0;

        while(event_index_1 != e1.number_of_events() || 
              event_index_2 != e2.number_of_events()) {

            CGAL_assertion(event_index_1 != e1.number_of_events() || 
                           event_index_2 != e2.number_of_events());
            if(event_index_1==e1.number_of_events()) {
                curr_lowest_arc=CGAL::CGALi::SECOND_CURVE;
            } 
            else if(event_index_2==e2.number_of_events()) {
                curr_lowest_arc=CGAL::CGALi::FIRST_CURVE;
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
                bool e1_smaller = split_compare(e1,event_index_1,e2,event_index_2);
                curr_lowest_arc = e1_smaller ? CGAL::CGALi::FIRST_CURVE : CGAL::CGALi::SECOND_CURVE;
            }

            curr_multiplicity = -1;

            // Move the iterators
            size_type arcs_of_other_curve_left=0, arcs_of_other_curve_right=0;
            if(curr_lowest_arc==CGAL::CGALi::FIRST_CURVE) {
                size_type j=0;
                while(j<e1.number_of_incident_branches(event_index_1).first) {
                    if(*left_it==CGAL::CGALi::FIRST_CURVE) {
                        j++;
                    } else {
                        CGAL_assertion(event_index_2 < e2.number_of_events());
                        arcs_of_other_curve_left++;
                    }
                    left_it++;
                }

                j=0;
                while(j<e1.number_of_incident_branches(event_index_1).second) {
                    if(*right_it==CGAL::CGALi::FIRST_CURVE) {
                        j++;
                    } else {
                        CGAL_assertion(event_index_2 < e2.number_of_events());
                        arcs_of_other_curve_right++;
                    }
                    right_it++;
                }
                event_index_1++;
                if(arcs_of_other_curve_left+arcs_of_other_curve_right>0) {
                    // There was an intersection! Iterate over the remaining arcs
                    // on both sides belonging to this intersection
                    for(size_type j=arcs_of_other_curve_left;
                        j<e2.number_of_incident_branches(event_index_2).first;
                        j++) {
                        CGAL_assertion(*left_it==CGAL::CGALi::SECOND_CURVE);
                        left_it++;
                    }
                    for(size_type j=arcs_of_other_curve_right;
                        j<e2.number_of_incident_branches(event_index_2).second;
                        j++) {
                        CGAL_assertion(*right_it==CGAL::CGALi::SECOND_CURVE);
                        right_it++;
                    }
                    event_index_2++;
                    curr_lowest_arc=CGAL::CGALi::INTERSECTION;
                    curr_multiplicity=1;
                }
            } else { // curr_lowest_arc=CGAL::CGALi::SECOND_CURVE
                size_type j=0;
                while(j<e2.number_of_incident_branches(event_index_2).first) {
                    if(*left_it==CGAL::CGALi::SECOND_CURVE) {
                        j++;
                    } else {
                        CGAL_assertion(event_index_1 < e1.number_of_events());
                        arcs_of_other_curve_left++;
                    }
                    left_it++;
                }
              
                j=0;
                while(j<e2.number_of_incident_branches(event_index_2).second) {
                    if(*right_it==CGAL::CGALi::SECOND_CURVE) {
                        j++;
                    } else {
                        CGAL_assertion(event_index_1 < e1.number_of_events());
                        arcs_of_other_curve_right++;
                    }
                    right_it++;
                }
                event_index_2++;
                if(arcs_of_other_curve_left+arcs_of_other_curve_right>0) {
                    // There was an intersection! Iterate over the remaining arcs
                    // on both sides belonging to this intersection
                    for(size_type j=arcs_of_other_curve_left;
                        j<e1.number_of_incident_branches(event_index_1).first;
                        j++) {
                        CGAL_assertion(*left_it==CGAL::CGALi::FIRST_CURVE);
                        left_it++;
                    }
                    for(size_type j=arcs_of_other_curve_right;
                        j<e1.number_of_incident_branches(event_index_1).second;
                        j++) {
                        CGAL_assertion(*right_it==CGAL::CGALi::FIRST_CURVE);
                        right_it++;
                    }
                    event_index_1++;
                    curr_lowest_arc=CGAL::CGALi::INTERSECTION;
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
        

    // Helping functor
    struct Xval_of_status_line_CA_1 {
        typedef Status_line_CA_1 argument_type;
        typedef X_coordinate_1 result_type;
        X_coordinate_1 operator() (const Status_line_CA_1& status_line) const {
            return status_line.x();
        }
    };

    // Prepares the Status_line_CPA_1 objects for lazy evaluation
    void compute_event_and_intermediate_values() {
        
	Xval_of_status_line_CA_1 xval;
        Curve_analysis_2& c1=this->ptr()->c1_, c2=this->ptr()->c2_;
	
	size_type number_of_one_curve_events 
            = c1.num_events() + c2.num_events();

	std::vector<X_coordinate_1> one_curve_events;
	one_curve_events.reserve(number_of_one_curve_events);

        std::vector<CGAL::CGALi::Three_valued> one_curve_events_type;
	one_curve_events_type.reserve(number_of_one_curve_events);

        typename CGAL::Real_embeddable_traits<X_coordinate_1>::Compare compare;

	
        CGAL::CGALi::set_union_with_source
            (::boost::make_transform_iterator(c1.event_begin(),xval),
             ::boost::make_transform_iterator(c1.event_end(),xval),
             ::boost::make_transform_iterator(c2.event_begin(),xval),
             ::boost::make_transform_iterator(c2.event_end(),xval),
             std::back_inserter(one_curve_events),
             std::back_inserter(one_curve_events_type),
             compare);

	size_type number_of_events = number_of_one_curve_events +
            static_cast<size_type>(this->ptr()->resultant_roots.size());

        std::vector<X_coordinate_1>& events = this->ptr()->event_x_coordinates;
	events.reserve(number_of_events);
        events.clear();
        std::vector<CGAL::CGALi::Three_valued> events_type;
        CGAL::CGALi::set_union_with_source
            (one_curve_events.begin(),
             one_curve_events.end(),
             this->ptr()->resultant_roots.begin(),
             this->ptr()->resultant_roots.end(),
             std::back_inserter(events),
             std::back_inserter(events_type),
             compare);

        typename std::vector<CGAL::CGALi::Three_valued>::iterator one_curve_it
            =one_curve_events_type.begin();
        size_type inter_count=0, f_count=0,g_count=0;
        std::vector<Index_triple>& index_triples=this->ptr()->index_triples;

	this->ptr()->event_slices.reserve(number_of_events);
	this->ptr()->index_triples.reserve(number_of_events);

        for(size_type i=0;i<static_cast<size_type>(events.size());i++) {
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << CGAL::to_double(events[i]) << std::flush;
#endif
*/
            switch(events_type[i]) {
            case(CGAL::CGALi::ROOT_OF_FIRST_SET): {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << " one curve event" << std::endl;
#endif
*/
                this->ptr()->event_slices.push_back(Lazy_status_line_1());
                switch(*(one_curve_it++)) {
                case(CGAL::CGALi::ROOT_OF_FIRST_SET): {
                    index_triples.push_back(Index_triple(-1,f_count,-1));
                    f_count++;
                    break;
                }
                case(CGAL::CGALi::ROOT_OF_SECOND_SET): {
                    index_triples.push_back(Index_triple(-1,-1,g_count));
                    g_count++;
                    break;
                }
                case(CGAL::CGALi::ROOT_OF_BOTH_SETS): {
                    index_triples.push_back(Index_triple(-1,f_count,g_count));
                    f_count++;
                    g_count++;
                    break;
                }
                }
                break;
            }
            case(CGAL::CGALi::ROOT_OF_SECOND_SET): {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << " two curve event" << std::endl;
#endif
*/
                this->ptr()->
                    event_slices.push_back(Lazy_status_line_1());

                index_triples.push_back
                    (Index_triple(inter_count,-1,-1));
                inter_count++;
                break;
            }
            case(CGAL::CGALi::ROOT_OF_BOTH_SETS): {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << " one and two curve event" 
                                     << std::endl;
#endif
*/
                this->ptr()->
                    event_slices.push_back(Lazy_status_line_1());


                switch(*(one_curve_it++)) {
                case(CGAL::CGALi::ROOT_OF_FIRST_SET): {
                    index_triples.push_back
                        (Index_triple(inter_count,f_count,-1));
                    f_count++;
                    break;
                }
                case(CGAL::CGALi::ROOT_OF_SECOND_SET): {
                    index_triples.push_back
                        (Index_triple(inter_count,-1,g_count));
                    g_count++;
                    break;
                }
                case(CGAL::CGALi::ROOT_OF_BOTH_SETS): {
                    index_triples.push_back
                        (Index_triple(inter_count,f_count,g_count));
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
                           (this->ptr()->resultant_roots.size()));
        CGAL_assertion(one_curve_it==one_curve_events_type.end());
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "Prepare intermediate slices.." << std::flush;
#endif

	this->ptr()->intermediate_values.reserve(number_of_events+1);
	this->ptr()->intermediate_slices.reserve(number_of_events+1);

        for(size_type i=0;
            i<=static_cast<size_type>(events.size());
            i++) {
            this->ptr()->intermediate_values.push_back(Lazy_boundary());
            //boundary_value_in_interval(i);
            this->ptr()->intermediate_slices.push_back(Lazy_status_line_1());
        }
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    }

    // Creates an intermediate slice at a rational value
    Status_line_CPA_1 create_intermediate_slice_at(int i) const {
        
        Boundary r = boundary_value_in_interval(i);


        typedef typename CGAL::Fraction_traits<Poly_coer_1> FT;
        
        CGAL_assertion(static_cast<bool>((boost::is_same
                                          < typename FT::Numerator_type,
                                          Polynomial_1 >::value)));

        typename FT::Numerator_type p1,p2;
        typename FT::Denominator_type denom;
        
        Poly_coer_1 f1_at_r_with_denom 
            = typename Polynomial_traits_2::Swap() 
                (this->ptr()->c1_.f(), 0, 1).evaluate(r);

        Poly_coer_1 f2_at_r_with_denom 
            = typename Polynomial_traits_2::Swap() 
                (this->ptr()->c2_.f(), 0, 1).evaluate(r);

        typename FT::Decompose() (f1_at_r_with_denom, 
                                  p1, denom); 
        typename FT::Decompose() (f2_at_r_with_denom, 
                                  p2, denom); 
        
        Solve_1 solve_1;
        std::vector<X_coordinate_1> p1_roots,p2_roots;
        solve_1(p1,std::back_inserter(p1_roots));
        solve_1(p2,std::back_inserter(p2_roots));
        
	size_type number_of_roots 
            = static_cast<size_type>(p1_roots.size() + p2_roots.size());
        std::vector<X_coordinate_1> p12_roots;
	p12_roots.reserve(number_of_roots);
	std::vector<CGAL::CGALi::Three_valued> p12_order;
	p12_order.reserve(number_of_roots);

        CGAL::CGALi::Distinct_compare<X_coordinate_1> distinct_compare;
        set_union_with_source(p1_roots.begin(),
                              p1_roots.end(),
                              p2_roots.begin(),
                              p2_roots.end(),
                              std::back_inserter(p12_roots),
                              std::back_inserter(p12_order),
                              distinct_compare);

        Slice_info slice_info;
        
        for(typename std::vector<CGAL::CGALi::Three_valued>::const_iterator 
                it = p12_order.begin();
            it!=p12_order.end();
            it++) {
            switch(*it){
            case(CGAL::CGALi::ROOT_OF_FIRST_SET): {
                slice_info.push_back
                    (std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
                break;
            }
            case(CGAL::CGALi::ROOT_OF_SECOND_SET): {
                slice_info.push_back
                    (std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
                break;
            }
            case(CGAL::CGALi::ROOT_OF_BOTH_SETS): {  
                CGAL_assertion(false);
                break;
            }
            }
        }

        Status_line_CPA_1 new_slice 
            = create_slice_from_slice_info(i,slice_info,false);

        return new_slice;          
    }
    
    // Create a slice with id \c id from the Slice_info object
    Status_line_CPA_1 create_slice_from_slice_info(size_type id,
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
            CGAL_assertion(it->first != CGAL::CGALi::CANDIDATE);
            switch(it->first) {
            case(CGAL::CGALi::FIRST_CURVE): {
                if(event_flag) {
                    arc_container.push_back(std::make_pair(0,it->second));
                } else {
                    int_container.push_back(0);
                }
                break;
            }
            case(CGAL::CGALi::SECOND_CURVE): {
                if(event_flag) {
                    arc_container.push_back(std::make_pair(1,it->second));
                } else {
                    int_container.push_back(1);
                }
                break;
            }
            case(CGAL::CGALi::INTERSECTION): {
                CGAL_assertion(event_flag);
                arc_container.push_back(std::make_pair(2,it->second));
                break;
            }
            case(CGAL::CGALi::CANDIDATE): {
                CGAL_assertion(false);
                break;
            }
            }
        }
        
        return event_flag 
            ? Status_line_CPA_1(id,arc_container,*this)
            : Status_line_CPA_1(id,int_container,*this);  
    }

    // Computes the subresultant coefficients of \c f and \c g
    void compute_subresultant_coefficients() const {
        typedef std::vector<Polynomial_1> Polynomial_container;
        this->ptr()->principal_subresultants = Polynomial_container();
        this->ptr()->coprincipal_subresultants = Polynomial_container();
        const Polynomial_2& f = this->ptr()->f, g=this->ptr()->g;
        std::vector<Polynomial_2> subresultants;
        if(f.degree()<g.degree()) {
#if CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS 
            CGAL::CGALi::bezout_polynomial_subresultants
                (g,f,std::back_inserter(subresultants));
#else
            typename CGAL::Polynomial_traits_d<Polynomial_2>
                ::Polynomial_subresultants()(g,
                                             f,
                                         std::back_inserter(subresultants));
#endif
        } else {
#if CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS 
            CGAL::CGALi::bezout_polynomial_subresultants(f,
                                                 g,
                                                 std::back_inserter(subresultants));
#else
            typename CGAL::Polynomial_traits_d<Polynomial_2>
                ::Polynomial_subresultants()(f,
                                             g,
                                         std::back_inserter(subresultants)); 
#endif
        }

        size_type n = static_cast<size_type>(subresultants.size());

        for(size_type i=0;i<n;i++) {
            if(subresultants[i].degree() < i) {
                this->ptr()->principal_subresultants->
                    push_back(Polynomial_1(0));
            }
            else {
                this->ptr()->principal_subresultants->
                    push_back(subresultants[i][i]);
            }
        }
        for(size_type i=1;i<n;i++) {
            if(subresultants[i].degree() < i-1) {
                this->ptr()->coprincipal_subresultants->
                    push_back(Polynomial_1(0));
            }
            else {
                this->ptr()->coprincipal_subresultants->
                    push_back(subresultants[i][i-1]);
            }
        }
        // This must be corrected, if f and g have same degree:
        if(f.degree() == g.degree()) {
            if(n>=1) {
                this->ptr()->principal_subresultants.get()[n-1] 
                    = Polynomial_1(g.lcoeff());
            }
            if(n>=2) {
                this->ptr()->coprincipal_subresultants.get()[n-2] 
                    = Polynomial_1(g[g.degree()-1]);
            }
        }


    }

    // Computes a slice_info object at X_coordinate_1 \c alpha
    Slice_info construct_slice_info(X_coordinate_1 alpha) const
        throw(CGAL::CGALi::Non_generic_position_exception) {
        
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
                    (std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
                i2++;
                continue;
            }
            if(i2==n2) {
                slice_info.push_back
                    (std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
                i1++;
                continue;
            }
            if(match!=matchings.end() && 
               i1==match->first && 
               i2==match->second) {
                slice_info.push_back(std::make_pair(CGAL::CGALi::CANDIDATE,1));
                i1++;
                i2++;
                match++;
                continue;
            }
            CGAL_assertion(!overlap(e1,i1,e2,i2));
            if(e1.lower_boundary(i1) < e2.lower_boundary(i2)) {
                slice_info.push_back
                    (std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
                i1++;
                continue;
            } else {
                slice_info.push_back
                    (std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
                i2++;
                continue;
            }
        }
        CGAL_assertion(match==matchings.end());
        return slice_info;
    }
      
      
    Status_line_CPA_1 construct_generic_case(size_type i) const 
        throw(CGAL::CGALi::Non_generic_position_exception) {
        
        X_coordinate_1 alpha = this->ptr()->event_x_coordinates[i];
        
        Slice_info slice_info;

        size_type index_of_fg = this->ptr()->index_triples[i].fg;
        size_type index_of_ffy =this->ptr()->index_triples[i].ffy;
        size_type index_of_ggy =this->ptr()->index_triples[i].ggy;
        if(index_of_fg>=0) {
            if(alpha.is_root_of(this->ptr()->c1_.f().lcoeff())
               || alpha.is_root_of(this->ptr()->c2_.f().lcoeff())) {
                throw CGAL::CGALi::Non_generic_position_exception();
            }
            size_type k = -1; // not yet computed
            if(index_of_ffy==-1 && index_of_ggy==-1) {
                // this means, we need the multiplicity of the intersections
                if(! this->ptr()->principal_subresultants) {
                    compute_subresultant_coefficients();
                }
                if(alpha.is_root_of
                   (this->ptr()->principal_subresultants.get()[1])) {
                    // multiplicity cannot be determined, throw exception
                    throw CGAL::CGALi::Non_generic_position_exception();
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
                while(slice_it->first!=CGAL::CGALi::CANDIDATE) {
                    if(slice_it->first==CGAL::CGALi::FIRST_CURVE) {
                        i1++;
                    }
                    if(slice_it->first==CGAL::CGALi::SECOND_CURVE) {
                        i2++;
                    }
                    if(slice_it->first==CGAL::CGALi::INTERSECTION) {
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
      
    bool check_candidate_by_arc_pattern(size_type index,
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
        
        CGAL::CGALi::Slice_type curr;
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
       

    /*
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
                         size_type root_index) const {
        
        X_coordinate_1 xval = e1.x();
        CGAL_assertion(xval==e2.x());
        bool is_intersection=false;
        size_type index_in_fg = this->ptr()->index_triples[root_index].fg;
        size_type mult_of_resultant
            = this->ptr()->multiplicities_of_resultant_roots[index_in_fg];

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
                              std::make_pair(CGAL::CGALi::INTERSECTION,
                                             mult_of_intersection));
        }
        else {
            bool e1_smaller=split_compare(e1,i1,e2,i2);

            slice_it=slice_info.erase(slice_it);
            if(e1_smaller) {
                slice_it = slice_info.insert
                    (slice_it,std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
                slice_it++;
                slice_it = slice_info.insert
                    (slice_it,std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
            } else {
                slice_it = slice_info.insert
                    (slice_it,std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
                slice_it++;
                slice_it = slice_info.insert
                    (slice_it,std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
            }
        }
    }

    /*
     * Checks intersection with symbolic methods
     */
    bool check_candidate_symbolically(Status_line_CA_1& e1,size_type i1,
                                      Status_line_CA_1& e2,size_type i2,
                                      size_type k) const {
        Polynomial_1 p = -this->ptr()->coprincipal_subresultants.get()[k-1];
        Polynomial_1 q 
            = this->ptr()->principal_subresultants.get()[k]*Coefficient(k);
        X_coordinate_1 alpha = e1.x();
        CGAL_assertion(alpha==e2.x());
        if(zero_test_bivariate(alpha,this->ptr()->f,p,q) && 
           zero_test_bivariate(alpha,this->ptr()->g,p,q)) {
            return true;
        }
        else {
            throw CGAL::CGALi::Non_generic_position_exception();
        }
        return false; // never happens
    }

    /*
     * Checks whether the isolting intervals for the point on \c e1 with
     * index \c index1, and for the point on \c e2 with index \c index2
     * overlap
     */
    bool overlap(Status_line_CA_1& e1, size_type index1,Status_line_CA_1& e2, size_type index2) const {
        if(e1.lower_boundary(index1) > e2.upper_boundary(index2)) {
            return false;
        }
        else if(e1.upper_boundary(index1) < e2.lower_boundary(index2)) {
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
                                     Status_line_CA_1& e2) const {

        std::vector<size_type> possible_overlaps;
        for(size_type i=0;i<e2.number_of_events();i++) {
            if(overlap(e1,index1,e2,i)) {
                possible_overlaps.push_back(i);
            }
        }
        while(possible_overlaps.size()>1) {
            if(possible_overlaps.size()==2) {
                // Prevent that both intervals touch in a boundary
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


    size_type degree_of_local_gcd(size_type index_of_fg,
                            X_coordinate_1 alpha) const {
        
        if(this->ptr()->multiplicities_of_resultant_roots[index_of_fg] == 1) {
            return 1;
        } else {
            size_type k=1;
            if(! this->ptr()->principal_subresultants) {
                compute_subresultant_coefficients();

            }
            while(alpha.is_root_of(this->ptr()->principal_subresultants.get()[k])) {
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
        Index_triple& triple = this->ptr()->index_triples[i];
        return c ? triple.ggy : triple.ffy;
    }

    size_type event_of_curve_analysis(size_type i, 
                                      const Curve_analysis_2& c) const {
        CGAL_assertion(c.id()==curve_analysis(false).id() ||
                       c.id()==curve_analysis(true).id());
        Index_triple& triple = this->ptr()->index_triples[i];
        return (c.id()==curve_analysis(false).id()) ? triple.ffy : triple.ggy;
    }

    /*! 
     * \brief Returns the number of event slices
     *
     * Precisely, this is the number of points which are either root of
     * the resultant of the two curves, or root of discriminant of one
     * of the curves
     */
    size_type number_of_status_lines_with_event() const {
        return static_cast<size_type>(this->ptr()->event_x_coordinates.size());
    }

    //! For convenience
    size_type num_events() const {
        return number_of_status_lines_with_event();
    }

    //! Returns the x-coordinate of the <tt>i</tt>th event
    X_coordinate_1 event_x(size_type i) const {
        return this->ptr()->event_x_coordinates[i];
    }

    /*!
     * \brief The index of the x-coordinate
     *
     * For x-value \c x, the index of the suitable slice is computed. For
     * event value, the \c event flag is set to true, otherwise to false
     * and the slice of the interval to which \c x belongs is returned
     */
    void x_to_index(X_coordinate_1 x, 
                    size_type& idx, bool& event) const {
        const std::vector<X_coordinate_1>& sl 
            = this->ptr()->event_x_coordinates;
        idx = std::lower_bound(sl.begin(),
                               sl.end(),
                               x) - sl.begin();
        event = (idx < static_cast<size_type>(sl.size()) && (sl[idx] == x));

    }

    /*
     * \brief Returns the indices of the <tt>i</tt>th event value
     *
     * Returns a Index_triple <tt>(fg,ffy,ggy)</tt> such that
     * the <tt>i</tt>th event root is the <tt>fg</tt>th root of the 
     * resultant of \c f and \c g, the <tt>ffy</tt>th root of the 
     * discriminant of \c f, and  the <tt>ggy</tt>th root of the 
     * discriminant of \c g.
     */
    Index_triple event_indices(size_type i) const {
        return this->ptr()->index_triples[i];
    }

    Status_line_CPA_1 status_line_for_x(X_coordinate_1 x, 
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
        

    Status_line_CPA_1 status_line_at_exact_x(X_coordinate_1 x) {
        return status_line_for_x(x);
    }


private:
      
    const Status_line_CPA_1& _status_line_at_event(size_type i,
                                             bool use_shear) const {
        
        if(! this->ptr()->event_slices[i]) {
            this->ptr()->event_slices[i] = create_event_slice(i,use_shear);
        }
        
        return *(this->ptr()->event_slices[i]);
    }

public:
      
    //! Returns the Status_line_CPA_1 at the <tt>i</tt>th event
    const Status_line_CPA_1& status_line_at_event(size_type i) const {
        return this->_status_line_at_event(i,true);
    }
      


    //! Returns the Status_line_CPA_1 at the <tt>i</tt>th interval
    const Status_line_CPA_1& status_line_of_interval(size_type i) const {

        if(! this->ptr()->intermediate_slices[i]) {

            this->ptr()->intermediate_slices[i] 
                = create_intermediate_slice_at(i);
          
        }

        return *(this->ptr()->intermediate_slices[i]);
    }
        
    //!  Returns boundary representative value at the <tt>i</tt>th interval
    const Boundary boundary_value_in_interval(size_type i) const {

        const std::vector<X_coordinate_1>& events 
	    = this->ptr()->event_x_coordinates;
        /*	  
	  if(i==0) {
          this->ptr()->intermediate_values[i] = events[0].low() - 1;
	  } else if(i==static_cast<size_type>(events.size())) {
          this->ptr()->intermediate_values[i] = events[i-1].high() + 1;
	  } else {
          const Curve_analysis_2& c1 = this->ptr()->c1_, c2=this->ptr()->c2_;
          Polynomial_1 p1 = events[i-1].polynomial(),
          p2 = events[i].polynomial();
          Index_triple tr1 = this->ptr()->index_triples[i-1],
          tr2 = this->ptr()->index_triples[i];
          Boundary left1, right1, left2, right2;
          if(tr1.fg>=0) {
          left1=this->ptr()->stripe_values[tr1.fg];
          right1=this->ptr()->stripe_values[tr1.fg+1];
          } else if(tr1.ffy>=0) {
          left1=c1.boundary_value_in_interval(tr1.ffy);
          right1=c1.boundary_value_in_interval(tr1.ffy+1);
          } else {
          left1=c2.boundary_value_in_interval(tr1.ggy);
          right1=c2.boundary_value_in_interval(tr1.ggy+1);
          }
          if(tr2.fg>=0) {
          left2=this->ptr()->stripe_values[tr2.fg];
          right2=this->ptr()->stripe_values[tr2.fg+1];
          } else if(tr2.ffy>=0) {
          left2=c1.boundary_value_in_interval(tr2.ffy);
          right2=c1.boundary_value_in_interval(tr2.ffy+1);
          } else {
          left2=c2.boundary_value_in_interval(tr2.ggy);
          right2=c2.boundary_value_in_interval(tr2.ggy+1);
          }
          X_coordinate_1 x1(p1,left1,right1), x2(p2,left2,right2);
          this->ptr()->intermediate_values[i] = x1.rational_between(x2);
	  }
        */
	  
        if(! this->ptr()->intermediate_values[i]) {
            // Create the intermediate x-coordinate first
            if(events.size()==0) {
                CGAL_assertion(i==0);
                this->ptr()->intermediate_values[0]=Boundary(0);
            } else {
                if(i==0) {
                    this->ptr()->intermediate_values[i] 
                        = simple_rational_left_of(events[i]);
                } else if(i == static_cast<size_type>(events.size())) {
                    this->ptr()->intermediate_values[i]
                        = simple_rational_right_of(events[i-1]);

                } else {
                    this->ptr()->intermediate_values[i]
                        = simple_rational_between(events[i-1],events[i]);
                    //= events[i-1].rational_between(events[i]);
                }
            }
        }
     
        return *(this->ptr()->intermediate_values[i]);

    }
        
private:
      
    struct Boundary_to_coercion_functor {
        
        typedef Boundary argument_type;
        typedef Coercion_type result_type;

        result_type operator() (argument_type x) const {
            typename CGAL::Coercion_traits<Boundary,Coefficient>::Cast cast;
            return cast(x);
        }
    };

    struct Coefficient_to_coercion_functor {
        
        typedef Coefficient argument_type;
        typedef Coercion_type result_type;

        result_type operator() (argument_type x) const {
            typename CGAL::Coercion_traits<Boundary,Coefficient>::Cast cast;
            return cast(x);
        }
    };

    /* 
     * \brief  Symbolic zero test.
     *
     * Checks whether <tt>h(x,y(x))=0</tt>, where <tt>y(x)</tt> is a rational 
     * expression in terms of \c x, i.e. <tt>y=p/q</tt> with <tt>p,q</tt> 
     * univariate polynomials
     */
    bool zero_test_bivariate(const X_coordinate_1& alpha,
                             const Polynomial_2& h,
                             const Polynomial_1& p,
                             const Polynomial_1& q) const {

        // Compute h_0(x)=q(x)^n*h(x,p(x)/q(x))
      
        bool result;
#if !CGAL_ACK_USE_NO_REDUCTION_MODULO_RESULTANT
        bool general = ! alpha.is_rational();
      
        if(general) {

            Poly_coer_1 p_rat
                (boost::make_transform_iterator
                 (p.begin(),Coefficient_to_coercion_functor()),
                 boost::make_transform_iterator
                 (p.end(),Coefficient_to_coercion_functor())),
            q_rat
                (boost::make_transform_iterator
                 (q.begin(),Coefficient_to_coercion_functor()),
                 boost::make_transform_iterator
                 (q.end(),Coefficient_to_coercion_functor()));
            Poly_coer_1 modulus
                (boost::make_transform_iterator
                 (alpha.polynomial().begin(),Coefficient_to_coercion_functor()),
                 boost::make_transform_iterator
                 (alpha.polynomial().end(),Coefficient_to_coercion_functor()));

/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Mod: " << modulus << std::endl;
#endif 
*/
            p_rat=this->mod(p_rat,modulus);
            q_rat=this->mod(q_rat,modulus);

            size_type n = h.degree();
            // Create the powers of p and q mod modulus
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "precomp powers.." << std::flush;
#endif
*/
            std::vector<Poly_coer_1> p_powers(n+1),q_powers(n+1);
            p_powers[0]=Poly_coer_1(Boundary(1));
            q_powers[0]=Poly_coer_1(Boundary(1));
            Poly_coer_1 intermediate;
            for(size_type i=1;i<=n;i++) {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << i << ": mult.." << std::flush;
#endif
*/
                intermediate=p_powers[i-1]*p_rat;
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "mod.." << std::flush;
#endif
*/
                p_powers[i]=this->mod(intermediate,modulus);
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "simpl.." << std::flush;
#endif
*/
                p_powers[i].simplify_coefficients();
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "mult.." << std::flush;
#endif
*/
                intermediate=q_powers[i-1]*q_rat;
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "mod.." << std::flush;
#endif
*/
                q_powers[i]=this->mod(intermediate,modulus);
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "simpl.." << std::flush;
#endif
*/
                q_powers[i].simplify_coefficients();
            }
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done\ncomp rat pol.." << std::flush;
#endif
*/
	
            Poly_coer_1 curr_coeff,curr_fac;
            Poly_coer_1 h_0_rat(Coercion_type(0));
            for(size_type i=0;i<=n;i++) {
                Poly_coer_1 tmp_pol
                    (boost::make_transform_iterator
                     (h[i].begin(),Coefficient_to_coercion_functor()),
                     boost::make_transform_iterator
                     (h[i].end(),Coefficient_to_coercion_functor()));
                curr_fac=this->mod
                    (tmp_pol*p_powers[i]*q_powers[n-i], modulus);
                h_0_rat+=curr_fac;
            }
            typedef typename CGAL::Fraction_traits<Poly_coer_1> FT;
            
            CGAL_assertion
                (static_cast<bool>((boost::is_same
                                    < typename FT::Numerator_type,
                                      Polynomial_1 >::value)));

            typename FT::Numerator_type integralized_pol;
            typename FT::Denominator_type denom;

            typename FT::Decompose() (h_0_rat, integralized_pol, denom); 

            return CGAL::CGALi::is_root_of(alpha,integralized_pol);
        }
        else {
            typename Coercion::Cast cast;
            Coercion_type b = cast(alpha.rational()),
                p_b=p.evaluate(b),q_b=q.evaluate(b);
            size_type n = h.degree();
            Coercion_type eval(0);
            for(size_type i=0;i<=n;i++) {
                eval+=h[i].evaluate(b)*CGAL::ipower(p_b,i)*CGAL::ipower(q_b,n-i);
            }
            result=(CGAL::sign(eval)==CGAL::ZERO);
        }
#else
#warning Uses no reduction modulo resultant!
        Polynomial_1 h_0=h.evaluate_homogeneous(p,q);
        result=AcX::is_root_of(alpha,h_0);
#endif      
	
        return result;

    }

    Poly_coer_1 mod(Poly_coer_1 a,Poly_coer_1 b) const {
        Poly_coer_1 ret=CGAL::mod(a,b);
        return ret;
    }


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
        CGAL_assertion(n == static_cast<size_type>( new_info_container.size()));
        //iterate through the vector and update 
        // (-1 stands for "multiplicity unknown")
        for(size_type i=0;i<n;i++) {
            size_type m = old_info_container[i].size();
            CGAL_assertion(m == static_cast<size_type>(new_info_container[i].size()));
            for(size_type j=0;j<m;j++) {
                old_info_container[i][j].mult
                    = std::max(new_info_container[i][j].mult,
                               old_info_container[i][j].mult);
            }
          
        }
    }

    void new_shear_for_intersection_info
    (Intersection_info_container& info_container) const {
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Use shear for intersections.." << std::endl;
#endif
        bool good_direction_found=false;
        Self sh_pair;
        Integer s;

        while(! good_direction_found) {
            try {
                info_container.clear();
                info_container.resize(this->ptr()->resultant_roots.size());                
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
                sh_pair=Self(sh1,sh2);

#if CGAL_ACK_DEBUG_FLAG 
                CGAL_ACK_DEBUG_PRINT << "Shear back intersection points..." 
                                     << std::flush;
#endif
                for(size_type i=0;
                    i<static_cast<size_type>
                        (sh_pair.ptr()->event_x_coordinates.size());
                    i++) {
                    if(sh_pair.ptr()->index_triples[i].fg==-1) {
                        continue;
                    }
                    Status_line_CPA_1 slice 
                        = sh_pair._status_line_at_event(i,false);
                    Curves_at_event_functor functor(slice);
                    for(size_type j=0;j<slice.number_of_events();j++) {
                        if(functor(j) == CGAL::CGALi::INTERSECTION) {
                            this->update_intersection_info(info_container,
                                                           sh_pair,
                                                           slice,
                                                           i,j,s);
                        }
                    }
                }
                good_direction_found=true;
            }
            catch(CGAL::CGALi::Non_generic_position_exception ex) {
                this->ptr()->shear_controller.report_failure(s);
            }
        }

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
        return;
    }


    Status_line_CPA_1 create_event_slice_with_shear(size_type i) const {
        while(true) { // we know that it works at some point
            try {
                if(! this->ptr()->intersection_info_container) {
                    // We shear for the first time, so compute the stripe values now
                    find_intermediate_values
                        (this->ptr()->resultant_roots.begin(),
                         this->ptr()->resultant_roots.end(),
                         std::back_inserter(this->ptr()->stripe_values));
                    Intersection_info_container info_container;
                    new_shear_for_intersection_info(info_container);
                    merge_new_intersection_info(info_container);
                }		  
                Status_line_CPA_1 slice 
                    = create_event_slice_from_current_intersection_info(i);

                return slice;
            } catch(CGAL::CGALi::Non_generic_position_exception ex) {
                // just try the next one
                Intersection_info_container info_container;
                new_shear_for_intersection_info(info_container);
                merge_new_intersection_info(info_container);
            }
        }
    }

      
    Status_line_CPA_1 
        create_event_slice_from_current_intersection_info (size_type i) 
        const throw(CGAL::CGALi::Non_generic_position_exception){
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Reduce the candidates.." << std::flush;
#endif
        Index_triple triple = this->ptr()->index_triples[i];
        size_type index_of_fg = triple.fg;
        Intersection_info_container& intersection_info_container
            = *(this->ptr()->intersection_info_container);
        X_coordinate_1 alpha = this->ptr()->event_x_coordinates[i];
        CGAL_assertion(index_of_fg <
                       static_cast<size_type>
                           (intersection_info_container.size()));
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << i << "th slice has " 
                             << intersection_info_container[index_of_fg].size()
                             << " intersections" << std::endl;
#endif
        Status_line_CA_1 e1=this->ptr()->c1_.
            status_line_at_exact_x(this->ptr()->resultant_roots[index_of_fg]),
            e2=this->ptr()->c2_.
            status_line_at_exact_x(this->ptr()->resultant_roots[index_of_fg]);
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
            if(slice[j].first==CGAL::CGALi::INTERSECTION) {
                inter_info_it++;
            }
            if(slice[j].first==CGAL::CGALi::CANDIDATE) {
                slice[j].first=CGAL::CGALi::INTERSECTION;
                if(triple.ffy==-1 && triple.ggy==-1 && inter_info_it->mult==-1) {
                    // Multiplicity unknown for case where we need it
                    throw CGAL::CGALi::Non_generic_position_exception();
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

    // THINK COPY AND PASTE FROM Shear_transformation.h!
    Boundary x_sheared(Boundary x, Boundary y,Integer sh) const {
        return x-sh*y;
    }

    void update_intersection_info(Intersection_info_container& 
                                  info_container,
                                  Self& sh_pair,
                                  Status_line_CPA_1 slice,
                                  size_type i,
                                  size_type j,
                                  Integer s) const {
        typedef typename Rep::Intersection_info Intersection_info;
        X_coordinate_1& xval = sh_pair.ptr()->event_x_coordinates[i];
        CGAL_assertion(Curves_at_event_functor(slice)(j)
                       ==CGAL::CGALi::INTERSECTION);
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
        const std::vector<Boundary>& stripe_values=this->ptr()->stripe_values;
        size_type left_index = -1, 
            right_index = static_cast<size_type>(stripe_values.size()-1);
        X_coordinate_1 xv = ev.x();
        Boundary lx = xv.low(), rx=xv.high(),
            x_iv_size = rx-lx;
        Boundary ly = ev.lower_boundary(index),
            ry = ev.upper_boundary(index);;
        while(left_index < right_index) {
            if(x_iv_size > ry-ly) {
                xv.refine();
                lx = xv.low();
                rx=xv.high();
                x_iv_size=rx-lx;
                continue;
            }
            ev.refine(index);
            ly = ev.lower_boundary(index);
            ry = ev.upper_boundary(index);
            Boundary left=(s<0) ? x_sheared(lx,ry,-s): x_sheared(lx,ly,-s);
            Boundary right = (s<0) ? x_sheared(rx,ly,-s) : x_sheared(rx,ry,-s);
            CGAL_assertion(left<right);
            while(left_index<right_index && stripe_values[left_index+1]<left) {
                ++left_index;
            }
            while(left_index<right_index && right<stripe_values[right_index]) {
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
            bool index_smaller = split_compare(ev,index,comp_ev,comp_index);
            if(!index_smaller) {
                info_it++;
            } else {
                break;
            }
        }
        info_vec.insert(info_it,intersection_info);        
    }
      
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
    size_type reduce_number_of_candidates_and_intersections_to(size_type n,
                                                         Status_line_CA_1& e1,
                                                         Status_line_CA_1& e2,
                                                         Slice_info& slice,
                                                         size_type k=-1) const {
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
            if(slice[i].first==CGAL::CGALi::CANDIDATE) {
                number_of_candidates++;
            }
            if(slice[i].first==CGAL::CGALi::INTERSECTION) {
                number_of_intersections++;
            }
        }
        CGAL_assertion(number_of_intersections<=n);
  
        typename Slice_info::iterator slice_it=slice.begin();
        size_type i1=0,i2=0;
        size_type max_candidate_mult=0;
        while(n<number_of_candidates+number_of_intersections) {
            if(slice_it==slice.end()) {
                CGAL_assertion(e1.number_of_events()==i1 && e2.number_of_events()==i2);
                if(max_candidate_mult<k) {
                    throw CGAL::CGALi::Non_generic_position_exception();
                } else {
                    slice_it=slice.begin();
                    max_candidate_mult=0;
                    i1=i2=0;
                }
            }
            switch(slice_it->first) {
            case(CGAL::CGALi::FIRST_CURVE): {
                i1++;
                break;
            }
            case(CGAL::CGALi::SECOND_CURVE): {
                i2++;
                break;
            }
            case(CGAL::CGALi::CANDIDATE): {
                if(e1.interval_length(i1)<e2.interval_length(i2)) {
                    e2.refine(i2);
                }
                else {
                    e1.refine(i1);
                }
                if(! overlap(e1,i1,e2,i2)) {
                    number_of_candidates--;
                    slice_it=slice.erase(slice_it);
                    if(e1.lower_boundary(i1)<e2.lower_boundary(i2)) {
                        slice_it=slice.
                            insert(slice_it,std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
                        slice_it++;
                        slice_it=slice.
                            insert(slice_it,std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
                    } else {
                        slice_it=slice.
                            insert(slice_it,std::make_pair(CGAL::CGALi::SECOND_CURVE,-1));
                        slice_it++;
                        slice_it=slice.
                            insert(slice_it,std::make_pair(CGAL::CGALi::FIRST_CURVE,-1));
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
            case(CGAL::CGALi::INTERSECTION): {
                i1++;
                i2++;
                break;
            }
            }
            slice_it++;
        }
        return number_of_intersections+number_of_candidates;
    }      
        

    void debug_print() const {
        for(size_type i=0;i<static_cast<size_type>(this->ptr()->slice_info_container.size());i++) {
            for(size_type j=0;j<static_cast<size_type>(this->ptr()->slice_info_container[i].get().size());j++) {
                switch(this->ptr()->slice_info_container[i].get()[j].first) {
                case(CGAL::CGALi::FIRST_CURVE): {
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "1" << std::flush;
#endif
                    break;
                }
                case(CGAL::CGALi::SECOND_CURVE): {
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "2" << std::flush;
#endif
                    break;
                }
                case(CGAL::CGALi::INTERSECTION): {
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT 
                        << "(S," 
                        << *(this->ptr()->slice_info_container[i])[j].second 
                        << ")" << std::flush;
#endif
                    break;
                }
                case(CGAL::CGALi::CANDIDATE): {
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "C" << std::flush;
#endif
                    break;
                }
                }
            }
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << std::endl;
#endif
        }
    }

    // Handle provides
    // .id()
    // .is_identical
};

//! \brief Prints the objects.
template<typename AlgebraicKernel_2>
std::ostream& operator<< 
    (std::ostream& out, 
     const Curve_pair_analysis_2<AlgebraicKernel_2>& curve_pair) {
    typedef Curve_pair_analysis_2<AlgebraicKernel_2> Curve_pair_analysis_2;
    typedef typename Curve_pair_analysis_2::size_type size_type;
    typedef typename Curve_pair_analysis_2::Index_triple Index_triple;
    typedef typename Curve_pair_analysis_2::Status_line_CPA_1 Slice;
    out << "--------------- Analysis results ---------------" << std::endl;
    out << "Number of constructed event lines: " << curve_pair.num_events() 
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
    for(size_type j = 0; j<curve_pair.num_events();j++) {

        out << "Event line at " << CGAL::to_double(curve_pair.event_x(j))
            << ": " << std::endl;
        out << "Indices: ";
        Index_triple triple = curve_pair.event_indices(j);
        out << "fg: " << triple.fg << ", ffy: " 
            << triple.ffy <<", ggy: " << triple.ggy 
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


CGAL_END_NAMESPACE

#endif
