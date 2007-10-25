// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_ALGEBRAIC_REAL_TRAITS_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_ALGEBRAIC_REAL_TRAITS_H

/*! \file Algebraic_real_traits.h
 *  \brief defines a class \c Algebraic_real_traits
 *  
 *  This class implements a model of \c AlgebraicRealTraits_1 to handle 
 *  \c AlgebraicCurveKernel_2::Algebraic_real_2 separately depending on 
 *  underlying curve pair analysis. 
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

#include <AcX/Algebraic_curve_pair_2.h>

#define CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS \
    typedef typename Algebraic_curve_pair_2::Algebraic_curve_2 Curve_2; \
    typedef typename Algebraic_real_2::X_coordinate_1 X_coordinate_1; \
    typedef Algebraic_real_2 Type; \
    typedef typename Curve_2::Boundary Boundary;

CGAL_BEGIN_NAMESPACE
    
namespace CGALi {

// one needs AlgebraicReal_2 and AlgebraicCurvePair_2 in order to
// instantiate Algebraic_real_traits

//! algebraic real traits template
template <class AlgebraicReal_2, class AlgebraicCurvePair_2>
struct Algebraic_real_traits {
    
    //! this instance's first template argument
    typedef AlgebraicReal_2 Algebraic_real_2;  
    
    //! this instance's second template argument
    typedef AlgebraicCurvePair_2 Algebraic_curve_pair_2; 

    CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS
    
    // instance of algebraic curve kernel
    typedef Null_functor Boundary_between;
    typedef Null_functor Lower_boundary;
    typedef Null_functor Upper_boundary;
    typedef Null_functor Refine;
};

//! specialization for AlciX     
template <class AlgebraicReal_2, class Curve_>
struct Algebraic_real_traits<AlgebraicReal_2, 
    AcX::Algebraic_curve_pair_2<Curve_> > {
    
    //! this instance's first template argument
    typedef AlgebraicReal_2 Algebraic_real_2;  
    
    //! this instance's second template argument
    typedef AcX::Algebraic_curve_pair_2<Curve_> Algebraic_curve_pair_2; 

    CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS
    
    //! type of curve vertical line
    typedef typename Curve_2::Curve_vertical_line Event_line;

    //! computes boundary in between y-coordinates of two algebraic reals
    //! defined over the same vertical line
    //!
    //! \pre r1.arcno() != r2.arcno() 
    //! \pre r1.x() == r2.x() and r1.curve() == r2.curve()
    struct Boundary_between 
            : public Binary_function< Type, Type, Boundary > {
        
        Boundary operator()(const Type& r1, const Type& r2) const {
            
            CGAL_precondition(r1.arcno() != r2.arcno());
            CGAL_precondition(r1.x() == r2.x() && r1.curve() == r2.curve());
            
            Boundary res;
            Event_line vline = r1.curve().event_info_at_x(r1.x());
            if(r1.arcno() < r2.arcno()) {
                res = (vline.upper_boundary(r1.arcno()) +
                    vline.lower_boundary(r2.arcno())) / Boundary(2);
            } else {
                res = (vline.lower_boundary(r1.arcno()) +
                    vline.upper_boundary(r2.arcno())) / Boundary(2);
            }
            CGAL::simplify(res);
            return res;
        }
    };
    
    //! returns current lower boundary of an isolating interval defining 
    //! y-coordinate of an algebraic real
    struct Lower_boundary
            : public Unary_function<Type, Boundary> {
        
        Boundary operator()(const Type& r) const {
            Event_line vline = r.curve().event_info_at_x(r.x());
            return vline.lower_boundary(r.arcno());
        }
    };

    //! returns current upper boundary of an isolating interval defining 
    //! y-coordinate of an algebraic real
    struct Upper_boundary
            : public Unary_function<Type, Boundary> {
         
        Boundary operator()(const Type& r) const {
            Event_line vline = r.curve().event_info_at_x(r.x());
            return vline.upper_boundary(r.arcno());
        }
    };
                
    struct Refine
            : public Unary_function<Type, void> {
        
        //! \brief refines isolating interval of an y-coorinate of algebraic 
        //! real to make it at least half of the original interval
        //!
        //! note that an interval may also degenerate to a single point
        void operator()(const Type& r) const {
        
            r.curve().event_info_at_x(r.x()).refine(r.arcno());
        }
        
        //! \brief refines isolating interval of an y-coorinate of algebraic 
        //! real w.r.t. the given relative precision
        //!
        //! the resulting interval is: 
        //! <tt>|lower - upper|/|r.y()| <= 2^(-rel_prec)</tt> 
        void operator()(Type& r, int rel_prec) const {
            
            Event_line vline = r.curve().event_info_at_x(r.x());
            Boundary prec = (vline.interval_length(r.arcno())) / 
                CGAL::POLYNOMIAL::ipower(Boundary(2), rel_prec);
                
            vline.refine_to(r.arcno(), prec);
        }
    };                
}; // class Algebraic_real_traits<AcX>

} // namespace CGALi 
            
CGAL_END_NAMESPACE
            
#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_ALGEBRAIC_REAL_TRAITS_H
