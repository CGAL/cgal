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

#if !CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE
#include <AcX/Algebraic_curve_pair_2.h>
#endif

#define CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS \
    typedef typename Algebraic_curve_pair_2::Algebraic_curve_2 Curve_2; \
    typedef Algebraic_real_2 Type; \
    typedef typename Curve_2::Boundary Boundary; \
    typedef typename Curve_2::Coefficient Coefficient;

CGAL_BEGIN_NAMESPACE
    
namespace CGALi {

//! algebraic real traits template
template <class AlgebraicReal_2, class AlgebraicCurvePair_2>
struct Algebraic_real_traits_for_y {
    
    //! this instance's first template argument
    typedef AlgebraicReal_2 Algebraic_real_2;
    //! this instance's second template argument
    typedef AlgebraicCurvePair_2 Algebraic_curve_pair_2;
    
    CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS
    
    typedef Null_functor Boundary_between;
    typedef Null_functor Lower_boundary;
    typedef Null_functor Upper_boundary;
    typedef Null_functor Refine;
};

//! algebraic real traits template
template <class AlgebraicReal>
struct Algebraic_real_traits {
    
    //! this instance's first template argument
    typedef AlgebraicReal Algebraic_real;
    
    //CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS
    
    typedef Null_functor Boundary_between;
    typedef Null_functor Lower_boundary;
    typedef Null_functor Upper_boundary;
    typedef Null_functor Refine;
};

#if !CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE

template <class AlgebraicCurveKernel_2, class Curve_>
struct Algebraic_real_traits_for_y<Xy_coordinate_2<
    AlgebraicCurveKernel_2>, AcX::Algebraic_curve_pair_2<Curve_> > {
    
     //! this instance's first template argument
    typedef Xy_coordinate_2<AlgebraicCurveKernel_2> Algebraic_real_2;
    
    //! this instance's second template argument
    typedef AcX::Algebraic_curve_pair_2<Curve_> Algebraic_curve_pair_2; 

    CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS

    //! myself
    typedef Algebraic_real_traits_for_y<Algebraic_real_2,
        Algebraic_curve_pair_2> Self;
    
    //! type of curve vertical line
    typedef typename Curve_2::Status_line Event_line;

    //! computes boundary between y-coordinates of two algebraic reals
    //! defined over the same vertical line
    //!
    //! \pre r1.arcno() != r2.arcno() 
    //! \pre r1.x() == r2.x() and r1.curve() == r2.curve()
    struct Boundary_between 
            : public Binary_function< Type, Type, Boundary > {
        
        Boundary operator()(const Type& r1, const Type& r2) const {
            
            CGAL_precondition(r1.arcno() != r2.arcno());
            CGAL_precondition(r1.x() == r2.x() &&
                r1.curve().id() == r2.curve().id());
            
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
        //! resulting interval is:
        //! <tt>|lower - upper|/|r.y()| <= 2^(-rel_prec)</tt> 
        void operator()(const Type& r, int rel_prec) const {
            
            Event_line vline = r.curve().event_info_at_x(r.x());
            int arcno = r.arcno();
            Boundary prec = (vline.interval_length(arcno)) /
                CGAL::POLYNOMIAL::ipower(Boundary(2), rel_prec);

            /////////// attention!! need to test for exact zero !!
               
            // Refine until both boundaries have the same sign
            while(CGAL::sign(vline.lower_boundary(arcno)) !=
                    CGAL::sign(vline.upper_boundary(arcno))) {
                vline.refine(arcno);
            }
            
            CGAL_assertion(
                CGAL::sign(vline.lower_boundary(arcno)) != CGAL::ZERO &&
                CGAL::sign(vline.upper_boundary(arcno)) != CGAL::ZERO);

            // Refine until precision is reached
            while((vline.upper_boundary(arcno) - vline.lower_boundary(arcno)) /
                   CGAL::max(CGAL::abs(vline.upper_boundary(arcno)),
                        CGAL::abs(vline.lower_boundary(arcno))) > prec ) {
                vline.refine(arcno);
            }
        }
    };
};

#endif // CGAL_ACK_2_NO_ALG_REAL_TRAITS_FOR_XY_COORDINATE
    
template <class Kernel_2>
struct Algebraic_real_traits<Xy_coordinate_2<Kernel_2> > :
    public Algebraic_real_traits_for_y< Xy_coordinate_2<Kernel_2>,
        typename Kernel_2::Curve_pair_2 > {
};
  
template <class Coefficient_, class FieldWithSqrt, class Rational_,
          class HandlePolicy, class AlgebraicRealRep >
struct Algebraic_real_traits<NiX::Algebraic_real<Coefficient_, FieldWithSqrt,
    Rational_, HandlePolicy, AlgebraicRealRep> > {

    //! this instances first template argument
    typedef NiX::Algebraic_real<Coefficient_, FieldWithSqrt, Rational_,
        HandlePolicy, AlgebraicRealRep> Algebraic_real_1;

    //! just a Type ?
    typedef Algebraic_real_1 Type;

    //! boundary type
    typedef typename Algebraic_real_1::Rational Boundary;

    //! computes rational boundary between two algebraic reals
    struct Boundary_between 
        : public Binary_function< Type, Type, Boundary > {
        
        Boundary operator()(const Type& t1, const Type& t2 ) const {
            return t1.rational_between(t2);
        }
    };

    //! returns current lower boundary of an algebraic real
    struct Lower_boundary
        : public Unary_function< Type, Boundary > {
                    
        Boundary operator()(const Type& t) const {
            return t.low();
        }
    };

    //! returns current upper boundary of an algebraic real
    struct Upper_boundary
        : public Unary_function< Type, Boundary > {
        
        Boundary operator()(const Type& t) const {
            return t.high();
        }
    };

    struct Refine
        : public Unary_function< Type, void > {

        //! \brief refines isolating interval defining algebraic real to make
        //! it at least half of the original interval
        //!
        //! note that an interval may also degenerate to a single point
        void operator()(const Type& t) const {
            t.refine();
        }

        //! \brief refines isolating interval of an algebraic real w.r.t. the
        //! given relative precision
        //!
        //! resulting interval is:
        //! <tt>|lower - upper|/|t| <= 2^(-rel_prec)</tt>
        void operator()(Type& t, int rel_prec) const {
            
            if(CGAL::is_zero(t)) {
                t = Type(0);
                return;
            } 
            Boundary len = t.high() - t.low(), prec = len /
                CGAL::POLYNOMIAL::ipower(Boundary(2), rel_prec);
            while(len > prec) {
                t.refine();
                len = t.high() - t.low();
            }
        }
    };
};

} // namespace CGALi 
            
CGAL_END_NAMESPACE
            
#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_ALGEBRAIC_REAL_TRAITS_H
