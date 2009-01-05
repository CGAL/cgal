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

#ifndef CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
#define CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION 0
#endif

#if CGAL_ACK_USE_EXACUS
#if !CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
#include <AcX/Algebraic_curve_pair_2.h>
#endif
#endif

#define CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS \
    typedef typename Curve_pair_2::Algebraic_curve_2 Curve_2; \
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
    typedef AlgebraicCurvePair_2 Curve_pair_2;
    
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
    
    typedef Null_functor Boundary_between;
    typedef Null_functor Lower_boundary;
    typedef Null_functor Upper_boundary;
    typedef Null_functor Refine;
};

template <class AlgebraicCurveKernel_2, class CurvePair_2>
struct Algebraic_real_traits_for_y<Xy_coordinate_2<
    AlgebraicCurveKernel_2>, CurvePair_2 > {
    
     //! this instance's first template argument
    typedef Xy_coordinate_2<AlgebraicCurveKernel_2> Algebraic_real_2;
    
    //! this instance's second template argument
    typedef CurvePair_2 Curve_pair_2; 

    CGAL_SNAP_ALGEBRAIC_REAL_TRAITS_2_TYPEDEFS

    //! myself
    typedef Algebraic_real_traits_for_y< Algebraic_real_2, Curve_pair_2 > Self;
    

#if !CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
    //! type of curve vertical line
    typedef typename Curve_2::Status_line Event_line;
#endif

    //! computes boundary between y-coordinates of two algebraic reals
    //! defined over the same vertical line
    struct Boundary_between 
            : public std::binary_function< Type, Type, Boundary > {
        
        Boundary operator()(const Type& r1, const Type& r2) const {

            CGAL_precondition(r1.y() != r2.y());

            Boundary res(0);
#if !CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
            Event_line vline1 =
                r1.curve()._internal_curve().event_info_at_x(r1.x());

            Event_line vline2 =
                r2.curve()._internal_curve().event_info_at_x(r2.x());
#endif            

            
            Boundary low1, low2, high1, high2;
            
            while (true) {
#if CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
                low1 = r1.y().low();
                high1 = r1.y().high();

                low2 = r2.y().low();
                high2 = r2.y().high();
#else
                low1 = vline1.lower_boundary(r1.arcno());
                high1 = vline1.upper_boundary(r1.arcno());
                
                low2 = vline2.lower_boundary(r2.arcno());
                high2 = vline2.upper_boundary(r2.arcno());
#endif
                
                if (low1 > high2) {
                    res = ((low1 + high2)/Boundary(2));
                    break;
                }
                if (low2 > high1) {
                    res = ((low2 + high1)/Boundary(2));
                    break;
                }
                
                // else
#if CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
                r1.y().refine();
                r2.y().refine();
#else
                vline1.refine(r1.arcno());
                vline2.refine(r2.arcno());
#endif
            }

            CGAL::simplify(res);

            CGAL_postcondition_code(
                    CGAL::Comparison_result exp = CGAL::SMALLER
            );
            CGAL_postcondition_code(
                    if (r1.y() > r2.y()) {
                        exp = CGAL::LARGER;
                    }
            );
            CGAL_postcondition(r1.y().compare(res) == exp);
            CGAL_postcondition(r2.y().compare(res) == -exp);

            return res;
        }
    };
    
    //! returns current lower boundary of an isolating interval defining 
    //! y-coordinate of an algebraic real
    struct Lower_boundary
            : public std::unary_function<Type, Boundary> {
        
        Boundary operator()(const Type& r) const {
#if CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
            return r.y().low();
#else
            Event_line vline = r.curve()._internal_curve().
                event_info_at_x(r.x());
            return vline.lower_boundary(r.arcno());
#endif
        }
    };

    //! returns current upper boundary of an isolating interval defining 
    //! y-coordinate of an algebraic real
    struct Upper_boundary
            : public std::unary_function<Type, Boundary> {
         
        Boundary operator()(const Type& r) const {
#if CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
            return r.y().high();
#else
            Event_line vline = r.curve()._internal_curve().
                event_info_at_x(r.x());
            return vline.upper_boundary(r.arcno());
#endif
        }
    };
                
    struct Refine
            : public std::unary_function<Type, void> {
        
        //! \brief refines isolating interval of an y-coorinate of algebraic 
        //! real to make it at least half of the original interval
        //!
        //! note that an interval may also degenerate to a single point
        void operator()(const Type& r) const {
#if CGAL_ACK_2_USE_EXPENSIVE_Y_MEMBER_FOR_APPROXIMATION
            return r.y().refine();
#else
            r.curve()._internal_curve().event_info_at_x
                (r.x()).refine(r.arcno());
#endif
        }
        
        //! \brief refines isolating interval of an y-coorinate of algebraic 
        //! real w.r.t. the given relative precision
        //!
        //! resulting interval is:
        //! <tt>|lower - upper|/|r.y()| <= 2^(-rel_prec)</tt> 
        void operator()(const Type& r, int rel_prec) const {
            Boundary low = Lower_boundary()(r);
            Boundary high = Upper_boundary()(r);
            
            Boundary prec = (high - low) /
                CGAL::ipower(Boundary(2), rel_prec);
            
            /////////// attention!! need to test for exact zero !!
               
            // Refine until both boundaries have the same sign
            while (CGAL::sign(low) != CGAL::sign(high)) {
                Refine()(r);
                low = Lower_boundary()(r);
                high = Upper_boundary()(r);
            }
            
            CGAL_assertion(
                    CGAL::sign(low) != CGAL::ZERO &&
                    CGAL::sign(high) != CGAL::ZERO);
            
            // Refine until precision is reached
            while (((high - low) / CGAL::max(CGAL::abs(high), CGAL::abs(low)))
                   > prec ) {
                Refine()(r);
                low = Lower_boundary()(r);
                high = Upper_boundary()(r);
            }
        }
    };
};

template <class Kernel_2>
struct Algebraic_real_traits<Xy_coordinate_2<Kernel_2> > :
    public Algebraic_real_traits_for_y< Xy_coordinate_2<Kernel_2>,
        typename Kernel_2::Internal_curve_pair_2 > {
};
  
template <class Coefficient_, class Rational_,
          class HandlePolicy, class AlgebraicRealRep >
struct Algebraic_real_traits<CGAL::CGALi::Algebraic_real_pure
    < Coefficient_, Rational_, HandlePolicy, AlgebraicRealRep > > {

    //! this instances first template argument
    typedef CGAL::CGALi::Algebraic_real_pure<Coefficient_, Rational_,
        HandlePolicy, AlgebraicRealRep> Algebraic_real_1;

    //! just a Type ?
    typedef Algebraic_real_1 Type;

    //! boundary type
    typedef typename Algebraic_real_1::Rational Boundary;

    //! computes rational boundary between two algebraic reals
    struct Boundary_between 
        : public std::binary_function< Type, Type, Boundary > {
        
        Boundary operator()(const Type& t1, const Type& t2 ) const {
            return t1.rational_between(t2);
        }
    };

    //! returns current lower boundary of an algebraic real
    struct Lower_boundary
        : public std::unary_function< Type, Boundary > {
                    
        Boundary operator()(const Type& t) const {
            return t.low();
        }
    };

    //! returns current upper boundary of an algebraic real
    struct Upper_boundary
        : public std::unary_function< Type, Boundary > {
        
        Boundary operator()(const Type& t) const {
            return t.high();
        }
    };

    struct Refine
        : public std::unary_function< Type, void > {

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
                CGAL::ipower(Boundary(2), rel_prec);
            while(len > prec) {
                t.refine();
                len = t.high() - t.low();
            }
        }
    };
};

#if !CGAL_ACK_USE_EXACUS
// !!!! TODO: That is preliminary!


template <class AlgebraicKernel_2>
struct Algebraic_real_traits_for_y<Xy_coordinate_2<
    AlgebraicKernel_2>, CGAL::Null_functor > {
    
    typedef AlgebraicKernel_2 Algebraic_kernel_2;

     //! this instance's first template argument
    typedef Xy_coordinate_2<Algebraic_kernel_2> Algebraic_real_2;
    
    //! myself
    typedef Algebraic_real_traits_for_y< Algebraic_real_2, CGAL::Null_functor > Self;
    
    //! For convenience
    typedef Algebraic_real_2 Type;

    typedef typename Algebraic_real_2::Boundary Boundary;

    typedef typename Algebraic_kernel_2::Curve_analysis_2 Curve_analysis_2;
    
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    typedef typename Status_line_1::Bitstream_descartes Isolator;

    //! computes boundary between y-coordinates of two algebraic reals
    //! defined over the same vertical line
    struct Boundary_between 
            : public std::binary_function< Type, Type, Boundary > {
        
        Boundary operator()(const Type& r1, const Type& r2) const {

            CGAL_precondition(r1.y() != r2.y());

            Boundary res(0);

            Isolator isol1 =
                r1.curve().status_line_at_exact_x(r1.x()).isolator();

            Isolator isol2 =
                r2.curve().status_line_at_exact_x(r2.x()).isolator();
            
            Boundary low1, low2, high1, high2;
            
            while (true) {
                low1 = isol1.left_boundary(r1.arcno());
                high1 = isol1.right_boundary(r1.arcno());
                
                low2 = isol2.left_boundary(r2.arcno());
                high2 = isol2.right_boundary(r2.arcno());
                
                if (low1 > high2) {
                    res = ((low1 + high2)/Boundary(2));
                    break;
                }
                if (low2 > high1) {
                    res = ((low2 + high1)/Boundary(2));
                    break;
                }
                
                // else
                isol1.refine_interval(r1.arcno());
                isol2.refine_interval(r2.arcno());
            }

            CGAL::simplify(res);

            CGAL_postcondition_code(
                    CGAL::Comparison_result exp = CGAL::SMALLER
            );
            CGAL_postcondition_code(
                    if (r1.y() > r2.y()) {
                        exp = CGAL::LARGER;
                    }
            );
            CGAL_postcondition(r1.y().compare(res) == exp);
            CGAL_postcondition(r2.y().compare(res) == -exp);

            return res;
        }
    };
    
    //! returns current lower boundary of an isolating interval defining 
    //! y-coordinate of an algebraic real
    struct Lower_boundary
            : public std::unary_function<Type, Boundary> {
        
        Boundary operator()(const Type& r) const {
            return r.curve().status_line_at_exact_x(r.x()).
                lower_boundary(r.arcno());
        }
    };

    //! returns current upper boundary of an isolating interval defining 
    //! y-coordinate of an algebraic real
    struct Upper_boundary
            : public std::unary_function<Type, Boundary> {
         
        Boundary operator()(const Type& r) const {
            return r.curve().status_line_at_exact_x(r.x()).
                upper_boundary(r.arcno());
        }
    };
                
    struct Refine
            : public std::unary_function<Type, void> {
        
        //! \brief refines isolating interval of an y-coorinate of algebraic 
        //! real to make it at least half of the original interval
        //!
        //! note that an interval may also degenerate to a single point
        void operator()(const Type& r) const {
            r.curve().status_line_at_exact_x(r.x()).refine(r.arcno());
        }
        
        //! \brief refines isolating interval of an y-coorinate of algebraic 
        //! real w.r.t. the given relative precision
        //!
        //! resulting interval is:
        //! <tt>|lower - upper|/|r.y()| <= 2^(-rel_prec)</tt> 
        void operator()(const Type& r, int rel_prec) const {
            Boundary low = Lower_boundary()(r);
            Boundary high = Upper_boundary()(r);
            
            Boundary prec = (high - low) /
                CGAL::ipower(Boundary(2), rel_prec);
            
            /////////// attention!! need to test for exact zero !!
               
            // Refine until both boundaries have the same sign
            while (CGAL::sign(low) != CGAL::sign(high)) {
                Refine()(r);
                low = Lower_boundary()(r);
                high = Upper_boundary()(r);
            }
            
            CGAL_assertion(
                    CGAL::sign(low) != CGAL::ZERO &&
                    CGAL::sign(high) != CGAL::ZERO);
            
            // Refine until precision is reached
            while (((high - low) / CGAL::max(CGAL::abs(high), CGAL::abs(low)))
                   > prec ) {
                Refine()(r);
                low = Lower_boundary()(r);
                high = Upper_boundary()(r);
            }
        }
    };
};

#endif

} // namespace CGALi 
            
CGAL_END_NAMESPACE
            
#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_ALGEBRAIC_REAL_TRAITS_H
