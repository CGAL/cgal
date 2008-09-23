// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file Filtered_algebraic_curve_kernel_2.h
 *  \brief defines class \c Filtered_algebraic_curve_kernel_2
 *  
 *  Algebraic_curve_kernel_2 with some filter techniques
 */

#ifndef CGAL_FILTERED_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_FILTERED_ALGEBRAIC_CURVE_KERNEL_2_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_curve_kernel_2.h>

CGAL_BEGIN_NAMESPACE

/*! 
 * \b Filtered_algebraic_curve_kernel_2 is an extension of
 * \c Algebraic_curve_kernel_2, thus also a model of
 * \c AlgebraicKernelWithAnalysis_d_2 and \c CurveKernel_2.
 * It redefines some functors in a way that numeric filters are used before
 * the exact computation is triggered.
 *
 * The basic idea behind all numeric computations in this kernel is that
 * coordinates are refined up to a certain threshold, and the exact method
 * is called if the approximate functor did not succeed until the threshold
 * was reached. This threshold is a double value, defined by the flag
 * CGAL_ACK_THRESHOLD_FOR_FILTERED_KERNEL, its default value is set
 * inside the file Algebraic_curve_kernel/flags.h
 */
#if CGAL_ACK_USE_EXACUS
template < class AlgebraicCurvePair_2, class AlgebraicKernel_1 >
#else
template < class AlgebraicKernel_1 >
#endif
class Filtered_algebraic_curve_kernel_2 
    : public CGAL::Algebraic_curve_kernel_2
#if CGAL_ACK_USE_EXACUS
    < AlgebraicCurvePair_2, AlgebraicKernel_1 > {
#else
    < AlgebraicKernel_1 > {
#endif

// for each predicate functor defines a member function returning an instance
// of this predicate
#define CGAL_Algebraic_Kernel_pred(Y,Z) \
    Y Z() const { return Y(); }

// the same for construction functors
#define CGAL_Algebraic_Kernel_cons(Y,Z) CGAL_Algebraic_Kernel_pred(Y,Z)

public:

    typedef AlgebraicKernel_1 Algebraic_kernel_1;
    
#if CGAL_ACK_USE_EXACUS
    typedef AlgebraicCurvePair_2 Algebraic_curve_pair_2;

    typedef CGAL::Algebraic_curve_kernel_2
        < Algebraic_curve_pair_2, Algebraic_kernel_1 > 
        Algebraic_curve_kernel_2;
#else
    typedef CGAL::Algebraic_curve_kernel_2 < Algebraic_kernel_1 > 
        Algebraic_curve_kernel_2;
#endif

    typedef Algebraic_curve_kernel_2 Base;

#if CGAL_ACK_USE_EXACUS
protected:
    // type of an internal curve pair
    typedef typename Base::Internal_curve_pair_2 Internal_curve_pair_2;
#endif

public:
    //! \name types and functors for \c GPA_2< >
    //!@{
    
    //! myself
#if CGAL_ACK_USE_EXACUS
    typedef Filtered_algebraic_curve_kernel_2
        < AlgebraicCurvePair_2, Algebraic_kernel_1 > Self;
#else
    typedef Filtered_algebraic_curve_kernel_2 < Algebraic_kernel_1 > Self;
#endif
    
    //! type of coefficient
    typedef typename Base::Coefficient Coefficient;

    //! type of x-coordinate
    typedef typename Base::X_coordinate_1 X_coordinate_1;

    typedef typename Algebraic_kernel_1::Boundary Boundary;
        
    //!@}
                
    //! \brief default constructor
    Filtered_algebraic_curve_kernel_2() {  
    }
    
    //! type of a curve point 
    typedef typename Base::Xy_coordinate_2 Xy_coordinate_2;

    //! type of 1-curve analysis
    typedef typename Base::Curve_analysis_2 Curve_analysis_2;

    // Polynomial type
    typedef typename Base::Polynomial_2 Polynomial_2;

public:
    
    typedef typename Xy_coordinate_2::Bbox_2 Bbox_2;

    static double& threshold() {
        static boost::optional<double> _b;
        if(! _b) {
            _b = CGAL_ACK_THRESHOLD_FOR_FILTERED_KERNEL; 
        }
        return _b.get();
    }
    
protected:

    struct Approximate_compare_y_2 :
        public std::binary_function< Xy_coordinate_2, Xy_coordinate_2, 
                Comparison_result > {
        
        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                     const Xy_coordinate_2& xy2) const {

            Bbox_2 bbox1 = xy1.approximation_box_2(threshold()),
                bbox2 = xy2.approximation_box_2(threshold());
            if(bbox1.ymin() > bbox2.ymax()) {
                return CGAL::LARGER;
            }
            if(bbox1.ymax() < bbox2.ymin()) {
                return CGAL::SMALLER;
            }
            return CGAL::EQUAL;
        }
    };
    CGAL_Algebraic_Kernel_pred(Approximate_compare_y_2, 
                               approximate_compare_y_2_object);

public:

    //! Filtered comparison of y-coordinates of two points
    struct Compare_y_2 :
        public Base::Compare_y_2 {

        typedef typename Base::Compare_y_2 Base;
        
        Compare_y_2(Self *kernel) :
            Base(kernel) {
        }

        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                     const Xy_coordinate_2& xy2) const {

            CGAL::Comparison_result res = Approximate_compare_y_2()(xy1,xy2);
            if(res != CGAL::EQUAL) 
                return res;
            
            return Base::operator()(xy1, xy2);
        }
    };
    Compare_y_2 compare_y_2_object() const {
        return Compare_y_2((Self *)this);
    }

    //! Filtered lexicographic comparison of two points
    struct Compare_xy_2 :
          public Base::Compare_xy_2 {

        typedef typename Base::Compare_xy_2 Base;

        Compare_xy_2(Self *kernel) :
            Base(kernel) {
        }
        
        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                     const Xy_coordinate_2& xy2, 
                                     bool equal_x = false) const {

            CGAL::Comparison_result res = 
                Base::_m_kernel->compare_x_2_object()(xy1.x(),xy2.x());
            if(res != CGAL::EQUAL) 
                return res;
            
            res = Approximate_compare_y_2()(xy1,xy2);
            if(res != CGAL::EQUAL) 
                return res;
            
            return Base::operator()(xy1, xy2, true);
        }
    };
    //CGAL_Algebraic_Kernel_pred(Compare_xy_2, compare_xy_2_object);

    Compare_xy_2 compare_xy_2_object() const {
        return Compare_xy_2((Self *)this);
    }

    //! Filtered sign computation of polynomial and point
    struct Sign_at_2 :
        public Base::Sign_at_2 {

        typedef typename Base::Sign_at_2 Base;

        typedef typename Xy_coordinate_2::Boundary_interval Interval;
        
        typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
            ::template Rebind<Boundary,1>::Other::Type
            Poly_rat_1;
        typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
            ::template Rebind<Boundary,2>::Other::Type
            Poly_rat_2;
        
        Sign operator()(const Curve_analysis_2& ca,
                const Xy_coordinate_2& r) const
        {
            if(ca.is_identical(r.curve())) // point lies on the same curve
                return CGAL::ZERO;
            
            r.approximation_box_2(threshold()); 
                // makes sure refined boundaries 
            Interval iv = r.interval_evaluate_2(ca.polynomial_2());
            CGAL::Sign s_lower = CGAL::sign(iv.lower());
            if( s_lower == CGAL::sign(iv.upper()) ) {
                return s_lower;
            }
            return Base::operator()(ca, r);
        }

        Sign operator()(Polynomial_2& f,
                        const Xy_coordinate_2& r) const
        {
            return (*this)(typename Base::Construct_curve_2()(f),r);
        }
        
    };
    CGAL_Algebraic_Kernel_pred(Sign_at_2, sign_at_2_object);

#undef CGAL_Algebraic_Kernel_pred    
#undef CGAL_Algebraic_Kernel_cons 
    
    //!@}
 
}; // class Filtered_algebraic_curve_kernel_2

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_ALGEBRAIC_CURVE_KERNEL_2_H
