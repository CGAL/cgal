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
 *  Curve and curve pair analysis for algebraic plane curves
 */

#ifndef CGAL_FILTERED_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_FILTERED_ALGEBRAIC_CURVE_KERNEL_2_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_curve_kernel_2.h>

CGAL_BEGIN_NAMESPACE

template < class AlgebraicCurvePair_2, class AlgebraicKernel_1 >
class Filtered_algebraic_curve_kernel_2 
    : public CGAL::Algebraic_curve_kernel_2
          < AlgebraicCurvePair_2, AlgebraicKernel_1 > {

// for each predicate functor defines a member function returning an instance
// of this predicate
#define CGAL_Algebraic_Kernel_pred(Y,Z) \
    Y Z() const { return Y(); }

// the same for construction functors
#define CGAL_Algebraic_Kernel_cons(Y,Z) CGAL_Algebraic_Kernel_pred(Y,Z)

public:

    typedef AlgebraicCurvePair_2 Algebraic_curve_pair_2;
    typedef AlgebraicKernel_1 Algebraic_kernel_1;
    typedef CGAL::Algebraic_curve_kernel_2
        < Algebraic_curve_pair_2, Algebraic_kernel_1 > 
        Algebraic_curve_kernel_2;
    typedef Algebraic_curve_kernel_2 Base;

protected:

    // type of an internal curve
    typedef typename Base::Internal_curve_2 Internal_curve_2;

    // type of an internal curve pair
    typedef typename Base::Internal_curve_pair_2 Internal_curve_pair_2;

public:
    //! \name types and functors for \c GPA_2< >
    //!@{
    
    //! myself
    typedef Filtered_algebraic_curve_kernel_2
        < Algebraic_curve_pair_2, Algebraic_kernel_1 > Self;
    
    //! type of coefficient
    typedef typename Base::Coefficient Coefficient;

    //! type of x-coordinate
    typedef typename Base::X_coordinate_1 X_coordinate_1;

    // TODO remove when deriving from AK_1
    typedef typename Algebraic_kernel_1::Boundary Boundary;
        
    //! new CGAL univariate polynomial type (_CGAL postfix is temporary to
    //! avoid type clashes with \c Polynomial_2 type defined later
    typedef typename Base::Polynomial_1 Polynomial_1;
    //! new CGAL bivariate polynomial type
    typedef typename  Base::Polynomial_2 Polynomial_2;
    //! bivariate polynomial traits
    typedef typename Base::Polynomial_traits_2 Polynomial_traits_2;
    
    //!@}
                
    //! \brief default constructor
    Filtered_algebraic_curve_kernel_2() //: _m_curve_cache() 
    {  }
    
    typedef typename Base::Construct_curve_2 Construct_curve_2;

    CGAL_Algebraic_Kernel_cons(Construct_curve_2, construct_curve_2_object);
    
    //! type of a curve point 
    typedef typename Base::Xy_coordinate_2 Xy_coordinate_2;

    //! type of 1-curve analysis
    typedef typename Base::Curve_analysis_2 Curve_analysis_2;

    //! type of 2-curve analysis
    typedef typename Base::Curve_pair_analysis_2 Curve_pair_analysis_2;
    
    //! traits class for \c X_coordinate
    typedef typename Base::X_real_traits_1 X_real_traits_1;

    //! traits class for \c Xy_coorinate_2
    typedef typename Base::Y_real_traits_1 Y_real_traits_1;

public:
    
    typedef typename Xy_coordinate_2::Bbox_2 Bbox_2;

    static double& threshold() {
        static boost::optional<double> _b;
        if(! _b) {
            _b = .01; 
        }
        return _b.get();
    }
    
    //! returns the first coordinate of \c Xy_coordinate_2
    typedef typename Base::Get_x_2 Get_x_2;

    CGAL_Algebraic_Kernel_cons(Get_x_2, Get_x_2_object);

    
    //! returns the second coordinate of \c Xy_coordinate_2
    typedef typename Base::Get_y_2 Get_y_2;

    CGAL_Algebraic_Kernel_cons(Get_y_2, Get_y_2_object);


    typedef typename Base::Refine_x_2 Refine_x_2;

    CGAL_Algebraic_Kernel_pred(Refine_x_2, refine_x_2_object);


    typedef typename Base::Refine_y_2 Refine_y_2;

    CGAL_Algebraic_Kernel_pred(Refine_y_2, refine_y_2_object);


    //! computes the current lower boundary of the first coordinate of \c r
    typedef typename Base::Lower_boundary_x_2 Lower_boundary_x_2;
       
    CGAL_Algebraic_Kernel_cons(Lower_boundary_x_2, lower_boundary_x_2_object);

    
    //! computes the current upper boundary of the first coordinate of \c r
    typedef typename Base::Upper_boundary_x_2 Upper_boundary_x_2;

    CGAL_Algebraic_Kernel_cons(Upper_boundary_x_2, upper_boundary_x_2_object);

    
    //! computes the current lower boundary of the second coordinate of \c r
    typedef typename Base::Lower_boundary_y_2 Lower_boundary_y_2;
       
    CGAL_Algebraic_Kernel_cons(Lower_boundary_y_2, lower_boundary_y_2_object);
    
    //! computes the current lower boundary of the second coordinate of \c r
    typedef typename Base::Upper_boundary_y_2 Upper_boundary_y_2;

    CGAL_Algebraic_Kernel_cons(Upper_boundary_y_2, upper_boundary_y_2_object);
    
    //! returns the number of boundary type in-between x-coordinates of two
    //! Xy_coordinate_2 objects
    typedef typename Base::Boundary_between_x_2 Boundary_between_x_2;

    CGAL_Algebraic_Kernel_cons(Boundary_between_x_2, 
                               boundary_between_x_2_object);
            
    typedef typename Base::Boundary_between_y_2 Boundary_between_y_2;

    CGAL_Algebraic_Kernel_cons(Boundary_between_y_2, 
            boundary_between_y_2_object);
    
    //! \brief comparison of x-coordinates 
    struct Compare_x_2 :
        public Binary_function<X_coordinate_1, X_coordinate_1, 
                Comparison_result > {

        Comparison_result operator()(const X_coordinate_1& x1, 
                                     const X_coordinate_1& x2) const {
            return typename Base::Compare_x_2()(x1,x2);
        }
        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                     const Xy_coordinate_2& xy2) const {
            return ((*this)(xy1.x(), xy2.x()));
        }
    };
    CGAL_Algebraic_Kernel_pred(Compare_x_2, compare_x_2_object);

protected:

    struct Approximate_compare_y_2 :
        public Binary_function< Xy_coordinate_2, Xy_coordinate_2, 
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

    //! \brief comparison of y-coordinates of two points
    struct Compare_y_2 :
        public Binary_function< Xy_coordinate_2, Xy_coordinate_2, 
                Comparison_result > {
        
        Compare_y_2(Algebraic_curve_kernel_2 *kernel) :
            _m_kernel(kernel) {
        }

        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                     const Xy_coordinate_2& xy2) const {
            CGAL::Comparison_result approx_compare
                = Approximate_compare_y_2()(xy1,xy2);
            if(approx_compare!=CGAL::EQUAL) {
                return approx_compare;
            }
            return typename Base::Compare_y_2(_m_kernel)(xy1, xy2);
        }
        
    private:
        Algebraic_curve_kernel_2 *_m_kernel; 
    };
    //CGAL_Algebraic_Kernel_pred(Compare_y_2, compare_y_2_object);
    
    Compare_y_2 compare_y_2_object() const {
        return Compare_y_2((Algebraic_curve_kernel_2 *)this);
    }

    //! lexicographical comparison of two objects of type \c Xy_coordinate_2
    //!
    //! \c equal_x specifies that only y-coordinates need to be compared
    struct Compare_xy_2 :
          public Binary_function<Xy_coordinate_2, Xy_coordinate_2, 
                Comparison_result > 
    {

        Compare_xy_2(Algebraic_curve_kernel_2 *kernel) :
            _m_kernel(kernel) {
        }
        
        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                     const Xy_coordinate_2& xy2, 
                                     bool equal_x = false) const {

            CGAL::Comparison_result x_compare = Compare_x_2()(xy1.x(),xy2.x());
            if( x_compare != CGAL::EQUAL ) {
                return x_compare;
            }
            CGAL::Comparison_result approx_compare
                = Approximate_compare_y_2()(xy1,xy2);
            if(approx_compare!=CGAL::EQUAL) {
                return approx_compare;
            }
            return typename Base::Compare_xy_2(_m_kernel)(xy1,xy2,true);
        }

    private:
        Algebraic_curve_kernel_2 *_m_kernel;    

    };
    //CGAL_Algebraic_Kernel_pred(Compare_xy_2, compare_xy_2_object);

    Compare_xy_2 compare_xy_2_object() const {
        return Compare_xy_2((Algebraic_curve_kernel_2 *)this);
    }


    //! \brief checks whether curve has only finitely many self-intersection
    //! points, i.e., it has no self-overlapped continuous parts
    //!
    //! for algerbaic curves this means that supporting polynomial is 
    //! square-free
    typedef typename Base::Has_finite_number_of_self_intersections_2
        Has_finite_number_of_self_intersections_2;
    
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_self_intersections_2, 
            has_finite_number_of_self_intersections_2_object);
            
    //! \brief checks whether a curve pair has finitely many intersections,
    //! in other words, whether two curves have no continuous common part
    //!
    //! in case of algerbaic curves: checks whether supporting polynomials are
    //! coprime
    typedef typename Base::Has_finite_number_of_intersections_2
        Has_finite_number_of_intersections_2;
    
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_intersections_2, 
            has_finite_number_of_intersections_2_object);
    
    //! set of various curve and curve pair decomposition functions
    typedef typename Base::Decompose_2 Decompose_2;
    
    CGAL_Algebraic_Kernel_cons(Decompose_2, decompose_2_object);
    
    //!@}
public:
    //! \name types and functors for \c GPA_2<Algebraic_curve_kernel_2>
    //!@{
    
    typedef typename Base::Algebraic_real_1 Algebraic_real_1;
    typedef typename Base::Algebraic_real_2 Algebraic_real_2;
    
    typedef typename Base::Is_square_free_2 Is_square_free_2;
    typedef typename Base::Is_coprime_2 Is_coprime_2;

    typedef typename Base::Make_square_free_2 Make_squre_free_2;
    typedef typename Base::Square_free_factorization Square_free_factorization;
    typedef typename Base::Make_coprime_2 Make_coprime_2;
    
    //! \brief computes the derivative w.r.t. the first (innermost) variable
    typedef typename Base::Derivative_x_2 Derivative_x_2;

    CGAL_Algebraic_Kernel_cons(Derivative_x_2, derivative_x_2_object);

    //! \brief computes the derivative w.r.t. the first (outermost) variable
    typedef typename Base::Derivative_y_2 Derivative_y_2;

    CGAL_Algebraic_Kernel_cons(Derivative_y_2, derivative_y_2_object);

    typedef typename Base::X_critical_points_2 X_critical_points_2;

    CGAL_Algebraic_Kernel_cons(X_critical_points_2,
        x_critical_points_2_object);
    
    typedef typename Base::Y_critical_points_2 Y_critical_points_2;

    CGAL_Algebraic_Kernel_cons(Y_critical_points_2,
        y_critical_points_2_object);

    /*!\brief 
     * computes the sign of a bivariate polynomial \c p evaluated at the root 
     * \c r of a system of two bivariate polynomial equations
     *
     * returns a value convertible to \c CGAL::Sign
     */
    struct Sign_at_2 :
        public Binary_function< Curve_analysis_2, Xy_coordinate_2, Sign > {

        typedef typename Xy_coordinate_2::Boundary_interval Interval;
        
        typedef CGAL::Polynomial<Boundary> Poly_rat_1;
        typedef CGAL::Polynomial<Poly_rat_1> Poly_rat_2;
        
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
            return typename Base::Sign_at_2()(ca, r);
        }
        
    };

    CGAL_Algebraic_Kernel_pred(Sign_at_2, sign_at_2_object);

    /*!\brief
     * copies in the output iterator \c roots the common roots of polynomials
     * \c p1 and \c p2 and copies in the output iterator \c mults their 
     * respective multiplicity as intergers, in the same order
     *
     * template argument type of \c roots is \c Xy_coordinate_2 , returns the
     * pair of respective past-the-end iterators
     *
     * \pre p1 and p2 are square-free and the set of solutions of the system
     * is 0-dimensional
     */  
    typedef typename Base::Solve_2 Solve_2;
    
    CGAL_Algebraic_Kernel_cons(Solve_2, solve_2_object);

#undef CGAL_Algebraic_Kernel_pred    
#undef CGAL_Algebraic_Kernel_cons 
    
    //!@}
 
}; // class Filtered_algebraic_curve_kernel_2

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_ALGEBRAIC_CURVE_KERNEL_2_H
