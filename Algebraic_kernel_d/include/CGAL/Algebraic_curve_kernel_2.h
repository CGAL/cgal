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
//				   Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

// SweepX/test/include/simple_models.h

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_2_H

#include <CGAL/basic.h>

#include <CGAL/Algebraic_kernel_1.h>

#include <CGAL/Algebraic_curve_kernel_2/Xy_coordinate_2.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_vertical_line_1.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_analysis_2.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_pair_vertical_line_1.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_pair_analysis_2.h>

#include <algorithm>

CGAL_BEGIN_NAMESPACE

template < class AlgebraicCurvePair_2, class AlgebraicKernel_1 >
class Algebraic_curve_kernel_2 {

#define CGAL_Algebraic_Kernel_pred(Y,Z) \
    Y Z() const { return Y(); }

private:
    // WRAPPING TYPES
    //! type of an internal curve pair
    typedef AlgebraicCurvePair_2 Internal_curve_pair_2;

    //! type of an internal curve
    typedef typename AlgebraicCurvePair_2::Algebraic_curve_2 Internal_curve_2;

    //! type of internal x_coordinate
    typedef typename Internal_curve_2::X_coordinate Internal_x_coordinate;
    
    //! type of internal coefficient
    typedef typename Internal_curve_2::Coefficient Internal_coefficient;

	//! type of internal polynomial
	typedef typename Internal_curve_2::Poly_d	Internal_polynomial_2;
	
	typedef typename NiX::Polynomial_traits<Internal_polynomial_2>::
		Innermost_coefficient Innermost_coefficient;

public:

    // TYPES AND FUNCTORS for Curved_kernel_2< CurveAnalyser_2 >

	//! type of 1D algebraic kernel 
	typedef AlgebraicKernel_1 Algebraic_kernel_1;

	//! myself
	typedef Algebraic_curve_kernel_2<AlgebraicCurvePair_2, AlgebraicKernel_1>
		Self;
    
    //! type of coefficient
    typedef Internal_coefficient Coefficient;

	//! type of curve pair
	typedef Internal_curve_pair_2 Curve_pair_2;

    //! type of single curve
    typedef Internal_curve_2 Curve_2;

	//! type of x-coordinate
    typedef Internal_x_coordinate X_coordinate_1;
		
	//! new CGAL univariate polynomial type
	typedef ::CGAL::Polynomial<Innermost_coefficient> Polynomial_1;
	//! new CGAL bivariate polynomial type
    typedef ::CGAL::Polynomial<Polynomial_1> Polynomial_2;
	//! bivariate polynomial traits
	typedef ::CGAL::Polynomial_traits_d< Polynomial_2 > Polynomial_traits_2;
	
	// temporary functor providing conversion from Poly_in type to
	// Poly_out type
	template <class Poly_2_from, class Poly_2_to>
	struct Polynomial_converter 
	{
		typedef typename Poly_2_from::NT Poly_1_from;
		typedef typename Poly_2_to::NT Poly_1_to;
		
		// needed for make_transform_iterator
		typedef Poly_1_to result_type;
		
		Poly_2_to operator()(const Poly_2_from& p) const
		{
			//Polynomial_converter<Poly_in_2, Poly_out_2> convert;
			return Poly_2_to(
				::boost::make_transform_iterator(p.begin(), *this),
				::boost::make_transform_iterator(p.end(), *this));
		}
		
		Poly_1_to operator()(const Poly_1_from& p) const
		{
			return Poly_1_to(p.begin(), p.end());
		}
	};
	
	typedef Polynomial_converter<Internal_polynomial_2, Polynomial_2>
				NiX2CGAL_converter;
	
	typedef Polynomial_converter<Polynomial_2, Internal_polynomial_2>
				CGAL2NiX_converter;
		
	struct Construct_curve_2 :
			public Unary_function< Internal_polynomial_2, Curve_2 >
	{
		//CC_2(ACK_2& ack)

		Curve_2 operator()(const Internal_polynomial_2& f) const
		{
			return Curve_2(f);
		}
		
		Curve_2 operator()(const Polynomial_2& f) const
		{
			CGAL2NiX_converter cvt;
			return Curve_2(cvt(f));
		}
		
	};
	CGAL_Algebraic_Kernel_pred(Construct_curve_2, construct_curve_2_object);

	/*CC_2 cc_2_object() {
		return CC_2(*this);
	}*/
	
	//! type of a curve point 
	typedef CGALi::Xy_coordinate_2<Self> Xy_coordinate_2;
      
	struct Compare_x_2 :
		 public Binary_function<X_coordinate_1, X_coordinate_1, 
				Comparison_result >
	{
        Comparison_result operator()(const X_coordinate_1& x1, 
                                         const X_coordinate_1& x2) const {
		// not yet implemented in Algebraic_kernel_1
//             Algebraic_kernel_1::C ak;
// 			return (ak.compare_x_2_object()(x1, x2));
			return x1.compare(x2);
        }

        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                         const Xy_coordinate_2& xy2) const {
			return ((*this)(xy1.x(), xy2.x()));
        }
    };
	CGAL_Algebraic_Kernel_pred(Compare_x_2, compare_x_2_object);

    struct Compare_y_2 {
        // later!
    };
	CGAL_Algebraic_Kernel_pred(Compare_y_2, compare_y_2_object);
    
    //! compares two points lexicographically
	struct Compare_xy_2 :
          public Binary_function<Xy_coordinate_2, Xy_coordinate_2, 
				Comparison_result > 
	{
		Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                          const Xy_coordinate_2& xy2) const {
             return xy1.compare_xy(xy2);
        }
    };
	CGAL_Algebraic_Kernel_pred(Compare_xy_2, compare_xy_2_object);
    
	//! \brief checks whether curve has only finitely many self-intersection 
	//! points, i.e., it has no self-overlapped continuous parts
	//!
	//! for algerbaic curves this means that supporting polynomial is 
	//! square-free
	struct Is_not_self_overlapped_2 : // or Has_no_overlaps_2 ?, Is_simple_2 ?
			public Unary_function< Curve_2, bool > {
		
		bool operator()(const Curve_2& c) const {
			typename Polynomial_traits_2::Is_square_free is_square_free;
			NiX2CGAL_converter cvt;
			Polynomial_2 res = cvt(c.f());
			return is_square_free(res);
        }
    };
	CGAL_Algebraic_Kernel_pred(Is_not_self_overlapped_2, 
			is_not_self_overlapped_2_object);
            
    //! \brief checks whether a pair of curves has common part
	//!
	//! in case of algerbaic curves: checks whether supporting polynomials are
	//! coprime
	struct Is_not_decomposable_2 :  // No_common_one_dimensional_part_2 ?
			public Binary_function< Curve_2, Curve_2, bool > { 
				
		bool operator()(const Curve_2& c1, const Curve_2& c2) {
			// if curve ids are the same - non-decomposable
			if(c1.id() == c2.id()) 
				return true;
        
			typename Polynomial_traits_2::Gcd_up_to_constant_factor gcd_utcf;
			typename Polynomial_traits_2::Total_degree total_degree;
			NiX2CGAL_converter cvt;
			Polynomial_2 p1 = cvt(c1.f()), p2 = cvt(c2.f());
		/*	std::cout << "p1 = " << c1.f() << "; p2 = " << c2.f() << "\n" <<
				cvt_back(gcd_utcf(p1, p2)) << "\n";
		*/	
			return (total_degree(gcd_utcf(p1, p2)) == 0);  
        }
    };
	CGAL_Algebraic_Kernel_pred(Is_not_decomposable_2, 
			is_not_decomposable_2_object);
    
    struct Decompose_2 {
        
		//! \brief returns a curve without self-overlapping parts 
		//!
		//! in case of algebraic curves computes square-free part of supporting
		//! polynomial
        Curve_2 operator()(const Curve_2& c) {
            typename Polynomial_traits_2::Make_square_free make_square_free;
			NiX2CGAL_converter cvt;
			CGAL2NiX_converter cvt_back;
			
			return Curve_2(cvt_back(make_square_free(cvt(c.f()))));
        }
        
		//! \brief computes a square-free factorization of a curve \c c, 
		//! retunrs the number of pairwise coprime square-free factors
		//! 
		//! returns square-free pairwise coprime factors in \c fit and
		//! multiplicities in \c mit. Template argument type of \c fit is
		//! \c Curve_2, and \c mit is \c int
		template< class OutputIterator1, class OutputIterator2 >
        int operator()( const Curve_2& c, OutputIterator1 fit, 
				OutputIterator2 mit ) const {
			typename Polynomial_traits_2::
				Square_free_factorization_up_to_constant_factor factorize;
			NiX2CGAL_converter cvt;
			CGAL2NiX_converter cvt_back;
			std::vector<Polynomial_2> factors;
			int n_factors = factorize(cvt(c.f()), std::back_inserter(factors),
					mit); 
	        for(int i = 0; i < (int)factors.size(); i++) {
                *fit++=Curve_2(cvt_back(factors[i]));
            }
			return n_factors;
        }
        
		//! \brief computes for a given pair of curves \c c1 and \c c2 their 
		//! common  part \c oib and coprime parts \c oi1 and \c oi2 
		//! respectively; returns \c true if the curves were decomposed
		//!
		//! returns true if \c c1 and \c c2 are coprime. Template argument
		//! type of \c oi* is \c Curve_2
        template < class OutputIterator > 
        bool operator()(const Curve_2& c1, const Curve_2& c2,
            OutputIterator oi1, OutputIterator oi2, OutputIterator oib) {
			
			typedef std::vector<Curve_2> Curves;
			Curves parts_f, parts_g;
            if(Curve_2::decompose(c1, c2, 
				std::back_inserter(parts_f), std::back_inserter(parts_g))) {
				// copy the common part returned through both iterators
				// we need to move a common part from oi1/oi2 to oib
				*oib++ = parts_f[0];
				if(parts_f.size() > 1)
					std::copy(parts_f.begin() + 1, parts_f.end(), oi1);
				if(parts_g.size() > 1)
					std::copy(parts_g.begin() + 1, parts_g.end(), oi2);
				return true;
			}
			return false;
        }
    };
	CGAL_Algebraic_Kernel_pred(Decompose_2, decompose_2_object);
    
    
    // TYPES AND FUNCTORS for Curved_kernel_2< Algebraic_kernel_2 >
    
    typedef Curve_2 Polynomial_2;
    
    typedef Construct_curve_2 Construct_polynomial_2;

    typedef X_coordinate_1 Algebraic_real_1;
    typedef Xy_coordinate_2 Algebraic_real_2;
    
    typedef Is_not_self_overlapped_2 Is_square_free_2;
    typedef Is_not_decomposable_2 Is_coprime_2;

    typedef Decompose_2 Make_square_free_2;
    typedef Decompose_2 Square_free_factorization;
    typedef Decompose_2 Make_coprime_2;
    
    struct Derivative_x_2 {
        // later
    };

    struct Derivative_y_2 {
        // later
    };

    struct X_critical_points_2 {
        // later
    };

    struct Y_critical_points_2 {
        // later
     };

    struct Sign_at_2 {
        // can be implemented using Curve_pair_analysis -> later
    };

    struct Solve_2 {
        // can be implemented using Curve_pair_analysis -> later
    };

    // TYPES AND FUNCTORS for Curved_kernel_2< both >
    
	//! type of 1-curve analysis
	typedef CGALi::Curve_analysis_2<Self> Curve_analysis_2; 

	//! type of a curve pair analysis
    typedef CGALi::Curve_pair_analysis_2<Self> Curve_pair_analysis_2; 
      
}; // class Algebraic_curve_kernel_2

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_2_H
