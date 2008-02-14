// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: 4515 rschindl$
//
//
// Author(s)     : Ralf Schindlbeck <rschindl@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


/*! \file CGAL/Float_traits.h
   \brief Defines class CGAL::Float_traits. 
 
  Provides dependent types and function objects for all the functions
  beyond operators of the different float types.
*/
#ifndef CGAL_FLOAT_TRAITS_H
#define CGAL_FLOAT_TRAITS_H

#include <CGAL/basic.h>

// #include <NiX/enums.h>

#include <stack>

// namespace NiX {
CGAL_BEGIN_NAMESPACE

//compensation for the #include <NiX/enums.h>

//Enumerator indicating the desired rounding mode for an operation
//and the global rounding mode, if exists.
enum Rounding_mode {
    TO_NEAREST = 0, //!< round to nearest number
    TO_ZERO    = 1, //!< round towards zero
    TO_P_INF   = 2, //!< round up (towards positive infinity)
    TO_N_INF   = 3, //!< round down (towards negative infinity)
    TO_INF     = 4  //!< round away from zero
};
//Enumerator indicating a special value represented by a floating point value.
enum Special_value {
    SV_NONE   =  0, //!< indicates non-specialness of a value
    SV_N_ZERO = -1, //!< zero with a negative sign
    SV_P_ZERO =  1, //!< zero with a positive sign
    SV_N_INF  = -2, //!< negative infinity
    SV_P_INF  =  2, //!< positive infinity
    SV_NAN    =  3  //!< "Not Any Number (not even infinity)": error indicator
};

typedef long Precision_type;

// The tags for Float type corresponding to the float type concepts
// ================================================================
  
//! corresponds to the \c FloatNumber concept.
struct Float_number_tag {};
  
//! corresponds to the \c ExactFloatNumber concept.
struct Exact_float_number_tag {};
  
//! corresponds to the \c FixedPointNumber concept.
struct Fixed_point_number_tag {};

// The float traits kernel class 
// =========================================================================
	
// This class provides implementation dependent access to the precision and
// rounding mode settings along with information about the float type, the
// corresponding integer type and the float concept.
	
template< class FT_ >
class Float_traits_kernel {
public:
    typedef FT_							FT;
    typedef CGAL::Null_tag			IT;
    typedef CGAL::Null_tag Float_type;
				
    typedef CGAL::Null_tag			Precision_type;
    typedef CGAL::Null_functor		Get_precision;
    typedef CGAL::Null_functor		Set_precision;
    typedef CGAL::Null_functor		Get_rounding_mode;
    typedef CGAL::Null_functor		Set_rounding_mode;
};

// The float traits template
// =========================

// The general float traits is only declared and not defined.
// That leads to an immediate error message at the first attempt to use
// a non-supported float type.
// Documentation for the concept of being a valid specialization
// of Float_traits is provided in NumeriX/manual/AlgebraicKernel/AK_Float.dox

template< class FT_ >
class Float_traits;
 

// The float traits base class
// =========================================================================


//! The template specialization that can be used for types that are not any
//! of the float type concepts. All functors are set to \c CGAL::Null_functor.
//! See also \link NiX_Float_traits_functors concept Float_traits \endlink .
//! \ingroup NiX_Float_traits_bases
//
template< class FT_kernel_ > 
class Float_traits_base {
private:
    typedef FT_kernel_	FT_kernel;

public:  		   
    typedef CGAL::Null_functor  Construct;
    typedef CGAL::Null_functor  Get_mantissa;
    typedef CGAL::Null_functor  Get_exponent;
    typedef CGAL::Null_functor  Get_mantissa_and_exponent;
    typedef CGAL::Null_functor  Round_relative;
    typedef CGAL::Null_functor  Round_absolute;

    typedef CGAL::Null_functor  Add;
    typedef CGAL::Null_functor  Sub;
    typedef CGAL::Null_functor  Mul;
    typedef CGAL::Null_functor  Div;
    typedef CGAL::Null_functor	Sqrt;
    typedef CGAL::Null_functor  Mul_by_pow_of_2;
    typedef CGAL::Null_functor  Div_by_pow_of_2;
    	
    typedef CGAL::Null_functor	Floor;
    typedef CGAL::Null_functor	Ceil;

    typedef CGAL::Null_functor	Construct_special_value;

    typedef CGAL::Null_functor	Is_special_value;
    typedef CGAL::Null_functor	Get_special_value;
    typedef CGAL::Null_functor	Is_zero;
    typedef CGAL::Null_functor	Is_inf;
					
    typedef CGAL::Null_functor  Set_exponent;

    class Rounding_mode_controller {
    public:
			
        inline Rounding_mode_controller( const Rounding_mode& mode ) {
            // store current rounding mode
            rounding_modes.push( get_rounding_mode() );

#ifndef NDEBUG
            // store scope id
            scope_id = rounding_modes.size();
#endif					
            set_rounding_mode( mode );
        };
			
        inline ~Rounding_mode_controller() {
            // restore old rounding mode
            set_rounding_mode( rounding_modes.top() );
					
            rounding_modes.pop();
        };
			
        inline void set_rounding_mode( 
                const Rounding_mode& mode ) const {
#ifndef NDEBUG
            if( scope_id != rounding_modes.size() )
                CGAL_error_msg( "Rounding mode controller out of scope" );
#endif
            typename FT_kernel::Set_rounding_mode set_mode;
            set_mode( mode );				
        };
			
        inline static Rounding_mode get_rounding_mode() {
            typename FT_kernel::Get_rounding_mode get_mode;
            return get_mode();
        };
			
    private:				
        static std::stack<Rounding_mode> rounding_modes;	
#ifndef NDEBUG
        unsigned scope_id;
#endif
    };

    class Precision_controller {            
    public:
        typedef typename FT_kernel::Precision_type PT;
			
        Precision_controller( const PT& precision ) {
            // store current precision
            precisions.push( get_precision() );
					
#ifndef NDEBUG
            // store scope id
            scope_id = precisions.size();
#endif
            set_precision( precision );					
        };
			
        inline ~Precision_controller() {
            // restore previous precision
            set_precision( precisions.top() );

            precisions.pop();
        };
			
        inline void set_precision( const PT& precision ) const {
#ifndef NDEBUG
            if( scope_id != precisions.size() )
                CGAL_error_msg( "Precision controller out of scope" );
#endif					
            typename FT_kernel::Set_precision set_prec;
            set_prec( precision );				
        };
			
        inline static PT get_precision() {
            typename FT_kernel::Get_precision get_prec;
            return get_prec();
        };
															
    private:			
        static std::stack<PT> precisions;				
#ifndef NDEBUG
        unsigned scope_id;
#endif				
    };				
};

template< class FT_kernel >
std::stack< typename FT_kernel::Precision_type > 
CGAL::Float_traits_base< FT_kernel >::Precision_controller::precisions = 
    std::stack< typename FT_kernel::Precision_type >();

template< class FT_kernel >
std::stack<Rounding_mode> 
CGAL::Float_traits_base< FT_kernel >::Rounding_mode_controller::rounding_modes = 
    std::stack<Rounding_mode>();


// namespace Intern {
namespace CGALi {
// IT_float_traits provide functionality for arbitary integer values
// such as Division (without and with rounding) and Shifting
// which are not explizit definied for integer value, but needed 
// for the float_type
template< class IT_  >
class IT_float_traits {
public:
    typedef IT_ IT;

    //! returns the length (in bits) of an integer value $v$.
    typedef CGAL::Null_functor Length;

    //! calculates for an integer value $e$, $2^e$
    class Pow_2 {
    public:
	IT operator() ( const IT& e ) const {
            if( e<0 ) return 0;
            if( e==0 ) return 1;
            IT a = 1;
            IT p = 2;
            IT n = e;
            while( n > 0 ) {
                if( n%2 == 1 )
                    a = p*a;
                n = n/2;
                p = p*p;
            }
            return a;
	}
    };
      

    class Ilog_2 { 
    public:
	IT operator()( const IT& a) const {
            IT log=0;
            IT pow = 1;
            while(pow < a)
                {
                    pow = pow * 2;
                    log = log + 1;
                }
            return log;
	}
    };

    class Floor_Ilog_2 { 
    public:
	IT operator()( const IT& a) const {
            IT log=0;
            IT pow = 1;
            IT floor_log;
            while(pow <= a)
                {
                    pow = pow * 2;
                    floor_log = log;
                    log = log + 1;
                }
            return floor_log;
	}
    };
      

    //! count the number of zeros at the end in the binary reprasentation of an integer
    class Zeros {
    public:
	IT operator() ( const IT& value ) const {
            if ( value == IT(0) ) return 0;
            IT zeros = -1;
            IT d = value;
            IT f = 1;
            do {
                zeros++;
                d /= IT(2);
                f *= IT(2);
            } while ( d * f == value ) ;

            return zeros;
	}
    };


    //! Division of two integer values $a/b$
    //! return an integer div which satisfy the equation
    //! a = div * b + r (r < b )
    //! Division by zero is not definied.
    class Integer_Div {
    public:
	void help ( const IT&a, const IT& b, IT& sol, IT& sum ) {
            if ( a<b ){ 
                sol = 0;
                sum = 0;
                return;
            }

            if ( a==b ) {
                sol = 1;
                sum = b;
                return ;
            }

            sol = IT( 1 );
            sum = b;
            while ( (sum+sum) <= a ) {
                sum += sum;
                sol += sol;
            }

            return ;
	}

	IT operator() ( const IT& a, const IT& b ) {
            typename CGAL::Real_embeddable_traits<IT>::Abs abs;

            CGAL_precondition( b != 0 );
	  
            IT sol = 0; // solution of division
            IT sum_rem = abs( a );
            IT abs_b = abs( b );
            do {
                IT tmpsol, tmpsum;

                help( sum_rem, abs_b, tmpsol, tmpsum );
                sol += tmpsol;
                sum_rem -= tmpsum;
            } while ( sum_rem >= abs_b );
	    // false => terminate 


            if (  ( a < 0 && b > 0 )
                    || ( a > 0 && b < 0 ) ) {
                sol = (IT)0 - sol;
            }

            return sol;   
	}
    }; // class Integer_Div

    class Mod {
    public:
        IT operator() ( const IT& a, const IT&b ) {
            typename IT_float_traits<IT>::Integer_Div int_div;
            IT divisor = int_div( a, b );
            return (a-(b*divisor));
	}
    };
      


    //! shifts an integer value by shift to the left
    //! (multiplication by $2^shift$)
    class Shift_left_by {
    public:
	IT operator() ( const IT& value, const IT& shift ) const {
            typename IT_float_traits<IT>::Pow_2 pow_2;
            return value* pow_2( shift );
	}
    }; // class Shift_left_by

    //! shifts an integer value by shift to the right
    //! (division by $2^shift$)
    class Shift_right_by {
    public:
	IT operator() ( const IT& value, const IT& shift ) const {
            typename IT_float_traits<IT>::Pow_2 pow_2;
            typename IT_float_traits<IT>::Integer_Div Integer_div;
            return Integer_div( value, pow_2( shift ) );
	}
    }; // class Shift_right_by

}; // class IT_float_traits



}  // namespace CGALi


//! For the exact binary number types of EXACUS,
//! this base class provides implementations of
//! \c NiX::NT_traits::Floor_log2_abs and \c NiX::NT_traits::Ceil_log2_abs
//! based on \c NiX::FloatTraits::Get_mantissa and
//! \c NiX::FloatTraits::Get_exponent together with the respective
//! \c NiX::NT_traits functor for the mantissa.
//! \ingroup NiX_NT_traits_bases
template <class NT>
class NT_traits_log2_float_traits_base {
//TODO
//porting Floor_log2_abs and Ceil_log2_abs from EXACUS2CGAL
public:
    //! <tt>NiX::NT_traits::Floor_log2_abs()(NT x)</tt> implemented as
    //! <tt>NiX::NT_traits<Get_mantissa::result_type>::Floor_log2_abs()(Get_mantissa()(x)) + Get_exponent()(x)</tt>
//     class Floor_log2_abs {
//     public:
//         //! type for the \c AdaptableUnaryFunction concept. 
//         typedef typename Float_traits<NT>::Get_exponent::result_type
//         result_type;
//         //! type for the \c AdaptableUnaryFunction concept. 
//         typedef NT argument_type;
//         //! the function call.
//         result_type operator() (argument_type x) {
//             typedef typename Float_traits<NT>::Get_mantissa Get_mantissa;
//             typedef typename Float_traits<NT>::Get_exponent Get_exponent;
//             typename NT_traits<typename Get_mantissa::result_type>
//                 ::Floor_log2_abs log;
//             return log(Get_mantissa()(x)) + Get_exponent()(x);
//         }
//     };
    //! <tt>NiX::NT_traits::Ceil_log2_abs()(NT x)</tt> implemented as
    //! <tt>NiX::NT_traits<Get_mantissa::result_type>::Ceil_log2_abs()(Get_mantissa()(x)) + Get_exponent()(x)</tt>
//     class Ceil_log2_abs {
//     public:
//         //! type for the \c AdaptableUnaryFunction concept. 
//         typedef typename Float_traits<NT>::Get_exponent::result_type
//         result_type;
//         //! type for the \c AdaptableUnaryFunction concept. 
//         typedef NT argument_type;
//         //! the function call.
//         result_type operator() (argument_type x) {
//             typedef typename Float_traits<NT>::Get_mantissa Get_mantissa;
//             typedef typename Float_traits<NT>::Get_exponent Get_exponent;
//             typename NT_traits<typename Get_mantissa::result_type>
//                 ::Ceil_log2_abs log;
//             return log(Get_mantissa()(x)) + Get_exponent()(x);
//         }
//     };
};

// } // namespace NiX
CGAL_END_NAMESPACE
#endif // CGAL_FT_TRAITS_H
// EOF
