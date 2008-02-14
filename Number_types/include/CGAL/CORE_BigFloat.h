// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_CORE_BIGFLOAT_H
#define CGAL_CORE_BIGFLOAT_H

#include <CGAL/number_type_basic.h>
#include <CGAL/CORE_coercion_traits.h>

CGAL_BEGIN_NAMESPACE

//template<typename BFI> class  Bigfloat_interval_traits;

// forward declarations of interval support
// namespace CGALi {
CORE::BigFloat 
inline 
round(const CORE::BigFloat& x, long rel_prec );

//TODO
// functions: upper and lower are defined twice the time, 
// please correct this unaesthetic thing
CORE::BigFloat 
inline
upper(CORE::BigFloat x){
    CORE::BigFloat result = ::CORE::BigFloat(x.m()+x.err(),0,x.exp());
    CGAL_postcondition(result >= x);
    return result; 
}

CORE::BigFloat 
inline 
lower(CORE::BigFloat x){
    CORE::BigFloat result = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
    CGAL_postcondition(result <= x);
    return result; 
}
// CORE::BigFloat 
// inline
// upper(CORE::BigFloat x);
// 
// CORE::BigFloat 
// inline 
// lower(CORE::BigFloat x);

// } // namespace CGALi

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigFloat >
  : public Algebraic_structure_traits_base< CORE::BigFloat,
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false          Is_exact;
    typedef Tag_true           Is_numerical_sensitive;

    class Sqrt
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            // What I want is a sqrt computed with ::CORE::defRelPrec bits.
            // And not ::CORE::defBFsqrtAbsPrec as CORE does. 
            
            CGAL_precondition(::CORE::defRelPrec.toLong() > 0);
            CGAL_precondition(x > 0);
            
            Type a = CGAL::round(x, ::CORE::defRelPrec.toLong()*2);
            CGAL_postcondition(a > 0); 

            Type tmp1 = 
                CORE::BigFloat(a.m(),0,0).sqrt(::CORE::defRelPrec.toLong());
            Type err  =  
                Type(0,long(std::sqrt(double(a.err()))),0) 
                * CORE::BigFloat::exp2(a.exp()*7);
            Type result = tmp1*CORE::BigFloat::exp2(a.exp()*7) + err;
           
            CGAL_postcondition(result >= 0);
//#ifndef NDEBUG
//            Type tmp = result * result; 
//            typedef Bigfloat_interval_traits<CORE::BigFloat> bfi_traits;
            //bfi_traits::Upper upper;
            //bfi_traits::Lower lower;
//            CGAL_postcondition(NiX::lower(tmp) <= NiX::lower(x));
//            CGAL_postcondition(NiX::upper(tmp) >= NiX::upper(x));
//#endif
            return result;
        }
    };

    class Kth_root
      : public Binary_function<int, Type, Type> {
      public:
        Type operator()( int k,
                                        const Type& x) const {
            CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
            // CORE::radical isn't implemented for negative values of x, so we
            //  have to handle this case separately
            if( x < 0 && k%2 != 0) {
              return Type(-CORE::radical( -x, k ) );
            }


            return Type( CORE::radical( x, k ) );
        }
    };
};

//
// Real embeddable traits
//
template <> class Real_embeddable_traits< CORE::BigFloat >
  : public Real_embeddable_traits_base< CORE::BigFloat > {
  public:

    class Abs
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            Type result; 
          
            if(x.isZeroIn()){
                CORE::BigInt m; 
                if(x.m() < 0 ){
                    m = -(x.m()-x.err());
                }else{
                    m =  x.m()+x.err();
                }
                if(m % 2 == 1) m += 1;
                
                Type upper(m,0,x.exp());
                result = CORE::centerize(CORE::BigFloat(0),upper);
                
                CGAL_postcondition(result.m()-result.err() <= 0); 
                if(result.m()-result.err() != 0){
                    result = this->operator()(result);
                }
                CGAL_postcondition(result.m()-result.err() == 0); 
            }else{
                result = CORE::abs(x);
            }
            CGAL_postcondition(result.m()-result.err() >= 0); 
            CGAL_postcondition(Type(result.m()+result.err(),0,result.exp()) 
                         >= Type(x.m()+x.err(),0,x.exp()));       
            return result;
        }
    };

    class Sign
      : public Unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
            ::CGAL::Sign result =  sign( x.sign());
            return result; 
        }
    };

    class Compare
      : public Binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) sign( (x-y).sign());
        }
    };

    class To_double
      : public Unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          // this call is required to get reasonable values for the double
          // approximation
          return x.doubleValue();
        }
    };

    class To_interval
      : public Unary_function< Type, std::pair< double, double > > {
    public:
        std::pair<double, double> operator()( const Type& x ) const {
                        
            double lb,ub;
           
            Type x_lower = CGAL::lower(CGAL::round(CGAL::lower(x),52));
            Type x_upper = CGAL::upper(CGAL::round(CGAL::upper(x),52));
            
            // since matissa has 52 bits only, conversion to double is exact 
            lb = x_lower.doubleValue();
            CGAL_postcondition(lb == x_lower);
            ub = x_upper.doubleValue();
            CGAL_postcondition(ub == x_upper);             
            
            std::pair<double, double> result(lb,ub);
            CGAL_postcondition( result.first <=  CORE::Expr(CGAL::lower(x)));
            CGAL_postcondition( result.second >=  CORE::Expr(CGAL::upper(x)));
            return result;      
        }
    };
};

//**********************NEW**********************

// Specialization of CGAL::Float_traits_kernel for CORE::BigFloat

template <> class Float_traits_kernel< ::CORE::BigFloat > {
public:
	typedef ::CORE::BigFloat		FT;
	typedef ::CORE::BigInt		  	IT;
	typedef ::CGAL::Float_number_tag	Float_type;
	typedef CGAL::Precision_type     Precision_type;
	
	class Get_precision {
	public:
		typedef CGAL::Precision_type	result_type;
		
		CGAL::Precision_type operator()() const {
			return CGAL::Precision_type( ::CORE::defRelPrec.toLong());
		}
		
	};
	
	class Set_precision {							
	public:
		typedef void 				result_type;
		typedef CGAL::Precision_type first_argument_type;
		
		void operator()( const CGAL::Precision_type& precision ) const {
			::CORE::defRelPrec = precision;
			::CORE::defBFdivRelPrec = precision;
		}
		
	};			

	class Get_rounding_mode {
	public:
		typedef CGAL::Rounding_mode	result_type;
		
		/*??
		CGAL::Rounding_mode operator()() const {
			return CGAL::Rounding_mode( CORE::BigFloat::get_rounding_mode() );
		}
		*/
	};

	class Set_rounding_mode {			
	public:
		typedef void				result_type;
		typedef CGAL::Rounding_mode 	first_argument_type;
	
		/* ??
		void operator()( const CGAL::Rounding_mode& rounding_mode ) const {
			CORE::BigFloat::set_rounding_mode( (CORE::rounding_modes)rounding_mode );
		}
		*/
	};
	
};


// Specialization of CGAL::Float_traits for CORE::BigFloat

template <> class Float_traits< ::CORE::BigFloat > 
	: public Float_traits_base< Float_traits_kernel< ::CORE::BigFloat> >  {
		
private:
	/* ??
	static Rounding_mode get_rounding_mode() {
		return Rounding_mode( CORE::BigFloat::get_rounding_mode() );
	};
	*/		
	
	static Precision_type get_precision() {
		return Precision_type( ::CORE::defRelPrec.toLong() );
	};
	
public:			
	typedef Float_traits_kernel< CORE::BigFloat > 	FT_kernel;
	typedef FT_kernel::FT							FT;
	typedef FT_kernel::IT							IT;
	typedef FT_kernel::Float_type					Float_type;

	class Construct 
		: public Float_traits_base< FT_kernel >::Construct {
	public:
		FT operator() ( const IT& m, const IT& e = IT(0) ) {
		return FT(m,0,0)*FT::exp2(e.intValue());
		}
	}; // class Constructor
	
	class Set_precision {
	public:
		// type for the \c AdaptableUnaryFunction concept.
		typedef long  argument_type;
		// type for the \c AdaptableUnaryFunction concept.
		typedef long  result_type;  
		
		long operator() ( long prec ) const {
			FT_kernel::Set_precision set_precision;  
			FT_kernel::Get_precision get_precision;  
			long old_precision = get_precision();
			set_precision(prec);
			return old_precision;
		}
	};
	
	class Get_precision {
	public:
		// type for the \c AdaptableGenerator concept.
		typedef long  result_type;  
	
		long operator() () const {
			FT_kernel::Get_precision get_precision;  
			return get_precision();
		}
	};



		
	/*
	class Construct_special_value {
	public:
		FT operator()( const Special_value& sv ) {
			switch( sv ) {
			case SV_N_ZERO: return leda::bigfloat::nZero;
			case SV_P_ZERO: return leda::bigfloat::pZero;
			case SV_N_INF: return leda::bigfloat::nInf;
			case SV_P_INF: return leda::bigfloat::pInf;
			case SV_NAN: return leda::bigfloat::NaN; 
			default: return FT(0);
			}
		}
	};
	*/

	/* ??
	class Is_special_value {
	public:				
		bool operator()( const FT& bf ) const {
			return leda::isSpecial( bf );
		}
	};
	*/
	
	/* ??
	class Get_special_value {
	public:
		Special_value operator()( const FT& bf ) const {
			if( leda::isNaN( bf ) )
				return SV_NAN;
			else if( leda::isnZero( bf ) )
				return SV_N_ZERO;
			else if( leda::ispZero( bf ) )
				return SV_P_ZERO;						
			else if( leda::isnInf( bf ) )
				return SV_N_INF;					
			else if( leda::ispInf( bf ) )
				return SV_P_INF;
			else	
				return SV_NONE;
		}
	};
	*/
	
	class Is_zero {
	public:
		bool operator()( const FT& bf ) const {
			CGAL_precondition(bf.isExact());
			return bf == FT(0);
		}
	};

	/* ??
	class Is_inf {
	public:
		bool operator()( const FT& bf ) const {
			return leda::isInf( bf );
		}
	};
	*/
			
	class Get_mantissa 
		: public Float_traits_base< FT_kernel >::Get_mantissa {
	public:
		//! type for the \c AdaptableUnaryFunction concept.
		typedef IT result_type;
		//! type for the \c AdaptableUnaryFunction concept.
		typedef FT argument_type;
		IT operator()( const FT& a ) const {
			return a.m();
		}
	}; // class GetMantissa

	class Get_exponent 
		: public Float_traits_base< FT_kernel >::Get_exponent {
	public:
		//! type for the \c AdaptableUnaryFunction concept.
		typedef long result_type;
		//! type for the \c AdaptableUnaryFunction concept.
		typedef FT argument_type;
		
		result_type operator()( const FT& a ) const {
			return 14*a.exp(); // The basis is 8092 
		}
	}; // class GetExponent
	
	class Get_mantissa_and_exponent 
		: public Float_traits_base< FT_kernel >::Get_mantissa_and_exponent {
	public:
		//! type for the \c AdaptableUnaryFunction concept.
		typedef std::pair< IT, long > result_type;
		//! type for the \c AdaptableUnaryFunction concept.
		typedef FT argument_type;
		
		result_type operator()( const FT& a ) const {
			IT m = a.m();
			long e = 14*a.exp();
			result_type p(  m, e );
			
			return p;
		}
	}; // class GetMantissaAndExponent
	
	/* ??
	class Round_relative
		: public Float_traits_base< FT_kernel >::Round_relative {
	public:
		//! type for the \c AdaptableTernaryFunction concept.
		typedef FT result_type;
		//! type for the \c AdaptableTernaryFunction concept.
		typedef FT first_argument_type;
		//! type for the \c AdaptableTernaryFunction concept.
		typedef CGAL::Precision_type second_argument_type;
		//! type for the \c AdaptableTernaryFunction concept.
		typedef CGAL::Rounding_mode third_argument_type;
		
		FT operator()( const FT& a, 
						const CGAL::Precision_type& p = get_precision(), 
						const CGAL::Rounding_mode& m = get_rounding_mode() ) const {
			return leda::round(a,(leda::sz_t)p,(leda::rounding_modes)m ) ;
		}
	}; // class Round
	
	*/

	/* ??
		class Round_absolute : public Float_traits_base< FT_kernel >::Round_absolute {
		public:
			//! type for the \c AdaptableTernaryFunction concept.
			typedef FT result_type;
			//! type for the \c AdaptableTernaryFunction concept.
			typedef FT first_argument_type;
			//! type for the \c AdaptableTernaryFunction concept.
			typedef IT second_argument_type;
			//! type for the \c AdaptableTernaryFunction concept.
			typedef CGAL::Rounding_mode third_argument_type;
			
			FT operator()( const FT& a,
							const IT& new_exponent,
							const CGAL::Rounding_mode& r_mode = get_rounding_mode() ) const {
				typedef CGAL::Float_number<IT> FN;
				FN fn( a.get_significant(), a.get_exponent() );
				fn = fn.round_absolute( new_exponent, r_mode );
				return FT( fn.get_mantissa(), fn.get_exponent() );			    				    	
			}
		};
	*/
		// special functors ( only for Big_float )
		
		/* ??
		class Add
		: public Float_traits_base< FT_kernel >::Add {
		public:
			typedef FT result_type;
			typedef FT first_argument_type;
			typedef FT second_argument_type;
			typedef CGAL::Rounding_mode third_argument_type;
			typedef CGAL::Precision_type fourth_argument_type;
			
			FT operator()( const FT& a, const FT& b, 
						const Rounding_mode& m = get_rounding_mode(),	
							const Precision_type& p = get_precision() ) const {
					return leda::add(a,b,(leda::sz_t)p,(leda::rounding_modes)m);
				}
		}; // class Add

		class Sub 
			: public Float_traits_base< FT_kernel >::Sub{
			public:
				typedef FT result_type;
				typedef FT first_argument_type;
				typedef FT second_argument_type;
				typedef CGAL::Rounding_mode third_argument_type;
				typedef CGAL::Precision_type fourth_argument_type;

				FT operator()( const FT& a, const FT& b , 
							const Rounding_mode& m = get_rounding_mode(),	
							const Precision_type& p = get_precision() ) const {
					return leda::sub(a,b,(leda::sz_t)p,(leda::rounding_modes)m);
				}
		}; // class Sub


		class Mul 
			: public Float_traits_base< FT_kernel >::Mul {
			public:
				typedef FT result_type;
				typedef FT first_argument_type;
				typedef FT second_argument_type;
				typedef CGAL::Rounding_mode third_argument_type;
				typedef CGAL::Precision_type fourth_argument_type;

				FT operator()( const FT& a, const FT& b , 
							const Rounding_mode& m = get_rounding_mode(),	
							const Precision_type& p = get_precision() ) const {
					return leda::mul(a,b,(leda::sz_t)p,(leda::rounding_modes)m);
				}
		}; // class Mul

		class Div
			: public Float_traits_base< FT_kernel >::Div{
			public:
				typedef FT result_type;
				typedef FT first_argument_type;
				typedef FT second_argument_type;
				typedef CGAL::Rounding_mode third_argument_type;
				typedef CGAL::Precision_type fourth_argument_type;

				FT operator()( const FT& a, const FT& b , 
							const Rounding_mode& m = get_rounding_mode(),	
							const Precision_type& p = get_precision() ) const {
					return leda::div(a,b,(leda::sz_t)p,(leda::rounding_modes)m);
				}
		}; // class Div
	
		class Sqrt
			: public Float_traits_base< FT_kernel >::Sqrt{
			public:
				typedef FT result_type;
				typedef FT first_argument_type;
				typedef CGAL::Rounding_mode second_argument_type;
				typedef CGAL::Precision_type third_argument_type;

				FT operator()( const FT& bf,
							const Rounding_mode& r_mode = get_rounding_mode(),	
							const Precision_type& prec = get_precision() ) const {
					return leda::sqrt( bf, (leda::sz_t)prec, (leda::rounding_modes)r_mode);
				}
		}; // class Sqrt

		*/

	class Mul_by_pow_of_2 
		: public Float_traits_base< FT_kernel >::Mul_by_pow_of_2{
	public:
		//! type for the \c AdaptableBinaryFunction concept.
		typedef FT result_type;
		//! type for the \c AdaptableBinaryFunction concept.
		typedef FT first_argument_type;
		//! type for the \c AdaptableBinaryFunction concept.
		typedef IT second_argument_type;
	
		FT operator()( const FT& a, const IT& e ) const {
			return a*FT::exp2(e.intValue());
		}
	}; // class Mul_by_pow_of_2

	class Div_by_pow_of_2 
		: public Float_traits_base< FT_kernel >::Div_by_pow_of_2{
	public:
		//! type for the \c AdaptableBinaryFunction concept.
		typedef FT result_type;
		//! type for the \c AdaptableBinaryFunction concept.
		typedef FT first_argument_type;
		//! type for the \c AdaptableBinaryFunction concept.
		typedef IT second_argument_type;
		
		FT operator()( const FT& a, const IT& e ) const {
			return a*FT::exp2(-e.intValue());
		}
	}; // class Div_by_pow_of_2

	class Floor {
	public:
		typedef IT result_type;
		typedef FT argument_type;
		
		IT operator()( const FT& fn ) const {
			IT fl = fn.BigIntValue();
			if(fn.sign() < 0 && fn.cmp(fl)!=0) {
				fl--;
			}
			return fl;
		}
	};
	
	class Ceil {
	public:
		typedef IT result_type;
		typedef FT argument_type;
		
		IT operator()( const FT& fn ) const {
			IT fl = fn.BigIntValue();
			if(fn.sign() >0 && fn.cmp(fl)!=0) {
				fl++;
			}
			return fl;
		}
	};
		
	
};

//**********************NEW**********************

CGAL_END_NAMESPACE

// #include <CGAL/Number_types/core_interval_support.h>
#include <CGAL/core_interval_support.h>

//since types are included by CORE_coercion_traits.h:
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>

#endif // CGAL_CORE_BIGFLOAT_H
