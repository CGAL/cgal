// Copyright (c) 1999,2007  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Stefan Schirra, Michael Hemmer

#ifndef CGAL_LEDA_BIGFLOAT_H
#define CGAL_LEDA_BIGFLOAT_H

#include <CGAL/number_type_basic.h>

#ifdef CGAL_USE_LEDA

#include <utility>

#include <CGAL/leda_coercion_traits.h>
#include <CGAL/Interval_nt.h>



#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/bigfloat.h>
#else
#include <LEDA/numbers/bigfloat.h>
#endif

#include <CGAL/Float_traits.h>

CGAL_BEGIN_NAMESPACE

template <> class Algebraic_structure_traits< leda_bigfloat >
  : public Algebraic_structure_traits_base< leda_bigfloat,
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false           Is_exact;
    typedef Tag_true            Is_numerical_sensitive;

    class Sqrt
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return CGAL_LEDA_SCOPE::sqrt( x );
        }
    };

    class Kth_root
      : public Binary_function<int, Type, Type> {
      public:
        Type operator()( int k,
                                        const Type& x) const {
            CGAL_precondition_msg(k > 0, "'k' must be positive for k-th roots");
            // heuristic: we ask for as many precision as the argument has
            long d = x.get_significant_length();
            if ( d < 53) // O.K. we want at least double precision
                d = 53;
            return CGAL_LEDA_SCOPE::sqrt_d( x, d, k);
        }
    };

};

template <> class Real_embeddable_traits< leda_bigfloat >
  : public Real_embeddable_traits_base< leda_bigfloat > {
  public:

    class Abs
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return CGAL_LEDA_SCOPE::abs( x );
        }
    };

    class Sign
      : public Unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
          return (::CGAL::Sign) CGAL_LEDA_SCOPE::sign( x );
        }
    };

    class Compare
      : public Binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) CGAL_LEDA_SCOPE::compare( x, y );
        }
    };

    class To_double
      : public Unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return x.to_double();
        }
    };

    class To_interval
      : public Unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {

          // assuming leda_bigfloat guarantee 1 bit error max
          Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
          Interval_nt_advanced approx (CGAL_LEDA_SCOPE::to_double(x));
          FPU_set_cw(CGAL_FE_UPWARD);
          approx += Interval_nt<false>::smallest();
          return approx.pair();
        }
    };

    class Is_finite
      : public Unary_function< Type, bool > {
      public:
        bool operator()( const Type& x )  const {
          return !( CGAL_LEDA_SCOPE::isInf(x) || CGAL_LEDA_SCOPE::isNaN(x) );
        }
    };
};

template<>
class Is_valid< leda_bigfloat >
  : public Unary_function< leda_bigfloat, bool > {
  public :
    bool operator()( const leda_bigfloat& x ) const {
      return !( CGAL_LEDA_SCOPE::isNaN(x) );
    }
};

//**********************NEW**********************

// Specialization of CGAL::Float_traits_kernel for leda::bigfloat

template <> class Float_traits_kernel< ::leda::bigfloat > {
	public:
		typedef ::leda::bigfloat		FT;
		typedef ::leda::integer		  	IT;
		typedef ::CGAL::Float_number_tag	Float_type;
		typedef CGAL::Precision_type     Precision_type;
				
		class Get_precision {
			public:
				typedef CGAL::Precision_type	result_type;
					
				CGAL::Precision_type operator()() const {
					return CGAL::Precision_type( leda::bigfloat::get_precision() );
				}
		};

		class Set_precision {							
			public:
				typedef void 				result_type;
				typedef CGAL::Precision_type first_argument_type;
					
				void operator()( const CGAL::Precision_type& precision ) const {
					leda::bigfloat::set_precision( (leda::sz_t)precision );
				}
		};			

		class Get_rounding_mode {
			public:
				typedef CGAL::Rounding_mode	result_type;
					
				CGAL::Rounding_mode operator()() const {
					return CGAL::Rounding_mode( leda::bigfloat::get_rounding_mode() );
				}
		};

		class Set_rounding_mode {			
			public:
				typedef void				result_type;
				typedef CGAL::Rounding_mode 	first_argument_type;
					
				void operator()( const CGAL::Rounding_mode& rounding_mode ) const {
					leda::bigfloat::set_rounding_mode( (leda::rounding_modes)rounding_mode );
				}
		};
		
	};


// Specialization of CGAL::Float_traits for leda::bigfloat

template <> class Float_traits< ::leda::bigfloat > 
	: public Float_traits_base< Float_traits_kernel< ::leda::bigfloat> >  {
	
	private:
		static Rounding_mode get_rounding_mode() {
			return Rounding_mode( leda::bigfloat::get_rounding_mode() );
		};
		
		static Precision_type get_precision() {
			return Precision_type( leda::bigfloat::get_precision() );
		};
	
	public:			
		typedef Float_traits_kernel< leda::bigfloat > 	FT_kernel;
		typedef FT_kernel::FT							FT;
		typedef FT_kernel::IT							IT;
		typedef FT_kernel::Float_type					Float_type;


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
					


	class Construct 
		: public Float_traits_base< FT_kernel >::Construct {
		public:
			FT operator() ( const IT& m, const IT& e = IT(0) ) {
				return FT( m, e);
			}
	}; // class Constructor
	
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

	class Is_special_value {
		public:				
			bool operator()( const FT& bf ) const {
				return leda::isSpecial( bf );
			}
	};
	
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
	
	class Is_zero {
		public:
			bool operator()( const FT& bf ) const {
				return leda::isZero( bf ) || bf == FT(0);
			}
	};

	class Is_inf {
		public:
			bool operator()( const FT& bf ) const {
				return leda::isInf( bf );
			}
	};
			
	class Get_mantissa 
		: public Float_traits_base< FT_kernel >::Get_mantissa {
		public:
			//! type for the \c AdaptableUnaryFunction concept.
			typedef IT result_type;
			//! type for the \c AdaptableUnaryFunction concept.
			typedef FT argument_type;

			IT operator()( const FT& a ) const {
				return a.get_significant();
			}
	}; // class GetMantissa

	class Get_exponent 
		: public Float_traits_base< FT_kernel >::Get_exponent {
		public:
			//! type for the \c AdaptableUnaryFunction concept.
			typedef long result_type;
			//! type for the \c AdaptableUnaryFunction concept.
			typedef FT argument_type;

			long operator()( const FT& a ) const {
								return a.get_exponent().to_long();
			}
	}; // class GetMantissa

	class Get_mantissa_and_exponent 
		: public Float_traits_base< FT_kernel >::Get_mantissa_and_exponent {
		public:
			//! type for the \c AdaptableUnaryFunction concept.
			typedef std::pair< IT, IT > result_type;
			//! type for the \c AdaptableUnaryFunction concept.
			typedef FT argument_type;

			std::pair< IT, IT > operator()( const FT& a ) const {
				IT m = a.get_significant();
				IT e = a.get_exponent();
				std::pair< IT, IT > p(  m, e );

				return p;
			}
	}; // class GetMantissaAndExponent
	
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


//TODO
//porting Round_absolute and NiX::Float_number from EXACUS2CGAL

// 	class Round_absolute : public Float_traits_base< FT_kernel >::Round_absolute {
// 		public:
// 			//! type for the \c AdaptableTernaryFunction concept.
// 			typedef FT result_type;
// 			//! type for the \c AdaptableTernaryFunction concept.
// 			typedef FT first_argument_type;
// 			//! type for the \c AdaptableTernaryFunction concept.
// 			typedef IT second_argument_type;
// 			//! type for the \c AdaptableTernaryFunction concept.
// 			typedef CGAL::Rounding_mode third_argument_type;
// 			
// 			FT operator()( const FT& a,
// 							const IT& new_exponent,
// 							const CGAL::Rounding_mode& r_mode = get_rounding_mode() ) const {
// 				typedef NiX::Float_number<IT> FN;
// 				FN fn( a.get_significant(), a.get_exponent() );
// 				fn = fn.round_absolute( new_exponent, r_mode );
// 				return FT( fn.get_mantissa(), fn.get_exponent() );
// 			}
// 	};

	// special functors ( only for Big_float )

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
				return FT(a.get_significant(), a.get_exponent()+e);
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
				return FT(a.get_significant(), a.get_exponent()-e);
			}
	}; // class Div_by_pow_of_2
	
	class Floor {
			public:
				typedef IT result_type;
				typedef FT argument_type;
				
				IT operator()( const FT& fn ) const {
					return leda::to_integer( fn, leda::TO_N_INF );
				}
		};
		
	class Ceil {
			public:
				typedef IT result_type;
				typedef FT argument_type;
				
				IT operator()( const FT& fn ) const {
					return leda::to_integer( fn, leda::TO_P_INF );
				}
		};
	

};

//**********************NEW**********************

CGAL_END_NAMESPACE

// Unary + is missing for leda::bigfloat
namespace leda {
    inline bigfloat operator+( const bigfloat& i) { return i; }
} // namespace leda

//since types are included by leda_coercion_traits.h:
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>

#endif // CGAL_USE_LEDA

#endif // CGAL_LEDA_BIGFLOAT_H
