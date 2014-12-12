// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//============================================================================

#ifndef CGAL_CORE_BIGFLOAT_H
#define CGAL_CORE_BIGFLOAT_H

#include <CGAL/basic.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/CORE/BigFloat.h>
#include <CGAL/CORE_coercion_traits.h>
#include <CGAL/Interval_traits.h> 
#include <CGAL/Bigfloat_interval_traits.h> 

namespace CGAL {

// ######### Interval_traits 

template<> 
class Interval_traits<CORE::BigFloat> 
    : public internal::Interval_traits_base<CORE::BigFloat>{
    typedef CORE::BigFloat Interval;
public: 
    typedef Interval_traits<CORE::BigFloat> Self; 
    typedef CORE::BigFloat Type;
    typedef CORE::BigFloat Bound;
    typedef CGAL::Tag_true Is_interval; 
    typedef CGAL::Tag_true Is_bigfloat_interval; 
  
 
    struct Lower :public std::unary_function<Interval,Bound>{
        Bound operator() ( Interval x ) const {   
            CORE::BigFloat result = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
            CGAL_postcondition(result <= x);
            return result; 
        }
    };
    
    struct Upper :public std::unary_function<Interval,Bound>{
        Bound operator() ( Interval x ) const {     
            CORE::BigFloat result = ::CORE::BigFloat(x.m()+x.err(),0,x.exp());
            CGAL_postcondition(result >= x);
            return result; 
        }
    };

    struct Width :public std::unary_function<Interval,Bound>{
         
        Bound operator() ( Interval x ) const {    
            unsigned long err = 2*x.err();
            return Bound(CORE::BigInt(err),0,x.exp());
        }
    };

    struct Median :public std::unary_function<Interval,Bound>{
         
        Bound operator() ( Interval x ) const {   
            return Bound(x.m(),0,x.exp());
        }
    };

    struct Norm :public std::unary_function<Interval,Bound>{
        Bound operator() ( Interval x ) const {
          BOOST_USING_STD_MAX();
          return max BOOST_PREVENT_MACRO_SUBSTITUTION (Upper()(x).abs(),Lower()(x).abs());
        }
    };
    
    struct Zero_in :public std::unary_function<Interval,bool>{
        bool operator() ( Interval x ) const {      
            return x.isZeroIn(); 
        }
    };

    struct In :public std::binary_function<Bound,Interval,bool>{  
        bool operator()( Bound x, const Interval& a ) const {    
            CGAL_precondition(CGAL::singleton(x));
            return (Lower()(a) <= x && x <= Upper()(a));
        }
    };

    struct Equal :public std::binary_function<Interval,Interval,bool>{  
        bool operator()( const Interval& a, const Interval& b ) const { 
            return (Upper()(a) == Upper()(b) &&  Lower()(a) == Lower()(b));
        }
    };
    
    struct Subset :public std::binary_function<Interval,Interval,bool>{  
        bool operator()( const Interval& a, const Interval& b ) const {   
            return Lower()(b) <= Lower()(a) && Upper()(a) <= Upper()(b);
        }
    };
    
    struct Proper_subset :public std::binary_function<Interval,Interval,bool>{ 
        bool operator()( const Interval& a, const Interval& b ) const { 
            return Subset()(a,b) && (!Equal()(a,b));
        }
    };
    
    struct Intersection :public std::binary_function<Interval,Interval,Interval>{ 
      Interval operator()( const Interval& a, const Interval& b ) const {
            BOOST_USING_STD_MAX();
            BOOST_USING_STD_MIN();
            // std::cout <<"a= (" << a.m() << "+-" << a.err() << ")*2^" << a.exp() << std::endl;
            Bound l(max BOOST_PREVENT_MACRO_SUBSTITUTION (Lower()(a),Lower()(b)));
            Bound u(min BOOST_PREVENT_MACRO_SUBSTITUTION (Upper()(a),Upper()(b)));

            if(u < l ) throw Exception_intersection_is_empty();
            return Construct()(l,u);
        }
    };
 

    struct Overlap :public std::binary_function<Interval,Interval,bool>{
        bool operator() ( Interval x, Interval y ) const {       
            Self::Zero_in Zero_in;
            bool result = Zero_in(x-y);
            return result;
        }
    };
   
    struct Hull :public std::binary_function<Interval,Interval,Interval>{

      // for debugging
/*      void print_bf(CORE::BigFloat bf, std::string s) const {
        
        std::cout << s << ".m()=" << bf.m() << ","
                  << s << ".err()=" << bf.err() << ","
                  << s << ".exp()=" << bf.exp() << ","
                  << "td=" << bf << std::endl;
      }
*/

        Interval operator() ( Interval x, Interval y ) const {
            BOOST_USING_STD_MAX();
            BOOST_USING_STD_MIN();
#if 0
            // this is not possible since CORE::centerize has a bug.
            Interval result = CORE::centerize(x,y);
#else 

            //print_bf(x,"x");
            //print_bf(y,"y");
            
             CORE::BigFloat result;
             
            // Unfortunately, CORE::centerize(x,y) has bugs. 
            if ((x.m() == y.m()) && (x.err() == y.err()) && (x.exp() == y.exp())) { 
                return x;
            }
                         
            CORE::BigFloat lower = min BOOST_PREVENT_MACRO_SUBSTITUTION (CGAL::lower(x), CGAL::lower(y));
            CORE::BigFloat upper = max BOOST_PREVENT_MACRO_SUBSTITUTION (CGAL::upper(x), CGAL::upper(y));

            CORE::BigFloat mid = (lower + upper)/2;
             
            //print_bf(lower,"lower");
            //print_bf(upper,"upper");
            //print_bf(mid,"mid");

            // Now we have to compute the error. The problem is that .err() is just a long
            CORE::BigFloat err = (upper - lower)/CORE::BigFloat(2);
                   
            //print_bf(err,"err");

            //std::cout << "lower    " << lower << std::endl;
            //std::cout << "upper    " << upper << std::endl;
            //std::cout << "mid      " << mid << std::endl;
            //std::cout << "err I    " << err << std::endl;
            
            // shift such that err.m()+err.err() fits into long 
            int digits_long = std::numeric_limits<long>::digits;
            if(::CORE::bitLength(err.m()+err.err()) >= digits_long){
                long shift = ::CORE::bitLength(err.m()) - digits_long + 1 ; 
                //std::cout << "shift " << shift<< std::endl;
                long new_err = ((err.m()+err.err()) >> shift).longValue()+1; 
                err = CORE::BigFloat(0,new_err,0) * CORE::BigFloat::exp2(err.exp()*CORE::CHUNK_BIT+shift);
            }else{           
                err = CORE::BigFloat(0,err.m().longValue()+err.err(),err.exp());
            }
            //print_bf(err,"new_err");

            // TODO: This is a workaround for a bug in operator+ 
            // of CORE::Bigfloat. If the exponent difference is too big,
            // this might cause problems, since the error is a long
            if(mid.exp() > err.exp()) {
                long mid_err = mid.err();
                CORE::BigInt mid_m = mid.m();
                mid_err = mid_err << (mid.exp()-err.exp())*CORE::CHUNK_BIT;
                mid_m = mid_m << (mid.exp()-err.exp())*CORE::CHUNK_BIT;
                mid = CORE::BigFloat(mid_m,mid_err,err.exp());
                //print_bf(mid,"corr_mid");
            }
            
            //print_bf(result,"result");        

            result = mid + err;  
             
#endif 

            CGAL_postcondition( 
                    CGAL::lower(result) 
                    <=  min BOOST_PREVENT_MACRO_SUBSTITUTION (CGAL::lower(x), CGAL::lower(y)));
            CGAL_postcondition( 
                    CGAL::upper(result) 
                    >= max BOOST_PREVENT_MACRO_SUBSTITUTION (CGAL::upper(x), CGAL::upper(y)));

            

            return result ;
        }
    };

    struct Singleton :public std::unary_function<Interval,bool> {
        bool operator() ( Interval x ) const {       
            return (x.err() == 0); 
        }
    };

    struct Construct :public std::binary_function<Bound,Bound,Interval>{
        Interval operator()( const Bound& l,const Bound& r) const {
            CGAL_precondition( l < r ); 
            return Hull()(l,r);
        }
    };
};


// ########### Bigfloat_interval_traits 


// template<typename BFI> long relative_precision(BFI bfi);
namespace internal{

CORE::BigFloat 
inline 
round(const CORE::BigFloat& x, long rel_prec = CORE::defRelPrec.toLong() ){
    CGAL_postcondition(rel_prec >= 0);   

    // since there is not rel prec defined if Zero_in(x)
    if (x.isZeroIn()) return x; 
    // if (CGAL::get_significant_bits(x) <= rel_prec) return x;
   
// if 1 
//    CORE::BigFloat xr;
//    xr.approx(x,rel_prec,1024);
//    typedef CORE::BigFloat BF; 
// else       
    typedef CORE::BigFloat BF; 
    BF xr;
   
    CORE::BigInt m = x.m();
    long         err = x.err();
    long         exp = x.exp(); 
   

//    std::cout <<"(" << m << "+-" <<err << ")*2^"<<(CORE::CHUNK_BIT*exp) << std::endl; 
//    if (err != 0) 
//      std::cout <<"current prec: " <<  CGAL::relative_precision(x) << std::endl;
//    else 
//      std::cout <<"current prec: " << " SINGLETON " << std::endl;
//    std::cout <<"desired prec: " << rel_prec << std::endl; 
//    std::cout <<"bitLength: " << CORE::bitLength(m) << std::endl; 
//    long shift = ::CORE::bitLength(m) - rel_prec - 1;
   
    long shift ;
    if (err == 0)
      shift = ::CORE::bitLength(m) - rel_prec - 3;
    else      
      shift = CGAL::relative_precision(x) - rel_prec -1; 
    
    if( shift > 0 ){    
      m   >>= shift ; 
      err >>= shift; 
      xr = BF(m,err+1,0)*BF::exp2(exp*CORE::CHUNK_BIT+shift);     
    }else{    // noting to do
        xr = x; 
    }

//    std::cout <<"(" <<m << "+-" <<err+1 << ")*2^"<<(CORE::CHUNK_BIT*exp) << std::endl; 
//    if (xr.err() != 0) 
//      std::cout <<"current prec: " <<  CGAL::relative_precision(xr) << std::endl;
//    else 
//      std::cout <<"current prec: " << " SINGLETON "<< std::endl;
//    std::cout <<"desired prec: " << rel_prec << std::endl; 
    
// endif     
    CGAL_postcondition(singleton(xr) || CGAL::relative_precision(xr) - rel_prec >= 0); 
    CGAL_postcondition(singleton(xr) || CGAL::relative_precision(xr) - rel_prec <= 32);   
    CGAL_postcondition(BF(xr.m()-xr.err(),0,xr.exp()) <= BF(x.m()-x.err(),0,x.exp()));
    CGAL_postcondition(BF(xr.m()+xr.err(),0,xr.exp()) >= BF(x.m()+x.err(),0,x.exp()));
    return xr;     
}
}

template<> class Bigfloat_interval_traits<CORE::BigFloat> 
:public Interval_traits<CORE::BigFloat>
{

    typedef CORE::BigFloat NT;
    typedef CORE::BigFloat BF;
public:
  typedef Bigfloat_interval_traits<NT> Self;
  
   struct Relative_precision {
        // type for the \c AdaptableUnaryFunction concept.
        typedef NT  argument_type;
        // type for the \c AdaptableUnaryFunction concept.
        typedef long  result_type;

        long operator()( NT x) const { 
          CGAL_precondition(!Singleton()(x));
          CGAL_precondition(!CGAL::zero_in(x));
          
          x = x.abs();
          NT w = Width()(x);
          w /= ::CORE::BigFloat(x.m()-x.err(),0,x.exp());    
          w = w.abs();
          return -(CORE::ceilLg(w.m()+w.err())+w.exp()*CORE::CHUNK_BIT);
        }
    };
       
    struct Set_precision {
        // type for the \c AdaptableUnaryFunction concept.
        typedef long  argument_type;
        // type for the \c AdaptableUnaryFunction concept.
        typedef long  result_type;  
     
        long operator() ( long prec ) const {    
            long result =  ::CORE::defRelPrec.toLong();
            ::CORE::defRelPrec = prec; 
            ::CORE::defBFdivRelPrec = prec;
            return result; 
        }
    };
     
    struct Get_precision {
        // type for the \c AdaptableGenerator concept.
        typedef long  result_type;  
     
        long operator() () const {
            return  ::CORE::defRelPrec.toLong(); 
        }
    };
};




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
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            // What I want is a sqrt computed with ::CORE::defRelPrec bits.
            // And not ::CORE::defBFsqrtAbsPrec as CORE does. 
            
            CGAL_precondition(::CORE::defRelPrec.toLong() > 0);
            CGAL_precondition(x > 0);
            
            Type a = CGAL::internal::round(x, ::CORE::defRelPrec.toLong()*2);
            CGAL_postcondition(a > 0); 

            Type tmp1 = 
                CORE::BigFloat(a.m(),0,0).sqrt(::CORE::defRelPrec.toLong());
            Type err  =  
                Type(0,long(std::sqrt(double(a.err()))),0) 
                * CORE::BigFloat::exp2(a.exp()*7);
            Type result = tmp1*CORE::BigFloat::exp2(a.exp()*7) + err;
           
            CGAL_postcondition(result >= 0);
            CGAL_postcondition(CGAL::lower(result*result) <= CGAL::lower(x));
            CGAL_postcondition(CGAL::upper(result*result) >= CGAL::upper(x));

            return result;
        }
    };

    class Kth_root
      : public std::binary_function<int, Type, Type> {
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
  : public INTERN_RET::Real_embeddable_traits_base< CORE::BigFloat , CGAL::Tag_true  > {
  public:
    class Abs
      : public std::unary_function< Type, Type > {
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

    class Sgn
      : public std::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
            ::CGAL::Sign result =  sign( x.sign());
            return result; 
        }
    };

    class Compare
      : public std::binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) sign( (x-y).sign());
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type, 
                Comparison_result )
    };

    class To_double
      : public std::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          // this call is required to get reasonable values for the double
          // approximation
          return x.doubleValue();
        }
    };

    class To_interval
      : public std::unary_function< Type, std::pair< double, double > > {
    public:
        std::pair<double, double> operator()( const Type& x ) const {
                        
            double lb,ub;
           
            Type x_lower = CGAL::lower(CGAL::internal::round(CGAL::lower(x),50));
            Type x_upper = CGAL::upper(CGAL::internal::round(CGAL::upper(x),50));
            
            // since matissa has 50 bits only, conversion to double is exact 
            lb = x_lower.doubleValue();
            CGAL_postcondition(lb == x_lower);
            ub = x_upper.doubleValue();
            CGAL_postcondition(ub == x_upper);             
            
            std::pair<double, double> result(lb,ub);
            CGAL_postcondition( result.first  <=  CORE::Expr(CGAL::lower(x)));
            CGAL_postcondition( result.second >=  CORE::Expr(CGAL::upper(x)));
            return result;      
        }
    };
};

} //namespace CGAL

//since types are included by CORE_coercion_traits.h:
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_arithmetic_kernel.h>

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CORE::BigFloat>
  {
    typedef CORE::BigFloat Real;
    typedef CORE::BigFloat NonInteger;
    typedef CORE::BigFloat Nested;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 60,
      MulCost = 60
    };
  };
}

#endif // CGAL_CORE_BIGFLOAT_H
