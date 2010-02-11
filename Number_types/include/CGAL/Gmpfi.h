// Copyright (c) 2007-2009 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@loria.fr>

#ifndef CGAL_GMPFI_H
#define CGAL_GMPFI_H

#include <CGAL/number_type_basic.h>
#include <CGAL/mpfi_coercion_traits.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/Bigfloat_interval_traits.h>

namespace CGAL{

template <>
struct Algebraic_structure_traits<Gmpfi>:
public Algebraic_structure_traits_base<Gmpfi,Field_with_kth_root_tag>{

        typedef Tag_false       Is_exact;
        typedef Tag_true        Is_numerical_sensitive;
        typedef Uncertain<bool> Boolean;

        struct Is_zero:
        public std::unary_function<Type,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_zero();
                }
        };

        struct Is_one:
        public std::unary_function<Type,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_one();
                }
        };

        struct Square:
        public std::unary_function<Type,Type>{
                Type operator()(const Type &x)const{
                        return x.square();
                };
        };

        struct Is_square:
        public std::binary_function<Type,Type&,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_square();
                };
                Boolean operator()(const Type &x,Type &y)const{
                        return x.is_square(y);
                };
        };

        struct Sqrt:
        public std::unary_function<Type,Type>{
                Type operator()(const Type &x)const{
                        return x.sqrt();
                };
        };

        struct Kth_Root:
        public std::binary_function<int,Type,Type>{
                Type operator()(int k,const Type &x)const{
                        return (k==3?x.cbrt():x.kthroot(k));
                };
        };

        struct Divides:
        public std::binary_function<Type,Type,Boolean>{
                Boolean operator()(const Type &d,const Type &n)const{
                        return !(d.is_zero());
                };
                Boolean operator()(const Type &d,const Type &n,Type &c)const{
                        return d.divides(n,c);
                };
        };
};

template <>
class Real_embeddable_traits<Gmpfi>:
public INTERN_RET::Real_embeddable_traits_base<Gmpfi,CGAL::Tag_true>{

        typedef Algebraic_structure_traits<Type>        AST;

        public:

        typedef Tag_true                        Is_real_embeddable;
        typedef Uncertain<bool>                 Boolean;
        typedef Uncertain<Comparison_result>    Comparison_result;
        typedef Uncertain<CGAL::Sign>           Sign;
 
        typedef AST::Is_zero    Is_zero;

        struct Is_finite:
        public std::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return(x.is_number());
                };
        };

        struct Abs:
        public std::unary_function<Type,Type>{
                inline Type operator()(const Type &x)const{
                        return x.abs();
                };
        };

        struct Sgn:
        public std::unary_function<Type,Sign>{
                inline Sign operator()(const Type &x)const{
                        return x.sign();
                };
        };

        struct Is_positive:
        public std::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return x.is_positive();
                };
        };

        struct Is_negative:
        public std::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return x.is_negative();
                };
        };

        struct Compare:
        public std::binary_function<Type,Type,Comparison_result>{
                inline Comparison_result operator()
                        (const Type &x,const Type &y)const{
                                return x.compare(y);
                        };
        };

        struct To_double:
        public std::unary_function<Type,double>{
                inline double operator()(const Type &x)const{
                        return x.to_double();
                };
        };

        struct To_interval:
        public std::unary_function<Type,std::pair<double,double> >{
                inline std::pair<double,double> operator()(const Type &x)const{
                                return x.to_interval();
                        };
                };
};



template<>
class Interval_traits<Gmpfi>
{
public: 
  typedef Interval_traits<Gmpfi> Self; 
  typedef Gmpfi Interval; 
  typedef CGAL::Gmpfr Bound; 
  typedef CGAL::Tag_false With_empty_interval; 
  typedef CGAL::Tag_true  Is_interval; 

  struct Construct :public std::binary_function<Bound,Bound,Interval>{
    Interval operator()( const Bound& l,const Bound& r) const {
      CGAL_precondition( l < r ); 
      return Interval(std::make_pair(l,r));
    }
  };

  struct Lower :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.inf();
    }
  };

  struct Upper :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.sup();
    }
  };

  struct Width :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return Gmpfr::sub(a.sup(),a.inf(),std::round_toward_infinity);
    }
  };

  struct Median :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return (a.inf()+a.sup())/2;
    }
  };
    
  struct Norm :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.abs().sup();
    }
  };

  struct Singleton :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return a.inf() == a.sup();
    }
  };

  struct Zero_in :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return a.inf() <= 0  &&  0 <= a.sup();
    }
  };

  struct In :public std::binary_function<Bound,Interval,bool>{
    bool operator()( Bound x, const Interval& a ) const {
      return a.inf() <= x && x <= a.sup();
    }
  };

  struct Equal :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.is_same(b);
    }
  };
    
  struct Overlap :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.do_overlap(b);
    }
  };
    
  struct Subset :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return b.inf() <= a.inf() && a.sup() <= b.sup() ;  
    }
  };
    
  struct Proper_subset :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return Subset()(a,b) && ! Equal()(a,b); 
    }
  };
    
  struct Hull :public std::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      return Interval( std::make_pair(min(a.inf(),b.inf()) , max(a.sup(),b.sup())) ) ;
    }
  };
    
  
  struct Empty :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return false;
    }
  };
  
  struct Intersection :public std::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      CGAL_precondition(a.do_overlap(b));
      return Interval( std::make_pair(max(a.inf(),b.inf()) , min(a.sup(),b.sup()))) ;
    }
  };
};

template<>
class Bigfloat_interval_traits<Gmpfi>:
  public Interval_traits<Gmpfi>
{
public:
  typedef Gmpfi NT;

  typedef CGAL::Gmpfr BF;

  struct Get_significant_bits: public std::unary_function<NT,long>{

    long operator()( NT x) const {
      CGAL_precondition(x.inf() <= x.sup());
      if(CGAL::zero_in(x)) return -1;

      BF w(CGAL::width(x));
      // BF labs(CGAL::lower(CGAL::abs(rrx))) ;
//       std::cout << "get sig " << std::endl; 
//       std::cout << x.inf() << std::endl; 
//       std::cout << x.sup() << std::endl; 
//       std::cout << labs << std::endl; 
//       std::cout << w << std::endl;
//       std::cout << err << std::endl;  

      mpfr_div(w.fr(), w.fr(), CGAL::lower(CGAL::abs(x)).fr(), GMP_RNDU);
      mpfr_log2(w.fr(), w.fr(), GMP_RNDD);
      return -mpfr_get_si(w.fr(), GMP_RNDU);
    }
  };
  
  struct Set_precision {
    // type for the \c AdaptableUnaryFunction concept.
    typedef long  argument_type;
    // type for the \c AdaptableUnaryFunction concept.
    typedef long  result_type;  
     
    long operator()( long prec ) const {
      return Gmpfi::set_default_precision(prec); 
    }
  };
  
  struct Get_precision {
    // type for the \c AdaptableGenerator concept.
    typedef long  result_type;  
    long operator()() const {
      return Gmpfi::get_default_precision(); 
    }
  };

};




} // namespace CGAL

#include <CGAL/GMP/Gmpfi_type.h>

#endif  // CGAL_GMPFI_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
