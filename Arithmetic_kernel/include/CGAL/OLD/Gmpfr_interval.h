// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany), 
// National University of Athens (Greece).   
// Copyright (c) 2009 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL:$
// $Id:$
//
// Author(s)     : George Tzoumas <geotz@di.uoa.gr>,
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================
//
//    \brief provide CGAL support for class CGAL::Gmpfr_interval. 
//

#ifndef CGAL_GMPFR_INTERVAL_H
#define CGAL_GMPFR_INTERVAL_H

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/GMP/Gmpfr_interval_type.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/Bigfloat_interval_traits.h>

namespace CGAL {

template <> class Algebraic_structure_traits< Gmpfr_interval >
  : public Algebraic_structure_traits_base < Gmpfr_interval, 
                                               Field_with_sqrt_tag >{
    
public:
  typedef Tag_false           Is_exact;
  typedef Tag_true            Is_numerical_sensitive;
    
  class Sqrt
    : public std::unary_function< Type, Type > {
  public:
    Type operator()( const Type& x ) const {
      return ::boost::numeric::sqrt(x);
    }
  };
};

template <> class Real_embeddable_traits< Gmpfr_interval >
  : public INTERN_RET::Real_embeddable_traits_base<Gmpfr_interval, CGAL::Tag_true> {
public:
        
  class Abs
    : public std::unary_function< Type, Type > {
  public:
    Type operator()( const Type& x ) const {
      return ::boost::numeric::abs(x);
    }
  };
    
  class To_double
    : public std::unary_function< Type, double > {
  public:
    double operator()( const Type& x ) const {
      return ::boost::numeric::median(x).to_double();
    }
  };

  class To_interval
    : public std::unary_function< Type, std::pair< double, double > > {
  public:
    std::pair<double, double> operator()( const Type& x ) const {            
      std::pair<double, double> lower_I(x.lower().to_interval());
      std::pair<double, double> upper_I(x.upper().to_interval());
      return std::pair< double, double >(
          (CGAL::min)(lower_I.first , upper_I.first ),
          (CGAL::max)(lower_I.second, upper_I.second));
    }
  };
};

template<>
class Interval_traits<Gmpfr_interval>
{
public: 
  typedef Interval_traits<Gmpfr_interval> Self; 
  typedef Gmpfr_interval Interval; 
  typedef CGAL::Gmpfr Bound; 
  typedef CGAL::Tag_true With_empty_interval; 
  typedef CGAL::Tag_true Is_interval; 

  struct Construct :public std::binary_function<Bound,Bound,Interval>{
    Interval operator()( const Bound& l,const Bound& r) const {
      CGAL_precondition( l < r ); 
      return Interval(l,r);
    }
  };

  struct Lower :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.lower();
    }
  };

  struct Upper :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.upper();
    }
  };

  struct Width :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return ::boost::numeric::width(a);
    }
  };

  struct Median :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return ::boost::numeric::median(a);
    }
  };
    
  struct Norm :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return ::boost::numeric::norm(a);
    }
  };

  struct Empty :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return ::boost::numeric::empty(a);
    }
  };

  struct Singleton :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return ::boost::numeric::singleton(a);
    }
  };

  struct Zero_in :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return ::boost::numeric::in_zero(a);
    }
  };

  struct In :public std::binary_function<Bound,Interval,bool>{
    bool operator()( Bound x, const Interval& a ) const {
      return ::boost::numeric::in(x,a);
    }
  };

  struct Equal :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return ::boost::numeric::equal(a,b);
    }
  };
    
  struct Overlap :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return ::boost::numeric::overlap(a,b);
    }
  };
    
  struct Subset :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return ::boost::numeric::subset(a,b);
    }
  };
    
  struct Proper_subset :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return ::boost::numeric::proper_subset(a,b);
    }
  };
    
  struct Hull :public std::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      return ::boost::numeric::hull(a,b);
    }
  };
    
  struct Intersection :public std::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      Interval r = ::boost::numeric::intersect(a,b);
      return r;
    }
  };
};

template<>
class Bigfloat_interval_traits<Gmpfr_interval>:
  public Interval_traits<Gmpfr_interval>
{
public:
  typedef Gmpfr_interval NT;

  typedef CGAL::Gmpfr BF;

  struct Get_significant_bits: public std::unary_function<NT,long>{

    long operator()( NT x) const {
      if(CGAL::zero_in(x)) return -1;
      BF labs = CGAL::lower(CGAL::abs(x)) ;
      BF w = CGAL::width(x);
      BF err;
      mpfr_div(err.fr(), w.fr(), labs.fr(), GMP_RNDU);
      mpfr_log2(err.fr(), err.fr(), GMP_RNDD);
      return -mpfr_get_si(err.fr(), GMP_RNDU);
    }
  };
  
  struct Set_precision {
    // type for the \c AdaptableUnaryFunction concept.
    typedef long  argument_type;
    // type for the \c AdaptableUnaryFunction concept.
    typedef long  result_type;  
     
    long operator()( long prec ) const {
      long old_prec = mpfr_get_default_prec();
//            std::cerr << "precision set to " << prec << " from " << old_prec << std::endl;
      mpfr_set_default_prec(prec); 
      return old_prec;
    }
  };
     
  struct Get_precision {
    // type for the \c AdaptableGenerator concept.
    typedef long  result_type;  
    long operator()() const {
      return mpfr_get_default_prec(); 
    }
  };

};

//Gmp internal coercions:
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpfr_interval)

// The following definitions reflect the interaction of the Gmpfr

// built in types :
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,Gmpfr_interval)
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,Gmpfr_interval)
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,Gmpfr_interval)
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float    ,Gmpfr_interval)
  CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double   ,Gmpfr_interval)


template <> 
struct Coercion_traits<CGAL::Gmpfr_interval, CGAL::Gmpz>{
  typedef Tag_true  Are_explicit_interoperable;
  typedef Tag_false Are_implicit_interoperable;
  typedef CGAL::Gmpfr_interval Type;
    
  struct Cast{
    typedef Type result_type;
    Type operator()(const CGAL::Gmpfr_interval& x)  const { return x;}
    Type operator()(const CGAL::Gmpz x) const {
      CGAL::Gmpfr lower, upper;
      mpfr_set_z (lower.fr(), x.mpz(), GMP_RNDD);
      mpfr_set_z (upper.fr(), x.mpz(), GMP_RNDU);
      Type bfi(lower, upper);
      CGAL_postcondition( bfi.lower() <= x );
      CGAL_postcondition( bfi.upper() >= x );
      return bfi; 
    }
  };
};

template <> // mirror
struct Coercion_traits<CGAL::Gmpz,CGAL::Gmpfr_interval>
  :public Coercion_traits<CGAL::Gmpfr_interval, CGAL::Gmpz>{};

template <> 
struct Coercion_traits<CGAL::Gmpfr_interval, CGAL::Gmpq>{
  typedef Tag_true  Are_explicit_interoperable;
  typedef Tag_false Are_implicit_interoperable;
  typedef CGAL::Gmpfr_interval Type;
  typedef Coercion_traits<CGAL::Gmpfr_interval, CGAL::Gmpz> CTZ;
    
  struct Cast{
    typedef Type result_type;
    Type operator()(const CGAL::Gmpfr_interval& x)  const { return x;}
    Type operator()(const CGAL::Gmpq x) const {
      // early exits  
      if (CGAL::is_zero(x)) return Type(0,0);
      if (CGAL::is_one(x.denominator())){
        return CTZ::Cast()(x.numerator());
      }
      // TODO: ensure that prec is reached for resulting interval ?
      Gmpfr lower, upper, nf, df;
      CGAL::Gmpz num = x.numerator();
      CGAL::Gmpz den = x.denominator();
      mp_prec_t prec = mpfr_get_default_prec();
      CGAL_assertion( mpfr_get_prec(lower.fr()) == prec);
      CGAL_assertion( mpfr_get_prec(upper.fr()) == prec );
      
      mpfr_set_z (nf.fr(), num.mpz(), GMP_RNDD);
      mpfr_set_z (df.fr(), den.mpz(), 
          (CGAL::sign(num) == CGAL::NEGATIVE)? GMP_RNDD: GMP_RNDU);
      mpfr_div(lower.fr(), nf.fr(), df.fr(), GMP_RNDD);

      mpfr_set_z (nf.fr(), num.mpz(), GMP_RNDU);
      mpfr_set_z (df.fr(), den.mpz(), 
          (CGAL::sign(num) == CGAL::NEGATIVE)? GMP_RNDU: GMP_RNDD);
      mpfr_div(upper.fr(), nf.fr(), df.fr(), GMP_RNDU);

      Type bfi(lower, upper);
      
      CGAL_postcondition( bfi.lower() <= x );
      CGAL_postcondition( bfi.upper() >= x );
      return bfi; 
    }
  };
};

template <> // mirror
struct Coercion_traits<CGAL::Gmpq,CGAL::Gmpfr_interval>
  :public Coercion_traits<CGAL::Gmpfr_interval, CGAL::Gmpq>{};


// lower GMP types:
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpfr,Gmpfr_interval)

} //namespace CGAL
#endif //  CGAL_GMPFR_INTERVAL_H
