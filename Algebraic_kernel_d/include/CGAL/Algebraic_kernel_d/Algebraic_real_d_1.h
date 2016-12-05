// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================

// The comments are all original EXACUS comments and aren't adapted.

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/tss.h>

#include <iterator>
#include <list>
#include <vector>
#include <queue>

namespace CGAL {
namespace internal {
template <class Coefficient_,  class Rational_, class HandlePolicy , class AlgebraicRealRep_d_1 >
class Algebraic_real_d_1;
}

template <
class Coefficient_, class Rational_, class HandlePolicy, class AlgebraicRealRep_d_1 >
typename Get_arithmetic_kernel<Coefficient_>::Arithmetic_kernel::Bigfloat_interval
inline convert_to_bfi(
    const internal::Algebraic_real_d_1< Coefficient_, Rational_, HandlePolicy, AlgebraicRealRep_d_1 >& x);

namespace internal {
/*! \ingroup NiX_Algebraic_real
  \brief  An Algebraic_real_d_1 \a x is represented by a polynomial and an
  isolating interval. It is guaranteed that the polynomial is square free.
  The isolating interval is an open interval and it is guaranteed that the
  polynomial is non-zero at the endpoints of the interval.

  The algebraic real are reference counted by default.
  The template parameters are:
  - \b Coefficient_: a model of the \c IntegralDomain concept.
  - \b Rational_: a model of the \c Field concept.
  Must be assignable to \c Field_with_sqrt_.
  - \b HandlePolicy: a model of the \c HandlePolicy concept or the
  \c Handle_policy_in_place class template that selects a specialized
  implementation without reference counting. Has the
  default \c Handle_policy_union.

  THIS CLASS IS CONSIDERED AS EXPERIMENTAL !
*/
template < class Coefficient_, 
           class Rational_,
           class HandlePolicy = ::CGAL::Handle_policy_no_union,
           class AlgebraicRealRep_d_1 = internal::Algebraic_real_rep< Coefficient_, Rational_ > >
class Algebraic_real_d_1 :
    public ::CGAL::Handle_with_policy< AlgebraicRealRep_d_1, HandlePolicy > {

  // currently Rational is the only supported Bound type. 
  CGAL_static_assertion(
      (   ::boost::is_same <Rational_, 
          typename Get_arithmetic_kernel<Coefficient_>::Arithmetic_kernel::Rational>::value));
    


public :
  typedef ::CGAL::Handle_with_policy<AlgebraicRealRep_d_1,HandlePolicy>                   Base;
  typedef Algebraic_real_d_1<Coefficient_,Rational_, HandlePolicy, AlgebraicRealRep_d_1>  Self;

  typedef Coefficient_                                  Coefficient;
  typedef Rational_                                     Bound; 
  typedef HandlePolicy                                  Handle_policy;
  typedef AlgebraicRealRep_d_1                          Algebraic_real_rep_d_1; 
  typedef typename Algebraic_real_rep_d_1::Polynomial_1 Polynomial_1; 

  // These public typedefs should be removed in the long run
  // typedef typename Algebraic_real_rep_d_1::Polynomial_1 Polynomial; 
  typedef Rational_                                     Rational;
  
private:
  typedef CGAL::Fraction_traits<Rational> FT_rational;
  typedef typename FT_rational::Numerator_type Integer;
  
private: 
  static inline Self& get_default_instance(){
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(Self, x,0); 
    return x; 
  }
public:
  
  //! Default constructor
  Algebraic_real_d_1() : Base(static_cast<const Base&>(get_default_instance())) {}

  //! copy constructor: copy existing Algebraic_real_d_1 (shares rep)
  Algebraic_real_d_1(const Self& p) : Base(static_cast<const Base&>(p)) {}
  
  //! creates the algebraic real from \a i.
  Algebraic_real_d_1(int i ) : Base(Algebraic_real_rep_d_1(i)) { }

  //! creates the algebraic real from \a x.
  explicit Algebraic_real_d_1(const Rational& x) : Base(Algebraic_real_rep_d_1(x)) { }

  /*! \brief creates the algebraic real as the unique root of \a P
   *   in the open interval <var>]low,high[</var>.
   *   \pre P is square free.
   *   \pre P is not zero at low
   *   \pre P is not zero at high
   */
  Algebraic_real_d_1(const Polynomial_1& P, Rational low, Rational high)
    : Base (Algebraic_real_rep_d_1(P,low,high)) {}

  //! returns the polynomial defining \a x
  const Polynomial_1& polynomial() const { return this->ptr()->polynomial(); }

  /*! \brief returns the degree of the polynomial defining \a x
   *  This is not necessarily the algebraic degree of \a x, since
   *  the polynomial may be reducible.
   */
  int degree() const {
    return CGAL::degree(this->ptr()->polynomial());
  }
  
  //! returns the lower endpoint of the isolating interval
  Rational low() const { return this->ptr()->low(); }
  Rational lower() const { return this->ptr()->low(); }

  //! returns the upper endpoint of the isolating interval
  Rational high() const { return this->ptr()->high(); }
  Rational upper() const { return this->ptr()->high(); }

  /*! \brief returns the sign of the defining polynomial
   *  at the lower endpoint of the isolating interval
   */
  int sign_at_low() const { return this->ptr()->sign_at_low(); }

  /*! \brief returns whether a rational representation
   * (as Rational) is known.
   */
  bool is_rational() const { return this->ptr()->is_rational(); }

  //! returns a rational representation of \a x
  /*! \pre: type() == NiX::IS_Rational_IONAL */
  Rational rational() const {
    CGAL_precondition(is_rational());
    return this->ptr()->rational();
  }

  //! compute a \c double approximation (without guarantees)
  double to_double() const {
    return (to_interval().first +to_interval().second)/2;
  }

  /*! \brief returns a double Interval approximation,
   *  it is guaranteed that \a x is contained in the Interval.
   *
   *  This is not the isolating interval.
   */
  std::pair<double, double> to_interval() const {
    if( this->ptr()->interval_option ) {
      return *(this->ptr()->interval_option);
    } else {
      typedef typename Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Bigfloat_interval BFI;
      long old_precision = get_precision( BFI() );
      set_precision( BFI(), 53 );
      std::pair<double, double> interval = CGAL::to_interval( convert_to_bfi( (*this)));
      this->ptr()->interval_option = boost::optional< std::pair<double, double> >(interval);
      set_precision( BFI(), old_precision );
      return *(this->ptr()->interval_option);
    }
  }

  /*! \brief Refines the isolating interval. */
  void refine() const{ this->ptr()->refine(); }

  /*! \brief Bisects the isolating interval. */
  void bisect() const{ this->ptr()->bisect(); }

  /*! \brief Refines the isolating interval until \a m is outside
   *  the \c closed interval
   */
  template < class NTX >
  void strong_refine(const NTX& m) const{ compare(m); }

  /*! \brief compares \a x with respect to \a y
   *  It returns CGAL::SMALLER if \a x is smaller than <var>y</var>,
   *  CGAL::EQUAL if they are equal and CGAL::LARGER otherwise.
   */
  template < class NTX >
  CGAL::Comparison_result
  compare(const NTX& y) const { return intern_compare(y,false); }


  /*! \brief compares \a x with respect to \a y
   *  \pre \a x != \a y
   */
  template < class NTX >
  CGAL::Comparison_result
  compare_distinct(const NTX& y) const { return intern_compare(y,true); }

  template <class Algebraic_real_iterator>
  static void conjugate(Algebraic_real_iterator begin,
      Algebraic_real_iterator end){
    if(begin == end) return;

    CGAL_precondition(begin->ptr()->next==begin->ptr());
    CGAL_precondition(begin->ptr()->prev==begin->ptr());

    for(Algebraic_real_iterator it = begin; it != end; ++it){
      CGAL_precondition(it->ptr()->next==it->ptr());
      CGAL_precondition(it->ptr()->prev==it->ptr());
      CGAL_precondition(it->polynomial()==begin->polynomial());
      it->ptr()->next         =begin->ptr();
      it->ptr()->prev         =begin->ptr()->prev;
      begin->ptr()->prev->next=it->ptr();
      begin->ptr()->prev      =it->ptr();
    }
  }


private:
  CGAL::Comparison_result
  intern_compare(const Self& y, bool distinct) const {
    if(distinct)
      CGAL_precondition(this->ptr()!=y.ptr());
    else if(this->ptr()==y.ptr()) {
      return CGAL::EQUAL;
    }
    CGAL::Comparison_result result =
      this->ptr()->compare(*(y.ptr()), distinct);
    if (result == CGAL::EQUAL) {
      this->unify(y);
    }
    return result;
  }

  CGAL::Comparison_result
  intern_compare(const Rational& y, bool distinct) const {
    return this->ptr()->compare(y,distinct);
  }

  template < class NTX >
  CGAL::Comparison_result
  intern_compare(const NTX& y, bool distinct) const {
    if(y           < NTX(low())) return CGAL::LARGER;
    if(NTX(high()) < y         ) return CGAL::SMALLER;

    if(!distinct){
      typename Real_embeddable_traits<NTX>::To_interval to_interval;
      if( CGAL::Interval_nt<true>( this->to_interval() ).do_overlap(
              CGAL::Interval_nt<true>(to_interval(y)) ) ) {
        if(polynomial().sign_at(y)==CGAL::ZERO)
          return CGAL::EQUAL;
      }
    }

    while(NTX(low()) <= y && y <= NTX(high())) {
      refine();
    }
    if(y            <  NTX(low())) return CGAL::LARGER;
    CGAL_assertion(NTX(high())  <  y); return CGAL::SMALLER;
  }

public:
  //! returns if y is contained in the \c closed isolating interval
  template < class NTX >
  bool contains(const NTX& y) const {
    return ((NTX(low()) <= y) && (y <= NTX(high())));
  }

  //! return if \a x is a root of Q
  bool is_root_of(const Polynomial_1& Q) const {return this->ptr()->is_root_of(Q); }

  /*! \brief returns a rational (Rational) between this number and \c y.
   *  \pre x != y
   */
  Rational rational_between (const Self& y) const{
    CGAL::Comparison_result s = compare(y);
    CGAL_precondition(s != CGAL::EQUAL);
    if(s == CGAL::SMALLER){
      Rational r((high()+y.low())/Rational(2));
      CGAL::simplify(r);
      return r;
    }else{
      Rational r((y.high()+low())/Rational(2));
      CGAL::simplify(r);
      return r;
    }
  }

public:
  /*! \brief refines the isolating interval to ]<var>lo</var>,
   *  <var>hi</var>[.
   *
   *  This function can be used to inform an Algebraic_real_d_1 \a x of an
   *  externally refined isolating interval. Its arguments must be
   *  the boundaries of an isolating interval that contains \a x in its
   *  interior. (Use other functions like \c .strong_refine() in case you
   *  want to communicate an explicit value.)
   */
  void refine_to(const Rational& lo, const Rational& hi) const {
    // test whether lo < x < hi
    // and refines isolating interval until in ]lo,hi[
    CGAL_assertion_code(CGAL::Comparison_result s =) compare_distinct(lo);
    CGAL_assertion(CGAL::LARGER == s);
    CGAL_assertion_code(s =) compare_distinct(hi);
    CGAL_assertion(CGAL::SMALLER == s);
  }


public:
  template <class NTX>
  bool operator==( const NTX& y) const {return compare(y)==CGAL::EQUAL;}
  template <class NTX>
  bool operator!=( const NTX& y) const {return compare(y)!=CGAL::EQUAL;}
  template <class NTX>
  bool operator< ( const NTX& y) const {return compare(y)==CGAL::SMALLER;}
  template <class NTX>
  bool operator> ( const NTX& y) const {return compare(y)==CGAL::LARGER;}
  template <class NTX>
  bool operator<=( const NTX& y) const {return compare(y)!=CGAL::LARGER;}
  template <class NTX>
  bool operator>=( const NTX& y) const {return compare(y)!=CGAL::SMALLER;}

  //! unary operator +
  const Self& operator+() const { return *this; }

  //! unary operator -
  Self operator-() const {
    Polynomial_1 P(polynomial());
    P.scale_up(Coefficient(-1));
    Rational high_(-low());
    Rational low_ (-high());
    return Self(P,low_,high_);
  }

  //! Simplifies the algebraic number
  void simplify() const {
    this->ptr()->simplify();
  }

}; // class Algebraic_real_d_1

} // namespace internal


//----------------------------------------------------------

/*! \ingroup NiX_Algebraic_real_d_1
 *  \ingroup NiX_NT_traits_spec
 *  \brief NT_traits class for NiX::Algebraic_real_d_1, which is a model of the
 *  RealComparable concept.
 *
 *  NiX::Algebraic_real_d_1 does not support any arithmetic operations, thus they
 *  are not even a model of the IntegralDomainWithoutDiv concept. \see NiX_NT_Concepts
 */
template< class Coefficient, class Rational, class HandlePolicy, class RepClass >
class Real_embeddable_traits< internal::Algebraic_real_d_1< Coefficient, Rational, HandlePolicy, RepClass > >
  : public INTERN_RET::Real_embeddable_traits_base< internal::Algebraic_real_d_1< Coefficient, Rational, HandlePolicy, RepClass > , CGAL::Tag_true > {

public:

  typedef internal::Algebraic_real_d_1< Coefficient, Rational, HandlePolicy, RepClass > Type;

  class Compare
    : public std::binary_function< Type, Type, CGAL::Comparison_result > {
  public:
    CGAL::Comparison_result operator()( const Type& a, const Type& b ) const
    { return a.compare( b ); }
    CGAL::Comparison_result operator()( int a, const Type& b ) const
    { return - b.compare( Rational(a) ); }
    CGAL::Comparison_result operator()( const Type& a, int b ) const
    { return a.compare( Rational(b) ); }
    CGAL::Comparison_result operator()( const Rational&  a, const Type& b ) const
    { return - b.compare( a ); }
    CGAL::Comparison_result operator()( const Type& a, const Rational&  b ) const
    { return a.compare( b ); }
    CGAL::Comparison_result operator()(
        const typename First_if_different<Coefficient,Rational>::Type&  a,
        const Type& b ) const
    {
      typename Coercion_traits<Coefficient,Type>::Cast cast;
      return cast(a).compare(b);
    }
    CGAL::Comparison_result operator()(
        const Type& a,
        const typename First_if_different<Coefficient,Rational>::Type&  b ) const
    {
      typename Coercion_traits<Coefficient,Type>::Cast cast;
      return a.compare(cast(b));
    }
  };

  class Sgn
    : public std::unary_function< Type, CGAL::Sign > {
  public:
    CGAL::Sign operator()( const Type& a ) const {
      return a.compare( Rational(0) );
    }
  };

  class To_double
    : public std::unary_function< Type, double > {
  public:
    double operator()(const Type& a) const {
      return a.to_double();
    }
  };

  class To_interval
    : public std::unary_function< Type, std::pair<double, double> > {
  public:
    typename std::pair<double, double> operator()(const Type& a) const {
      return a.to_interval();
    }
  };
};



/*! \relates NiX::Algebraic_real_d_1
 *  \brief outputs \c x to \c os
 */
template<class Coefficient, class Rational, class HandlePolicy, class RepClass >
std::ostream&
operator << (std::ostream& os,
    const CGAL::internal::Algebraic_real_d_1<Coefficient, Rational, HandlePolicy, RepClass >& x){
  os << "["   << x.polynomial() 
     << ",["  << oformat(x.low()) 
     << " , " << oformat(x.high()) << " ]]";
  return os;
}

/*! \relates NiX::Algebraic_real_d_1
 *  \brief read an NiX::Algebraic_real_d_1 from \c is into \c x.
 */
template<class Coefficient, class Rational, class HandlePolicy, class RepClass>
std::istream&
operator >> (std::istream& is,
    CGAL::internal::Algebraic_real_d_1<Coefficient, Rational, HandlePolicy, RepClass>& x){
  
  typedef  CGAL::internal::Algebraic_real_d_1<Coefficient, Rational, HandlePolicy, RepClass > ALGNUM;
 
  Rational low, high;
  typename CGAL::Polynomial_type_generator<Coefficient,1>::Type poly;
  
  swallow(is, '[');// read the "["
  is >> poly;
  swallow(is, ',');// read the ","
  swallow(is, '[');// read the ","
  is >> iformat(low);
  swallow(is, ',');// read the ","
  is >> iformat(high);
  swallow(is, ']');// read the "]"
  swallow(is, ']');// read the "]"
  x = ALGNUM(poly, low, high);
  return is;
}

template <class Coefficient_, class Rational_, class HandlePolicy, class AlgebraicRealRep_d_1 >
typename Get_arithmetic_kernel<Coefficient_>::Arithmetic_kernel::Bigfloat_interval
inline
convert_to_bfi(const internal::Algebraic_real_d_1< Coefficient_, Rational_, HandlePolicy, AlgebraicRealRep_d_1 >& x){
  typedef typename Get_arithmetic_kernel<Coefficient_>::Arithmetic_kernel AT;
  typedef typename AT::Bigfloat_interval BFI;

  if (x.is_rational()) return convert_to_bfi(x.rational());

  if(CGAL::sign(x) == CGAL::ZERO) return (BFI(0));

  CGAL_postcondition(CGAL::sign(x.low()) == CGAL::sign(x.high()));
  long final_prec = set_precision( BFI(),get_precision( BFI())+4);

  BFI bfi = CGAL::hull(convert_to_bfi(x.low()), convert_to_bfi(x.high()));

  while( !singleton(bfi) &&  get_significant_bits(bfi) < final_prec  ){
    x.refine();
    bfi = CGAL::hull(
        convert_to_bfi(x.low()),
        convert_to_bfi(x.high()));
  }

  set_precision(BFI(),final_prec);
  return bfi;
}


template <class Coefficient,class Rational,class Handle_policy,class Rep_class>
struct Get_arithmetic_kernel<CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> >{

  typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel
  Arithmetic_kernel;

};

template <class Coefficient,class Rational,class Handle_policy,class Rep_class>
inline CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class>
min BOOST_PREVENT_MACRO_SUBSTITUTION(
    const CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class>& x,
    const CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class>& y){
  return (x<=y)?x:y;
}
template <class Coefficient,class Rational,class Handle_policy,class Rep_class>
inline CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class>
max BOOST_PREVENT_MACRO_SUBSTITUTION(
    const CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class>& x,
    const CGAL::internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class>& y){
  return (x>=y)?x:y;
}

template <class Coefficient, class Rational, class Handle_policy, class Rep_class>
struct Coercion_traits<int, internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> >{
  typedef internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> Type;
  typedef Tag_true Are_explicit_interoperable;
  typedef Tag_false Are_implicit_interoperable;
  struct Cast{
    typedef Type result_type;
    Type operator()(const Type& x) const { return x; }
    Type operator()(const int& x) const { return Type(x); }
  };
};

template <class A, class B, class C, class D>
struct Coercion_traits<internal::Algebraic_real_d_1<A,B,C,D>,int >
  :public Coercion_traits<int,internal::Algebraic_real_d_1<A,B,C,D> >{};


template <class Coefficient, class Rational, class Handle_policy, class Rep_class>
struct Coercion_traits<Rational, internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> >{
  typedef internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> Type;
  typedef Tag_true Are_explicit_interoperable;
  typedef Tag_false Are_implicit_interoperable;
  struct Cast{
    typedef Type result_type;
    Type operator()(const Type& x) const { return x; }
    Type operator()(const Rational& x) const { return Type(x); }
  };
};

template <class A, class Rational, class C, class D>
struct Coercion_traits<internal::Algebraic_real_d_1<A,Rational,C,D>, Rational >
  :public Coercion_traits<Rational,internal::Algebraic_real_d_1<A,Rational,C,D> >{};


template <class Coefficient, class Rational, class Handle_policy, class Rep_class>
struct Coercion_traits<
  typename First_if_different<Coefficient,Rational>::Type,
  internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> >
{
  typedef internal::Algebraic_real_d_1<Coefficient,Rational,Handle_policy,Rep_class> Type;
  typedef Tag_true Are_explicit_interoperable;
  typedef Tag_false Are_implicit_interoperable;
  typedef Coercion_traits<Coefficient,Rational> CTCR;
  struct Cast{
  private:
    Type operator()(const Coefficient& a, Tag_true) const {
      return Type(typename CTCR::Cast()(a));
    }
    Type operator()(const Coefficient& a, Tag_false) const {
      typedef typename Type::Polynomial_1                  Poly;
      typedef CGAL::Polynomial_traits_d<Poly>            PT;
      Poly p = typename PT::Shift()(Poly(1),1) - a ; // (x-a)
      Rational b(2);
      Coefficient aa = CGAL::abs(a);
      while( CGAL::compare(b,aa) != LARGER ){b*=2;};
      return Type(p,-b,b);
    }
  public:
    typedef Type result_type;
    Type operator()(const Type& a) const { return a; }
    Type operator()(const Coefficient& a) const {
      static const bool b = boost::is_same<Rational,typename CTCR::Type>::value;
      return (*this)(a,Boolean_tag<b>());
    }
  };
};

template <class Coefficient, class Rational, class C, class D>
struct Coercion_traits<
  internal::Algebraic_real_d_1<Coefficient,Rational,C,D>,
  typename First_if_different<Coefficient,Rational>::Type>
  :public Coercion_traits<Coefficient,internal::Algebraic_real_d_1<Coefficient,Rational,C,D> >{};

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H

// EOF
