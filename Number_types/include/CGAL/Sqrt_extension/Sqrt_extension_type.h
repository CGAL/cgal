// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>
//                 Ron Wein         <wein@post.tau.ac.il>

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


/*! \file CGAL/Sqrt_extension.h
    \brief Defines class CGAL::Sqrt_extension.

    Provides the Number type Sqrt_extension that extends a given number type
    \c NT by a square root. One can add several square roots recursively.

*/

#ifndef CGAL_SQRT_EXTENSION_TYPE_H
#define CGAL_SQRT_EXTENSION_TYPE_H

#include <CGAL/disable_warnings.h>

#include <CGAL/number_type_basic.h>
#include <boost/operators.hpp>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Sqrt_extension_fwd.h>
#include <boost/optional.hpp>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/NT_converter.h>

#define CGAL_int(T)    typename First_if_different<int,    T>::Type

namespace CGAL {


template<class NT, class ROOT, class ACDE_TAG, class FP_TAG >
void output_maple(std::ostream&, const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>&);

/*! \ingroup CGAL_Sqrt_extension
\brief represents an extension of a number type by one square root.

 An instance of this class
represents  an extension of the type NT by a square root of the
type ROOT. In case NT and ROOT do not coincide,
NT must be constructible from ROOT.  The number type NT
must be at least a model of the IntegralDomainWithoutDiv concept.

An Sqrt_extension is a model of RealComparable if NT is RealComparable.\n
The <B>algebraic type</B> of CGAL::Sqrt_extension depends on the algebraic type
of NT:
- IntegralDomainWithoutDiv -> IntegralDomainWithoutDiv
- IntegralDomain           -> IntegralDomain
- UFDomain                 -> IntegralDomain
- EuclideanRing            -> IntegralDomain
- Field                    -> Field
- FieldWithSqrt            -> Field (not recommended)


Note that NT and ROOT can themselves be an instance of
CGAL::Sqrt_extension, yielding a nested extension.\n
Note that the extension of an UFDomain or EuclideanRing is just an
IntegralDomain, since the extension in general destroys the unique
factorization property.
*/

namespace internal {
template <class Tag>
class Interval_optional_caching;

//empty class
template <>
class Interval_optional_caching< ::CGAL::Tag_false >
{
protected:
  void invalidate_interval() {}
  bool is_cached() const {return false;}
  std::pair<double,double> cached_value() const {return std::pair<double,double>();}
  void update_cached_value(const std::pair<double,double>&) const {}
};


template <>
class Interval_optional_caching< ::CGAL::Tag_true >
{
protected:
  typedef boost::optional< std::pair<double,double> > Cached_interval;
  mutable Cached_interval interval_;
  void invalidate_interval() {interval_=Cached_interval();}
  bool is_cached() const {return (interval_?true:false);}
  std::pair<double,double> cached_value() const {return *interval_;}
  void update_cached_value(const std::pair<double,double>& i) const {interval_=i;}
};

} //namespace internal

template <class NT_,class ROOT_, class ACDE_TAG_, class FP_TAG>
class Sqrt_extension : public internal::Interval_optional_caching<FP_TAG>,
    boost::ordered_field_operators1< Sqrt_extension<NT_,ROOT_,ACDE_TAG_,FP_TAG>,
    boost::ordered_field_operators2< Sqrt_extension<NT_,ROOT_,ACDE_TAG_,FP_TAG>, NT_ ,
    boost::ordered_field_operators2< Sqrt_extension<NT_,ROOT_,ACDE_TAG_,FP_TAG>, CGAL_int(NT_)
> > >{
public:
    typedef NT_ NT;
    typedef ROOT_ ROOT;
    typedef ACDE_TAG_ ACDE_TAG;
    typedef Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> Self;
private:
    NT a0_;
    NT a1_;
    ROOT root_;
    bool is_extended_;

    typedef CGAL::Algebraic_structure_traits<NT> Algebraic_structure_traits_nt;
    typedef CGAL::Real_embeddable_traits<NT> Real_embeddable_traits_nt;

    typedef typename CGAL::Coercion_traits< ROOT, NT >::Cast Root_nt_cast;

public:
    //! Default constructor of \c Sqrt_extension
    Sqrt_extension()
        : a0_( NT(0)), a1_( NT(0)), root_(ROOT(0)), is_extended_(false) {}

    Sqrt_extension(CGAL_int(NT) i)
        : a0_(NT(i)), a1_(NT(0)), root_(ROOT(0)), is_extended_(false) {}

    Sqrt_extension(const NT& i)
        : a0_(NT(i)), a1_(NT(0)), root_(ROOT(0)), is_extended_(false) {}

    //non-documented: used for Make_sqrt
    Sqrt_extension(const ROOT& root,bool)
          :a0_(NT(0)), a1_(NT(1)), root_(root), is_extended_(true) {}

    /*!\brief Explicit constructor of Sqrt_extension, from any type NTX.
     * \pre NT must constructible from NTX */
    template <class NTX>
    explicit Sqrt_extension(const NTX& i)
        : a0_(NT(i)), a1_(NT(0)), root_(ROOT(0)), is_extended_(false) {}



    /*! \brief Expicite constructor of Sqrt_extension, from
     *  \c Sqrt_extension<NTX,ROOTX>.
     *  \pre NT must constructible from NTX
     *  \pre ROOT must constructible from ROOTX */
//     template <class NTX,class ROOTX>
//     explicit Sqrt_extension(const Sqrt_extension<NTX,ROOTX>& x)
//         : a0_(x.a0()),
//           a1_(x.a1()),
//           root_(x.root()),
//           is_extended_(x.is_extended()){}

    /*! \brief Constructor from some type NTX and ROOTX
     *  NT must be constructible from NTX and NTY\\
     *  ROOT must be construcible from ROOTX\\
     */
    template <class NTX,class NTY,class ROOTX>
    explicit Sqrt_extension(const NTX& a0, const NTY& a1, const ROOTX& root)
        : a0_(a0),
          a1_(a1),
          root_(root),
          is_extended_( true ){
      CGAL_precondition( ACDE_TAG_::value || !CGAL_NTS is_zero(root) ); //if ACDE_TAG_ is Tag_true, we allow root to be 0
    }

    /*! \brief Constructor of Sqrt_extension, from
     *  \c polynomial. The bool
     *  \pre NT must constructible from (NTX,NTX)
     *  \pre ROOT must constructible from (NTX,NTX)
     */
  template <class NTX>
  explicit Sqrt_extension(const NTX& a, const NTX& b, const NTX& c, const bool is_smaller,
    typename boost::enable_if< boost::mpl::and_<
      boost::is_same< typename Fraction_traits<NT>::Numerator_type,NTX >,
      boost::is_same< typename Fraction_traits<ROOT>::Numerator_type,NTX >
    > >::type* = 0  )
  {
    typename Fraction_traits<NT>::Compose compose_nt;
    typename Fraction_traits<ROOT>::Compose compose_root;
    if ( a != 0 ) {
      a0_ = compose_nt(-b,2*a);
      root_ = CGAL_NTS square(a0_) - compose_root(c,a);
      if(CGAL_NTS is_zero(root_)) {
        is_extended_ = false;
      } else {
        a1_ = (is_smaller ? -1 : 1);
        is_extended_ = true;
      }
    }
    else {
      CGAL_assertion( b != 0 );
      is_extended_ = false;
      a0_ = compose_nt(-c,b);
      a1_ = 0;
      root_ = 0;
    }
  }


  Self conjugate() const
  {
    if(!is_extended_) return *this;
    return Self(a0_,-a1_,root_);
  }

    //! Access operator for a0_, \c const
    inline const NT& a0() const { return a0_; }
    inline const NT& alpha() const { return a0_; } //for backward compatibility
    //! Access operator for a0_
    NT&        a0()       { return a0_; }
    NT&        alpha()       { return a0_; } //for backward compatibility
    //! Access operator for a1_, \c const
    inline const NT& a1() const { return a1_; }
    inline const NT& beta() const { return a1_; }  //for backward compatibility
    //! Access operator for a1_
    NT&        a1()       { return a1_; }
    NT&        beta()       { return a1_; }  //for backward compatibility
    //! Access operator for root_, \c const
    inline const ROOT& root() const { return root_; }
    inline const ROOT& gamma() const { return root_; }  //for backward compatibility
    //! Access operator for is_extended_, \c const
    inline const bool& is_extended() const { return is_extended_; }
    inline bool is_rational() const {
      CGAL_precondition( (boost::is_same<NT,ROOT>::value) || !"NT and ROOT should be identical and rational");
      return !is_extended_;} //for backward compatibility

    //!check if the number is an extension (test the values) and update the internal flag
    inline bool check_if_is_extended() {
      if(CGAL_NTS is_zero(a1_) || CGAL_NTS is_zero(root_))
        is_extended_ = false;
      return is_extended_;
    }

    //! Access operator for root_
    //ROOT& root() { return root_; }

    // TODO: add to documentation
    static void check_roots(const Self& a, const Self& b ) {
      if(a.is_extended() && b.is_extended())
        CGAL_precondition (a.root() == b.root());
    }

    // output_maple function for EXACUS compatibility
    // TODO: remove if no longer needed
    inline void output_maple(std::ostream& os ) const {
        CGAL::output_maple( os, *this );
    }

    std::pair<double, double> compute_to_interval() const
    {
      if (! is_extended_)
          return CGAL_NTS to_interval(a0_);

      const CGAL::Interval_nt<false>&  a0_int = CGAL_NTS to_interval(a0_);
      const CGAL::Interval_nt<false>&  a1_int = CGAL_NTS to_interval(a1_);
      const CGAL::Interval_nt<false>&  root_int = CGAL_NTS to_interval(root_);

      CGAL::Interval_nt<false>::Protector p;
      const CGAL::Interval_nt<false>&  x_int =
          a0_int + (a1_int * CGAL::sqrt(root_int));

      std::pair<double, double> res(x_int.inf(), x_int.sup());
      this->update_cached_value(res);
      return res;
    }

    std::pair<double,double> to_interval() const
    {
      std::pair<double,double> ret=this->is_cached()?this->cached_value():compute_to_interval();
      return ret;
    }

    //! propagates the simplify command to the members of xx
    void simplify(){
      if(is_extended_){
        if(CGAL_NTS is_zero(a1_)){
          is_extended_ = false;
        }else{
          CGAL_NTS simplify(a1_);
          CGAL_NTS simplify(root_);
        }
      }
      CGAL_NTS simplify(a0_);
    }

    //! determines the sign of xx by repeated squaring (with no filtering).
    ::CGAL::Sign sign_() const
    {
        ::CGAL::Sign s0,s1;

        s0 = CGAL_NTS sign(a0_);
        s1 = CGAL_NTS sign(a1_);

        if (s0 == s1) return s0;
        if (s0 == CGAL::ZERO) return s1;
        if (s1 == CGAL::ZERO) return s0;

        // s0*s1=-1
        NT r = a1_*a1_*Root_nt_cast()(root_) - a0_*a0_;
        // if(r>0) return s1 else s0
        if (s1 == CGAL::POSITIVE)
            return CGAL_NTS sign(r);
        else
            return CGAL::opposite (CGAL_NTS sign(r));
    }

    //! determines the sign of xx by repeated squaring. (with filtering).
    ::CGAL::Sign sign() const {
        if (! is_extended_)
            return CGAL_NTS sign(a0());

        if (FP_TAG::value){
          const std::pair<double, double>&  x_in = this->to_interval();

          if (x_in.first > 0)
            return (CGAL::POSITIVE);
          else if (x_in.second < 0)
            return (CGAL::NEGATIVE);
        }

        return (this->sign_());
    }

    template < class BOOL_TAG >
    bool is_zero_(const BOOL_TAG&) const{
        // Is_real_comparable == CGAL::Tag_true
        return sign()==CGAL::ZERO;
    }
    bool is_zero_(const CGAL::Tag_false) const {
        // Is_real_comparable == CGAL::Tag_false
        // i.e. CGAL::Modular
        if(is_extended()){
            if(a0() == (NT)0 && a1()== (NT)0) {
                return true;
            }else{
                return (a0()*a0()-a1()*a1()*(NT)root() == NT(0));
            }
        }else{
            return a0() == NT(0);
        }
    }

    //! returns \a true if xx is zero
    bool is_zero() const {
//        typedef typename Algebraic_structure_traits_nt::Is_real_comparable Is_real_comparable;
      typedef typename Real_embeddable_traits_nt::Is_real_embeddable
                                                             Is_real_embeddable;
        return is_zero_(Is_real_embeddable());
    }

    //! returns the absolute value of xx
    Self abs() const {
        if (sign() == CGAL::NEGATIVE)
            return -*this;
        return *this;
    }

    Self& operator += (const Self& p) {
      this->invalidate_interval();
        if(is_extended_){
            if (p.is_extended_){
                #ifndef NDEBUG
                Self::check_roots(*this, p);
                #endif
                return *this = Self (a0_+p.a0_, a1_+p.a1_, root_);
            }else{
                return *this = Self (a0_+p.a0_, a1_, root_);
            }
        }else{
            if (p.is_extended_)
                return *this = Self (a0_+p.a0_, p.a1_, p.root_);
            else
                return *this = Self (a0_+p.a0_);
        }
    }

    Self& operator -= (const Self& p) {
      this->invalidate_interval();
        if(is_extended_){
            if (p.is_extended_){
                #ifndef NDEBUG
                Self::check_roots(*this, p);
                #endif
                return *this = Self (a0_-p.a0_, a1_-p.a1_, root_);
            }else{
                return *this = Self (a0_-p.a0_, a1_, root_);
            }
        }else{
            if (p.is_extended_)
                return *this = Self (a0_-p.a0_, -p.a1_, p.root_);
            else
                return *this = Self (a0_-p.a0_);
        }
    }

    Self& operator *= (const Self& p){
      this->invalidate_interval();
        if(is_extended_){
            if (p.is_extended_){
                #ifndef NDEBUG
                Self::check_roots(*this, p);
                #endif
                return *this = Self (
                        a0_*p.a0_+a1_*p.a1_*Root_nt_cast()(root_),
                        a0_*p.a1_+a1_*p.a0_,
                        root_);
            }else{
                return *this = Self (a0_*p.a0_, a1_*p.a0_, root_);
            }
        }else{
            if (p.is_extended_)
                return *this = Self (a0_*p.a0_, a0_*p.a1_, p.root_);
            else
                return *this = Self (a0_*p.a0_);
        }
    }
    Self& operator /= (const Self& p) {
      this->invalidate_interval();

        CGAL_assertion( ! CGAL_NTS is_zero(p) );
        typename CGAL::Algebraic_structure_traits<NT>::Integral_division idiv;

        if(p.is_extended_){
            NT denom = p.a0_*p.a0_ - p.a1_*p.a1_ * Root_nt_cast()(p.root_);
            if ( CGAL_NTS is_zero(denom) ) {
                // this is for the rare case in which root is a square
                // and the (pseudo) algebraic conjugate of p becomes zero
                *this = Self(idiv(a0_, NT(2)*p.a0_) + idiv( a1_, NT(2)*p.a1_));
            }else{
                *this *= Self(p.a0_,-p.a1_,p.root_);
                *this =  Self( idiv(a0_,denom), idiv(a1_,denom), root_);
            }
        }else{
            if(is_extended_){
                *this = Self(idiv(a0_,p.a0_), idiv(a1_,p.a0_), root_);
            }else{
                *this = Self(idiv(a0_,p.a0_));
            }
        }
        return *this;
    }

//------------------------------------------------------------------
// SPECIALIZE_MEMBERS NT

    Self& operator += (const NT& num) {
      this->invalidate_interval();
        a0() += NT(num);
        return *this;
    }
    Self& operator -= (const NT& num) {
      this->invalidate_interval();
        a0() -= NT(num);
        return *this;
    }
    Self& operator *= (const NT& num) {
      this->invalidate_interval();
        a0() *= NT(num);
        a1() *= NT(num);
        return *this;
    }
    Self& operator /= (const NT& num) {
      this->invalidate_interval();
        CGAL_assertion(! CGAL_NTS is_zero(num));
        a0() /= NT(num);
        a1() /= NT(num);
        return *this;
    }

   Self& operator += (CGAL_int(NT) num) {
      this->invalidate_interval();
        a0() += NT(num);
        return *this;
    }
    Self& operator -= (CGAL_int(NT) num) {
      this->invalidate_interval();
        a0() -= NT(num);
        return *this;
    }
    Self& operator *= (CGAL_int(NT) num) {
      this->invalidate_interval();
        a0() *= NT(num);
        a1() *= NT(num);
        return *this;
    }
    Self& operator /= (CGAL_int(NT) num) {
      this->invalidate_interval();
        CGAL_assertion(! CGAL_NTS is_zero(num));
        a0() /= NT(num);
        a1() /= NT(num);
        return *this;
    }


    CGAL::Comparison_result
        compare (const CGAL_int(NT)& num) const {
        return this->compare(NT(num));
    }

// compare with NT
CGAL::Comparison_result
compare (const NT& num) const {

    if (! is_extended_)
        return (CGAL::compare (a0_, num));

    if (FP_TAG::value){
      const std::pair<double, double>&  x_in = this->to_interval();
      const std::pair<double, double>&  y_in = CGAL::to_interval (num);

      if (x_in.second < y_in.first)
        return (SMALLER);
      else if (x_in.first > y_in.second)
        return (LARGER);
    }

    return ((Self(a0_ - num, a1_, root_)).sign_());
}

// compare of two values that may be in different extension
// However, the default is, that the the numbers are defined over the same
// extension.
CGAL::Comparison_result
  compare(const Self& y, bool in_same_extension = !ACDE_TAG::value ) const
{
    if (! is_extended_)
        return (CGAL::opposite (y.compare (a0_)));
    else if (! y.is_extended_)
        return (this->compare (y.a0_));

    if (in_same_extension)
        return ((*this) - y).sign();

    if (FP_TAG::value)
    {
      const std::pair<double, double>&  x_in = this->to_interval();
      const std::pair<double, double>&  y_in = y.to_interval();

      if (x_in.second < y_in.first)
        return (SMALLER);
      else if (x_in.first > y_in.second)
        return (LARGER);
    }

  // Perform the exact comparison:
  // Note that the comparison of (a1 + b1*sqrt(c1)) and (a2 + b2*sqrt(c2))
  // is equivalent to comparing (a1 - a2) and (b2*sqrt(c2) -  b1*sqrt(c1)).
  // We first determine the signs of these terms.
  const NT          diff_a0 = a0_ - y.a0_;
  const CGAL::Sign  sign_left = CGAL::sign (diff_a0);
  const NT          x_sqr = a1_*a1_ * Root_nt_cast()(root_);
  const NT          y_sqr = y.a1_*y.a1_ * Root_nt_cast()(y.root_);
  Comparison_result right_res = CGAL::compare (y_sqr, x_sqr);
  CGAL::Sign        sign_right = ZERO;

  if (right_res == LARGER)
  {
    // Take the sign of b2:
    sign_right = CGAL::sign (y.a1_);
  }
  else if (right_res == SMALLER)
  {
    // Take the opposite sign of b1:
    switch (CGAL::sign (a1_))
    {
    case POSITIVE :
      sign_right = NEGATIVE;
      break;
    case NEGATIVE:
      sign_right = POSITIVE;
      break;
    case ZERO:
      sign_right = ZERO;
      break;
    default:
      // We should never reach here.
      CGAL_error();
    }
  }
  else
  {
    // We take the sign of (b2*sqrt(c2) -  b1*sqrt(c1)), where both terms
    // have the same absolute value. The sign is equal to the sign of b2,
    // unless both terms have the same sign, so the whole expression is 0.
    sign_right = CGAL::sign (y.a1_);
    if (sign_right == CGAL::sign (a1_))
      sign_right = ZERO;
  }

  // Check whether on of the terms is zero. In this case, the comparsion
  // result is simpler:
  if (sign_left == ZERO)
  {
    if (sign_right == POSITIVE)
      return (SMALLER);
    else if (sign_right == NEGATIVE)
      return (LARGER);
    else
      return (EQUAL);
  }
  else if (sign_right == ZERO)
  {
    if (sign_left == POSITIVE)
      return (LARGER);
    else if (sign_left == NEGATIVE)
      return (SMALLER);
    else
      return (EQUAL);
  }

  // If the signs are not equal, we can determine the comparison result:
  if (sign_left != sign_right)
  {
    if (sign_left == POSITIVE)
      return (LARGER);
    else
      return (SMALLER);
  }

  // We now square both terms and look at the sign of the one-root number:
  //   ((a1 - a2)^2 - (b12*c1 + b22*c2)) + 2*b1*b2*sqrt(c1*c2)
  //
  // If both signs are negative, we should swap the comparsion result
  // we eventually compute.
  const NT          A = diff_a0*diff_a0 - (x_sqr + y_sqr);
  const NT          B = 2 * a1_ * y.a1_;
  const ROOT        C = root_ * y.root_;
  const CGAL::Sign  sgn = (Self(A, B, C)).sign_();
  const bool        swap_res = (sign_left == NEGATIVE);

  if (sgn == POSITIVE)
    return (swap_res ? SMALLER : LARGER);
  else if (sgn == NEGATIVE)
    return (swap_res ? LARGER : SMALLER);

  return (EQUAL);
}

  private:
  template <class A,class B,class C> static bool
  is_equal (const Sqrt_extension<A,B,Tag_false,C>& p1, const Sqrt_extension<A,B,Tag_false,C>& p2)
  { return (p1-p2).is_zero(); }
  template <class A,class B,class C> static bool
  is_equal (const Sqrt_extension<A,B,Tag_true,C>& p1, const Sqrt_extension<A,B,Tag_true,C>& p2)
  { return ( p1.compare(p2) == CGAL::ZERO ); }

  public:
  // BINARY operators
  friend bool operator == (const Sqrt_extension& p1, const Sqrt_extension& p2)
    { return Sqrt_extension::is_equal(p1, p2); } // if constexpr (ACDE_TAG::value) ...
  friend bool operator <  (const Sqrt_extension& p1, const Sqrt_extension& p2)
    { return ( p1.compare(p2) == CGAL::SMALLER ); }

  // NT
  friend bool operator == (const Sqrt_extension& p, const NT& num)
    { return (p-num).is_zero();}
  friend bool operator <  (const Sqrt_extension& p, const NT& num)
    { return ( p.compare(num) == CGAL::SMALLER ); }
  friend bool operator >  (const Sqrt_extension& p, const NT& num)
    { return ( p.compare(num) == CGAL::LARGER ); }

  //CGAL_int(NT)
  friend bool operator == (const Sqrt_extension& p, CGAL_int(NT) num)
    { return (p-num).is_zero();}
  friend bool operator <  (const Sqrt_extension& p, CGAL_int(NT) num)
    { return ( p.compare(num) == CGAL::SMALLER ); }
  friend bool operator >  (const Sqrt_extension& p, CGAL_int(NT) num)
    { return ( p.compare(num) == CGAL::LARGER ); }
};

/*!
 * Compute the square of a one-root number.
 */
template <class NT,class ROOT,class ACDE_TAG,class FP_TAG>
Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> square (const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& x)
{
  if (!x.is_extended())
    return Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> (x.alpha() * x.alpha());

  // Use the equation:
  //
  //   (a + b*sqrt(c))^2 = (a^2 + b^2*c) + 2ab*sqrt(c)
  //
  return (Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> (x.alpha()*x.alpha() + x.beta()*x.beta() * NT(x.gamma()),
                                2 * x.alpha() * x.beta(),
                                x.gamma()));
}

//NT_converter specializations
template <class NT1,class ROOT1,class NT2,class ROOT2,class ACDE_TAG,class FP_TAG>
struct NT_converter < Sqrt_extension<NT1,ROOT1,ACDE_TAG,FP_TAG> , Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG> >
  : public CGAL::cpp98::unary_function< NT1, NT2 >
{
    Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG>
    operator()(const Sqrt_extension<NT1,ROOT1,ACDE_TAG,FP_TAG> &a) const
    {
      if(!a.is_extended()) {
        return Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG>(NT_converter<NT1,NT2>()(a.a0()));
      } else {
        return Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG>
            (NT_converter<NT1,NT2>()(a.a0()),
             NT_converter<NT1,NT2>()(a.a1()),
             NT_converter<ROOT1,ROOT2>()(a.root()));
      }
    }
};

template <class NT1,class NT2,class ROOT2,class ACDE_TAG,class FP_TAG>
struct NT_converter < NT1 , Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG> >
  : public CGAL::cpp98::unary_function< NT1, NT2 >
{
    Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG>
    operator()(const NT1 &a) const
    {
        return Sqrt_extension<NT2,ROOT2,ACDE_TAG,FP_TAG>(NT_converter<NT1,NT2>()(a));
    }
};

//needed because it's a better match than the specialization <NT1,NT1>
template <class NT1,class ROOT1,class ACDE_TAG,class FP_TAG>
struct NT_converter < Sqrt_extension<NT1,ROOT1,ACDE_TAG,FP_TAG>, Sqrt_extension<NT1,ROOT1,ACDE_TAG,FP_TAG> >
  : public CGAL::cpp98::unary_function< NT1, NT1 >
{
    const Sqrt_extension<NT1,ROOT1,ACDE_TAG,FP_TAG>&
    operator()(const Sqrt_extension<NT1,ROOT1,ACDE_TAG,FP_TAG> &a) const
    {
        return a;
    }
};

// UNARY
template <class NT,class ROOT,class ACDE_TAG,class FP_TAG> Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>
  operator + (const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& p) { return p; }

template <class NT,class ROOT,class ACDE_TAG,class FP_TAG> Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>
operator - (const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& p){
    if(p.is_extended())
      return  Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>(-p.a0(),-p.a1(),p.root());
    else
      return  Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>(-p.a0());
}


// BINARY

template<class NT, class ROOT,class ACDE_TAG,class FP_TAG> inline
Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> min BOOST_PREVENT_MACRO_SUBSTITUTION(
const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> & x,
const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> & y){
  return (std::min)(x,y);
}
template<class NT, class ROOT,class ACDE_TAG,class FP_TAG> inline
Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> max BOOST_PREVENT_MACRO_SUBSTITUTION(
const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> & x,
const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG> & y){
  return (std::max)(x,y);
}

template <class NT, class ROOT,class ACDE_TAG,class FP_TAG> inline
CGAL::Comparison_result compare (
    const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& x,
    const Sqrt_extension<NT,ROOT,ACDE_TAG,FP_TAG>& y,
    bool in_same_extension = !ACDE_TAG::value ){
  return x.compare(y,in_same_extension);
}


template <class NT_,class ROOT_, class ACDE_TAG_, class FP_TAG>
void
print(std::ostream &os, const Sqrt_extension<NT_,ROOT_,ACDE_TAG_,FP_TAG> &r)
{
  if(!r.is_extended()) {
    os << "(" << r.a0() << ")";
  } else {
    os << "(" << r.a0() << " + " << r.a1() <<
          "*sqrt(" << r.root() << ")"<< ")";
  }
}


} //namespace CGAL

#undef CGAL_int

#include <CGAL/enable_warnings.h>

#endif  // CGAL_SQRT_EXTENSION_TYPE_H

// EOF
