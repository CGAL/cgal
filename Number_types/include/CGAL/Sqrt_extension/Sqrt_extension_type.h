// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/Sqrt_extension.h $
// $Id: Sqrt_extension.h 38725 2007-05-15 14:39:12Z hemmer $
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

#include <CGAL/number_type_basic.h>
#include <boost/operators.hpp>

//#define SQRT_EXT_USE_FILTER 1
#ifdef SQRT_EXT_USE_FILTER
  #include <CGAL/Interval_arithmetic.h> 
#endif


#define CGAL_int(T)    typename First_if_different<int,    T>::Type

CGAL_BEGIN_NAMESPACE

template <class NT_,class ROOT_> class Sqrt_extension;
template<class NT, class ROOT>
void output_maple(std::ostream&, const Sqrt_extension<NT,ROOT>&);

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
template <class NT_,class ROOT_>
class Sqrt_extension : 
    boost::ordered_field_operators1< Sqrt_extension<NT_,ROOT_>,
    boost::ordered_field_operators2< Sqrt_extension<NT_,ROOT_>, NT_ ,
    boost::ordered_field_operators2< Sqrt_extension<NT_,ROOT_>, CGAL_int(NT_)
> > >{
public:
    typedef NT_ NT;
    typedef ROOT_ ROOT;
    typedef Sqrt_extension<NT,ROOT> Self;
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

    /*!\brief Explicit constructor of Sqrt_extension, from any type NTX.
     * \pre NT must constructible from NTX */
    template <class NTX>
    explicit Sqrt_extension(const NTX& i)
        : a0_(NT(i)), a1_(NT(0)), root_(ROOT(0)), is_extended_(false) {}

    /*! \brief copy constructor  */
    Sqrt_extension(const Sqrt_extension<NT,ROOT>& x)
        : a0_(x.a0()),
          a1_(x.a1()),
          root_(x.root()),
          is_extended_(x.is_extended()){}


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
     *  NT must be constructible from NTX\\
     *  ROOT must be construcible from ROOTX\\
     */
    template <class NTX,class ROOTX>
    explicit Sqrt_extension(const NTX& a0, const NTX& a1, const ROOTX& root)
        : a0_(a0),
          a1_(a1),
          root_(root),
          is_extended_(true){
        CGAL_precondition( !CGAL_NTS is_zero(root));
    }

    //! Access operator for a0_, \c const
    inline const NT& a0() const { return a0_; }
    //! Access operator for a0_
    NT&        a0()       { return a0_; }
    //! Access operator for a1_, \c const
    inline const NT& a1() const { return a1_; }
    //! Access operator for a1_
    NT&        a1()       { return a1_; }
    //! Access operator for root_, \c const
    inline const ROOT& root() const { return root_; }
    //! Access operator for is_extended_, \c const
    inline const bool& is_extended() const { return is_extended_; }
    //! Access operator for root_
    //ROOT& root() { return root_; }

    // TODO: add to documentation
    static void check_roots(const Self& a, const Self b ) {
        if(a.is_extended() && b.is_extended()){
            if (a.root() != b.root()){
                //CGAL_error_msg("Interoperation of incompatible algebraic extensions");
                throw std::range_error("Interoperation of incompatible algebraic extensions");
            }
        }
    }

    // output_maple function for EXACUS compatibility
    // TODO: remove if no longer needed
    inline void output_maple(std::ostream& os ) const {
        CGAL::output_maple( os, *this );
    }

    std::pair<double, double> to_interval() const{
   
    if (! is_extended_)
        return CGAL_NTS to_interval(a0_);
    
    const CGAL::Interval_nt<true>&  a0_int = CGAL_NTS to_interval(a0_);
    const CGAL::Interval_nt<true>&  a1_int = CGAL_NTS to_interval(a1_);
    const CGAL::Interval_nt<true>&  root_int = CGAL_NTS to_interval(root_);
    const CGAL::Interval_nt<true>&  x_int = 
        a0_int + (a1_int * CGAL::sqrt(root_int));

    return (std::make_pair (x_int.inf(), x_int.sup()));
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

#ifdef SQRT_EXT_USE_FILTER
        const std::pair<double, double>&  x_in = this->to_interval(); 

        if (x_in.first > 0)
          return (CGAL::POSITIVE);
        else if (x_in.second < 0)
          return (CGAL::NEGATIVE);
#endif

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
        if(is_extended_){
            if (p.is_extended_){
                Self::check_roots(*this, p);
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

        if(is_extended_){
            if (p.is_extended_){
                Self::check_roots(*this, p);
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
        if(is_extended_){
            if (p.is_extended_){
                Self::check_roots(*this, p);
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
        a0() += NT(num);
        return *this;
    }
    Self& operator -= (const NT& num) {
        a0() -= NT(num);
        return *this;
    }
    Self& operator *= (const NT& num) {
        a0() *= NT(num);
        a1() *= NT(num);
        return *this;
    }
    Self& operator /= (const NT& num) {
        CGAL_assertion(! CGAL_NTS is_zero(num));
        a0() /= NT(num);
        a1() /= NT(num);
        return *this;
    }

   Self& operator += (CGAL_int(NT) num) {
        a0() += NT(num);
        return *this;
    }
    Self& operator -= (CGAL_int(NT) num) {
        a0() -= NT(num);
        return *this;
    }
    Self& operator *= (CGAL_int(NT) num) {
        a0() *= NT(num);
        a1() *= NT(num);
        return *this;
    }
    Self& operator /= (CGAL_int(NT) num) {
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

#ifdef SQRT_EXT_USE_FILTER
    const std::pair<double, double>&  x_in = this->to_interval(); 
    const std::pair<double, double>&  y_in = CGAL::to_interval (num); 
    
    if (x_in.second < y_in.first)
      return (SMALLER);
    else if (x_in.first > y_in.second)
      return (LARGER);
#endif

    return ((Self(a0_ - num, a1_, root_)).sign_());
}

// compare of two values that may be in different extension
// However, the default is, that the the numbers are defined over the same
// extension. 
CGAL::Comparison_result
compare(const Self& y, bool in_same_extension = true ) const
{
    if (! is_extended_)
        return (CGAL::opposite (y.compare (a0_)));
    else if (! y.is_extended_)
        return (this->compare (y.a0_));

    if (in_same_extension)
        return ((*this) - y).sign();

#ifdef SQRT_EXT_USE_FILTER
    const std::pair<double, double>&  x_in = this->to_interval(); 
    const std::pair<double, double>&  y_in = y.to_interval(); 
    
    if (x_in.second < y_in.first)
      return (SMALLER);
    else if (x_in.first > y_in.second)
      return (LARGER);
#endif

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
};

// UNARY 
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator + (const Sqrt_extension<NT,ROOT>& p) { return p; }

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>& p){
    if(p.is_extended())
        return  Sqrt_extension<NT,ROOT>(-p.a0(),-p.a1(),p.root());
    else
        return  Sqrt_extension<NT,ROOT>(-p.a0());
}


// BINARY 
template <class NT,class ROOT> bool
operator == (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return (p1-p2).is_zero() ; }
template <class NT,class ROOT> bool
operator < (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( p1.compare(p2) == CGAL::SMALLER ); }


// NT
// righthand side
template <class NT,class ROOT>   bool operator ==
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return (p-num).is_zero();}
template <class NT,class ROOT>    bool operator <
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( p.compare(num) == CGAL::SMALLER ); }
template <class NT,class ROOT>    bool operator >
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( p.compare(num) == CGAL::LARGER ); }

// lefthand side
template <class NT,class ROOT> bool operator ==
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return  ( p == num );}
template <class NT,class ROOT>  bool operator <
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( p.compare(num) == CGAL::LARGER ); }
template <class NT,class ROOT>  bool operator >
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( p.compare(num) == CGAL::SMALLER ); }


//CGAL_int(NT)
template <class NT,class ROOT>    bool operator ==
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return (p-num).is_zero();}
template <class NT,class ROOT>    bool operator <
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( p.compare(num) == CGAL::SMALLER ); }
template <class NT,class ROOT>    bool operator >
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( p.compare(num) == CGAL::LARGER ); }

template <class NT,class ROOT> bool operator ==
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return  ( p == num );}
template <class NT,class ROOT>  bool operator <
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( p.compare(num) == CGAL::LARGER ); }
template <class NT,class ROOT>  bool operator >
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( p.compare(num) == CGAL::SMALLER ); }

template<typename NT, typename ROOT> inline 
Sqrt_extension<NT,ROOT> min BOOST_PREVENT_MACRO_SUBSTITUTION(
const Sqrt_extension<NT,ROOT> & x,
const Sqrt_extension<NT,ROOT> & y){
  return (std::min)(x,y);
}
template<typename NT, typename ROOT> inline 
Sqrt_extension<NT,ROOT> max BOOST_PREVENT_MACRO_SUBSTITUTION(
const Sqrt_extension<NT,ROOT> & x,
const Sqrt_extension<NT,ROOT> & y){
  return (std::max)(x,y);
}

template <class NT, class ROOT> inline 
CGAL::Comparison_result compare (
    const Sqrt_extension<NT,ROOT>& x, 
    const Sqrt_extension<NT,ROOT>& y,
    bool in_same_extension = true ){
  return x.compare(y,in_same_extension);
}

CGAL_END_NAMESPACE

#undef CGAL_int

#endif  // CGAL_SQRT_EXTENSION_TYPE_H

// EOF
