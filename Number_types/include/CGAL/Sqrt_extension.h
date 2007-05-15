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
//                 Ron Wein         <wein@post.tau.ac.il>

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


/*! \file CGAL/Sqrt_extension.h
    \brief Defines class CGAL::Sqrt_extension.

    Provides the Number type Sqrt_extension that extends a given number type
    \c NT by a square root. One can add several square roots recursively.

*/

#ifndef CGAL_SQRT_EXTENSION_H
#define CGAL_SQRT_EXTENSION_H

#include <CGAL/number_type_basic.h>


// #define SQRT_EXT_USE_FILTER 1
#ifdef SQRT_EXT_USE_FILTER
  #include <CGAL/Interval_arithmetic.h> 
#endif

#define CGAL_int(T)    typename First_if_different<int,    T>::Type

CGAL_BEGIN_NAMESPACE

template <class NT,class ROOT> class Sqrt_extension;

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator + (const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator + (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator * (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> inline Sqrt_extension<NT,ROOT>
operator / (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template<class NT,class ROOT>
std::ostream& operator << (std::ostream& os, const Sqrt_extension<NT,ROOT>& p);
template <class NT,class ROOT>
std::istream& operator >> (std::istream& is, Sqrt_extension<NT,ROOT>& p);

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
class Sqrt_extension {
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
        NT r = a1_*a1_*NT(root_) - a0_*a0_;
        // if(r>0) return s1 else s0
        if (s1 == CGAL::POSITIVE)
            return CGAL_NTS sign(r);
        else
            return CGAL::opposite (CGAL_NTS sign(r));
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

public:

    //! propagates the simplify command to the members of xx
    void simplify(){
        CGAL_NTS simplify(a0_);
        CGAL_NTS simplify(a1_);
        CGAL_NTS simplify(root_);
    }

    //! determines the sign of xx by repeated squaring.
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

    friend Self operator + <>(const Self&);
    friend Self operator - <>(const Self&);
    friend Self operator + <>(const Self&, const Self&);
    friend Self operator - <>(const Self&, const Self&);
    friend Self operator * <>(const Self&, const Self&);
    friend Self operator / <>(const Self& p1, const Self& p2);

    Self& operator += (const Self& p) {
        (*this)=(*this)+p;
        return (*this);
    }
    Self& operator -= (const Self& p){
        (*this)=(*this)-p;
        return (*this);
    }
    Self& operator *= (const Self& p){
        (*this)=(*this)*p;
        return (*this);
    }
    Self& operator /= (const Self& p) {
        (*this)=(*this)/p;
        return (*this);
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
        typename Real_embeddable_traits_nt::Sign sign;
        CGAL_assert(sign(num) != 0);
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
        typename Real_embeddable_traits_nt::Sign sign;
        CGAL_assert(sign(num) != 0);
        a0() /= NT(num);
        a1() /= NT(num);
        return *this;
    }

public:
    void output_maple(std::ostream& os) const;
public:
    /*! \brief write Sqrt_extension to \c os in a format readable
         by \c input_ascii()

         The output format is:
         <TT><B>EXT[</B></TT><I>a0</I><TT><B>,</B></TT>
         <I>a1</I><TT><B>,</B></TT><I>root</I><TT><B>]</B></TT>
     */
   void output_ascii(std::ostream& os) const{
       os<<"EXT["<<a0()<<","<<a1()<<","<<root()<<"]";
   }
    /*! \relates CGAL::Sqrt_extension
     *  \brief read a Sqrt_extension from \c is
     *
     * The expected format is:
     * <TT><B>EXT[</B></TT><I>a0</I><TT><B>,</B></TT>
     * <I>a1</I><TT><B>,</B></TT><I>root</I><TT><B>]</B></TT>
     *
     * The format of the coefficients must be understandable for
     * <TT> is >> iformat(ai) </TT>.
     *
     * Example: A \c CGAL::Sqrt_extension<int,Root<1,int> > with a value of
     * 4-2*sqrt(5) has to be written as
     * \c <TT>Self[4,-2,5]</TT> for this function.
     */
    static Sqrt_extension<NT,ROOT> input_ascii(std::istream& is);

public:

// compare with NT
CGAL::Comparison_result
compare (const NT& num) const
{
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

// compare of two values with different extension
CGAL::Comparison_result
compare(const Self& y, bool in_same_extension = false ) const
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
  const NT          x_sqr = a1_*a1_ * NT(root_);
  const NT          y_sqr = y.a1_*y.a1_ * NT(y.root_);
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
      CGAL_assertion (false);
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

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator + (const Sqrt_extension<NT,ROOT>& p) { return p; }

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>& p){
    if(p.is_extended())
        return  Sqrt_extension<NT,ROOT>(-p.a0(),-p.a1(),p.root());
    else
        return  Sqrt_extension<NT,ROOT>(-p.a0());
}

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator + (const Sqrt_extension<NT,ROOT>& p1,
            const Sqrt_extension<NT,ROOT>& p2){
    typedef Sqrt_extension<NT,ROOT> EXT;
    if(p1.is_extended()){
        if (p2.is_extended()){CGAL_precondition(p1.root()==p2.root());}
        return EXT (p1.a0()+p2.a0(), p1.a1()+p2.a1(), p1.root());
    }else{
        if (p2.is_extended())
            return EXT (p1.a0()+p2.a0(), p2.a1(), p2.root());
        else
            return EXT (p1.a0()+p2.a0());
    }
}

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2) {
    typedef Sqrt_extension<NT,ROOT> EXT;
    if(p1.is_extended()){
        if (p2.is_extended()){CGAL_precondition(p1.root()==p2.root());}
        return EXT (p1.a0()-p2.a0(), p1.a1()-p2.a1(), p1.root());
    }else{
        if (p2.is_extended())
            return EXT (p1.a0()-p2.a0(), -p2.a1(), p2.root());
        else
            return EXT (p1.a0()-p2.a0());
    }
}

template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator * (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2) {
    typedef Sqrt_extension<NT,ROOT> EXT;

    if(p1.is_extended()){
        if (p2.is_extended()){
            CGAL_precondition(p1.root()==p2.root());
            return EXT (p1.a0()*p2.a0()+p1.a1()*p2.a1()*NT(p1.root()),
                        p1.a0()*p2.a1()+p1.a1()*p2.a0(),
                        p1.root());
        }else{
            return EXT (p1.a0()*p2.a0(), p1.a1()*p2.a0(), p1.root());
        }
    }else{
        if (p2.is_extended())
            return EXT (p1.a0()*p2.a0(), p1.a0()*p2.a1(), p2.root());
        else
            return EXT (p1.a0()*p2.a0());
    }
}

template <class NT,class ROOT> inline Sqrt_extension<NT,ROOT>
operator / (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2) {

    typedef Sqrt_extension<NT,ROOT> EXT;
    CGAL_assertion(! p2.is_zero());
    typename CGAL::Algebraic_structure_traits<NT>::Integral_division Idiv;

    if(p1.is_extended()){
        if (p2.is_extended()){
            CGAL_precondition(p1.root()==p2.root());
            NT c = p2.a0()*p2.a0()-p2.a1()*p2.a1()*NT(p2.root());
            if (c == NT(0)) {                          //TR
                NT a0 = Idiv(p1.a0(), NT(2)*p2.a0())
                      + Idiv(p1.a1(), NT(2)*p2.a1());
                return EXT(a0);
            }
            NT a0 = Idiv(p1.a0()*p2.a0()-p1.a1()*p2.a1()*NT(p2.root()),c);
            NT a1 = Idiv(p1.a1()*p2.a0()-p1.a0()*p2.a1(),c);
            return EXT(a0,a1,p1.root());
        }else{
            NT a0 = Idiv(p1.a0(),p2.a0());
            NT a1 = Idiv(p1.a1(),p2.a0());
            return EXT(a0,a1,p1.root());
        }
    }else{
        if (p2.is_extended()){
            NT c = p2.a0()*p2.a0()-p2.a1()*p2.a1()*NT(p2.root());
            if(c == NT(0)){
                NT a0 = Idiv(p1.a0(), NT(2)*p2.a0());
                return EXT(a0);
            }
            NT a0 = Idiv(p1.a0()*p2.a0(),c);
            NT a1 = Idiv(-p1.a0()*p2.a1(),c);
            return EXT(a0,a1,p2.root());
        }else{
            return EXT (Idiv(p1.a0(),p2.a0()));
        }
    }
}

template <class NT,class ROOT> bool
operator == (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return (p1-p2).is_zero() ; }

template <class NT,class ROOT> bool
operator != (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return !(p1-p2).is_zero() ; }

template <class NT,class ROOT> bool
operator <  (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( (p1-p2).sign() < 0 ); }

template <class NT,class ROOT> bool
operator <= (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( (p1-p2).sign() <= 0 ); }

template <class NT,class ROOT> bool
operator >  (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( (p1-p2).sign() > 0 ); }

template <class NT,class ROOT> bool
operator >= (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
  { return ( (p1-p2).sign() >= 0 ); }


// lefthand side
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator +
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p2.is_extended())
        return EXT (num + p2.a0(), p2.a1(), p2.root());
    else
        return EXT (num + p2.a0());
}

template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator -
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p2.is_extended())
        return EXT (num - p2.a0(), -p2.a1(), p2.root());
    else
        return EXT (num - p2.a0());
}

template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator *
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p2.is_extended())
        return EXT (num * p2.a0(), num * p2.a1(), p2.root());
    else
        return EXT (num * p2.a0());
}

template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator /
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num)/p2); }

// righthand side
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator +
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p1.is_extended())
        return EXT (p1.a0() + num, p1.a1(), p1.root());
    else
        return EXT (p1.a0() + num);
}

template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator -
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p1.is_extended())
        return EXT (p1.a0() - num, p1.a1(), p1.root());
    else
        return EXT (p1.a0() - num);
}

template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator *
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p1.is_extended())
        return EXT (p1.a0() * num, p1.a1() * num, p1.root());
    else
        return EXT (p1.a0() * num);
}

template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator /
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{
    typedef Sqrt_extension<NT,ROOT> EXT;

    if (p1.is_extended())
        return EXT (p1.a0() / num, p1.a1() / num, p1.root());
    else
        return EXT (p1.a0() / num);
}

// lefthand side
template <class NT,class ROOT>    bool operator ==
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() == 0 );}
template <class NT,class ROOT>    bool operator !=
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() != 0 );}
template <class NT,class ROOT>    bool operator <
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() < 0 );}
template <class NT,class ROOT>    bool operator <=
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() <= 0 );}
template <class NT,class ROOT>    bool operator >
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() > 0 );}
template <class NT,class ROOT>    bool operator >=
(const NT& num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() >= 0 );}

// righthand side
template <class NT,class ROOT>    bool operator ==
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() == 0 );}
template <class NT,class ROOT>    bool operator !=
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() != 0 );}
template <class NT,class ROOT>    bool operator <
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() < 0 );}
template <class NT,class ROOT>    bool operator <=
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() <= 0 );}
template <class NT,class ROOT>    bool operator >
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() > 0 );}
template <class NT,class ROOT>    bool operator >=
(const Sqrt_extension<NT,ROOT>& p, const NT& num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() >= 0 );}

// Algebraic structure traits
template< class Type, class Algebraic_type >
class Sqrt_extension_algebraic_structure_traits_base;

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                        CGAL::Integral_domain_without_division_tag >
  : public Algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_without_division_tag Algebraic_category;

    class Simplify
      : public Unary_function< Type&, void > {
      public:
        typedef void result_type;
        typedef Type& argument_type;

        void operator()( Type& x ) const {
          x.simplify();
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Integral_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_tag Algebraic_category;

    class Integral_division
      : public Binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y ) const {
          return x/y;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Unique_factorization_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Euclidean_ring_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Field_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  public:
    typedef Field_tag Algebraic_category;

    class Unit_part
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return( x == Type(0) ? Type(1) : x );
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                    CGAL::Field_with_sqrt_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Field_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                CGAL::Field_with_kth_root_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type,
                                                CGAL::Field_with_root_of_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

template< class COEFF, class ROOT>
class Algebraic_structure_traits< Sqrt_extension< COEFF, ROOT > >
  : public Sqrt_extension_algebraic_structure_traits_base<
      Sqrt_extension< COEFF, ROOT >,
      typename Algebraic_structure_traits< COEFF >::Algebraic_category > {
  public:
    typedef Sqrt_extension< COEFF, ROOT > Type;

    // Tag_true if COEFF and ROOT are exact
    typedef typename ::boost::mpl::if_c<
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<ROOT >::Is_exact,::CGAL::Tag_true>::value )&&
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<COEFF>::Is_exact,::CGAL::Tag_true>::value )
           ,::CGAL::Tag_true,::CGAL::Tag_false>::type Is_exact;

    typedef typename Algebraic_structure_traits<COEFF>::Is_numerical_sensitive
    Is_numerical_sensitive;
};

//
// Real embeddable traits
//

template< class COEFF, class ROOT >
class Real_embeddable_traits< Sqrt_extension<COEFF, ROOT> >
  : public INTERN_RET::Real_embeddable_traits_base_selector<
                  Sqrt_extension<COEFF, ROOT>,
                  typename Real_embeddable_traits<COEFF>::Is_real_embeddable > {
  public:
    typedef Sqrt_extension<COEFF, ROOT> Type;

    class Sign
        : public Unary_function< Type, ::CGAL::Sign >{
    public:
        ::CGAL::Sign operator()( const Type& x ) const {
            return x.sign();
        }
    };
    
    class Compare
        : public Binary_function< Type, Type, Comparison_result > {
    public:
        Comparison_result operator()( const Type& x, const Type& y) const {
            // must be from the same extension 
            return x.compare(y);
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type, 
                Comparison_result )
    };

    class To_interval
      : public Unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double,double> operator()(const Type& x) const {
            return x.to_interval();
        }
    };

    class To_double
      : public Unary_function< Type, double > {
      public:
        // The main problem here is, that even tough the total
        // expression fits into double, one of the coefficients
        // or the root may not. ?? !
        double operator()(const Type& x) const {
            if(x.is_extended()){
                return CGAL_NTS to_double(x.a0())
                    +  int(CGAL_NTS sign(x.a1()))
                    * CGAL_NTS sqrt(CGAL_NTS to_double(x.a1()*x.a1() *
                                                    Type(x.root())));
            }else{
                return CGAL_NTS to_double(x.a0());
            }
        }
    };
};

// #################### IO BEGIN ############################################
// #################### INPUT


template<class NT, class ROOT>
Sqrt_extension<NT,ROOT>
Sqrt_extension<NT,ROOT>::input_ascii(std::istream& is){
    // expected input format: EXT[ext.a0(),ext.a1(),root()]
    typedef Sqrt_extension<NT,ROOT> EXT;
    char c;
    NT a0;
    NT a1;
    ROOT root;

    swallow(is, 'E');
    swallow(is, 'X');
    swallow(is, 'T');
    swallow(is, '[');
    is >> iformat(a0);
    do c = is.get(); while (isspace(c));
    // TODO: Replace CGAL_assertion_msg( false, ... ) with CGAL_error
    if (c != ',') CGAL_error( "input error: , expected" );

    is >> iformat(a1);
    do c = is.get(); while (isspace(c));
    if (c != ',') CGAL_error( "input error: , expected" );

    is >> iformat(root);
    do c = is.get(); while (isspace(c));
    if (c != ']') CGAL_error( "input error: ] expected" );

    if ( root  < ROOT(0)) CGAL_error("input error: non-negative root expected");
    if ( root == ROOT(0)) return EXT(a0);
    return EXT(a0,a1,root);
}

// ##################### OUTPUT
//! write Sqrt_extension to \c os in \c LiS::IO::PRETTY format
/*! The output is intended to be Maple-readable; see module
 *  \link CGAL_io CGAL I/O Support \endlink.
 */
template<class NT, class ROOT>
void
Sqrt_extension<NT,ROOT>::output_maple(std::ostream& os) const{
    CGAL::IO::Mode o_mode=::CGAL::get_mode(os);
    ::CGAL::set_mode(os,CGAL::IO::PRETTY);

    if ( a0() != NT(0)){
        if ( a1() != NT(0)){
            os << a0()
               << "+" << CGAL::oformat(a1(),CGAL::Parens_as_product_tag())
               << "*sqrt(" << root() << ")";
        }else{
            os << a0();
        }
    }
    else{
        if (a1() != NT(0)){
            os << CGAL::oformat(a1(),CGAL::Parens_as_product_tag())
               << "*sqrt(" << root() << ")";
        }else{
            os << 0;
        }
    }
    ::CGAL::set_mode(os,o_mode);
    return;
}

/*! \relates CGAL::Sqrt_extension
 *  \brief output \c ext to \c os
 *
 *  Output \c ext in a format as specified by
 *  \c LiS::get_mode(os), see \link LiS_io LiS I/O Support \endlink.
 *  Currently, the output for \c LiS::IO::BINARY happens to be
 *  identical to \c LiS::IO::ASCII.
 */
template <class NT,class ROOT>
std::ostream& operator << (std::ostream& os,
                           const Sqrt_extension<NT,ROOT>& ext){
    switch(CGAL::get_mode(os)) {
    case CGAL::IO::PRETTY:
        ext.output_maple(os); break;
    default:
        ext.output_ascii(os); break;
    }
    return os;
}

/*! \relates CGAL::Sqrt_extension
 *  \brief try to read a CGAL::Sqrt_extension from \c is into \c ext
 *
 *  \c is must be in a mode that supports input of CGAL::Sqrt_extension
 *  (\c LiS::IO::ASCII or \c LiS::IO::BINARY) and the input from
 *  \c is must have the format of output to a stream of the same mode.
 */
template <class NT,class ROOT>
std::istream& operator >> (std::istream& is, Sqrt_extension<NT,ROOT>& ext) {
    CGAL_precondition(!CGAL::is_pretty(is));
    ext = Sqrt_extension<NT,ROOT>::input_ascii(is);
    return is;
}
// ################################ IO END ###############################

//################################# CGAL::Fraction_traits ##################
// Select the right alternative as Fraction_traits
// The actual Type traits is Intern::Sqrt_ext_Ftr_base_2
// The selction is done in two steps:
// 1. Inter::Sqrt_ext_Ftr_base_1 selects by the BOOL_TAG whether the COEFF type
//    Is_fraction
// 2. Intern::Sqrt_ext_Ftr_base_2 checks whether the internal type of the ROOT
//    is still implicite convertible to the new COEFF type.
//    since the ROOT type it self can not be converted.
namespace Intern{
    template <class EXT, bool> class Sqrt_ext_Ftr_base_2;
    template <class EXT, class BOOL_TAG> class Sqrt_ext_Ftr_base_1;
}

/*! \ingroup CGAL_Sqrt_extension
    \ingroup CGAL_Fraction_traits_spec
    \brief Specialisation of CGAL::Fraction_traits for CGAL::Sqrt_extension.
 *
 *  Extensions provide suitable specializations of \c CGAL::Fraction_traits.
 *  They are decomposable iff their coefficient type is.
 *  The denominator \e d of a Extension \e ext is a low common multiple
 *  (see \c CGAL::Fraction_traits::Common_factor for details) of the
 *  denominators of its coefficients.  The numerator is the Extenion
 *  \e d*ext with a fraction-free coefficient type.
 *
 *  This works for nested Sqrt_extensions, too.
 */

template <class COEFF, class ROOT_NT >
class Fraction_traits< Sqrt_extension<COEFF,ROOT_NT > >
    : public Intern::Sqrt_ext_Ftr_base_1<
    Sqrt_extension<COEFF,ROOT_NT >,
    typename CGAL::Fraction_traits<COEFF>::Is_fraction >
{
    // nothing new
};

namespace Intern {

// Use this if the coefficients cannot be decomposed
// into numerator and denominator
template <class NT_ >
class Sqrt_ext_Ftr_base_2< NT_, false > {
public:
    typedef NT_ NT;
    typedef ::CGAL::Tag_false Is_fraction;
    typedef ::CGAL::Null_tag Numerator_type;
    typedef ::CGAL::Null_tag Denominator_type;
    typedef ::CGAL::Null_tag Common_factor;
    typedef ::CGAL::Null_tag Decompose;
    typedef ::CGAL::Null_tag Compose;
};

template <class COEFF, class ROOT_NT>
class Sqrt_ext_Ftr_base_2< Sqrt_extension<COEFF,ROOT_NT>, true > {
private:
    typedef Fraction_traits<COEFF> CFT;
public:
    typedef Sqrt_extension<COEFF,ROOT_NT> NT;
    typedef CGAL::Tag_true Is_fraction;
    typedef Sqrt_extension<typename CFT::Numerator_type,ROOT_NT> Numerator_type;
    typedef typename CFT::Denominator_type Denominator_type;
    typedef typename Algebraic_structure_traits<Denominator_type>::Gcd Common_factor;

    class Decompose {
    public:
        typedef NT first_argument_type;
        typedef Numerator_type second_argument_type;
        typedef Denominator_type& third_argument_type;
        void operator () (const NT& ext,
                          Numerator_type&   num,
                          Denominator_type& den){
            typename CFT::Decompose decompose;
            typename CFT::Common_factor common_factor;
            typedef typename CFT::Numerator_type NUM;
            typedef typename CFT::Denominator_type DEN;

            if(ext.is_extended()){
                NUM a0_num, a1_num;
                DEN a0_den, a1_den;
                DEN common_den;
                decompose(ext.a0(),a0_num,a0_den);
                decompose(ext.a1(),a1_num,a1_den);
                common_den=common_factor(a0_den,a1_den);

                a0_num = a0_num*CGAL::integral_division(a1_den,common_den);
                a1_num = a1_num*CGAL::integral_division(a0_den,common_den);
                den = CGAL::integral_division(a0_den,common_den)*a1_den;
                num = Numerator_type(a0_num,a1_num,ext.root());
            }else{
                NUM a0_num;
                decompose(ext.a0(),a0_num,den);
                num = Numerator_type(a0_num);
            }
        }
    };
    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef NT result_type;
        NT operator () (const Numerator_type&   num,
                        const Denominator_type& den){
            if(num.is_extended()){
                typename CFT::Compose compose;
                COEFF a0=compose(num.a0(),den);
                COEFF a1=compose(num.a1(),den);
                return NT(a0,a1,num.root());
            }else{
                typename CFT::Compose compose;
                COEFF a0=compose(num.a0(),den);
                return NT(a0);
            }
        }
    };
};

template <class EXT, class BOOL_TAG>
class Sqrt_ext_Ftr_base_1;

template <class COEFF, class ROOT_NT>
class Sqrt_ext_Ftr_base_1< Sqrt_extension<COEFF,ROOT_NT >, CGAL::Tag_true >
    : public Sqrt_ext_Ftr_base_2<
    Sqrt_extension<COEFF,ROOT_NT >,
    ::boost::is_same< typename CGAL::Coercion_traits<ROOT_NT,typename CGAL::Fraction_traits<COEFF>::Numerator_type>::Type,
                        typename CGAL::Fraction_traits<COEFF>::Numerator_type>::value >
{
    //nothing new
};

    template <class COEFF, class ROOT_NT>
    class Sqrt_ext_Ftr_base_1< Sqrt_extension<COEFF,ROOT_NT>, CGAL::Tag_false >
        : public Sqrt_ext_Ftr_base_2< Sqrt_extension<COEFF,ROOT_NT >, false>
    {
        //nothing new
    };
} // namespace Intern


/*
namespace Intern{
    template <class SqrtExt,class BoolTag> class Sqrt_ext_Coftr_base_1;
    template <class SqrtExt>
    class Sqrt_ext_Coftr_base_1< SqrtExt, CGAL::Tag_false >{
    public:
        typedef SqrtExt          Numerator_type;
        typedef ::CGAL::Tag_false Is_composable;
        typedef ::CGAL::Null_tag Denominator_type;
        typedef ::CGAL::Null_tag Type;
        typedef ::CGAL::Null_tag Compose;
    };
    template <class SqrtExt>
    class Sqrt_ext_Coftr_base_1< SqrtExt, CGAL::Tag_true >{
        typedef typename SqrtExt::NT Coeff;
        typedef typename SqrtExt::ROOT Root;
        typedef typename CGAL::Cofraction_traits<Coeff> CFT;
        typedef typename CFT::Type Type_coeff;

    public:
        typedef SqrtExt                                       Numerator_type;
        typedef ::CGAL::Tag_true                               Is_composable;
        typedef typename CFT::Denominator_type                Denominator;
        typedef CGAL::Sqrt_extension<Type_coeff,Root> Type;

        class Compose {
    public:
            //! first argument type
            typedef Numerator_type   first_argument_type;
            //! second argument type
            typedef Denominator_type second_argument_type;
            //! result type
            typedef Type    result_type;
            //! Compose fraction
            Type operator() (Numerator_type   num,
                                      Denominator_type den){
                if(num.is_extended()){
                    typename CFT::Compose compose_coeff;
                    Type_coeff a0_new(compose_coeff(num.a0(),den));
                    Type_coeff a1_new(compose_coeff(num.a1(),den));
                    return result_type(a0_new, a1_new, num.root());
                }else{
                    typename CFT::Compose compose_coeff;
                    return result_type(compose_coeff(num.a0(),den));
                }
            };
        };
    };
}

template <class Coeff, class Root>
class Cofraction_traits<Sqrt_extension<Coeff,Root> >
    :public Intern::Sqrt_ext_Coftr_base_1<
    Sqrt_extension<Coeff,Root>,
    typename CGAL::Cofraction_traits<Coeff>::Is_composable>{
    //nothing new;
};
*/


template <class COEFF, class ROOT>
class Needs_parens_as_product< Sqrt_extension<COEFF,ROOT> >{
public:
    typedef Sqrt_extension<COEFF,ROOT> NT;
    bool operator()(const NT& t){
        if( t.a0() != NT(0) && t.a1() != NT(0)){
            return true;
        }
        if( t.a1() == NT(0) ){
            Needs_parens_as_product<COEFF> npap;
            return npap(t.a0());
        }
        return false;
    }
};





/////////// COERCION_TRAITS BEGIN

// <EXT,EXT>
template <class A_coeff, class B_coeff, class Root>
struct Coercion_traits_for_level<Sqrt_extension<A_coeff, Root>,
                           Sqrt_extension<B_coeff, Root>,
                           CTL_SQRT_EXT>{
private:
    typedef Coercion_traits<A_coeff, B_coeff> CT;
    typedef Sqrt_extension<A_coeff,Root> A;
    typedef Sqrt_extension<B_coeff,Root> B;

public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Sqrt_extension<typename CT::Type, Root> Type;

    struct Cast{
    private:
        inline Type cast(const Type& x) const{ return x; }

        template <class T>
        inline Type cast(const T& x) const{
            typename CT::Cast cast;
            if (x.is_extended()) {
                return result_type(cast(x.a0()),cast(x.a1()),x.root());
            } else {
                return result_type(cast(x.a0()));
            }
        }
    public:
        typedef Type result_type;
        // this is in order to allow A and B only
        Type operator()(const A& x) const { return cast(x);}
        Type operator()(const B& x) const { return cast(x);}
    };
};

template <class Coeff, class Root_1, class Root_2>
struct Coercion_traits_for_level<Sqrt_extension<Sqrt_extension<Coeff,Root_1>,
Root_2>,
                           Sqrt_extension<Coeff,Root_1>,
                           CTL_SQRT_EXT>{
private:
    typedef Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2> A;
    typedef Sqrt_extension<Coeff,Root_1> B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    // Type = A
    typedef Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const { return x;}
        Type operator()(const B& x) const { return Type(x);}
    };
};

template <class Coeff, class Root_1, class Root_2>
struct Coercion_traits_for_level
<
            Sqrt_extension<Coeff,Root_1>,
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2>
            ,CTL_SQRT_EXT>
    :public Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2>,
            Sqrt_extension<Coeff,Root_1>
            ,CTL_SQRT_EXT>
{};

template <class Coeff, class Root_1>
struct Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1>,
            Sqrt_extension<Coeff,Root_1>
            ,CTL_SQRT_EXT>{
private:
    typedef  Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1> A;
    typedef  Sqrt_extension<Coeff,Root_1> B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef  Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const { return x;}
        Type operator()(const B& x) const { return Type(x);}
    };
};

template <class Coeff, class Root_1>
struct Coercion_traits_for_level
<
            Sqrt_extension<Coeff,Root_1>,
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1>
            ,CTL_SQRT_EXT>
    :public Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1>,
            Sqrt_extension<Coeff,Root_1>
            ,CTL_SQRT_EXT>
{};


namespace INTERN_CT{
// Coercion_traits for Sqrt_extenison to FieldWithSqrt
template <class A, class B> class CT_ext_to_fwsqrt;
// Coercion_traits for Sqrt_extenison not with FieldWithSqrt
template <class A, class B> class CT_ext_not_to_fwsqrt;
} // namespace INTERN_CT


//<EXT,ANY>
template <class Coeff, class Root, class B>
struct Coercion_traits_for_level<Sqrt_extension<Coeff, Root>, B , CTL_SQRT_EXT>
:public ::boost::mpl::if_c<
             // if B is fwsqrt
              ::boost::is_base_and_derived<
                  Field_with_sqrt_tag,
typename Algebraic_structure_traits<B>::Algebraic_category >::value ||
              ::boost::is_same<
                  Field_with_sqrt_tag,
typename Algebraic_structure_traits<B>::Algebraic_category >::value
            ,
            //then take Intern::Coercion_traits for fwsqrt
            INTERN_CT::CT_ext_to_fwsqrt<Sqrt_extension<Coeff,Root>, B>
            ,
            //else take Intern::Coercion_traits not for fwsqrt
            INTERN_CT::CT_ext_not_to_fwsqrt< Sqrt_extension<Coeff,Root> ,B>
              >::type
{};

// <ANY,EXT>
template <class Coeff, class Root, class B>
struct Coercion_traits_for_level
<B,Sqrt_extension<Coeff, Root>,CTL_SQRT_EXT >
    :public Coercion_traits_for_level<Sqrt_extension<Coeff,Root>,B,CTL_SQRT_EXT>
{};

namespace INTERN_CT{
// EXT coercion with FieldWithSqrt
template <class Coeff, class Root, class FieldWithSqrt>
struct CT_ext_to_fwsqrt<Sqrt_extension<Coeff,Root>,
                                         FieldWithSqrt>{
private:
    typedef Sqrt_extension<Coeff,Root> A;
    typedef FieldWithSqrt B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef FieldWithSqrt Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const {
            typedef Coercion_traits<Coeff,FieldWithSqrt> CT_coeff;
            typedef Coercion_traits<Root ,FieldWithSqrt> CT_root;
            typename CT_coeff::Cast coeff_cast;
            typename CT_root::Cast root_cast;
            if (x.is_extended()) {
                typename CGAL::Algebraic_structure_traits<
                typename CT_root::Type>::Sqrt sqrt;
                return // a0+a1*sqrt(root)
                    coeff_cast(x.a0())+
                    coeff_cast(x.a1())*
                    sqrt(root_cast(x.root()));
            } else {
                return coeff_cast(x.a0());
            }
        }
        Type operator()(const B& x) const { return x;}
    };
};

// EXT coercion not with FieldWithSqrt
template <class Coeff, class Root, class B_>
struct CT_ext_not_to_fwsqrt<Sqrt_extension<Coeff,Root>, B_>{
private:
    typedef Sqrt_extension<Coeff,Root> A;
    typedef B_ B;
    typedef Coercion_traits<Coeff,B> CT;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Sqrt_extension<typename CT::Type,Root> Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const A& x) const {
            typename CT::Cast cast;
            if (x.is_extended()) {
                return Type(cast(x.a0()),cast(x.a1()),x.root());
            } else {
                return Type(cast(x.a0()));
            }
        }
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        }
    };
};
} // namespace INTERN_CT

/////////// COERCION_TRAITS END

// lefthand side
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator +
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num) + p2); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator -
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num) - p2); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator *
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num) * p2); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator /
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num)/p2); }

// righthand side
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator +
(const Sqrt_extension<NT,ROOT>& p1, CGAL_int(NT) num)
{ return (p1 + Sqrt_extension<NT,ROOT>(num)); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator -
(const Sqrt_extension<NT,ROOT>& p1, CGAL_int(NT) num)
{ return (p1 - Sqrt_extension<NT,ROOT>(num)); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator *
(const Sqrt_extension<NT,ROOT>& p1, CGAL_int(NT) num)
{ return (p1 * Sqrt_extension<NT,ROOT>(num)); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator /
(const Sqrt_extension<NT,ROOT>& p1, CGAL_int(NT) num)
{ return (p1 / Sqrt_extension<NT,ROOT>(num)); }

// lefthand side
template <class NT,class ROOT>    bool operator ==
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() == 0 );}
template <class NT,class ROOT>    bool operator !=
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() != 0 );}
template <class NT,class ROOT>    bool operator <
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() < 0 );}
template <class NT,class ROOT>    bool operator <=
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() <= 0 );}
template <class NT,class ROOT>    bool operator >
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() > 0 );}
template <class NT,class ROOT>    bool operator >=
(CGAL_int(NT) num, const Sqrt_extension<NT,ROOT>& p)
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() >= 0 );}

// righthand side
template <class NT,class ROOT>    bool operator ==
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() == 0 );}
template <class NT,class ROOT>    bool operator !=
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() != 0 );}
template <class NT,class ROOT>    bool operator <
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() < 0 );}
template <class NT,class ROOT>    bool operator <=
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() <= 0 );}
template <class NT,class ROOT>    bool operator >
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() > 0 );}
template <class NT,class ROOT>    bool operator >=
(const Sqrt_extension<NT,ROOT>& p, CGAL_int(NT) num)
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() >= 0 );}

CGAL_END_NAMESPACE

#undef CGAL_int

#endif  // CGAL_SQRT_EXTENSION_H

// EOF
