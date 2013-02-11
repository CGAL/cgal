// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Seel <seel@mpi-inf.mpg.de>
//                 
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file CGAL/Polynomial.h
 *  \brief Defines class CGAL::Polynomial.
 *  
 *  Polynomials in one variable (or more, by recursion)
 */

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_TYPE_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_TYPE_H

#define CGAL_icoeff(T) typename CGAL::First_if_different<       \
typename CGAL::internal::Innermost_coefficient_type<T>::Type, T, 1>::Type  

#define CGAL_int(T) typename CGAL::First_if_different< int,   \
typename CGAL::internal::Innermost_coefficient_type<T>::Type , 2>::Type 


#include <CGAL/ipower.h>
#include <cstdio>
#include <sstream>
#include <CGAL/Polynomial/misc.h>

#include <CGAL/use.h>

#ifdef CGAL_HAS_THREADS
#  include <boost/thread/tss.hpp>
#endif

namespace CGAL {

template <class NT> class Polynomial;
template <class NT> class Scalar_factor_traits;
template <class NT> Polynomial<NT> operator - (const Polynomial<NT>& p);

namespace internal {

template <class NT> class Polynomial_rep;

// \brief tag type to distinguish a certain constructor of \c CGAL::Polynomial
class Creation_tag {};

//
// The internal representation class Polynomial_rep<NT>
//

// \brief  internal representation class for \c CGAL::Polynomial
template <class NT_> class Polynomial_rep 
{ 
  typedef NT_ NT;
  typedef std::vector<NT> Vector;
  typedef typename Vector::size_type      size_type;
  typedef typename Vector::iterator       iterator;
  typedef typename Vector::const_iterator const_iterator;
  Vector coeff;

  Polynomial_rep() : coeff() {}
  Polynomial_rep(Creation_tag, size_type s) : coeff(s,NT(0)) {}

  Polynomial_rep(size_type n, ...);

#ifdef CGAL_USE_LEDA
  Polynomial_rep(const LEDA::array<NT>& coeff_)
    : coeff(coeff_.size())
  {
    for (int i = 0; i < coeff_.size(); i++) {
      coeff[i] = coeff_[i + coeff_.low()];
    }
  }
#endif // CGAL_USE_LEDA

  template <class Forward_iterator>
  Polynomial_rep(Forward_iterator first, Forward_iterator last) 
    : coeff(first,last)
  {}

  void reduce() {
    while ( coeff.size()>1 && CGAL::is_zero(coeff.back())) coeff.pop_back();
  }

  void simplify_coefficients() {
    typename Algebraic_structure_traits<NT>::Simplify simplify;
    for (iterator i = coeff.begin(); i != coeff.end(); i++) {
      simplify(*i);
    }
  }

  friend class Polynomial<NT>;
};  // class Polynomial_rep<NT_>

template <class NT>
Polynomial_rep<NT>::Polynomial_rep(size_type n, ...)
  : coeff(n)
{
  // varargs, hence not inline, otherwise g++-3.1 -O2 makes trouble
  va_list ap; va_start(ap, n);
  for (size_type i = 0; i < n; i++) {
    coeff[i] = *(va_arg(ap, const NT*));
  }
  va_end(ap);
}

}// namespace internal

//
// The actual class Polynomial<NT>
//

/*! \ingroup CGAL_Polynomial
  \brief polynomials in one variable (or more, by recursion)

  An instance of the data type \c CGAL::Polynomial represents a
  polynomial <I>p = a<SUB>0</SUB> + a<SUB>1</SUB>*x + ...
  + a<SUB>d</SUB>*x<SUP>d</SUP></I> from the ring
  NT[x]. The data type offers standard ring operations, comparison
  operations (e.g. for symbolic computation with an infimaximal \e x ),
  and various algebraic operations (gcd, resultant).

  \c CGAL:Polynomial offers a full set of algebraic operators, i.e.
  binary +, -, *, / as well as +=, -=, *=, /=; not only for polynomials
  but also for a polynomial and a number of the coefficient type.
  (The / operator must only be used for integral divisions, i.e.
  those with remainder zero.)
  Unary + and - and (in)equality ==, != are also provided.
  If the member function \c sign() (see ibid.) is defined for
  a coefficient type, then so are the comparison operators
  <, >, <=, >= corresponding to the \c sign() of the difference.
  The \c sign() of a polynomial is the \c sign() of its leading
  coefficient, hence comparing by the \c sign() of the difference
  amounts to lexicographic comparison of the coefficient sequence,
  with the coefficient of the highest power taking precedence over
  those of lower powers.

  \c NT must be at least a model of \c IntegralDomainWithoutDiv.
  For all operations naturally involving division, a \c IntegralDomain
  is required. If more than a \c IntegralDomain is required, this is documented.

  \c NT can itself be an instance of \c CGAL::Polynomial, yielding a
  crude form of multivariate polynomials. They behave correctly from an
  algebraic point of view (in particular w.r.t. gcd and resultant
  computation), but not always as a user would expect. For example, the
  leading coefficient of a polynomial from (NT[x])[y] naturally is an
  element of NT[x]. Similarly, computing derivations, resultants, etc.
  always happens w.r.t. the outmost variable.

  Inexact and limited-precision types can be used as coefficients,
  but at the user's risk. The algorithms implemented were written with
  exact number types in mind.

  This data type is implemented as a handle type with value semantics
  (i.e. if you change it, it clones its representation by calling
  \c copy_on_write(), hence no
  aliasing occurs) using \c LiS::Handle (without unification).  The
  coefficients are stored in a vector of \c NT entries. Arithmetic
  operations are implemented naively: + and - need a number of NT
  operations which is linear in the degree while * is quadratic.
  The advanced algebraic operations are implemented with classical
  non-modular methods.

  The important invariant to be preserved by all methods is that
  the coefficient sequence does not contain leading zero coefficients
  (where leading means at the high-degree end), with the excpetion that
  the zero polynomial is represented by a single zero coefficient.
  An empty coefficient sequence denotes an undefined value.

  Many functions modifying a \c CGAL::Polynomial appear both as a member
  function (returning \c void ) which modifies the present object
  and as a non-member function returning a new \c CGAL::Polynomial
  while leaving their argument unchanged. The former is more efficient
  when the old value is no longer referenced elsewhere whereas the
  latter is more convenient.
  \b History: This data type has evolved out of \c RPolynomial 
  from Michael Seel's PhD thesis.  */ 

template <class NT_>
class Polynomial 
  : public Handle_with_policy< internal::Polynomial_rep<NT_> >,
    public boost::ordered_field_operators1< Polynomial<NT_> , 
           boost::ordered_field_operators2< Polynomial<NT_> , NT_ ,  
           boost::ordered_field_operators2< Polynomial<NT_> , CGAL_icoeff(NT_),
           boost::ordered_field_operators2< Polynomial<NT_> , CGAL_int(NT_)  > > > > 
{
  typedef typename internal::Innermost_coefficient_type<NT_>::Type Innermost_coefficient_type; 
public: 

  //! \name Typedefs 
  //@{ 
  //! coefficient type of this instance 
  typedef NT_ NT; 
  //! representation pointed to by this handle 
  typedef internal::Polynomial_rep<NT> Rep;
  //! base class  
  typedef Handle_with_policy< Rep > Base;
  //! container used to store coefficient sequence
  typedef typename Rep::Vector    Vector;
  //! container's size type
  typedef typename Rep::size_type size_type;
  //! container's iterator (random access)
  typedef typename Rep::iterator  iterator;
  //! container's const iterator (random access)
  typedef typename Rep::const_iterator const_iterator;
  //! the Self type
  typedef Polynomial<NT> Self; 
  //@}
  

protected:
  //! \name Protected methods
  //@{
  //! access to the internal coefficient sequence
  Vector& coeffs() { return this->ptr()->coeff; }
  //! const access to the internal coefficient sequence
  const Vector& coeffs() const { return this->ptr()->coeff; }
  //! create an empty polynomial with s coefficients (degree up to s-1)
  Polynomial(internal::Creation_tag f, size_type s)
    : Base(internal::Polynomial_rep<NT>(f,s) )
    {}
    //! non-const access to coefficient \c i
    /*! The polynomial's representation must not be shared between
     *  different handles when this function is called.
     *  This can be ensured by calling \c copy_on_write().
     *
     *  If assertions are enabled, the index \c i is range-checked.
     */
    NT& coeff(unsigned int i) {
      CGAL_precondition(!this->is_shared() && i<(this->ptr()->coeff.size()));
      return this->ptr()->coeff[i]; 
    }

    //! remove leading zero coefficients
    void reduce() { this->ptr()->reduce(); }
    //! remove leading zero coefficients and warn if there were any
    void reduce_warn() {
      CGAL_precondition( this->ptr()->coeff.size() > 0 );
      if (this->ptr()->coeff.back() == NT(0)) {
        CGAL_warning_msg(false, "unexpected degree loss (zero divisor?)");
        this->ptr()->reduce();
      }
    }
    //@}

//
// Constructors of Polynomial<NT>
//
private:
    static Self& get_default_instance(){
      #ifdef CGAL_HAS_THREADS  
        static boost::thread_specific_ptr< Self > safe_x_ptr;
          if (safe_x_ptr.get() == NULL) 
            safe_x_ptr.reset(new Self(0));
        return *safe_x_ptr.get();  
      #else
        static Self x = Self(0);
        return x;
      #endif        
    }
public:
    //! \name Constructors
    //! copy constructor: copy existing polynomial (shares rep)
    Polynomial(const Self& p = get_default_instance()) : Base(static_cast<const Base&>(p)) {}
        
    //! construct the constant polynomial a0 from any type convertible to NT
    template <class T>
      explicit Polynomial(const T& a0)
      : Base(Rep(internal::Creation_tag(), 1))
      { coeff(0) = NT(a0); reduce(); simplify_coefficients(); } 
          
    //! construct the constant polynomial a0
    explicit Polynomial(const NT& a0)
      : Base(Rep(1, &a0))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x
    Polynomial(const NT& a0, const NT& a1)
      : Base(Rep(2, &a0,&a1))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + a2*x^2
    Polynomial(const NT& a0,const NT& a1,const NT& a2)
      : Base(Rep(3, &a0,&a1,&a2))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x + ... + a3*x^3
    Polynomial(const NT& a0,const NT& a1,const NT& a2, const NT& a3)
      : Base(Rep(4, &a0,&a1,&a2,&a3))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + ... + a4*x^4
    Polynomial(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4)
      : Base(Rep(5, &a0,&a1,&a2,&a3,&a4))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x + ... + a5*x^5
    Polynomial(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5)
      : Base(Rep(6, &a0,&a1,&a2,&a3,&a4,&a5))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + ... + a6*x^6
    Polynomial(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5, const NT& a6)
      : Base(Rep(7, &a0,&a1,&a2,&a3,&a4,&a5,&a6))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + ... + a7*x^7
    Polynomial(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5, const NT& a6, const NT& a7)
      : Base(Rep(8, &a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x + ... + a8*x^8
    Polynomial(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5, const NT& a6, const NT& a7,
        const NT& a8)
      : Base(Rep(9, &a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8))
      { reduce(); simplify_coefficients(); }
      
    /*! \brief construct the polynomial whose coefficients are determined
     *  by the iterator range.
     *
     *  Let <TT>a0 = *first</TT>, <TT>a1 = *++first</TT>, ...,
     *  <TT>ad = *it</TT>, where <TT>++it == last</TT>.
     *  Then the polynomial constructed is a0 + a1*x + ... + ad*x<SUP>d</SUP>
     */
    template <class Forward_iterator>
    Polynomial(Forward_iterator first, Forward_iterator last)
      : Base(Rep(first,last)) 
      { reduce(); simplify_coefficients(); }
                              
#if defined(CGAL_USE_LEDA) || defined(DOXYGEN_RUNNING)
    /*! \brief construct a polynomial from a LEDA \c array
     *
     *  The coefficients are determined by \c c[c.low()] up to
     *  \c c[c.high()] of the array \c c. The lowest array element
     *  \c c[c.low()] is always used for the constant term and so on,
     *  irrespective of whether \c c.low() is zero or not.
     */
    Polynomial(const LEDA::array<NT>& c)
      : Base(Rep(c))
      { reduce(); simplify_coefficients(); }
#endif // defined(CGAL_USE_LEDA) || defined(DOXYGEN_RUNNING)
                                
//
// Public member functions
//

public:
    //! \name Public Methods
    //@{

    //! a random access iterator pointing to the constant coefficient
    const_iterator begin() const { return this->ptr()->coeff.begin(); }
    //! a random access iterator pointing beyond the leading coefficient
    const_iterator end()   const { return this->ptr()->coeff.end(); }

    //! the degree \e d of the polynomial
    /*! The degree of the zero polynomial is 0.
     *  The degree of an undefined polynomial is -1.
     */
  int degree() const { return static_cast<int>(this->ptr()->coeff.size())-1; } 

    //! const access to the coefficient of x^i
    const NT& operator[](unsigned int i) const {
      CGAL_precondition( i<(this->ptr()->coeff.size()) );
      return this->ptr()->coeff[i];
    }

    //! the number of non-zero terms of the polynomial
    /*! For an undefined polynomial, this is set to -1. */
    int number_of_terms() const {
      int terms = 0;
      if (degree() < 0) return -1;
      for (int i = 0; i <= degree(); i++) {
        if ((*this)[i] != NT(0)) terms++;
      }
      return terms;
    }

    //! the leading coefficient a<SUB>d</SUB> of the polynomial
    const NT& lcoeff() const {
      CGAL_precondition( this->ptr()->coeff.size() > 0 );
      return this->ptr()->coeff.back();
    }
    

    /*! \brief evaluate the polynomial at \c x
     *
     *  \c x can have another type \c NTX than the coefficient type \c NT.
     *  The result type is defined by CGAL::Coercion_traits<>
     */
    // Note that there is no need to provide a special version for intervals.
    // This was shown by some benchmarks of George Tzoumas, for the 
    // Interval Newton method used in the Voronoi Diagram of Ellipses
    template <class NTX>
      typename Coercion_traits<NTX,NT>::Type 
      evaluate(const NTX& x_) const {
      typedef Coercion_traits<NTX,NT> CT;
      typename CT::Cast cast;
        
      CGAL_precondition( degree() >= 0 );
      int d = degree();
      typename CT::Type x = cast(x_);
      typename CT::Type y=cast(this->ptr()->coeff[d]);
      while (--d >= 0){
        //    y = y*x + cast(this->ptr()->coeff[d]);
        y *= x;
        y += cast(this->ptr()->coeff[d]);
      }
      return y; 
    }
public:
    //! evaluates the polynomial as a homogeneous polynomial
    //! in fact returns evaluate(u/v)*v^degree()
   
    template <class NTX>
      typename Coercion_traits<NTX,NT>::Type 
      evaluate_homogeneous(const NTX& u_, 
          const NTX& v_,
          int hom_degree = -1) const {
      if(hom_degree == -1 ) hom_degree = degree();
      CGAL_precondition( hom_degree >= degree());
      CGAL_precondition( hom_degree >= 0 );
      typedef Coercion_traits<NTX,NT> CT;
      typedef typename CT::Type Type;
      typename CT::Cast cast;

      Type u = cast(u_);
      Type v = cast(v_);
      Type monom;
      Type y(0);
      for(int i = 0; i <= hom_degree; i++){
        monom = CGAL::ipower(v,hom_degree-i)*CGAL::ipower(u,i);
        if(i <= degree())
          y += monom * cast(this->ptr()->coeff[i]);  
      }
      return y;
    }

private:
    // NTX not decomposable
    template <class NTX, class TAG >
      CGAL::Sign sign_at_(const NTX& x, TAG) const{
      CGAL_precondition(degree()>=0);
      return CGAL::sign(evaluate(x));
    }
    // NTX decomposable
    
    template <class NTX>
      CGAL::Sign sign_at_(const NTX& x, CGAL::Tag_true) const{
      CGAL_precondition(degree()>=0);
      typedef Fraction_traits<NTX> FT;
      typedef typename FT::Numerator_type Numerator_type;
      typedef typename FT::Denominator_type Denominator_type;
      Numerator_type num;
      Denominator_type den;
      typename FT::Decompose decompose;
      decompose(x,num,den);
      CGAL_precondition(CGAL::sign(den) == CGAL::POSITIVE);

      typedef Coercion_traits< Numerator_type , Denominator_type > CT;
      typename CT::Cast cast;
      return CGAL::sign(evaluate_homogeneous(cast(num),cast(den)));
    }
public:
    //! evaluates the sign of the Polynomial at x
    template <class NTX>
      CGAL::Sign sign_at(const NTX& x) const{
      // the sign with evaluate_homogeneous is called
      // if NTX is decaomposable
      // and NT would be changed by NTX 
      typedef typename Fraction_traits<NTX>::Is_fraction Is_fraction;
      typedef typename Coercion_traits<NT,NTX>::Type Type;
      typedef typename ::boost::mpl::if_c<
      ::boost::is_same<Type,NT>::value, Is_fraction, CGAL::Tag_false
        >::type If_decomposable_AND_Type_equals_NT;
            
      return sign_at_(x,If_decomposable_AND_Type_equals_NT());
    }
    
    
    // for the benefit of mem_fun1 & friends who don't like const ref args
    template <class NTX>
      typename CGAL::Coercion_traits<NT,NTX>::Type 
      evaluate_arg_by_value(NTX x) const { return evaluate(x); } 

    /*!  \brief evaluate the polynomial with all coefficients replaced by
     *  their absolute values
     *
     *  That is, the function computes <I>|a<SUB>0</SUB>| +
     *  |a<SUB>1</SUB>|*x + ... + |a<SUB>d</SUB>|*x<SUP>d</SUP></I>.
     *  As with \c evaluate(), \c x can be of a type other than \c NT.
     *  \pre Requires \c CGAL::Algebraic_structure_traits::Abs for NT.
     */
   
    template <class NTX> 
      typename Coercion_traits<NTX,NT>::Type 
      evaluate_absolute(const NTX& x) const {
      typedef typename Coercion_traits<NTX,NT>::Type Type;
      typedef typename Coercion_traits<NTX,NT>::Cast Cast;
      Type xx(Cast()(x));
      CGAL_precondition( degree() >= 0 );
      typename Real_embeddable_traits<Type>::Abs abs;
      int d = degree();
      Type y(abs(Cast()(this->ptr()->coeff[d])));
      while (--d >= 0) y = y*xx + abs(Cast()(this->ptr()->coeff[d]));
      return y;
    } 

    /*! \brief evaluate the polynomial with all coefficients replaced by
     *  their absolute values
     *
     *  This is a specialization for \c x is of type CGAL::Interval.
     */
    // TODO: Interval isn't available either!!
/*    Interval evaluate_absolute(const Interval& x) const {
      CGAL_precondition( degree() >= 0 );
      typename Algebraic_structure_traits<Interval>::Abs abs;
      typename Algebraic_structure_traits<NT>::To_Interval to_Interval;
      int d = 0;
      Interval y(to_Interval(this->ptr()->coeff[d]));
      while (++d <= degree()) 
      y+=abs(pow(x,d)*to_Interval(this->ptr()->coeff[d]));
      return y;
      } */
    /*! \brief return the sign of the leading coefficient
     *
     *  For univariate real polynomials, this is the sign
     *  of the limit process \e x --> oo.
     *  Also available as non-member function.
     */
    CGAL::Sign sign() const {
//        CGAL_static_assertion( (boost::is_same< typename Real_embeddable_traits<NT>::Is_real_embeddable,
//                              CGAL::Tag_true>::value) );
      return CGAL::sign(lcoeff());
    }

    //! return sign of difference
    CGAL::Comparison_result compare(const Polynomial& p2) const {
      typename Real_embeddable_traits<NT>::Compare compare;
      typename Real_embeddable_traits<NT>::Sgn sign;
      CGAL_precondition(degree() >= 0);
      CGAL_precondition(p2.degree() >= 0);

      if (this->is_identical(p2)) return CGAL::EQUAL;

      int d1 = degree();
      int d2 = p2.degree();
      if (d1 > d2) {
        return sign((*this)[d1]);
      } else if (d1 < d2) {
        return -sign(p2[d2]);
      }

      for (int i = d1; i >= 0; i--) {
        CGAL::Comparison_result s = compare((*this)[i], p2[i]);
        if (s != CGAL::EQUAL) return s;
      }
      return CGAL::EQUAL;
    }

    //! return sign of difference with constant "polynomial"
    CGAL::Comparison_result compare(const NT& p2) const {
      typename Real_embeddable_traits<NT>::Compare compare;
      typename Real_embeddable_traits<NT>::Sgn sign;
      CGAL_precondition(degree() >= 0);

      if (degree() > 0) {
        return sign(lcoeff());
      } else {
        return compare((*this)[0], p2);
      }
    }

    //! return true iff this is the zero polynomial
    bool is_zero() const
    { return degree()==0 && this->ptr()->coeff[0]==NT(0); }

    //! return \c -p if \c p.sign()<0 and \c p otherwise
    Polynomial<NT> abs() const
    { if ( sign()<0 ) return -*this; return *this; }

    //! return the gcd of all coefficients
    /*! The content is defined as 1 for the zero polynomial. */
    NT content() const {
      CGAL_precondition(degree() >= 0);
      if (is_zero()) return NT(0);
        
      return content_( typename Algebraic_structure_traits< NT >::Algebraic_category() );
    }

private:
    NT content_( Unique_factorization_domain_tag ) const {
      typename Algebraic_structure_traits<NT>::Integral_division idiv;
      typename Algebraic_structure_traits<NT>::Unit_part    upart;
      typename Algebraic_structure_traits<NT>::Gcd          gcd;
      const_iterator it = this->ptr()->coeff.begin(), ite = this->ptr()->coeff.end();
      while (*it == NT(0)) it++;
      NT d = idiv(*it, upart(*it));
      for( ; it != ite; it++) {
        if (d == NT(1)) return d;
        if (*it != NT(0)) d = gcd(d, *it);
      }
      return d;
    }

    NT content_( Field_tag ) const {
      return NT(1);
    }

public:
    
    //! return the unit part of the polynomial
    /*! It is defined as the unit part of the leading coefficient. */
    NT unit_part() const {
      typename Algebraic_structure_traits<NT>::Unit_part upart;
      return upart(lcoeff());
    }

    //! turn p(x) into its derivative p'(x)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void diff() {
      CGAL_precondition( degree() >= 0 );
      if (is_zero()) return;
      this->copy_on_write();
      if (degree() == 0) { coeff(0) = NT(0); return; }
      coeff(0) = coeff(1); // avoid redundant multiplication by NT(1)
      for (int i = 2; i <= degree(); i++) coeff(i-1) = coeff(i) * NT(i);
      this->ptr()->coeff.pop_back();
      reduce(); // in case NT has positive characteristic
    }

    //! replace p(x) by p(a*x)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void scale_up(const NT& a) {
      CGAL_precondition( degree() >= 0 );
      if (degree() == 0) return;
      this->copy_on_write();
      NT a_to_j = a;
      for (int j = 1; j <= degree(); j++) {
        coeff(j) *= a_to_j; 
        a_to_j *= a;
      }
      reduce_warn();
    }

    //! replace p(x) by b<SUP>d</SUP> * p(x/b)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void scale_down(const NT& b)
    {
      CGAL_precondition( degree() >= 0 );
      if (degree() == 0) return;
      this->copy_on_write();
      NT b_to_n_minus_j = b;
      for (int j = degree() - 1; j >= 0; j--) {
        coeff(j) *= b_to_n_minus_j; 
        b_to_n_minus_j *= b;
      }
      reduce_warn();
    }

    //! replace p(x) by b<SUP>d</SUP> * p(a*x/b)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void scale(const NT& a, const NT& b) { scale_up(a); scale_down(b); }

    //! replace p(x) by p(x+1)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void translate_by_one()
    {   // implemented using Ruffini-Horner, see [Akritas, 1989]
      CGAL_precondition( degree() >= 0 );
      this->copy_on_write();
      const int n = degree();
      for (int j = n-1; j >= 0; j--) {
        for (int i = j; i < n; i++) coeff(i) += coeff(i+1); 
      }
    }

    //! replace p(x) by p(x+c)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void translate(const NT& c)
    {   // implemented using Ruffini-Horner, see [Akritas, 1989]
      CGAL_precondition( degree() >= 0 );
      this->copy_on_write();
      const int n = degree();
      for (int j = n-1; j >= 0; j--) {
        for (int i = j; i < n; i++) coeff(i) += c*coeff(i+1);
      }
    }

    //! replace p by b<SUP>d</SUP> * p(x+a/b)
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void translate(const NT& a, const NT& b) 
    {   // implemented using Mehlhorn's variant of Ruffini-Horner
      CGAL_precondition( degree() >= 0 );
      this->copy_on_write();
      const int n = degree();

      NT b_to_n_minus_j = b;
      for (int j = n-1; j >= 0; j--) {
        coeff(j) *= b_to_n_minus_j;
        b_to_n_minus_j *= b;
      }

      for (int j = n-1; j >= 0; j--) {
        coeff(j) += a*coeff(j+1);
        for (int i = j+1; i < n; i++) {
          coeff(i) = b*coeff(i) + a*coeff(i+1); 
        }
        coeff(n) *= b;
      }
      reduce_warn();
    }

    //! replace p by x<SUP>d</SUP> * p(1/x), i.e. reverse the coefficient sequence
    /*! Also available as non-member function which returns the result
     *  instead of overwriting the argument. */
    void reversal() {
      CGAL_precondition( degree() >= 0 );
      this->copy_on_write();
      NT t;
      for (int l = 0, r = degree(); l < r; l++, r--) {
        t = coeff(l); coeff(l) = coeff(r); coeff(r) = t;
      }
      reduce();
    }

    //! divide \e P(x) by \e x , discarding the remainder \e p(0)
    void divide_by_x() {
      CGAL_precondition(degree() >= 0);
      if (is_zero()) return;
      this->copy_on_write();
      for (int i = 0; i < degree(); i++) {
        coeff(i) = coeff(i+1);
      }
      coeffs().pop_back();
    }

    //! invoke \c CGAL::Algebraic_structure_traits::Simplify on all coefficients
    void simplify_coefficients() { this->ptr()->simplify_coefficients(); }

    //! write polynomial to \c os in \c LiS::IO::PRETTY format
    /*! The output is intended to be Maple-readable; see module
     *  \link CGAL_io CGAL I/O Support \endlink.
     *
     * Example: A \c CGAL::Polynomial<int> with a value of
     * 4<I>x</I><SUP>2</SUP> - 1 will be written as
     * <TT> 4*x^2 + (-1) </TT> by this function.
     */
    void output_maple(std::ostream& os) const;
    //! write polynomial to \c os in a format readable by \c input_ascii()
    void output_ascii(std::ostream& os) const;
    //! write polynomial to \c os in \c LiS::IO::BENCHMARK format
    void output_benchmark(std::ostream& os) const;

    //! implement \c CGAL::Scalar_factor_traits::Scalar_div for polynomials
    void scalar_div(const typename
        Scalar_factor_traits< Polynomial<NT> >::Scalar& b) {
      typename Scalar_factor_traits<NT>::Scalar_div sdiv;
      this->copy_on_write();
      for (int i = degree(); i >= 0; --i) {
        sdiv(coeff(i), b);
      }
    };

    //@}

//
// Static member functions
//

    //! \name Static member functions
    //@{

    /*! \brief division with remainder on polynomials
     *
     * Given \c f and \c g, compute quotient \c q and remainder \c r
     * such that <I>f = g*q + r</I> and deg(<I>r</I>) < deg(<I>g</I>).
     *
     * \pre \c g!=0. NT is a field, or \c f and \c g are such that
     * the division can be performed in NT anyway.
     */
    static void euclidean_division (const Polynomial<NT>& f,
        const Polynomial<NT>& g,
        Polynomial<NT>& q, Polynomial<NT>& r);

    /*! \brief pseudo division with remainder on polynomials
     *
     * Given \c f and \c g != 0, compute quotient \c q and remainder \c r
     * such that <I>D*f = g*q + r</I> and deg(<I>r</I>) < deg(<I>g</I>),
     * where \e D = lcoeff(<I>g</I>)^max(0, deg(<I>f</I>)-deg(<I>g</I>)+1)
     *
     * This is similar to \c euclidean_division() except that multiplying
     * by \e D makes sure that the division can be performed over any ring.
     */
    static void pseudo_division(const Polynomial<NT>& f,
        const Polynomial<NT>& g, 
        Polynomial<NT>& q, Polynomial<NT>& r, NT& D);


    /*! \brief read a polynomial from \c is
     *
     * The expected format is:
     * <TT><B>P[</B></TT><I>d</I> <TT><B>(</B></TT><I>i</I><TT><B>,</B></TT>
     * <I>ai</I><TT><B>)</B></TT> ... <TT><B>]</B></TT> 
     * with coefficients in arbitrary order (but without
     * repetitions). <I>d</I> is the degree and <I>ai</I> is the coefficient
     * of <I>x<SUP>i</SUP></I>. Missing coefficients are set to zero.
     * Whitespace is ignored.
     * The format of the coefficients must be understandable for
     * <TT> is >> iformat(ai) </TT>.
     *
     * Example: A \c CGAL::Polynomial<int> with a value of
     * 4<I>x</I><SUP>2</SUP> - 1 has to be written as
     * \c P[2(2,4)(0,-1)] or \c P[2(2,4)(1,0)(0,-1)]
     * or similarly with permuted coefficients.
     */
    static Polynomial<NT> input_ascii(std::istream& is);

    //@}

//
// Arithmetic Operations, Part II:
// implementing two-address arithmetic (incl. mixed-mode) by member functions
//

// ...for polynomials
    Polynomial<NT>& operator += (const Polynomial<NT>& p1) {
      this->copy_on_write();
      int d = (std::min)(degree(),p1.degree()), i;
      for(i=0; i<=d; ++i) coeff(i) += p1[i];
      while (i<=p1.degree()) this->ptr()->coeff.push_back(p1[i++]);
      reduce(); return (*this);
    }

    Polynomial<NT>& operator -= (const Polynomial<NT>& p1) 
      {
        this->copy_on_write();
        int d = (std::min)(degree(),p1.degree()), i;
        for(i=0; i<=d; ++i) coeff(i) -= p1[i];
        while (i<=p1.degree()) this->ptr()->coeff.push_back(-p1[i++]);
        reduce(); return (*this);
      }

    Polynomial<NT>& operator *= (const Polynomial<NT>& p2)
      { 
        // TODO: use copy on write 
        Polynomial<NT> p1 = (*this);
        typedef typename Polynomial<NT>::size_type size_type;
        CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
        internal::Creation_tag TAG;
        Polynomial<NT>  p(TAG, size_type(p1.degree()+p2.degree()+1) ); 
        // initialized with zeros
        for (int i=0; i <= p1.degree(); ++i)
          for (int j=0; j <= p2.degree(); ++j)
            p.coeff(i+j) += (p1[i]*p2[j]); 
        p.reduce();
        return (*this) = p ;
      }

    Polynomial<NT>& operator /= (const Polynomial<NT>& p2)
      { 
        // TODO: use copy on write 
        CGAL_precondition(!p2.is_zero());
        if ((*this).is_zero()) return (*this);

        Polynomial<NT> p1 = (*this);
        typedef Algebraic_structure_traits< Polynomial<NT> > AST; 
        // Precondition: q with p1 == p2 * q must exist within NT[x].
        // If this holds, we can perform Euclidean division even over a ring NT
        // Proof: The quotients of each division that occurs are precisely
        //   the terms of q and hence in NT.
        Polynomial<NT> q, r;
        Polynomial<NT>::euclidean_division(p1, p2, q, r);
        CGAL_USE_TYPE(AST);
        CGAL_postcondition( !AST::Is_exact::value || p2 * q == p1);
        return (*this) = q;
      }


    // ...and in mixed-mode arithmetic
    Polynomial<NT>& operator += (const NT& num)
      { this->copy_on_write(); coeff(0) += (NT)num; return *this; }

    Polynomial<NT>& operator -= (const NT& num)
      { this->copy_on_write(); coeff(0) -= (NT)num; return *this; }

    Polynomial<NT>& operator *= (const NT& num) {
      CGAL_precondition(degree() >= 0);
      this->copy_on_write();
      for(int i=0; i<=degree(); ++i) coeff(i) *= (NT)num; 
      reduce();
      return *this;
    }

    Polynomial<NT>& operator /= (const NT& num)
      {
        CGAL_precondition(num != NT(0));
        CGAL_precondition(degree() >= 0);
        if (is_zero()) return *this;
        this->copy_on_write();
        typename Algebraic_structure_traits<NT>::Integral_division idiv;
        for(int i = 0; i <= degree(); ++i) coeff(i) = idiv(coeff(i), num);
        reduce_warn();
        return *this;
      }// ...and in mixed-mode arithmetic
                              
    // TODO: avoid  NT(num)
    Polynomial<NT>& operator += (CGAL_int(NT) num)
      { return *this += NT(num); } 
    Polynomial<NT>& operator -= (CGAL_int(NT) num)
      { return *this -= NT(num); } 
    Polynomial<NT>& operator *= (CGAL_int(NT) num)
      { return *this *= NT(num); } 
    Polynomial<NT>& operator /= (CGAL_int(NT) num)
      { return *this /= NT(num); }  
                              
    // TODO: avoid  NT(num)
    Polynomial<NT>& operator += (const CGAL_icoeff(NT)& num)
      { return *this += NT(num); } 
    Polynomial<NT>& operator -= (const CGAL_icoeff(NT)& num)
      { return *this -= NT(num); } 
    Polynomial<NT>& operator *= (const CGAL_icoeff(NT)& num)
      { return *this *= NT(num); } 
    Polynomial<NT>& operator /= (const CGAL_icoeff(NT)& num)
      { return *this /= NT(num); }

    // special operation to implement (pseudo-)division and the like
    void minus_offsetmult(const Polynomial<NT>& p, const NT& b, int k)
    {
      CGAL_precondition(!this->is_shared());
      int pd = p.degree();
      CGAL_precondition(degree() >= pd+k);
      for (int i = 0; i <= pd; i++) coeff(i+k) -= b*p[i];
      reduce();
    }
    
    friend Polynomial<NT> operator - <> (const Polynomial<NT>&);   
}; // class Polynomial<NT_>

// Arithmetic Operators, Part III:
// implementation of unary operators and three-address arithmetic
// by friend functions
//

template <class NT> inline
Polynomial<NT> operator + (const Polynomial<NT>& p) {
  CGAL_precondition(p.degree() >= 0);
  return p;
}

template <class NT> inline
Polynomial<NT> operator - (const Polynomial<NT>& p) {
  CGAL_precondition(p.degree()>=0);
  Polynomial<NT> res(p.coeffs().begin(),p.coeffs().end());
  typename Polynomial<NT>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}




template <class NT> inline
Polynomial<NT> operator * (const Polynomial<NT>& p1, 
    const Polynomial<NT>& p2)
{
  typedef typename Polynomial<NT>::size_type size_type;
  CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
  internal::Creation_tag TAG;
  Polynomial<NT>  p(TAG, size_type(p1.degree()+p2.degree()+1) ); 
  // initialized with zeros
  for (int i=0; i <= p1.degree(); ++i)
    for (int j=0; j <= p2.degree(); ++j)
      p.coeff(i+j) += (p1[i]*p2[j]); 
  p.reduce();
  return p;
}


//
// Comparison Operators
//

// polynomials only
template <class NT> inline
bool operator == (const Polynomial<NT>& p1, const Polynomial<NT>& p2) {
  CGAL_precondition(p1.degree() >= 0);
  CGAL_precondition(p2.degree() >= 0);
  if (p1.is_identical(p2)) return true;
  if (p1.degree() != p2.degree()) return false;
  for (int i = p1.degree(); i >= 0; i--) if (p1[i] != p2[i]) return false;
  return true;
}
template <class NT> inline
bool operator < (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) < 0 ); } 
template <class NT> inline
bool operator > (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) > 0 ); } 

// operators NT
template <class NT> inline 
bool operator == (const NT& num, const Polynomial<NT>& p) {
  CGAL_precondition(p.degree() >= 0);
  return p.degree() == 0 && p[0] == num;
}
template <class NT> inline
bool operator == (const Polynomial<NT>& p, const NT& num)  {
  CGAL_precondition(p.degree() >= 0);
  return p.degree() == 0 && p[0] == num;
}
template <class NT> inline
bool operator < (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) > 0 );}
template <class NT> inline
bool operator < (const Polynomial<NT>& p,const NT& num) 
{ return ( p.compare(num) < 0 );}
template <class NT> inline
bool operator > (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) < 0 );}
template <class NT> inline
bool operator > (const Polynomial<NT>& p,const NT& num) 
{ return ( p.compare(num) > 0 );}


// compare int #################################
template <class NT> inline
bool operator == (const CGAL_int(NT)& num, const Polynomial<NT>& p)  {
  CGAL_precondition(p.degree() >= 0);
  return p.degree() == 0 && p[0] == NT(num);
}
template <class NT> inline
bool operator == (const Polynomial<NT>& p, const CGAL_int(NT)& num)  {
  CGAL_precondition(p.degree() >= 0);
  return p.degree() == 0 && p[0] == NT(num);
}
template <class NT> inline
bool operator < (const CGAL_int(NT)& num, const Polynomial<NT>& p) 
{ return ( p.compare(NT(num)) > 0 );}
template <class NT> inline
bool operator < (const Polynomial<NT>& p, const CGAL_int(NT)& num) 
{ return ( p.compare(NT(num)) < 0 );}
template <class NT> inline
bool operator > (const CGAL_int(NT)& num, const Polynomial<NT>& p) 
{ return ( p.compare(NT(num)) < 0 );}
template <class NT> inline
bool operator > (const Polynomial<NT>& p, const CGAL_int(NT)& num) 
{ return ( p.compare(NT(num)) > 0 );}

// compare icoeff ###################################
template <class NT> inline
bool operator == (const CGAL_icoeff(NT)& num, const Polynomial<NT>& p)  {
  CGAL_precondition(p.degree() >= 0);
  return p.degree() == 0 && p[0] == NT(num);
}
template <class NT> inline
bool operator == (const Polynomial<NT>& p, const CGAL_icoeff(NT)& num)  {
  CGAL_precondition(p.degree() >= 0);
  return p.degree() == 0 && p[0] == NT(num);
}
template <class NT> inline
bool operator < (const CGAL_icoeff(NT)& num, const Polynomial<NT>& p) 
{ return ( p.compare(NT(num)) > 0 );}
template <class NT> inline
bool operator < (const Polynomial<NT>& p, const CGAL_icoeff(NT)& num) 
{ return ( p.compare(NT(num)) < 0 );}


template <class NT> inline
bool operator > (const CGAL_icoeff(NT)& num, const Polynomial<NT>& p) 
{ return ( p.compare(NT(num)) < 0 );}
template <class NT> inline
bool operator > (const Polynomial<NT>& p, const CGAL_icoeff(NT)& num) 
{ return ( p.compare(NT(num)) > 0 );}

//
// Algebraically non-trivial operations
//

// 1) Euclidean and pseudo-division of polynomials
// (implementation of static member functions)

template <class NT>
void Polynomial<NT>::euclidean_division(
    const Polynomial<NT>& f, const Polynomial<NT>& g,
    Polynomial<NT>& q, Polynomial<NT>& r)
{
  typedef Algebraic_structure_traits<NT> AST;
  typename AST::Integral_division idiv;
  int fd = f.degree(), gd = g.degree();
  if ( fd < gd ) {
    q = Polynomial<NT>(NT(0)); r = f;

    CGAL_postcondition( !AST::Is_exact::value || f == q*g + r); 
    return;
  }
  // now we know fd >= gd 
  int qd = fd-gd, delta = qd+1, rd = fd;

  internal::Creation_tag TAG;    
  q = Polynomial<NT>(TAG, delta ); 
  r = f; r.copy_on_write();
  while ( qd >= 0 ) {
    NT Q = idiv(r[rd], g[gd]);
    q.coeff(qd) += Q;
    r.minus_offsetmult(g,Q,qd);
    r.simplify_coefficients();
    if (r.is_zero()) break;
    rd = r.degree();
    qd = rd - gd;
  }
  q.simplify_coefficients();

  CGAL_postcondition( !AST::Is_exact::value || f == q*g + r);
}

#ifndef CGAL_POLY_USE_OLD_PSEUDODIV

template <class NT>
void Polynomial<NT>::pseudo_division(
    const Polynomial<NT>& A, const Polynomial<NT>& B,
    Polynomial<NT>& Q, Polynomial<NT>& R, NT& D)
{
  typedef Algebraic_structure_traits<NT> AST;
  // pseudo-division with incremental multiplication by lcoeff(B)
  // see [Cohen, 1993], algorithm 3.1.2

  CGAL_precondition(!B.is_zero());
  int delta = A.degree() - B.degree();

  if (delta < 0 || A.is_zero()) {
    Q = Polynomial<NT>(NT(0)); R = A; D = NT(1);
       
    CGAL_USE_TYPE(AST);
    CGAL_postcondition( !AST::Is_exact::value || Polynomial<NT>(D)*A == Q*B + R);
    return;
  }
  const NT d = B.lcoeff();
  int e = delta + 1;
  D = CGAL::ipower(d, e);
  internal::Creation_tag TAG;
  Q = Polynomial<NT>(TAG, e);
  R = A; R.copy_on_write(); R.simplify_coefficients();

  // invariant: d^(deg(A)-deg(B)+1 - e) * A == Q*B + R
  do { // here we have delta == R.degree() - B.degree() >= 0 && R != 0
    NT lR = R.lcoeff();
    for (int i = delta+1; i <= Q.degree(); i++) Q.coeff(i) *= d;
    Q.coeff(delta) = lR;              // Q = d*Q + lR * X^delta
    for (int i = 0; i <= R.degree(); i++) R.coeff(i) *= d;
    R.minus_offsetmult(B, lR, delta); // R = d*R - lR * X^delta * B
    R.simplify_coefficients();
    e--;
    delta = R.degree() - B.degree();
  } while (delta > 0 || (delta == 0 && !R.is_zero()));
  // funny termination condition because deg(0) = 0, not -\infty

  NT q = CGAL::ipower(d, e);
  Q *= q; Q.simplify_coefficients();
  R *= q; R.simplify_coefficients();

  CGAL_postcondition( !AST::Is_exact::value || Polynomial<NT>(D)*A == Q*B + R);
}

#else

template <class NT>
void Polynomial<NT>::pseudo_division(
    const Polynomial<NT>& f, const Polynomial<NT>& g, 
    Polynomial<NT>& q, Polynomial<NT>& r, NT& D)
{
  typedef Algebraic_structure_traits<NT> AST;
  // pseudo-division with one big multiplication with lcoeff(g)^{...}
  typename Algebraic_structure_traits<NT>::Integral_division idiv;

  int fd=f.degree(), gd=g.degree();
  if ( fd < gd ) {
    q = Polynomial<NT>(NT(0)); r = f; D = NT(1); 

    CGAL_postcondition( !AST::Is_exact::value  || Polynomial<NT>(D)*f==q*g+r);
    return;
  }
  // now we know rd >= gd 
  int qd = fd-gd, delta = qd+1, rd = fd;
  internal::Creation_tag TAG;
  q = Polynomial<NT>(TAG, delta );
  NT G = g[gd]; // highest order coeff of g
  D = CGAL::ipower(G, delta);
  Polynomial<NT> res = D*f;
  res.simplify_coefficients();
  while ( qd >= 0 ) {
    NT F = res[rd];    // highest order coeff of res
    NT t = idiv(F, G); // sure to be integral by multiplication of D
    q.coeff(qd) = t;   // store q coeff
    res.minus_offsetmult(g,t,qd); 
    res.simplify_coefficients();
    if (res.is_zero()) break;
    rd = res.degree();
    qd = rd - gd;
  }
  r = res; // already simplified
  q.simplify_coefficients();

  CGAL_postcondition( !AST::Is_exact::value  || Polynomial<NT>(D)*f==q*g+r);
}

template <class NT> inline
Polynomial<NT> division(const Polynomial<NT>& p1, 
    const Polynomial<NT>& p2,
    Integral_domain_tag)
{
  typedef Algebraic_structure_traits<NT> AST;
  CGAL_precondition(!p2.is_zero());
  if ( p1.is_zero() ) return p1;
  Polynomial<NT> q,r; NT D;
  Polynomial<NT>::pseudo_division(p1,p2,q,r,D);
  q/=D;

  CGAL_postcondition( !AST::Is_exact::value || p2 * q == p1);
  return q;
}

template <class NT> inline
Polynomial<NT> division(const Polynomial<NT>& p1, 
    const Polynomial<NT>& p2,
    Field_tag)
{
  typedef Algebraic_structure_traits<NT> AST;
  CGAL_precondition(!p2.is_zero());
  if (p1.is_zero()) return p1;
  Polynomial<NT> q,r;
  Polynomial<NT>::euclidean_division(p1,p2,q,r);
  CGAL_postcondition( !AST::Is_exact::value  || p2 * q == p1);
  return q;
}

#endif // CGAL_POLY_USE_OLD_PSEUDODIV

//
// I/O Operations
//

/*! \ingroup CGAL_Polynomial
 *  \relates CGAL::Polynomial
 *  \brief output \c p to \c os
 *
 *  Output \c p in a format as specified by
 *  \c LiS::get_mode(os), see \link LiS_io LiS I/O Support \endlink.
 *  Currently, the output for \c LiS::IO::BINARY happens to be
 *  identical to \c LiS::IO::ASCII.
 */
template <class NT>
std::ostream& operator << (std::ostream& os, const Polynomial<NT>& p) {
  switch(CGAL::get_mode(os)) {
  case CGAL::IO::PRETTY:
    p.output_maple(os); break;
  default:
    p.output_ascii(os); break;
  }
  return os;
}

/*! \ingroup CGAL_Polynomial
 *  \relates CGAL::Polynomial
 *  \brief try to read a polynomial from \c is into \c p
 *
 *  \c is must be in a mode that supports input of polynomials
 *  (\c LiS::IO::ASCII or \c LiS::IO::BINARY) and the input from
 *  \c is must have the format of output to a stream of the same mode.
 */
template <class NT>
std::istream& operator >> (std::istream& is, Polynomial<NT>& p) {
  CGAL_precondition(!CGAL::is_pretty(is));
  p = Polynomial<NT>::input_ascii(is);
  return is;
}


template <class NT> inline
void print_maple_monomial(std::ostream& os, const NT& coeff,
    const char *var, int expn)
{
  if (expn == 0 || coeff != NT(1)) {
    os << CGAL::oformat(coeff, Parens_as_product_tag());
    if (expn >= 1) os << "*";
  }
  if (expn >= 1) {
    os << var;
    if (expn > 1) os << "^" << CGAL::oformat(expn);
  }
}

// fwd declaration of Polynomial_traits_d
template <typename Polynomial_d> class Polynomial_traits_d;

template <class NT>
void Polynomial<NT>::output_maple(std::ostream& os) const {
  const Polynomial<NT>& p = *this;
  const char *varname;
  char vnbuf[42];
    
  // use variable names x, y, z, w1, w2, w3, ...
  if (Polynomial_traits_d<NT>::d < 3) {
    static const char *varnames[] = { "x", "y", "z" };
    varname = varnames[Polynomial_traits_d<NT>::d];
  } else {
    std::sprintf(vnbuf, "w%d", Polynomial_traits_d<NT>::d - 2);
    varname = vnbuf;
  }
    
  int i = p.degree();
  CGAL::print_maple_monomial(os, p[i], varname, i);
  while (--i >= 0) {
    if (p[i] != NT(0)) {
      os << " + ";
      CGAL::print_maple_monomial(os, p[i], varname, i);
    }
  }
}

template <class NT>
void Polynomial<NT>::output_ascii(std::ostream &os) const {
  const Polynomial<NT> &p = *this;
  if (p.is_zero()) { os << "P[0 (0," << oformat(NT(0)) << ")]"; return; }

  os << "P[" << oformat(p.degree());
  for (int i = 0; i <= p.degree(); i++) {
    if (p[i] != NT(0)) os << "(" << CGAL::oformat(i) << ","
                          << CGAL::oformat(p[i]) << ")";
  }
  os << "]";
}

template <class NT>
void Polynomial<NT>::output_benchmark(std::ostream &os) const {
  typedef typename Polynomial_traits_d< Polynomial<NT> >::Innermost_coefficient_type 
    Innermost_coefficient_type;
  typedef std::pair< Exponent_vector, Innermost_coefficient_type >
    Exponents_coeff_pair;
  typedef typename Polynomial_traits_d< Polynomial<NT> >::Monomial_representation Gmr;
    
  std::vector< Exponents_coeff_pair > monom_rep;
  Gmr gmr;
  gmr( *this, std::back_inserter( monom_rep ) );
    
  os << Benchmark_rep< Polynomial< NT > >::get_benchmark_name() << "( ";
    
  for( typename std::vector< Exponents_coeff_pair >::iterator it = monom_rep.begin();
       it != monom_rep.end(); ++it ) {
    if( it != monom_rep.begin() )
      os << ", ";
    os << "( " << bmformat( it->second ) << ", ";
    it->first.output_benchmark(os);
    os << " )";        
  }
  os << " )";
}

// Benchmark_rep specialization 
template < class NT >
class Benchmark_rep< CGAL::Polynomial< NT > > {
  const CGAL::Polynomial< NT >& t;
public:
  //! initialize with a const reference to \a t.
  Benchmark_rep( const CGAL::Polynomial< NT >& tt) : t(tt) {}
  //! perform the output, calls \c operator\<\< by default.
  std::ostream& operator()( std::ostream& out) const { 
    t.output_benchmark( out );
    return out;
  }
    
  static std::string get_benchmark_name() {
    std::stringstream ss;
    ss << "Polynomial< " << Polynomial_traits_d< Polynomial< NT > >::d;
        
    std::string coeff_name = Benchmark_rep< NT >::get_benchmark_name();
        
    if( coeff_name != "" )
      ss << ", " << coeff_name;
        
    ss << " >";
    return ss.str();
  }
};

// Moved to internal namespace because of name clashes
// TODO: Is this OK?
namespace internal {

inline static void swallow(std::istream &is, char d) {
  char c;
  do c = is.get(); while (isspace(c));
  if (c != d) CGAL_error_msg( "input error: unexpected character in polynomial");
}
} // namespace internal

template <class NT>
Polynomial<NT> Polynomial<NT>::input_ascii(std::istream &is) {
  char c;
  int degr = -1, i=0;

  internal::swallow(is, 'P');
  internal::swallow(is, '[');
  is >> CGAL::iformat(degr);
  if (degr < 0) {
    CGAL_error_msg( "input error: negative degree of polynomial specified");
  }
  internal::Creation_tag TAG;
  Polynomial<NT> p(TAG, degr+1);

  do c = is.get(); while (isspace(c));
  do {
    if (c != '(') CGAL_error_msg( "input error: ( expected");
    is >> CGAL::iformat(i);
    if (!(i >= 0 && i <= degr && p[i] == NT(0))) {
      CGAL_error_msg( "input error: invalid exponent in polynomial");
    };
    internal::swallow(is, ',');
    is >> CGAL::iformat(p.coeff(i));
    internal::swallow(is, ')');
    do c = is.get(); while (isspace(c));
  } while (c != ']');

  p.reduce();
  p.simplify_coefficients();
  return p;
}

template <class COEFF>
struct Needs_parens_as_product<Polynomial<COEFF> >{
  typedef Polynomial<COEFF> Poly;
  bool operator()(const Poly& x){ return (x.degree() > 0); }
};



} //namespace CGAL

#undef CGAL_icoeff
#undef CGAL_int

#endif // CGAL_POLYNOMIAL_POLYNOMIAL_TYPE_H
