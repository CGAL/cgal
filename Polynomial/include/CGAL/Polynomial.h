// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Seel <seel@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Polynomial.h
 *  \brief Defines class NiX::Polynomial.
 *  
 *  Polynomials in one variable (or more, by recursion)
 */

#ifndef CGAL_POLYNOMIAL_H
#define CGAL_POLYNOMIAL_H

#include <cstdarg>
#include <cctype>
#include <vector>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/Flattening_iterator.h>
//#include <NiX/Modular.h>

#include <CGAL/Exponent_vector.h>

#include <boost/static_assert.hpp>

#ifdef CGAL_USE_LEDA
#include <LEDA/array.h>
#endif // CGAL_USE_LEDA

CGAL_BEGIN_NAMESPACE

// TODO: copied from EXACUS/NumeriX/include/NiX/number_type_utils.h
template <typename NT>
inline
NT ipower(const NT& base, int expn) {
    // compute base^expn using square-and-multiply
    CGAL_precondition(expn >= 0);

    // handle trivial cases efficiently
    if (expn == 0) return NT(1);
    if (expn == 1) return base;

    // find the most significant non-zero bit of expn
    int e = expn, msb = 0;
    while (e >>= 1) msb++;

    // computing base^expn by square-and-multiply
    NT res = base;
    int b = 1<<msb;
    while (b >>= 1) { // is there another bit right of what we saw so far?
        res *= res;
        if (expn & b) res *= base;
    }
    return res;
}

// TODO: END included from number_type_utils.h

template <class NT> class Polynomial;
template <class NT> class Scalar_factor_traits;

template <class NT>
Polynomial<NT> operator - (const Polynomial<NT>& p);

template <class NT> inline
Polynomial<NT> operator + (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2);

template <class NT> inline
Polynomial<NT> operator - (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2);

template <class NT> inline
Polynomial<NT> operator * (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2);
template <class NT> inline
Polynomial<NT> operator / (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2);

namespace INTERN_POLYNOMIAL {

template <class NT> class Polynomial_rep;

// \brief tag type to distinguish a certain constructor of \c NiX::Polynomial
class Creation_tag {};

//
// The internal representation class Polynomial_rep<NT>
//

// \brief  internal representation class for \c NiX::Polynomial
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
        while ( coeff.size()>1 && coeff.back()==NT(0) ) coeff.pop_back();
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

}// namespace INTERN_POLYNOMIAL

//
// The actual class Polynomial<NT>
//

/*! \ingroup NiX_Polynomial
    \brief polynomials in one variable (or more, by recursion)

An instance of the data type \c NiX::Polynomial represents a
polynomial <I>p = a<SUB>0</SUB> + a<SUB>1</SUB>*x + ...
+ a<SUB>d</SUB>*x<SUP>d</SUP></I> from the ring
NT[x]. The data type offers standard ring operations, comparison
operations (e.g. for symbolic computation with an infimaximal \e x ),
and various algebraic operations (gcd, resultant).

\c NiX:Polynomial offers a full set of algebraic operators, i.e.
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

\c NT can itself be an instance of \c NiX::Polynomial, yielding a
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

Many functions modifying a \c NiX::Polynomial appear both as a member
function (returning \c void ) which modifies the present object
and as a non-member function returning a new \c NiX::Polynomial
while leaving their argument unchanged. The former is more efficient
when the old value is no longer referenced elsewhere whereas the
latter is more convenient.

\b History: This data type has evolved out of \c RPolynomial 
from Michael Seel's PhD thesis.  */ 
template <class NT_>
class Polynomial 
  : public Handle_with_policy< INTERN_POLYNOMIAL::Polynomial_rep<NT_> >
{
public: 
    //! \name Typedefs 
    //@{ 
    //! coefficient type of this instance 
    typedef NT_ NT; 
    //! representation pointed to by this handle 
    typedef INTERN_POLYNOMIAL::Polynomial_rep<NT> Rep;
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
    //@}

protected:
    //! \name Protected methods
    //@{
    //! access to the internal coefficient sequence
    Vector& coeffs() { return this->ptr()->coeff; }
    //! const access to the internal coefficient sequence
    const Vector& coeffs() const { return this->ptr()->coeff; }
    //! create an empty polynomial with s coefficients (degree up to s-1)
    Polynomial(INTERN_POLYNOMIAL::Creation_tag f, size_type s)
        : Base(INTERN_POLYNOMIAL::Polynomial_rep<NT>(f,s) )
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
            CGAL_warning("unexpected degree loss (zero divisor?)");
            this->ptr()->reduce();
        }
    }
    //@}

//
// Constructors of Polynomial<NT>
//

public:
    //! \name Constructors
    //@{

    //! default constructor: a new polynomial of undefined value.
    Polynomial() : Base( Rep(INTERN_POLYNOMIAL::Creation_tag(), 1) ) 
    { coeff(0) = NT(0); }
    

    //! copy constructor: copy existing polynomial (shares rep)
    Polynomial(const Polynomial<NT>& p) : Base(static_cast<const Base&>(p)) {}

    //! construct the constant polynomial a0 from any type convertible to NT
    template <class T>
    explicit Polynomial(const T& a0)
        : Base(Rep(INTERN_POLYNOMIAL::Creation_tag(), 1))
    { coeff(0) = NT(a0); reduce(); simplify_coefficients(); }
   
    //! construct the constant polynomial a0
    Polynomial(const NT& a0)
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
    //@}


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
    int degree() const { return this->ptr()->coeff.size()-1; } 

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
     *  The result type is defined by NiX::Coercion_traits<>
     */

    template <class NTX>
    typename Coercion_traits<NTX,NT>::Type 
    evaluate(const NTX& x) const {
        typedef Coercion_traits<NTX,NT> CT;
        typename CT::Cast cast;
    
        CGAL_precondition( degree() >= 0 );
        int d = degree();
    
        typename CT::Type y=cast(this->ptr()->coeff[d]);
        while (--d >= 0) 
            y = y*cast(x) + cast(this->ptr()->coeff[d]);
        return y; 
    }
private:
    
    template <class ToInterval, class ANY >
    struct Interval_evaluation_traits{
        typedef typename ToInterval::result_type RET;
    };
    template <class ANY >
    struct Interval_evaluation_traits<CGAL::Null_functor,ANY>{
        typedef CGAL::Null_functor RET;
    }; 
    
public:
    /*! \brief evaluate the polynomial at \c x
     *
     *  This is a specialization for \c x is of type NiX::Interval.
     */

    // TODO: This is a specialization for "evaluate" which handles the special
    //       case where the coefficients are Intervals better. The class also works
    //       without this specialization, but maybe slower in some cases.
    //       But since it would be better to have a general case for all interval
    //       types and not only for the CGAL::Interval_nt it should be replaced
    //       by an evaluate function with an additional hint-parameter.
/*    typename Interval_evaluation_traits<
        typename Real_embeddable_traits<NT>::To_Interval,int>::RET
    evaluate(const Interval& x) const {        
        typedef typename Real_embeddable_traits<NT>::To_Interval To_Interval; 
        To_Interval to_Interval;
        typedef typename Interval_evaluation_traits<To_Interval,int>::RET RET;
        CGAL_precondition( degree() >= 0 );
        int d = 0;
        RET y=to_Interval(this->ptr()->coeff[d]);
        while (++d <= degree()) 
            y += pow(x,d)*to_Interval(this->ptr()->coeff[d]);
        return y;
    }*/
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
            monom = ipower(v,hom_degree-i)*ipower(u,i);
            if(i <= degree())
                y += monom * cast(this->ptr()->coeff[i]);  
        }
        return y;
    }

private:
    // NTX not decomposable
    template <class NTX, class TAG >
    CGAL::Sign sign_at_(const NTX& x, TAG tag) const{
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
        CGAL_precondition(CGAL::Sign(den) == CGAL::POSITIVE);

        typedef Coercion_traits< Numerator_type , Denominator_type > CT;
        typename CT::Cast cast;
        return CGAL::Sign(evaluate_homogeneous(cast(num),cast(den)));
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
    NTX evaluate_arg_by_value(NTX x) const { return evaluate(x); } 

    /*!  \brief evaluate the polynomial with all coefficients replaced by
     *  their absolute values
     *
     *  That is, the function computes <I>|a<SUB>0</SUB>| +
     *  |a<SUB>1</SUB>|*x + ... + |a<SUB>d</SUB>|*x<SUP>d</SUP></I>.
     *  As with \c evaluate(), \c x can be of a type other than \c NT.
     *  \pre Requires \c NiX::NT_traits::Abs for NT.
     */
   
    template <class NTX> 
    typename Coercion_traits<NTX,NT>::Type 
    evaluate_absolute(const NTX& x) const {
        typedef typename Coercion_traits<NTX,NT>::Type Type;
        Type xx(x);
        CGAL_precondition( degree() >= 0 );
        typename Real_embeddable_traits<Type>::Abs abs;
        int d = degree();
        Type y(abs(Type(this->ptr()->coeff[d])));
        while (--d >= 0) y = y*xx + abs(Type(this->ptr()->coeff[d]));
        return y;
    } 

    /*! \brief evaluate the polynomial with all coefficients replaced by
     *  their absolute values
     *
     *  This is a specialization for \c x is of type NiX::Interval.
     */
    // TODO: Interval isn't available either!!
/*    Interval evaluate_absolute(const Interval& x) const {
        CGAL_precondition( degree() >= 0 );
        typename NT_traits<Interval>::Abs abs;
        typename NT_traits<NT>::To_Interval to_Interval;
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
//        BOOST_STATIC_ASSERT( (boost::is_same< typename Real_embeddable_traits<NT>::Is_real_embeddable,
//                              CGAL::Tag_true>::value) );
        typename Real_embeddable_traits<NT>::Sign sign;
        return sign(lcoeff());
    }

    //! return sign of difference
    CGAL::Comparison_result compare(const Polynomial& p2) const {
        typename Real_embeddable_traits<NT>::Compare compare;
        typename Real_embeddable_traits<NT>::Sign sign;
        CGAL_precondition(degree() >= 0);
        CGAL_precondition(p2.degree() >= 0);

        if (is_identical(p2)) return CGAL::EQUAL;

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
        typename Real_embeddable_traits<NT>::Sign sign;
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
        if (is_zero()) return NT(1);

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

    //! invoke \c NiX::NT_traits::Simplify on all coefficients
    void simplify_coefficients() { this->ptr()->simplify_coefficients(); }

    //! write polynomial to \c os in \c LiS::IO::PRETTY format
    /*! The output is intended to be Maple-readable; see module
     *  \link NiX_io NiX I/O Support \endlink.
     *
     * Example: A \c NiX::Polynomial<int> with a value of
     * 4<I>x</I><SUP>2</SUP> - 1 will be written as
     * <TT> 4*x^2 + (-1) </TT> by this function.
     */
    void output_maple(std::ostream& os) const;
    //! write polynomial to \c os in a format readable by \c input_ascii()
    void output_ascii(std::ostream& os) const;
     //! write polynomial to \c os in \c LiS::IO::BENCHMARK format
    void output_benchmark(std::ostream& os) const;

    //! implement \c NiX::Scalar_factor_traits::Scalar_div for polynomials
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
// Arithmetic Operations, Part I:
// declaring three-address arithmetic operators friends
//

    friend Polynomial<NT> operator - <> (const Polynomial<NT>&);   
    friend Polynomial<NT> operator + <> (const Polynomial<NT>&,
                                         const Polynomial<NT>&);
    friend Polynomial<NT> operator - <> (const Polynomial<NT>&,
                                         const Polynomial<NT>&);
    friend Polynomial<NT> operator * <> (const Polynomial<NT>&,
                                         const Polynomial<NT>&);
    friend Polynomial<NT> operator / <> (const Polynomial<NT>& p1,
                                         const Polynomial<NT>& p2);

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
     * Example: A \c NiX::Polynomial<int> with a value of
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
        int d = std::min(degree(),p1.degree()), i;
        for(i=0; i<=d; ++i) coeff(i) += p1[i];
        while (i<=p1.degree()) this->ptr()->coeff.push_back(p1[i++]);
        reduce(); return (*this);
    }

    Polynomial<NT>& operator -= (const Polynomial<NT>& p1) {
        this->copy_on_write();
        int d = std::min(degree(),p1.degree()), i;
        for(i=0; i<=d; ++i) coeff(i) -= p1[i];
        while (i<=p1.degree()) this->ptr()->coeff.push_back(-p1[i++]);
        reduce(); return (*this);
    }

    Polynomial<NT>& operator *= (const Polynomial<NT>& p1)
    { (*this) = (*this) * p1; return (*this); }

    Polynomial<NT>& operator /= (const Polynomial<NT>& p1)
    { (*this) = (*this) / p1; return (*this); }


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
    }


    // special operation to implement (pseudo-)division and the like
    void minus_offsetmult(const Polynomial<NT>& p, const NT& b, int k)
    {
        CGAL_precondition(!this->is_shared());
        int pd = p.degree();
        CGAL_precondition(degree() >= pd+k);
        for (int i = 0; i <= pd; i++) coeff(i+k) -= b*p[i];
        reduce();
    }
}; // class Polynomial<NT_>

//
// Traits classes for polynomials, making them number types themselves
//

// Algebraic structure traits

template< class POLY, class Algebraic_type >
class Polynomial_algebraic_structure_traits_base;

// The most basic suite of algebraic operations, suitable for the
// most basic kind of coefficient range, viz. a IntegralDomainWithoutDiv.
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, 
                                             Integral_domain_without_division_tag >
  : public Algebraic_structure_traits_base< POLY, 
                                            Integral_domain_without_division_tag > {
  public:
    typedef Integral_domain_without_division_tag Algebraic_category;
    
    class Simplify 
      : public Unary_function< POLY&, void > {
      public:
        void operator()( POLY& p ) const {
          p.simplify_coefficients();
        }        
    };
    
    class Unit_part 
      : public Unary_function< POLY, POLY > {
      public:
        POLY operator()( const POLY& x ) const {
          return POLY( x.unit_part() );
        }
    };
};

// Extend to the case that the coefficient range is a IntegralDomain (with div)
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Integral_domain_tag >
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Integral_domain_without_division_tag > {
  public:
    typedef Integral_domain_tag Algebraic_category;
    
    class Integral_division 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()( const POLY& x, const POLY& y ) const {
          return x / y;
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
        
    };
};

// Extend to a UFDomain as coefficient range
// Forward declaration for <NiX/polynomial_gcd.h> for NT_traits<Poly...>::Gcd
namespace INTERN_POLYNOMIAL_GCD {
template <class NT> inline
Polynomial<NT> gcd(const Polynomial<NT>&, const Polynomial<NT>&);
} // namespace INTERN_POLYNOMIAL_GCD

template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Unique_factorization_domain_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Integral_domain_tag > {
  public:
    typedef Unique_factorization_domain_tag Algebraic_category;
    
    class Gcd 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()( const POLY& x, const POLY& y ) const {
          // First: the extreme cases and negative sign corrections.          
          if (x == POLY(0)) {
              if (y == POLY(0))  
                  return POLY(0);
              return y.abs(); // TODO: Is this correct?
          }
          if (y == POLY(0))
              return x.abs(); // TODO: Is this correct?
          return INTERN_POLYNOMIAL_GCD::gcd(x,y);
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
        
    };
};

// Clone this for a EuclideanRing
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Euclidean_ring_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Unique_factorization_domain_tag > {
  // nothing new
};

// Extend to a field as coefficient range
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Field_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                            Unique_factorization_domain_tag > {
  public:
    typedef Euclidean_ring_tag Algebraic_category;
    
    class Div_mod {
      public:
        typedef POLY first_argument_type;
        typedef POLY second_argument_type;
        typedef POLY& third_argument_type;
        typedef POLY& fourth_argument_type;
        typedef void result_type;
        typedef Arity_tag< 4 >  Arity;
        
        void operator()( const POLY& x, const POLY& y, 
                         POLY& q, POLY& r ) const {
          POLY::euclidean_division( x, y, q, r );
        }
        
        template < class NT1, class NT2 >
        void operator()( const NT1& x, const NT2& y,
                         POLY& q, POLY& r ) const {
          BOOST_STATIC_ASSERT((::boost::is_same<
                  typename Coercion_traits< NT1, NT2 >::Type, POLY
                                               >::value));
          
          typename Coercion_traits< NT1, NT2 >::Cast cast;
          operator()( cast(x), cast(y), q, r );          
        }
        
    };
    
    class Div 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator()(const POLY& a, const POLY& b) const {
          POLY q, r;
          Div_mod()(a, b, q, r);
          return q;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
    };
    
    class Mod 
      : public Binary_function< POLY, POLY, POLY > {
      public:
        POLY operator () (const POLY& a, const POLY& b) const {
          POLY q, r;
          Div_mod()(a, b, q, r);
          return r;
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( POLY )
    };
    
};

// Clone this for a FieldWithSqrt
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, Field_with_sqrt_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, Field_tag > {
  // nothing new
};

// Clone this for a FieldWithKthRoot
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, 
                                                  Field_with_kth_root_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                                       Field_with_sqrt_tag > {
  // nothing new
};

// Clone this for a FieldWithRootOf
template< class POLY >
class Polynomial_algebraic_structure_traits_base< POLY, 
                                                  Field_with_root_of_tag > 
  : public Polynomial_algebraic_structure_traits_base< POLY, 
                                                    Field_with_kth_root_tag > {
  // nothing new
};

// The actual algebraic structure traits class
template< class NT > class Algebraic_structure_traits< Polynomial< NT > >
  : public Polynomial_algebraic_structure_traits_base< Polynomial< NT >,
         typename Algebraic_structure_traits< NT >::Algebraic_category > {
  public:
    typedef Polynomial<NT> Algebraic_structure;
    typedef typename Algebraic_structure_traits< NT >::Is_exact Is_exact;
};

// Real embeddable traits
// TODO: Polynomials aren't Real_embeddable! But for debugging and testing
//       reasons, the real embeddable functors are provided.
template< class NT > class Real_embeddable_traits< Polynomial<NT> > 
  : public Real_embeddable_traits_base< Polynomial<NT> > {
  public:
      
//    typedef typename Real_embeddable_traits<NT>::Is_real_embeddable 
//                                                            Is_real_embeddable;
  
    typedef Tag_false Is_real_embeddable;
    
    class Abs {
      public:
        typedef Polynomial<NT> argument_type;
        typedef Polynomial<NT> result_type;
        Polynomial<NT> operator()( const Polynomial<NT>& x ) const {
          return x.abs(); 
        }
    };

    class Sign {
      public:
        typedef Polynomial<NT>              argument_type;
        typedef CGAL::Sign        result_type;
        CGAL::Sign operator()( const Polynomial<NT>& x ) const {
          return x.sign();
        }
    };
    
    class Compare {
      public:
        typedef Polynomial<NT>                    first_argument_type;
        typedef Polynomial<NT>                    second_argument_type;
        typedef CGAL::Comparison_result result_type;
        CGAL::Comparison_result operator()( const Polynomial<NT>& x, 
                                            const Polynomial<NT>& y ) const {
          return x.compare(y);
        }
        
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Polynomial<NT>,
                                                      CGAL::Comparison_result );
        
    };

    class To_double {
      public:
        typedef typename Real_embeddable_traits<NT>::To_double NT_to_double;
        typedef Polynomial<typename NT_to_double::result_type> result_type;
        typedef Polynomial<NT> argument_type;
        result_type operator()( const Polynomial<NT>& x ) const {
            CGAL_precondition(x.degree() >= 0);
            NT_to_double to_double;
            return RET(::boost::make_transform_iterator(x.begin(),to_double),
                       ::boost::make_transform_iterator(x.end()  ,to_double));
        }
    };

    class To_interval {
      public:
        typedef typename Real_embeddable_traits<NT>::To_interval NT_to_interval;
        typedef Polynomial<typename NT_to_interval::result_type> result_type;
        typedef Polynomial<NT> argument_type;
        result_type operator()( const Polynomial<NT>& x ) const {
            CGAL_precondition( x.degree() >= 0 );
            NT_to_interval to_interval;  
            return RET(::boost::make_transform_iterator(x.begin(),to_interval),
                       ::boost::make_transform_iterator(x.end()  ,to_interval));
        }
    };

};

// Now we wrap up all of this in the actual NT_traits
// specialization for Polynomial<NT>
/*! \ingroup NiX_Polynomial
    \brief \c NiX::NT_traits < \c NiX::Polynomial<NT> >
 *
 *  If \c NT is a model of a number type concept, then so is
 *  \c Polynomial<NT>. A specialization of \c NiX::NT_traits
 *  is provided automatically by NiX/Polynomial.h.
 *
 *  The number type concepts for the coefficient domain NT are
 *  mapped to those for the polynomials as follows:
 *  <PRE>
    IntegralDomainWithoutDiv --> IntegralDomainWithoutDiv
    IntegralDomain           --> IntegralDomain
    UFDomain       --> UFDomain
    EuclideanRing  --> UFDomain
    Field          --> EuclideanRing
    FieldWithSqrt  --> EuclideanRing
    </PRE>
 *
 *  \c Polynomial<NT> is \c RealComparable iff \c NT is.
 *  The ordering is determined by the \c sign() of differences, see ibid.
 *
 *  \c Polynomial<NT> offers a non-<TT>Null_tag</TT> \c To_double
 *  iff \c NT does. If non-<TT>Null_tag</TT>, it returns a coefficient-wise
 *  \c double approximation of the polynomial.
 */



// We need to play a similar game to provide Fraction_traits


template <class POLY, class TAG>
class Poly_Ftr_base;

template <class POLY>
POLY fractionalize_polynomial(
    const typename Fraction_traits<POLY>::Numerator_type& p,
    const typename Fraction_traits<POLY>::Denominator_type& c
);
template <class POLY>
typename Fraction_traits<POLY>::Numerator_type
integralize_polynomial(
    const POLY& p,
    typename Fraction_traits<POLY>::Denominator_type& c
);

// Use this if the coefficients cannot be decomposed
// into numerator and denominator
template <class NT_>
class Poly_Ftr_base< Polynomial<NT_>, CGAL::Tag_false > {
public:
    typedef Polynomial<NT_> Type;
    typedef CGAL::Tag_false Is_fraction;
    typedef CGAL::Null_tag Numerator;
    typedef CGAL::Null_tag Denominator_type;
    typedef CGAL::Null_functor Common_factor;
    typedef CGAL::Null_functor Decompose;
    typedef CGAL::Null_functor Compose;
};

// If they can, use this
template <class NT_>
class Poly_Ftr_base< Polynomial<NT_>, CGAL::Tag_true > {
public:
    typedef Polynomial<NT_> Type;
    typedef CGAL::Tag_true Is_fraction;
    typedef Polynomial<typename Fraction_traits<NT_>::Numerator_type>
        Numerator_type;
    typedef typename Fraction_traits<NT_>::Denominator_type Denominator_type;
    typedef typename Fraction_traits<NT_>::Common_factor Common_factor;
    class Decompose {
    public:
        typedef Type first_argument_type;
        typedef Numerator_type& second_argument_type;
        typedef Denominator_type& third_argument_type;
        inline void operator () (
                const Type& p,
                Numerator_type& n,
                Denominator_type& d){
            n = integralize_polynomial<Type>(p, d);
        }
    };
    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef Type result_type;
        inline Type operator () (const Numerator_type& n,
                                 const Denominator_type& d){
            return fractionalize_polynomial<Type>(n, d);
        };
    };
};


// Select the right alternative as Fraction_traits
/*! \ingroup NiX_Polynomial
    \brief \c NiX::Fraction_traits < \c NiX::Polynomial<NT> >
 *
 *  Polynomials provide suitable specializations of \c NiX::Fraction_traits.
 *  They are decomposable iff their coefficient type is.
 *  The denominator \e d of a polynomial \e p is a low common multiple
 *  (see \c NiX::Fraction_traits::Common_factor for details) of the
 *  denominators of its coefficients.  The numerator is the polynomial
 *  \e d*p with a fraction-free coefficient type.
 *
 *  This works for nested polynomials, too.
 */
template <class NT_>
class Fraction_traits< Polynomial<NT_> >
    : public Poly_Ftr_base< Polynomial<NT_>,
                 typename Fraction_traits<NT_>::Is_fraction >
{
    // nothing new
};

/*! \ingroup NiX_Polynomial
    \relates NiX::Polynomial
 *  \brief implement \c NiX::Fraction_traits::Decompose
 */
template <class POLY>
typename Fraction_traits<POLY>::Numerator_type
integralize_polynomial(
    const POLY& p,
    typename Fraction_traits<POLY>::Denominator_type& c
) {

    typedef Fraction_traits<POLY> PFTRAITS;
    typedef Fraction_traits<typename POLY::NT> CFTRAITS;
    typedef typename PFTRAITS::Numerator_type INTPOLY;
    typedef typename CFTRAITS::Numerator_type INTCOEFF;
    typedef typename PFTRAITS::Denominator_type DENOM;

    const int d = p.degree();
    std::vector<INTCOEFF> integ(d+1);
    std::vector<DENOM> denom(d+1);
  
    int i;

    // decompose each coefficient into integral part and denominator
    typename CFTRAITS::Decompose decomp_coeff;
    for (i = 0; i <= d; i++) {
        decomp_coeff(p[i], integ[i], denom[i]);
    }

    // c = lcm(denom[0], ..., denom[d])
    typename Algebraic_structure_traits<DENOM>::Integral_division idiv;
    typename CFTRAITS::Common_factor        gcd;  // not really `greatest'
    c = denom[0];
    for (i = 1; i <= d; i++) {
        c *= idiv(denom[i], gcd(c, denom[i]));
    }

    // expand each (integ, denom) pair to common denominator
    for (i = 0; i <= d; i++) {
        integ[i] *= INTCOEFF(idiv(c, denom[i]));
    }
    return INTPOLY(integ.begin(), integ.end());
    
    return INTPOLY(0);
}

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief implement \c NiX::Fraction_traits::Compose
 */
template <class POLY>
POLY fractionalize_polynomial(
    const typename Fraction_traits<POLY>::Numerator_type& p,
    const typename Fraction_traits<POLY>::Denominator_type& c
) {

    typename Fraction_traits<typename POLY::NT>::Compose comp_coeff; 
    (void)comp_coeff;
    
    std::vector<typename POLY::NT> coeffs(p.degree()+1);
    
    for (int i = 0; i <= p.degree(); i++) {
        coeffs[i] = comp_coeff(p[i], c);
    }
    
    return POLY(coeffs.begin(), coeffs.end());

}

// Internally, Polynomials also need this:

//! used internally for data exchanged between nesting levels of polynomials
// this traits-class provides
//   a) poly_nesting_depth: a counter for identifying a polynomials nesting level
//   b) Innermost_coefficient: a type which defines the numbertype of the
//        innermost nesting level, which is not a polynomial itself
//   c) Innermost_lcoeff: returns the leading coefficient of the polynomial in its
//       common sense
//   d) Innermost_coefficient_to_polynomial: transforms an innermost
//       numbertype into a polynomial like a typecast

template< class NT >
class Polynomial_traits {
  public: 
    typedef NT Polynomial;

  public:
    static const int poly_nesting_depth = 0;
    typedef NT Innermost_coefficient;

    // TODO: Without this I currently get compile errors in a part which is never
    //       called for Polynomial_traits< NT > where NT is not a polynomial.
/*    class Get_monomial_representation {
      public:
        template< class Back_insert_iterator >
        void operator()( const Polynomial&, Back_insert_iterator ) {
          CGAL_error( "This should never be called" );
        }
      
    };*/

    
    class Innermost_leading_coefficient 
      : public Unary_function< NT, Innermost_coefficient > {
      public:        
        //! return the innermost leading coefficient
        Innermost_coefficient operator()(const NT& a) { return a; }
    };
      
    class Innermost_coefficient_to_polynomial 
      : public Unary_function< NT, NT > {
      public:        
        //! typecast from an innermost coefficient to NT (no change here)
        NT operator()(const NT& a) { return a; }
    };
  
  public:
    class Total_degree 
      : public Unary_function< Polynomial, int > {
    public:
        //! returns the total degree of the polynomial
        int operator()(const Polynomial& constant){
            return 0;
        }
    };
    
    typedef CGAL::Null_functor Innermost_coefficient_iterator;
    typedef CGAL::Null_functor Innermost_coefficient_begin;
    typedef CGAL::Null_functor Innermost_coefficient_end;
};

template< class NT >
class Polynomial_traits< Polynomial<NT> > {
  public:
    typedef Polynomial<NT> Polynomial;
    typedef typename Polynomial_traits<NT>::Innermost_coefficient
                                                         Innermost_coefficient;
  public:
    static const int poly_nesting_depth =
        Polynomial_traits<NT>::poly_nesting_depth + 1;

    class Construct_polynomial {
      public:
        typedef Polynomial  result_type;
        typedef std::pair< Exponent_vector, Innermost_coefficient > Exponents_coeff_pair;
        typedef std::vector< Exponents_coeff_pair > Monomial_rep; 
        typedef typename Monomial_rep::iterator Monomial_rep_iterator;
        
        Polynomial operator()() {
          return Polynomial(0);
        }
        
        Polynomial operator()( int i ) {
          return Polynomial(i);
        }
        
        template< class Input_iterator >
        Polynomial operator()( Input_iterator begin, Input_iterator end ) {
          return Polynomial( begin, end );
        }
        
        Polynomial operator()( Monomial_rep_iterator begin,
                               Monomial_rep_iterator end ) {
          Monomial_rep_iterator it = begin;
          return Create_polynomial_from_monomial_rep< NT >()( it, end, 0 ); 
        }
        
      // TODO: This shouldn't be public, but I need public access for the recursion
        template< class Coeff_type >
        class Create_polynomial_from_monomial_rep {
          public:
            Polynomial operator()( Monomial_rep_iterator& it,
                                   Monomial_rep_iterator end,
                                   int exponent_pos ) {
              std::vector< NT > coefficients;
              int previous_exponent = (exponent_pos>0)? it->first[exponent_pos-1] :
                                                        -1;
              while( it != end && ((exponent_pos>0)? it->first[exponent_pos-1] : -1) ==
                     previous_exponent ) {
              
                while( it->first[exponent_pos] > (int)coefficients.size() )
                  coefficients.push_back( 0 );
              
                coefficients.push_back( it->second );  
                
                ++it;
              }
              return Polynomial( coefficients.begin(), coefficients.end() );            
            }            
        };
        
        template< class Coeff_type >
        class Create_polynomial_from_monomial_rep< CGAL::Polynomial< Coeff_type > > {
          public:
            Polynomial operator()( Monomial_rep_iterator& it,
                                   Monomial_rep_iterator end,
                                   int exponent_pos ) {                                                  
              std::vector< NT > coefficients;
              int previous_exponent = (exponent_pos>0)? it->first[exponent_pos-1] :
                                                      -1;
              while( it != end && ((exponent_pos>0)? it->first[exponent_pos-1] : -1) ==
                     previous_exponent ) {
              
                while( it->first[exponent_pos] > (int)coefficients.size() )
                  coefficients.push_back( 0 );
              
                typename Polynomial_traits< NT >::Construct_polynomial::template Create_polynomial_from_monomial_rep<Coeff_type> construct;
                coefficients.push_back( construct( it, end, exponent_pos+1 ) );
                
                // The iterator already gets increased in the construct-call.                    
              }

              return Polynomial( coefficients.begin(), coefficients.end() );
           }
        };
    };
    
    
    class Get_monomial_representation {      
      public:
        typedef std::pair< Exponent_vector, Innermost_coefficient > Exponents_coeff_pair;
        typedef std::vector< Exponents_coeff_pair > Monomial_rep; 
        typedef std::back_insert_iterator< Monomial_rep > Back_insert_iterator;

                  
        void operator()( const Polynomial& p, Back_insert_iterator back_insert_iterator ) {
          Create_monomial_representation< NT >()( p, back_insert_iterator );
        }
      
      private:
        template< class Coeff_type >
        class Create_monomial_representation {
          public:
            void operator()( const Polynomial& p, 
                             Back_insert_iterator back_insert_iterator ) {             
              for( int exponent = 0; exponent <= p.degree(); ++exponent ) {
                Exponent_vector exp_vec;
                exp_vec.push_front( exponent );
                *back_insert_iterator = Exponents_coeff_pair( exp_vec, p[exponent] );              
              } 
            }
        };
        
        template< class Coeff_type >
        class Create_monomial_representation< CGAL::Polynomial<Coeff_type> > {
          public:
            void operator()( const Polynomial& p, 
                             Back_insert_iterator back_insert_iterator ) {             
              for( int exponent = 0; exponent <= p.degree(); ++exponent ) {
                Monomial_rep monomial_rep;
                typename Polynomial_traits< NT >::Get_monomial_representation()( p[exponent], std::back_inserter( monomial_rep ) );
                for( typename Monomial_rep::iterator it = monomial_rep.begin();
                     it != monomial_rep.end(); ++it ) {
                  it->first.push_front( exponent );
                }
                copy( monomial_rep.begin(), monomial_rep.end(), back_insert_iterator );               
              }
            }
        };                             
    };
       
    /*class Degree 
      : public Binary_function< Polynomial, int, int > {
      public:
        int operator()( const Polynomial& p ) {
          return 0;
        }
        
        int operator()( const Polynomial& p, int i ) {
          
        }
    };*/

    class Swap {
      public:
        typedef Polynomial      result_type;
        typedef Polynomial      first_argument_type;
        typedef int             second_argument_type;
        typedef int             third_argument_type;
        typedef Arity_tag< 3 >         Arity;
      
        Polynomial operator()( const Polynomial& p, int x, int y ) {
          typedef std::pair< Exponent_vector, Innermost_coefficient > Exponents_coeff_pair;
          typedef std::vector< Exponents_coeff_pair > Monomial_rep; 
          Get_monomial_representation gmr;
          Construct_polynomial construct;
          Monomial_rep mon_rep;
          gmr( p, std::back_inserter( mon_rep ) );
          for( typename Monomial_rep::iterator it = mon_rep.begin(); it != mon_rep.end();
               ++it ) {
            it->first.swap( x, y );
          }
          
          std::sort( mon_rep.begin(), mon_rep.end() );
          
          return construct( mon_rep.begin(), mon_rep.end() );
        }
        
    };


    
// TODO: Made this public again for backward compatibility
//       Needs to be discussed  
//  private:
public:
    class Innermost_leading_coefficient 
      : public Unary_function< Polynomial, Innermost_coefficient > {
      public:        
        //! return the innermost leading coefficient
        Innermost_coefficient operator()( const Polynomial& p) {
          typename Polynomial_traits<NT>::Innermost_leading_coefficient ilc;
          return ilc(p.lcoeff());
        }
    };

    class Innermost_coefficient_to_polynomial 
      : public Unary_function< Innermost_coefficient, Polynomial > {
      public:
        //! typecast from an innermost coefficient to polynomial
        Polynomial operator()(const Innermost_coefficient& a) {
            typename Polynomial_traits<NT>
                ::Innermost_coefficient_to_polynomial ictp;
            return CGAL::Polynomial<NT>(ictp(a));
        }
    };

    typedef CGAL::Recursive_const_flattening<
        poly_nesting_depth-1, typename CGAL::Polynomial<NT>::const_iterator
    > Coefficient_flattening;

    typedef typename Coefficient_flattening::Recursive_flattening_iterator
        Innermost_coefficient_iterator;

public:
    class Innermost_coefficient_begin 
      : public Unary_function< Polynomial, Innermost_coefficient_iterator > {
      public:
        Innermost_coefficient_iterator
        operator () (const Polynomial& p) {
            return typename Coefficient_flattening::Flatten()(p.end(),p.begin());
        }
    };

public:
    class Innermost_coefficient_end 
      : public Unary_function< Polynomial, Innermost_coefficient_iterator > {
      public:
        Innermost_coefficient_iterator
        operator () (const Polynomial& p) {
            return typename Coefficient_flattening::Flatten()(p.end(),p.end());
        }
    };

public:      
    class Total_degree 
      : public Unary_function< Polynomial, int > {
      public:
        //! returns the total degree of the polynomial
        int operator()(const Polynomial& polynomial){
            
            typedef Polynomial_traits<NT> NT_POLY_TRAITS;
            typename NT_POLY_TRAITS::Total_degree total_degree;
            
            CGAL_precondition(polynomial.degree()>=0);
            
            int result = 0;
            for(int i = 0; i <= polynomial.degree(); i++){
                 if(polynomial[i] != NT(0))
                     result = std::max(result , total_degree(polynomial[i]) + i );
            } 
            return result;
        }
    };
};
 
//! returns the total degree of the polynomial
template <class Polynomial> 
int total_degree(const Polynomial& polynomial){
    typedef Polynomial_traits<Polynomial> PT;
    typename PT::Total_degree total_degree;
    return total_degree(polynomial);
}


// see <NiX/Scalar_factor_traits.h>
/*! \ingroup NiX_Polynomial
 *  \ingroup NiX_Scalar_factor_traits_spec
    \brief \c NiX::Scalar_factor_traits < \c NiX::Polynomial<NT> >
 *
 *  The \c NiX::Scalar_factor_traits::Scalar_factor of a polynomial,
 *  even of a nested polynomial, is essentially the gcd of its scalar
 *  coefficients. (This differs from the \c content() in that the content
 *  is the gcd of the coefficient sequence, which is again a polynomial,
 *  albeit with one variable less, for nested polynomials.) The scalar
 *  factor of the zero polynomial is 0.
 *
 *  This currently only works if the scalar type is a model of \c UFDomain.
 */
template <class Coeff>
class Scalar_factor_traits< Polynomial<Coeff> > {
public:
    typedef Polynomial<Coeff> NT;
    typedef typename Scalar_factor_traits<Coeff>::Scalar Scalar;
    class Scalar_factor {
    public:
        //! argument type
        typedef NT argument_type;
        //! first argument type
        typedef NT first_argument_type;
        //! second argument type
        typedef Scalar second_argument_type;
        //! result type
        typedef Scalar result_type;
        //! returns the gcd of a polynomials coefficients
        Scalar operator()(const NT& p, const Scalar& d_=Scalar(0)) const {
            
            typename Scalar_factor_traits<Coeff>::Scalar_factor sfac;
            const Scalar unity(1);
            
            Scalar d(d_);
            if (p.is_zero()) return d;
            
            int i = p.degree();
            while((d != unity) && (i >= 0)) {
                d = sfac(p[i--],d);
            }
            return d;
        }
    };
    class Scalar_div {
    public:
        //! first argument type
        typedef NT first_argument_type;
         //! second argument type
        typedef Scalar second_argument_type;
        //! divides a polynomial \c p by a scalar factor \c b
        void operator () (NT& p, const Scalar& b) const { 
            CGAL_precondition(b != Scalar(0));
            p.scalar_div(b); 
        }
    };
};

/*! \ingroup NiX_Polynomial
 *  \ingroup NiX_Modular_traits_spec
 *  \brief Specialization of Modular_traits for NiX::Polynomial.
 * 
 *  NiX::Modular_traits::Modular_image maps the coefficients of a polynomial
 *  to their Modular_image and returns the resulting polynomial.  
 */
/*template< class COEFF >
class Modular_traits< Polynomial<COEFF> > {
    
private:
    typedef Modular_traits<COEFF> Mtr;
public:
    typedef Polynomial<COEFF> NT;
    typedef Modular_traits<NT> Self;
    typedef typename Mtr::Is_convertible Is_convertible;
    typedef Polynomial<typename Mtr::Modular_NT> Modular_NT;
    
    struct Modular_image{
        Modular_NT operator()(const NT& p){ 
            std::vector<typename Mtr::Modular_NT> V;
            typename Mtr::Modular_image modular_image;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image(p[i]));
            return Modular_NT(V.begin(),V.end());           
        }
    };
};*/

namespace INTERN_POLYNOMIAL {

template <class NT>
Polynomial<NT> canonicalize_polynomial_(Polynomial<NT> p, CGAL::Tag_true)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits<POLY>::Innermost_coefficient_to_polynomial ictp;
    typename Algebraic_extension_traits<IC>::Normalization_factor nfac;
  
    IC tmp = nfac(ilcoeff(p));
    if(tmp != IC(1)){
        p *= ictp(tmp);
    }
    remove_scalar_factor(p);
    p /= p.unit_part();
    p.simplify_coefficients();

    CGAL_postcondition(nfac(ilcoeff(p)) == IC(1));
    return p;
};

template <class NT>
Polynomial<NT> canonicalize_polynomial_(Polynomial<NT> p, CGAL::Tag_false)
{  
    remove_scalar_factor(p);
    p /= p.unit_part();
    p.simplify_coefficients();
    return p;
};

} // namespace INTERN_POLYNOMIAL

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  
 *  \brief divide a polynomial \c p by its Scalar_factor and Unit_part
 *  
 *  ...making it a canonical representative of all its constant multiples.
 *  Depending on the number type of the innermost coefficient, this
 *  function does
 *    a) dividing \c p by the leading coefficient in fields
 *    b) dividing \c p by the gcd of all coefficients in UFDomains
 *    c) extending the leading coefficient in Sqrt_extensions, so it
 *        becomes integral, and dividing \c p by the gcd of all scalars
 *        \see NiX/Sqrt_extension.h
 *  The result is uniquely determined by setting the leading coefficient
 *  to the minimal integral rational.
 */
template <class NT> inline
  Polynomial<NT> canonicalize_polynomial(const Polynomial<NT>& p)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;

    if (p.is_zero()) return p;
    return canonicalize_polynomial_(p, Is_extended());
};

namespace INTERN_POLYNOMIAL {

// Polynomial<NT> / Polynomial<NT>  -  coefficient type is extended
template <class NT>
Polynomial<NT> div_utcf_(
    Polynomial<NT> f, const Polynomial<NT>& g, bool, CGAL::Tag_true)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits<POLY>::Innermost_coefficient_to_polynomial ictp;
    typename Polynomial_traits<POLY>::Innermost_coefficient_begin begin;
    typename Polynomial_traits<POLY>::Innermost_coefficient_end end;
    typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;

    IC tmp = ilcoeff(g);
    tmp *= dfai(begin(g), end(g));
    f *= ictp(tmp);
    return canonicalize_polynomial(f / g);
}

// Polynomial<NT> / Polynomial<NT>  -  coefficient type is NOT extended
template <class NT>
Polynomial<NT> div_utcf_(
    Polynomial<NT> f, const Polynomial<NT>& g, bool is_canonicalized, CGAL::Tag_false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits<POLY>::Innermost_coefficient_to_polynomial ictp;

    if (!is_canonicalized) {
        IC lcoeff = ilcoeff(g);
        f *= ictp(lcoeff);
    }
    return canonicalize_polynomial(f / g);
}

// Polynomial<NT> / NT  -  NT is already the coefficient type and is extended
template <class NT>
Polynomial<NT> div_utcf_NT_is_IC(
    Polynomial<NT> f, const NT& g, CGAL::Tag_false)
{
    return canonicalize_polynomial(f);
}

// Polynomial<NT> / NT  -  NT is again a polynomial  -   coefficient type is extended
template <class NT, class Is_nested>
Polynomial<NT> div_utcf_NT_is_IC(
    Polynomial<NT> f, const NT& g, Is_nested)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits<POLY>::Innermost_coefficient_to_polynomial ictp;
    typename Polynomial_traits<NT>::Innermost_coefficient_begin begin;
    typename Polynomial_traits<NT>::Innermost_coefficient_end end;
    typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;

    IC tmp = ilcoeff(g);
    tmp *= dfai(begin(g), end(g));
    f *= ictp(tmp);
    return canonicalize_polynomial(f / g);
}

// Polynomial<NT> / NT  -  coefficient type is extended
template <class NT> inline
Polynomial<NT> div_utcf_(
    const Polynomial<NT>& f, const NT& g, bool, CGAL::Tag_true)
{
    typedef CGAL::Boolean_tag< (Polynomial_traits<NT>::poly_nesting_depth >= 2) > Is_nested;
    return div_utcf_NT_is_IC(f, g, Is_nested() );
}

// Polynomial<NT> / NT  -  coefficient type is NOT extended
template <class NT>
Polynomial<NT> div_utcf_(
    Polynomial<NT> f, const NT& g, bool is_canonicalized, CGAL::Tag_false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits<POLY>::Innermost_coefficient_to_polynomial ictp;

    if (!is_canonicalized) {
        IC lcoeff = ilcoeff(g);
        f *= ictp(lcoeff);
    }
    return canonicalize_polynomial(f / g);
}

} // namespace INTERN_POLYNOMIAL

//! divide \c f by \c g with respect to constant factors
/*! This function provides a division of two polynomials, which takes
 *  no care of constant factors of the innermost scalar type.
 *  The boolean parameter decides whether the divisor has already been
 *  canonicalized due to running time optimisation.
 *  The result is made unique by canonicalizing it.
 */

template <class NT> inline
Polynomial<NT> div_utcf(
    const Polynomial<NT>& f, const Polynomial<NT>& g, bool is_canonicalized = false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;

    return div_utcf_(f, g, is_canonicalized, Is_extended());
}

//! overloaded version for divisors with a by one lower nesting level
template <class NT> inline
Polynomial<NT> div_utcf(
    const Polynomial<NT>& f, const NT& g, bool is_canonicalized = false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits<POLY>::Innermost_coefficient IC;
    typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;

    return div_utcf_(f, g, is_canonicalized, Is_extended());
}

//
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
Polynomial<NT> operator + (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{ 
    typedef typename Polynomial<NT>::size_type size_type;
    CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
    bool p1d_smaller_p2d = p1.degree() < p2.degree();
    int min,max,i;
    if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
    else                 { max = p1.degree(); min = p2.degree(); }
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT> p(TAG, size_type(max + 1));
    for (i = 0; i <= min; ++i ) p.coeff(i) = p1[i]+p2[i];
    if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)=p2[i];
    else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=p1[i];
    p.reduce();
    return p;
}

template <class NT> inline
Polynomial<NT> operator - (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{ 
    typedef typename Polynomial<NT>::size_type size_type;
    CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
    bool p1d_smaller_p2d = p1.degree() < p2.degree();
    int min,max,i;
    if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
    else                 { max = p1.degree(); min = p2.degree(); }
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT> p(TAG, size_type(max + 1));
    for (i = 0; i <= min; ++i ) p.coeff(i)=p1[i]-p2[i];
    if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)= -p2[i];
    else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=  p1[i];
    p.reduce();
    return p;
}

template <class NT> inline
Polynomial<NT> operator * (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{
    typedef typename Polynomial<NT>::size_type size_type;
    CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT>  p(TAG, size_type(p1.degree()+p2.degree()+1) ); 
    // initialized with zeros
    for (int i=0; i <= p1.degree(); ++i)
        for (int j=0; j <= p2.degree(); ++j)
            p.coeff(i+j) += (p1[i]*p2[j]); 
    p.reduce();
    return p;
}

#ifndef NiX_POLY_USE_SLOW_DIVISION

template <class NT> inline
Polynomial<NT> operator / (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{
    typedef Algebraic_structure_traits< Polynomial<NT> > AST; 
    // Precondition: q with p1 == p2 * q must exist within NT[x].
    // If this holds, we can perform Euclidean division even over a ring NT
    // Proof: The quotients of each division that occurs are precisely
    //   the terms of q and hence in NT.
    CGAL_precondition(!p2.is_zero());
    if (p1.is_zero()) return p1;
    Polynomial<NT> q, r;
    Polynomial<NT>::euclidean_division(p1, p2, q, r);
    // TODO: Replace by correct Makro
    CGAL_postcondition( !AST::Is_exact::value || p2 * q == p1);
    
    return q;
}

#else

template <class NT> inline
Polynomial<NT> operator / (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{
    typedef typename Algebraic_structure_traits<NT>::Algebric_structure_tag Algebra_type;
    return division(p1, p2, Algebra_type());
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

#endif // NiX_POLY_USE_SLOW_DIVISION

//
// Arithmetic Operators, Part IV:
// Mixed-mode three-address arithmetic on top of two-address operators
//

// lefthand side
template <class NT> inline
Polynomial<NT> operator + (const NT& num, Polynomial<NT> p2)
{ p2 += num; return p2; }
template <class NT> inline
Polynomial<NT> operator - (const NT& num, Polynomial<NT> p2)
{ p2 -= num; return -p2; }
template <class NT> inline
Polynomial<NT> operator * (const NT& num, Polynomial<NT> p2)
{ p2 *= num; return p2; }
template <class NT> inline
Polynomial<NT> operator / (const NT& num, const Polynomial<NT>& p2) {
    CGAL_precondition(p2.degree() == 0);
    CGAL_precondition(p2[0] != NT(0));
    typename Algebraic_structure_traits<NT>::Integral_division idiv;
    return Polynomial<NT>(idiv(num, p2[0]));
}

// righthand side
template <class NT> inline
Polynomial<NT> operator + (Polynomial<NT> p1, const NT& num)
{ p1 += num; return p1; }
template <class NT> inline
Polynomial<NT> operator - (Polynomial<NT> p1, const NT& num)
{ p1 -= num; return p1; }
template <class NT> inline
Polynomial<NT> operator * (Polynomial<NT> p1, const NT& num)
{ p1 *= num; return p1; }
template <class NT> inline
Polynomial<NT> operator / (Polynomial<NT> p1, const NT& num)
{ p1 /= num; return p1; }


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
bool operator != (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return !(p1 == p2); }    
template <class NT> inline
bool operator < (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) < 0 ); }    
template <class NT> inline
bool operator <= (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) <= 0 ); }    
template <class NT> inline
bool operator > (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) > 0 ); }    
template <class NT> inline
bool operator >= (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) >= 0 ); }    

// mixed-mode, lefthand side
template <class NT> inline
bool operator == (const NT& num, const Polynomial<NT>& p)  {
    CGAL_precondition(p.degree() >= 0);
    return p.degree() == 0 && p[0] == num;
}
template <class NT> inline
bool operator != (const NT& num, const Polynomial<NT>& p) 
{ return !(num == p);}
template <class NT> inline
bool operator < (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) > 0 );}
template <class NT> inline
bool operator <=  (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) >= 0 );}
template <class NT> inline
bool operator > (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) < 0 );}
template <class NT> inline
bool operator >= (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) <= 0 );}

// mixed-mode, righthand side
template <class NT> inline
bool operator == (const Polynomial<NT>& p, const NT& num) 
{ return num == p; }
template <class NT> inline
bool operator != (const Polynomial<NT>& p, const NT& num) 
{ return !(num == p); }
template <class NT> inline
bool operator < (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) < 0 );}
template <class NT> inline
bool operator <= (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) <= 0 );}
template <class NT> inline
bool operator > (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) > 0 );}
template <class NT> inline
bool operator >= (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) >= 0 );}


//
// I/O Operations
//

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
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

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
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

template <class NT>
void Polynomial<NT>::output_maple(std::ostream& os) const {
    const Polynomial<NT>& p = *this;
    const char *varname;
    char vnbuf[42];
    
    // use variable names x, y, z, w1, w2, w3, ...
    if (Polynomial_traits<NT>::poly_nesting_depth < 3) {
        static const char *varnames[] = { "x", "y", "z" };
        varname = varnames[Polynomial_traits<NT>::poly_nesting_depth];
    } else {
        sprintf(vnbuf, "w%d", Polynomial_traits<NT>::poly_nesting_depth - 2);
        varname = vnbuf;
    }
    
    int i = p.degree();
    print_maple_monomial(os, p[i], varname, i);
    while (--i >= 0) {
        if (p[i] != NT(0)) {
            os << " + ";
            print_maple_monomial(os, p[i], varname, i);
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
    const Polynomial<NT> &p = *this;
    if (p.is_zero()) { 
        os << "Polynomial_1(0)"; 
        return; 
    }
    os << "Polynomial_1(";
    for (int i = 0; i <= p.degree(); i++) {
        os << CGAL::oformat(p[i]);
        if (i != p.degree()) {
            os << ",";
        }
    }
    os << ")";
}

// Moved to internal namespace because of name clashes
// TODO: Is this OK?
namespace INTERN_POLYNOMIAL {

  inline static void swallow(std::istream &is, char d) {
      char c;
      do c = is.get(); while (isspace(c));
      if (c != d) CGAL_assertion_msg( false, "input error: unexpected character in polynomial");
  }
} // namespace INTERN_POLYNOMIAL

template <class NT>
Polynomial<NT> Polynomial<NT>::input_ascii(std::istream &is) {
    char c;
    int degr = -1, i;

    INTERN_POLYNOMIAL::swallow(is, 'P');
    INTERN_POLYNOMIAL::swallow(is, '[');
    is >> CGAL::iformat(degr);
    if (degr < 0) {
        CGAL_assertion_msg( false, "input error: negative degree of polynomial specified");
    }
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT> p(TAG, degr+1);

    do c = is.get(); while (isspace(c));
    do {
        if (c != '(') CGAL_assertion_msg( false, "input error: ( expected");
        is >> CGAL::iformat(i);
        if (!(i >= 0 && i <= degr && p[i] == NT(0))) {
            CGAL_assertion_msg( false, "input error: invalid exponent in polynomial");
        };
        INTERN_POLYNOMIAL::swallow(is, ',');
        is >> CGAL::iformat(p.coeff(i));
        INTERN_POLYNOMIAL::swallow(is, ')');
        do c = is.get(); while (isspace(c));
    } while (c != ']');

    p.reduce();
    p.simplify_coefficients();
    return p;
}



//
// Non-Member Functions
//

// sign() forwarded to the sign() member function
template <class NT> inline 
CGAL::Sign sign(const Polynomial<NT>& p) { return p.sign(); }

// the non-member variants of diff() etc.
template <class NT> inline
Polynomial<NT> diff(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.diff(); return q; }

template<class NT> inline
Polynomial<NT> scale_up(const Polynomial<NT>& p, const NT& a)
{ Polynomial<NT> q(p); q.scale_up(a); return q; }

template<class NT> inline
Polynomial<NT> scale_down(const Polynomial<NT>& p, const NT& b)
{ Polynomial<NT> q(p); q.scale_down(b); return q; }

template<class NT> inline
Polynomial<NT> scale(const Polynomial<NT>& p, const NT& a, const NT& b)
{ Polynomial<NT> q(p); q.scale(a, b); return q; }

template<class NT> inline
Polynomial<NT> translate_by_one(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.translate_by_one(); return q; }

template<class NT> inline
Polynomial<NT> translate(const Polynomial<NT>& p, const NT& c)
{ Polynomial<NT> q(p); q.translate(c); return q; }

template<class NT> inline
Polynomial<NT> translate(const Polynomial<NT>& p, const NT& a, const NT& b)
{ Polynomial<NT> q(p); q.translate(a, b); return q; }

template<class NT> inline
Polynomial<NT> reversal(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.reversal(); return q; }


//! return an upper bound on the absolute value of all real roots of \c P.
/*! The upper bound is a power of two. Only works for univariate polynomials.
 *  \pre \c NT must be \c RealComparable.
 *  \relates NiX::Polynomial
 */
template <class NT>
NT weak_upper_root_bound(const Polynomial<NT>& P) { 
    // code comes from Kurt Mehlhorn
    // see [Mignotte, 1992], p.144 for a proof
    CGAL_precondition(Polynomial_traits<NT>::poly_nesting_depth == 0);
    typename Real_embeddable_traits<NT>::Abs abs;
    const int n = P.degree();
    NT x(1); 
    NT val;
    for (;;) {
        val = -abs(P[n]);
        for (int i = n-1; i >= 0; i--) {
            val = val*x + abs(P[i]);
        }
        if (val < NT(0)) return x;
        x *= NT(2);
    }
}

//! return the number of sign variations in the coefficient sequence of \c P.
/*! This is the number of sign changes (+ to - or - to +) in the
 *  coefficient sequence of the polynomial, ignoring zeroes.
 *  Only meaningful for univariate polynomials.
 *  \pre \c NT must be \c RealComparable.
 *  \relates NiX::Polynomial
 */
template <class NT>
int sign_variations(const Polynomial<NT>& P) { 
    typename Real_embeddable_traits<NT>::Sign sign;
    const int n = P.degree();
    int variations = 0;
    int old_sign = sign(P[n]); // never zero unless P is zero
    for (int i = n-1; i >= 0; i--) {
        int s = sign(P[i]);
        if (s == 0) continue;
        if (old_sign != s) {
            old_sign = s;
            variations++;
        }
    }
    return variations;
}


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

    INTERN_POLYNOMIAL::Creation_tag TAG;    
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

#ifndef NiX_POLY_USE_OLD_PSEUDODIV

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
       
        CGAL_postcondition( !AST::Is_exact::value || Polynomial<NT>(D)*A == Q*B + R);
        return;
    }
    const NT d = B.lcoeff();
    int e = delta + 1;
    D = ipower(d, e);
    INTERN_POLYNOMIAL::Creation_tag TAG;
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
    } while (delta > 0 || delta == 0 && !R.is_zero());
    // funny termination condition because deg(0) = 0, not -\infty

    NT q = ipower(d, e);
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
    INTERN_POLYNOMIAL::Creation_tag TAG;
    q = Polynomial<NT>(TAG, delta );
    NT G = g[gd]; // highest order coeff of g
    D = ipower(G, delta);
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

#endif // NiX_POLY_USE_OLD_PSEUDODIV

// The following former parts of Polynomial.h have been moved to the new files
// <NiX/polynomial_gcd.h> (items 2,3,5) and <NiX/prs_resultant.h> (item 4):
//
// 2) gcd (basic form without cofactors)
// 3) extended gcd computation (with cofactors)
// 4) resultant computation from polynomial remainder sequences (PRS)
// 5) square-free factorization

// This subroutine has been retained here for use in both new files.
namespace INTERN_POLYNOMIAL {
template <class NT> inline
void hgdelta_update(NT& h, const NT& g, int delta) {
    typename Algebraic_structure_traits<NT>::Integral_division idiv;

    // compute h = h^(1-delta) * g^delta
    switch (delta) {
    case 0:
        // h = h;
        break;
    case 1:
        h = g;
        break;
    default:
        h = idiv(ipower(g, delta), ipower(h, delta-1));
        break;
    }
}
} // namespace INTERN_POLYNOMIAL

// } // namespace NiX
CGAL_END_NAMESPACE

#include <CGAL/polynomial_gcd.h> // used above for NT_traits<Poly...>::Gcd
#include <CGAL/prs_resultant.h>  // for compatibility

CGAL_BEGIN_NAMESPACE

// Cofraction_traits added by Michael Hemmer
/*namespace NiX{
namespace Intern{
    template <class NT, class TAG> class Cofraction_traits_base;

    template <class NT_, class TAG>
    class Cofraction_traits_base<NiX::Polynomial<NT_>, TAG > {
        typedef NT_ NT;
    public:
        typedef Polynomial<NT_>                        Numerator_type;
        typedef ::LiS::False_tag                       Is_composable;
        typedef ::LiS::Null_tag                       Denominator_type;
        typedef ::LiS::Null_tag                       Type;
        typedef ::LiS::Null_tag                       Compose;
    };
    
    template <class NT_>
    class Cofraction_traits_base<NiX::Polynomial<NT_>, LiS::True_tag > {
        typedef NT_ NT;
        typedef NiX::Cofraction_traits<NT_> CFT_NT;
    public:
        typedef Polynomial<NT>                             Numerator_type;
        typedef ::LiS::True_tag                            Is_composable;
        typedef typename CFT_NT::Denominator_type          Denominator_type;
        typedef Polynomial<typename CFT_NT::Type> Type;
        
        class Compose {
        public:
            //! first argument type
            typedef Numerator_type   first_argument_type;
            //! second argument type
            typedef Denominator_type second_argument_type;
            //! result type
            typedef Type    result_type;
            //! Compose fraction
            Type operator() (Numerator_type num, 
                                      Denominator_type den 
                                      = Denominator_type(1)){
                Type tmp1; NiX::convert_to(num,tmp1);
                Type tmp2; NiX::convert_to(den,tmp2);
                return tmp1/tmp2;
            }
        };
    };
} //namespace Intern
*/
/*! \ingroup NiX_Cofraction_traits_specs
 *  \brief Specialization of Cofraction_traits for NiX::Polynomial<NT>.
 */
//template<class NT>
/*class Cofraction_traits<Polynomial<NT> > :
    public Intern::Cofraction_traits_base<
      Polynomial<NT>,
      typename Cofraction_traits<NT>::Is_composable>{
    //nothing new
};*/ 

template <class COEFF>
struct Needs_parens_as_product<Polynomial<COEFF> >{
    typedef Polynomial<COEFF> Poly;
    bool operator()(const Poly& x){ return (x.degree() > 0); }
};


// COERCION_TRAITS BEGIN 

//Coercion_traits_polynomial-----------------------------------
// If there is a Polynomial_traits, valid for more than one Polynomial
// class this part should be adapted, using a Polynomial_traits 
// and the nesting_depth 
template <class A,class B>
class Coercion_traits_for_level<Polynomial<A>, Polynomial<B>, CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;            
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const { 
            typename CT::Cast cast; 
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const Polynomial<B>& poly) const {  
            typename CT::Cast cast;  
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
    }; 
};
        
template <class A,class B>
class Coercion_traits_for_level<Polynomial<A>,B ,CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const {
            typename CT::Cast cast;
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                       ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        } 
    };                                                        
}; 
template <class A,class B> 
class Coercion_traits_for_level<B,Polynomial<A>,CTL_POLYNOMIAL  >
    :public Coercion_traits_for_level<Polynomial<A>,B,CTL_POLYNOMIAL >
{};

// COERCION_TRAITS END










CGAL_END_NAMESPACE

//
// trailing documentation
//

// Literature reference
//
// [Akritas, 1989]
//      Alkiviadis G. Akritas
//      Elements of Computer Algebra With Applications
//      Wiley, New York, 1989.
//
// [Cohen, 1993]
//      Cohen, Henri
//      A Course in Computational Algebraic Number Theory
//      Springer GTM 138, 1993
//
// [Cox et al, 1997]
//      David Cox; John Little; Donal O'Shea
//      Ideals, Varieties, and Algorithms
//      2nd ed., Springer UTM, 1997
//
// [Geddes et al, 1992]
//      Geddes, Keith O. and Czapor, Stephen R. and Labahn, George
//      Algorithms for Computer Algebra
//      Kluwer, 1992
//
// [Mignotte, 1992]
//      Mignotte, Maurice
//      Mathematics for Computer Algebra
//      Springer, 1992
//
// [PARI]
//      (a computer algebra system by Henri Cohen and collaborators)
//      http://www.parigp-home.de/

#endif  // NiX_POLYNOMIAL_H

// EOF
