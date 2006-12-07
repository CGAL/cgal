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

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_TYPE_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_TYPE_H

CGAL_BEGIN_NAMESPACE

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


CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_POLYNOMIAL_TYPE_H
