namespace CGAL {

/*!
\ingroup PkgPolynomialClasses

An instance of the data type `Polynomial` represents a
polynomial \f$ p = a_0 + a_1*x + ...a_i*x^i\f$ from the ring \f$ \mathrm{Coeff}[x]\f$.
`Coeff` can itself be an instance of `Polynomial`, yielding a form of
multivariate polynomials.

The template argument `Coeff` must be at
least a model of `IntegralDomainWithoutDivision`.
For all operations naturally involving division, an `IntegralDomain`
is required.
`Polynomial` offers a full set of algebraic operators, i.e.
binary <TT>+</TT>, <TT>-</TT>, <TT>*</TT>, <TT>/</TT> as well as
<TT>+=</TT>, <TT>-=</TT>, <TT>*=</TT>, <TT>/=</TT>; not only for polynomials
but also for a polynomial and a number of the coefficient type.
(The <TT>/</TT> operator must only be used for integral divisions, i.e.
those with remainder zero.)
The operations are implemented naively: <TT>+</TT> and <TT>-</TT> need a number of `Coeff`
operations which is linear in the degree while * is quadratic.
Unary <TT>+</TT> and <TT>-</TT> and (in)equality <TT>==</TT>, <TT>!=</TT> are provided as well.

`Polynomial` is a model of `LessThanComparable` if `Coeff` is a
model of `LessThanComparable`. In this case `Polynomial` provides
comparison operators <TT><</TT>, <TT>></TT>, <TT><=</TT>, <TT>>=</TT>, where the comparison amounts to
lexicographic comparison of the coefficient sequence,
with the coefficient of the highest power taking precedence over
those of lower powers.

`Polynomial` is a model of `Fraction` if `Coeff` is a
model of `Fraction`. In this case Polynomial may be decomposed into a
(scalar) denominator and a compound numerator with a simpler coefficient type.
Often operations can be performed faster on these denominator-free multiples.

`Polynomial` is a model of `Modularizable` if `Coeff` is a
model of `Modularizable`, where the homomorphic map on the polynomials
is simply defined as the canonical extension of the homomorphic map which is
defined on the coefficient type.

\cgalHeading{Implementation}

Inexact and limited-precision types can be used as coefficients,
but at the user's risk. The algorithms implemented were written with
exact number types in mind.

This data type is implemented as a handle type with value semantics
using `CGAL::Handle_with_policy`, where `HandlePolicy`
is `Handle_policy_no_union`.
An important invariant to be preserved by all methods is that
the coefficient sequence does not contain leading zero coefficients
(where leading means at the high-degree end), with the exception that
the zero polynomial is represented by a single zero coefficient.

\cgalModelsBareBegin
\cgalModelsBare{`Polynomial_d`}
\cgalModelsBare{`Assignable`}
\cgalModelsBare{`CopyConstructible`}
\cgalModelsBare{`DefaultConstructible`}
\cgalModelsBare{`EqualityComparable`}
\cgalModelsBare{`ImplicitInteroperable` with `int`}
\cgalModelsBare{`ImplicitInteroperable` with `Coeff`}
\cgalModelsBare{`Fraction` if `Coeff` is model of `Fraction`}
\cgalModelsBare{`LessThanComparable` if `Coeff` is model of `LessThanComparable`}
\cgalModelsBare{`Modularizable` if `Coeff` is model of `Modularizable`}
\cgalModelsBareEnd
*/
template< typename Coeff >
class Polynomial {
public:

/// \name Creation
/// @{

/*!
Introduces an variable initialized with 0.
*/
Polynomial ();

/*!
copy constructor.
*/
Polynomial (const Polynomial& x);

/*!
Constructor from int.
*/
Polynomial (const int &i);

/*!
Constructor from type Coeff.
*/
Polynomial (const Coeff &x);

/*!
Constructor from iterator range with value type Coeff.
*/
template <class Forward_iterator>
Polynomial(Forward_iterator first, Forward_iterator last);

/// @}

/// \name Operations
/// @{

/*!
A const random access iterator pointing to the constant coefficient.
*/
const_iterator begin() const;

/*!
A const random access iterator pointing beyond the leading coefficient.
*/
const_iterator end() const;

/*!
The degree of the polynomial in \f$ x\f$. The degree of the zero polynomial is 0.
*/
int degree() const;

/*!
Const access to the coefficient of \f$ x^i\f$.
*/
const NT& operator[](unsigned int i) const;

/*!
Const access to the leading coefficient.
*/
const NT& lcoeff() const;

/// @}

}; /* end Polynomial */


/*!
Writes `poly` to ostream `os`.
The format depends on the `CGAL::IO::MODE` of `os`.
In case the mode is `CGAL::IO::ASCII` the format is \f$ P[d(0,a_0)(1,a_1)\dots(d,a_d)]\f$,
where \f$ d\f$ is the degree of the polynomial.
The format is output sensitive, that is, coefficients that are zero are not reported.
In case the mode is `CGAL::IO::PRETTY` the format is human readable.
\relates Polynomial
*/

std::ostream& operator<<(std::ostream& os, const Polynomial<Coeff> &poly);

/*!
Reads `poly` from istream `is` in format \f$ P[d(0,a_0)(1,a_1)\dots(d,a_d)]\f$,
the output format in mode `CGAL::IO::ASCII`.
\relates Polynomial
*/

std::istream& operator>>(std::istream& is, const Polynomial<Coeff> &poly);

} /* end namespace CGAL */
