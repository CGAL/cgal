
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `Functor` substitutes all variables of a given multivariate
`PolynomialTraits_d::Polynomial_d` \f$ p\f$ by the values given in the
iterator range, where begin refers the value for the innermost variable.
In contrast to `PolynomialTraits_d::Substitute` the given polynomial \f$ p\f$
is interpreted as a homogeneous polynomial.
Hence the iterator range is required to be of length `PolynomialTraits_d::d+1`.

For instance the polynomial \f$ p(x_0,x_1) = x_0^2x_1^3+x_1^4\f$ is interpreted as the homogeneous
polynomial \f$ p(x_0,x_1,w) = x_0^2x_1^3+x_1^4w^1\f$.

\cgalRefines{CopyConstructible,Assignable,DefaultConstructible}

\cgalHeading{Types}

Note that the `result_type` is the coercion type of the value type of the
given iterator range and `PolynomialTraits_d::Innermost_coefficient_type`.
In particular `std::iterator_traits<Input_iterator>::%value_type` must be
`ExplicitInteroperable` with `PolynomialTraits_d::Innermost_coefficient_type`.
Hence, it can not be provided as a public type in advance.

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::Substitute`
\sa  `CGAL::Coercion_traits`

*/

class PolynomialTraits_d::SubstituteHomogeneous {
public:

/// \name Operations
/// @{

/*!

Substitutes each variable of \f$ p\f$ by the values given in the iterator range,
where \f$ p\f$ is interpreted as a homogeneous polynomial in all variables.
The begin iterator refers to the innermost variable \f$ x_0\f$.
\pre (`end-begin` == `PolynomialTraits_d::d`)+1

*/
template<class Input_iterator>
result_type operator()(PolynomialTraits_d::Polynomial_d p,
Input_iterator begin, Input_iterator end);

/// @}

}; /* end PolynomialTraits_d::SubstituteHomogeneous */

