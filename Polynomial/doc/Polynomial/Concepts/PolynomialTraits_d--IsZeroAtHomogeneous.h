
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` returns whether a
`PolynomialTraits_d::Polynomial_d` \f$ p\f$ is zero at a given homogeneous point,
which is given by an iterator range.

The polynomial is interpreted as a homogeneous polynomial in all variables.

For instance the polynomial \f$ p(x_0,x_1) = x_0^2x_1^3+x_1^4\f$ is interpreted as the homogeneous
polynomial \f$ p(x_0,x_1,w) = x_0^2x_1^3+x_1^4w^1\f$.

\cgalRefines{AdaptableFunctor,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::IsZeroAtHomogeneous {
public:

/// \name Types
/// @{

/*!

*/
typedef bool result_type;

/// @}

/// \name Operations
/// @{

/*!

Computes whether \f$ p\f$ is zero at the homogeneous point given by the iterator range,
where `begin` is referring to the innermost variable.
\pre (end-begin==`PolynomialTraits_d::d`+1)
\pre `std::iterator_traits< InputIterator >::%value_type` is `ExplicitInteroperable` with `PolynomialTraits_d::Innermost_coefficient_type`.

*/
template <class InputIterator>
result_type operator()(PolynomialTraits_d::Polynomial_d p,
InputIterator begin,
InputIterator end );

/// @}

}; /* end PolynomialTraits_d::IsZeroAtHomogeneous */

