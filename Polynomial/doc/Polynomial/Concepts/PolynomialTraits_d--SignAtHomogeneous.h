
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` returns the sign of a
`PolynomialTraits_d::Polynomial_d` \f$ p\f$ at a given homogeneous point,
which is given by an iterator range.

The polynomial is interpreted as a homogeneous polynomial in all variables.

For instance the polynomial \f$ p(x_0,x_1) = x_0^2x_1^3+x_1^4\f$ is interpreted as the homogeneous
polynomial \f$ p(x_0,x_1,w) = x_0^2x_1^3+x_1^4w^1\f$.

This functor is well defined if `PolynomialTraits_d::Innermost_coefficient_type` is
`RealEmbeddable`.

\cgalRefines{AdaptableFunctor,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::SignAtHomogeneous {
public:

/// \name Types
/// @{

/*!

*/
typedef CGAL::Sign result_type;

/// @}

/// \name Operations
/// @{

/*!

Returns the sign of \f$ p\f$ at the given homogeneous point, where `begin` is
referring to the innermost variable.
\pre (`end-begin`==`PolynomialTraits_d::d`+1)
\pre `std::iterator_traits< InputIterator >::%value_type` is `ExplicitInteroperable` with `PolynomialTraits_d::Innermost_coefficient_type`.

*/
template <class InputIterator>
result_type operator()(PolynomialTraits_d::Polynomial_d p,
InputIterator begin,
InputIterator end );

/// @}

}; /* end PolynomialTraits_d::SignAtHomogeneous */

