
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` provides several operators
to construct objects of type `PolynomialTraits_d::Polynomial_d`.

\cgalRefines `AdaptableFunctor`
\cgalRefines `CopyConstructible`
\cgalRefines `DefaultConstructible`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/

class PolynomialTraits_d::ConstructPolynomial {
public:

/// \name Types
/// @{

/*!

*/
typedef PolynomialTraits_d::Polynomial_d result_type;

/// @}

/// \name Operations
/// @{

/*!
Construct the zero polynomial.
*/
result_type operator()();

/*!
Construct the constant polynomial equal to \f$ i\f$.
*/
result_type operator()(int i);

/*!
Construct the constant polynomial equal to \f$ a\f$.
*/
result_type operator()
(PolynomialTraits_d::Innermost_coefficient_type a);

/*!
Construct the polynomial equal to \f$ a\f$.
*/
result_type operator()
(PolynomialTraits_d::Coefficient_type a);

/*!
\pre The value type of `InputIterator` is `PolynomialTraits_d::Coefficient_type`.

The operator constructs the a polynomial from the iterator range,
with respect to the outermost variable, \f$ x_{d-1}\f$.

The range starts with the coefficient for \f$ x_{d-1}^0\f$.

In case the range is empty, the zero polynomial is constructed.

*/
template < class InputIterator >
result_type operator()(InputIterator begin, InputIterator end);

/*!

Constructs a `Polynomial_d` from a given iterator range of
`std::pair<CGAL::Exponent_vector, PolynomialTraits_d::Innermost_coefficient_type>`.

The optional parameter `is_sorted` indicates whether the given iterator range is
already sorted.

\pre The value type of `InputIterator` is `std::pair<CGAL::Exponent_vector, PolynomialTraits_d::Innermost_coefficient_type>`.
\pre Each `CGAL::Exponent_vector` must have size \f$ d\f$.
\pre All appearing `CGAL::Exponent_vector`s are different.

*/
template < class InputIterator >
result_type operator()(InputIterator begin, InputIterator end, bool is_sorted= false);

/// @}

}; /* end PolynomialTraits_d::ConstructPolynomial */

