
namespace CGAL {

/*!
\ingroup nt_rrational

The class `Rational_traits` can be used to determine the type of the numerator
and denominator of a rational number type as `Quotient`, `Gmpq`,
`mpq_class` or `leda_rational`.

*/
template< typename NT >
class Rational_traits {
public:

/// \name Types
/// @{

/*!
the type of the numerator and denominator.
*/
typedef unspecified_type RT;

/// @}

/// \name Operations
/// @{

/*!
returns the numerator of `r`.
*/
RT numerator (const NT & r) const;

/*!
returns the denominator of `r`.
*/
RT denominator (const NT & r) const;

/*!
returns self.
*/
NT make_rational(const NT & x) const;

/*!
constructs a rational number `p.first/p.second`.
*/
NT make_rational(const std::pair<RT,RT> & p) const;

/*!
constructs a rational number.
*/
NT make_rational(const RT & n, const RT & d) const;

/*!
constructs a rational number.
*/
NT make_rational(const NT & n, const NT & d) const;

/// @}

}; /* end Rational_traits */
} /* end namespace CGAL */
