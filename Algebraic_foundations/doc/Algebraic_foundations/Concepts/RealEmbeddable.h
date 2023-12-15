
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of this concepts represents numbers that are embeddable on the real
axis. The type obeys the algebraic structure and compares two values according
to the total order of the real numbers.

Moreover, `CGAL::Real_embeddable_traits< RealEmbeddable >` is a model of
`RealEmbeddableTraits`

with:

- \link RealEmbeddableTraits::Is_real_embeddable `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_real_embeddable` \endlink set to `Tag_true`

and functors :

- \link RealEmbeddableTraits::Is_zero `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_zero` \endlink which is a model of `RealEmbeddableTraits_::IsZero`

- \link RealEmbeddableTraits::Abs `CGAL::Real_embeddable_traits< RealEmbeddable >::Abs` \endlink which is a model of `RealEmbeddableTraits_::Abs`

- \link RealEmbeddableTraits::Sgn `CGAL::Real_embeddable_traits< RealEmbeddable >::Sgn` \endlink which is a model of `RealEmbeddableTraits_::Sgn`

- \link RealEmbeddableTraits::Is_positive `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_positive` \endlink which is a model of `RealEmbeddableTraits_::IsPositive`

- \link RealEmbeddableTraits::Is_negative `CGAL::Real_embeddable_traits< RealEmbeddable >::Is_negative` \endlink which is a model of `RealEmbeddableTraits_::IsNegative`

- \link RealEmbeddableTraits::Compare `CGAL::Real_embeddable_traits< RealEmbeddable >::Compare` \endlink which is a model of `RealEmbeddableTraits_::Compare`

- \link RealEmbeddableTraits::To_double `CGAL::Real_embeddable_traits< RealEmbeddable >::To_double` \endlink which is a model of `RealEmbeddableTraits_::ToDouble`

- \link RealEmbeddableTraits::To_interval `CGAL::Real_embeddable_traits< RealEmbeddable >::To_interval` \endlink which is a model of `RealEmbeddableTraits_::ToInterval`

Remark:

If a number type is a model of both `IntegralDomainWithoutDivision` and
`RealEmbeddable`, it follows that the ring represented by such a number type
is a sub-ring of the real numbers and hence has characteristic zero.

\cgalRefines{EqualityComparable,LessThanComparable}

\sa `RealEmbeddableTraits`

*/

class RealEmbeddable {
public:

/// \name Operations
/// @{

/*!

*/
bool operator==(const RealEmbeddable &a,
const RealEmbeddable &b);


/*!

*/
bool operator!=(const RealEmbeddable &a,
const RealEmbeddable &b);

/*!

*/
bool operator< (const RealEmbeddable &a,
const RealEmbeddable &b);

/*!

*/
bool operator<=(const RealEmbeddable &a,
const RealEmbeddable &b);

/*!


*/
bool operator> (const RealEmbeddable &a,
const RealEmbeddable &b);

/*!

\relates RealEmbeddable
*/
bool operator>=(const RealEmbeddable &a,
const RealEmbeddable &b);

/// @}

}; /* end RealEmbeddable */

