namespace CGAL {

/*!
\ingroup RangeQueryItemClasses

The class `Fuzzy_sphere` implements fuzzy `d`-dimensional spheres.
A fuzzy sphere with radius \f$ r\f$ and fuzziness value \f$ \epsilon\f$ has
as inner approximation a sphere with radius \f$ r-\epsilon\f$ and
as outer approximation a sphere with radius \f$ r+\epsilon\f$.

\attention The fuzziness of a `Fuzzy_sphere` is specified by a parameter \f$ \epsilon\f$
denoting a maximal allowed distance to the boundary of a sphere.
If the distance to the boundary is greater than \f$ \epsilon\f$, points inside the
object are always reported and points outside the sphere are never reported.
Points whose distance to the boundary is less than or equal to \f$ \epsilon\f$
may or may not be reported. Subsequently, points on the inner and outer spheres
may or may not be reported. Specifically when \f$ \epsilon = 0\f$, points
on the sphere of radius \f$ r\f$ may or may not be reported.

\tparam Traits must be a model of the concept
`SearchTraits`, for example `CGAL::Cartesian_d<double>`.

\cgalModels `FuzzyQueryItem`

\sa `FuzzyQueryItem`

*/
template< typename Traits >
class Fuzzy_sphere {
public:

/// \name Types
/// @{

/*!
Dimension Tag.
*/
typedef unspecified_type D;

/*!
Point type.
*/
typedef Traits::Point_d Point_d;

/*!
Number type.
*/
typedef Traits::FT FT;

/// @}

/// \name Creation
///
/// @{

/*!
Construct a fuzzy sphere
centered at `center` with radius `radius` and fuzziness value `epsilon`.
*/
Fuzzy_sphere(Point_d center, FT radius, FT epsilon=FT(0),Traits t=Traits());

/*!
Construct a fuzzy sphere centered at `center` with radius `radius` and fuzziness value `epsilon`.
\attention Only available in case `Traits` is `Search_traits_adapter<Key,PointPropertyMap,BaseTraits>`.
*/
Fuzzy_sphere(Traits::Base::Point_d center, FT radius, FT epsilon=FT(0), Traits t=Traits());

/// @}

/// \name Operations
/// @{

/*!
Test whether the fuzzy sphere contains `p`.

That is, whether the distance between the center of the fuzzy sphere and `p` is
less than \f$ r\f$.
*/
bool contains(const Point_d& p) const;

/*!
Test whether the inner sphere intersects the rectangle
associated with a node of a tree.

That is, whether the minimal distance between the center of the fuzzy sphere and
`rectangle` is less than \f$ r-\epsilon\f$.
*/
bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const;

/*!
Test whether the outer sphere encloses the rectangle associated with a node of a tree.

That is, whether the maximal distance between the center of the fuzzy sphere and
`rectangle` is less than \f$ r+\epsilon\f$.
*/
bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const;

/// @}

}; /* end Fuzzy_sphere */
} /* end namespace CGAL */
