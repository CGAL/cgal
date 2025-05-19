namespace CGAL {

/*!
\ingroup RangeQueryItemClasses

The class `Fuzzy_iso_box` implements fuzzy `d`-dimensional (closed) iso boxes.
A fuzzy iso box with fuzziness value \f$ \epsilon\f$ has as inner and outer approximations
a box respectively eroded and dilated by a `d`-dim square with side length \f$ \epsilon\f$.
Points are returned depending on their location respective to the approximations,
as follows:
- points in the inner approximation (boundary included) are always returned.
- points in the outer approximation (boundary included) and not in the previous
  category may or may not be returned.
- points that do not fit in the two previous categories are never returned.

\tparam Traits must be a model of the concept
`SearchTraits`, for example `CGAL::Search_traits_2<CGAL::Simple_cartesian<double> >`.

\cgalModels{FuzzyQueryItem}
*/
template< typename Traits >
class Fuzzy_iso_box {
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
/// @{

/*!
Construct a fuzzy iso box
specified by the minimal iso box containing `p` and `q` and fuzziness value `epsilon`.

\pre `p` must be lexicographically smaller than `q`.
*/
Fuzzy_iso_box(Point_d p, Point_d q, FT epsilon=FT(0),Traits t=Traits());

/*!
Construct a fuzzy iso box specified by the minimal iso box containing `p` and `q` and fuzziness value `epsilon`.

\attention Only available in case `Traits` is
`Search_traits_adapter<Key,PointPropertyMap,BaseTraits>`.

\pre `p` must be lexicographically smaller than `q`.
*/
Fuzzy_iso_box(Traits::Base::Point_d p, Traits::Base::Point_d q, FT epsilon=FT(0), Traits t=Traits());

/// @}

/// \name Operations
/// @{

/*!
Test whether the fuzzy iso box contains `p`.
*/
bool contains(Point_d p) const;

/*!
Test whether the fuzzy iso box contains the point whose Cartesian
coordinates are contained in the range [`begin`, `end`).
*/
template <typename Coord_iterator>
bool contains_point_given_as_coordinates(Coord_iterator begin, Coord_iterator end) const;

/*!
Test whether the inner box intersects the rectangle
associated with a node of a tree.
*/
bool inner_range_intersects(const Kd_tree_rectangle<FT,D>& rectangle) const;

/*!
Test whether the outer box encloses the rectangle associated with a node of a tree.
*/
bool outer_range_contains(const Kd_tree_rectangle<FT,D>& rectangle) const;

/// @}

}; /* end Fuzzy_iso_box */
} /* end namespace CGAL */
