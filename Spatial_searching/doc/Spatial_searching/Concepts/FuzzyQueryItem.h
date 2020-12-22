/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `FuzzyQueryItem` describes the requirements for fuzzy `d`-dimensional spatial objects.

\cgalHasModel `CGAL::Fuzzy_sphere<Traits>`
\cgalHasModel `CGAL::Fuzzy_iso_box<Traits>`

*/

class FuzzyQueryItem {
public:

/// \name Types
/// @{

/*!
Dimension Tag.
*/
typedef unspecified_type Dimension;

/*!
represents a `d`-dimensional point.
*/
typedef unspecified_type Point_d;

/*!
Number type.
*/
typedef unspecified_type FT;

/// @}

/// \name Operations
/// @{

/*!
tests whether the query item contains `p`.
*/
bool contains(Point_d p) const;

/*!
\note Optional: must be defined when used with a `Kd_tree` where `EnablePointsCache` is set to `Tag_true`.

tests whether the query item contains the point whose Cartesian coordinates
are contained in the range [`begin`, `end`).
*/
template <typename Coord_iterator>
bool contains_point_given_as_coordinates(Coord_iterator begin, Coord_iterator end) const;

/*!
tests whether the inner approximation of the spatial object intersects a rectangle
associated with a node of a tree.
*/
bool inner_range_intersects(const Kd_tree_rectangle<FT,Dimension>& rectangle) const;

/*!
tests whether the outer approximation of the spatial object encloses the rectangle
associated with a node of a tree.
*/
bool outer_range_contains(const Kd_tree_rectangle<FT,Dimension>& rectangle) const;

/// @}

}; /* end FuzzyQueryItem */
