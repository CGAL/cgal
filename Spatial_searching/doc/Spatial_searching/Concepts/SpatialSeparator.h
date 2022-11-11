/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

The concept `SpatialSeparator` defines the requirements for
a separator.

A separator is a `d-1`-dimensional subspace that
separates a `d`-dimensional space into two parts.  One part of
space is said to be on the negative side of the separator and the
other part of space is said to be on the positive side of the
separator.

\cgalHasModel `CGAL::Plane_separator<FT>`

*/

class SpatialSeparator {
public:

/// \name Types
/// @{

/*!
Number type.
*/
typedef unspecified_type FT;

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
Separator();

/// @}

/// \name Operations
/// @{

/*!
Sets the cutting dimension to `d`.
*/
void set_cutting_dimension(int d);

/*!
Sets the cutting value to `v`.
*/
void set_cutting_value(FT v);

/*!
Returns the number of the cutting dimension.
*/
int cutting_dimension();

/*!
Returns the cutting value.
*/
FT cutting_value();

/*!
Returns true if and only if the point `p` is on the negative side of the separator.
*/
template <class Point_d>
bool has_on_negative_side(Point_d p);

/// @}

}; /* end SpatialSeparator */
