/*!
\ingroup PkgSpatialSearchingDConcepts
\cgalConcept

\cgalAdvancedConcept
\cgalAdvancedBegin
The concept `Splitter` defines the requirements for a function object class implementing a splitting rule.
\cgalAdvancedEnd

\cgalHasModelsBegin
\cgalHasModels{CGAL::Fair<Traits, SpatialSeparator>}
\cgalHasModels{CGAL::Median_of_rectangle<Traits, SpatialSeparator>}
\cgalHasModels{CGAL::Median_of_max_spread<Traits, SpatialSeparator>}
\cgalHasModels{CGAL::Midpoint_of_rectangle<Traits, SpatialSeparator>}
\cgalHasModels{CGAL::Midpoint_of_max_spread<Traits, SpatialSeparator>}
\cgalHasModels{CGAL::Sliding_fair<Traits, SpatialSeparator>}
\cgalHasModels{CGAL::Sliding_midpoint<Traits, SpatialSeparator>}
\cgalHasModelsEnd

*/

class Splitter {
public:

/// \name Types
/// The parameters `aspect_ratio` and `bucket_size` define the way in which \f$k\f$ - \f$d\f$ tree is constructed.
/// @{

/*!
Number type.
*/
typedef unspecified_type FT;

/*!
Separator.
*/
typedef unspecified_type Separator;

/*!
Typedef to an instantiation of `CGAL::Point_container<Traits>`.
*/
typedef unspecified_type Container;

/// @}

/// \name Operations
/// @{

/*!
Returns the maximal ratio between the largest and smallest side
of a cell allowed for fair splitting.
*/
FT aspect_ratio() const;

/*!
Returns the bucket size of the leaf nodes.
*/
unsigned int bucket_size() const;

/*!

Sets up `sep` and splits points of `c0` into `c0` and `c1` using `sep`.
Container `c0` should contain at least two points and `c1` must be empty.

*/
void operator()(Separator& sep, Container& c0, Container& c1) const;

/// @}

}; /* end Splitter */
