/*!
\ingroup PkgMesh3SecondaryConcepts
\cgalConcept

The concept `MeshPolyline_3` implements a container of points designed to represent a polyline (i.e., a sequence of points).
Types and functions provided in this concept are such as standard template library containers
are natural models of this concept.

\cgalHasModelsBegin
\cgalHasModelsBare{`std::vector<Kernel::Point_3>` for any Kernel of \cgal is a natural model of this concept.}
\cgalHasModelsEnd

\sa `CGAL::Mesh_domain_with_polyline_features_3<MD>`

*/
class MeshPolyline_3 {
public:

/// \name Types
/// @{

/*!
Point type. Must match the type `MeshDomain_3::Point_3`.
*/
typedef unspecified_type value_type;

/*!
A constant iterator on points. Must be a model of Bidirectional iterator and have `value_type` as value type.
*/
typedef unspecified_type const_iterator;

/// @}

/// \name Operations
/// @{

/*!
Returns an iterator on the first point of the polyline.
*/
const_iterator begin();

/*!
Returns the past-the-end iterator for the above iterator.
*/
const_iterator end();

/// @}

}; /* end MeshPolyline_3 */
