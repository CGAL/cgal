/*!
\ingroup PkgMesh3SecondaryConcepts
\cgalConcept

The concept `MeshDomainField_3` describes a scalar field which could be queried
at any point of the space.

\cgalHasModel `CGAL::Mesh_constant_domain_field_3<Gt,%Index>`

\sa `MeshDomain_3`
\sa `MeshDomainWithFeatures_3`
\sa `CGAL::Mesh_edge_criteria_3<Tr>`
\sa `CGAL::Mesh_facet_criteria_3<Tr>`
\sa `CGAL::Mesh_cell_criteria_3<Tr>`

*/

class MeshDomainField_3 {
public:

/// \name Types
/// @{

/*!
Numerical type.
*/
typedef unspecified_type FT;

/*!
Point type.
*/
typedef unspecified_type Point_3;

/*!
%Index type for points. Must match the type `MeshDomain_3::Index`.
*/
typedef unspecified_type Index;

/// @}

/*! \name Operations
The `operator()` returns the field value at a query point.
The field value may depend on the query point location and/or
on the input feature including the query point.
*/
/// @{

/*!

returns the value of the sizing field at the point `p`,
assumed to be included in the input complex feature with dimension `dimension`
and mesh vertex index `index`.
*/
FT operator()(const Point_3& p, int dimension, const Index& index) const;

/// @}

}; /* end MeshDomainField_3 */
