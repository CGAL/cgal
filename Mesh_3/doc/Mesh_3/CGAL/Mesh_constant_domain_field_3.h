namespace CGAL {

/*!
\ingroup PkgMesh3Domains

The class `Mesh_constant_domain_field_3` is a model of concept `MeshDomainField_3`. It provides
a constant field accessible using queries on 3D-points.

The class `Mesh_constant_domain_field_3` can also be customized through `set_size()` operations to become
a piecewise constant field, i.e.\ a sizing field with a constant size on each subpart
of the domain.

\tparam Gt is the geometric traits class. It must match the type `Triangulation::Geom_traits`,
where `Triangulation` is the nested type of the model of `MeshComplex_3InTriangulation_3` used
in the meshing process.

\tparam Index is the type of index of the vertices of the triangulation.
It must match the type `%Index` of the model of `MeshDomain_3` used in the meshing process.

\cgalModels `MeshDomainField_3`

\sa `MeshDomainField_3`

*/
template< typename Gt, typename Index >
class Mesh_constant_domain_field_3 {
public:

/// \name Types
/// @{

/*!
Numerical type.
*/
typedef Gt::FT FT;

/*!
Point type.
*/
typedef Gt::Point_3 Point_3;

/*!
Type of index of the vertices of the triangulation.
*/
typedef Index Index;

/// @}

/// \name Creation
/// @{

/*!

Builds a constant domain field with size `size`.
*/
Mesh_constant_domain_field_3(FT size);

/// @}

/// \name Operations
/// @{

/*!

Sets the size such as `operator()` will return size `size`
at any query point of dimension `dimension` and index `index`.
*/
void set_size(FT size, int dimension, const Index& index);

/// @}

}; /* end Mesh_constant_domain_field_3 */
} /* end namespace CGAL */
