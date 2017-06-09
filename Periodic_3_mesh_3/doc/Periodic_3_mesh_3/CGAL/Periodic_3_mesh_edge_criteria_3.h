namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3MeshClasses

The function object class `Periodic_3_mesh_edge_criteria_3` is a model of `MeshEdgeCriteria_3`. It
provides a bound for the size criterion.

This class is in essence identical to the non-periodic `CGAL::Mesh_cell_criteria_3`.
The difference is that the fundamental domain in which the construction
of the mesh is done (see \ref PkgPeriodic_3_mesh_3) must be known
by the cell criteria to be able to consider other copies of points or vertices.
For example, an edge of the periodic mesh might intersect the boundary of the
fundamental domain; to properly compute its length, one of the vertices
must be considered in a neighboring fundamental domain.

\cgalModels `MeshEdgeCriteria_3`

\tparam Tr must be identical to the nested type
`Triangulation` of the instance used as model of
`MeshComplex_3InTriangulation_3`.

\sa `CGAL::Periodic_3_mesh_criteria_3<Tr>`
\sa `CGAL::Periodic_3_mesh_facet_criteria_3<Tr>`
\sa `CGAL::Periodic_3_mesh_cell_criteria_3<Tr>`
\sa `CGAL::make_periodic_3_mesh_3()`

\sa `CGAL::Mesh_criteria_3<Tr>`
\sa `MeshDomainField_3`
*/
template< typename Tr >
class Periodic_3_mesh_edge_criteria_3 {
public:

/// \name Types
/// @{

/*!
Numerical type.
*/
typedef Tr::Geom_traits::FT FT;

/// @}

/// \name Creation
/// @{

/*!
Returns an object to serve as criteria for edges.

\param periodic_domain the fundamental domain.
\param length_bound is an upper bound for the length of the edges which are used
                    to discretize the curve segments.
Note that if one parameter is set to 0, then its corresponding criteria is ignored.

\tparam MD must be either a cuboid of type `Tr::Iso_cuboid` or a periodic mesh domain
        that is a model of `Periodic_3MeshDomain_3`.
*/
template< class MD >
Periodic_3_mesh_edge_criteria_3(const MD& periodic_domain,
                                FT length_bound);

/*!
Returns an object to serve as criteria for edges. The type `SizingField`
must be a model of concept `MeshDomainField_3`. The behavior and semantic of the argument are the same
as above, except that the length
parameter is a functional instead of a constant.
*/
template <class MD, class SizingField>
Periodic_3_mesh_edge_criteria_3(const MD& periodic_domain,
                                SizingField length_bound);

/// @}

}; /* end Periodic_3_mesh_edge_criteria_3 */
} /* end namespace CGAL */
