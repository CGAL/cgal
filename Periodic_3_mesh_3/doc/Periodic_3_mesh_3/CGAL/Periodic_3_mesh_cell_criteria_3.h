namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3MeshClasses

The class `Periodic_3_mesh_cell_criteria_3` is a model of `MeshCellCriteria_3`. It provides,
for the mesh tetrahedra, a uniform shape criteria
and a sizing field which may be a uniform or variable field.

This class is in essence identical to the non-periodic `CGAL::Mesh_cell_criteria_3`.
The difference is that the fundamental domain in which the construction
of the mesh is done (see \ref PkgPeriodic_3_mesh_3) must be known
by the cell criteria to be able to consider other copies of points or vertices.
For example, a cell of the periodic mesh might intersect the boundary of the
fundamental domain; to properly compute its cell ratio, some of the vertices
must be considered in a neighboring fundamental domain.

\tparam Tr must be identical to the nested type
`Triangulation` of the instance used as model of
`MeshComplex_3InTriangulation_3`.

\cgalModels `MeshCellCriteria_3`

\sa `CGAL::Periodic_3_mesh_criteria_3<Tr>`
\sa `CGAL::Periodic_3_mesh_edge_criteria_3<Tr>`
\sa `CGAL::Periodic_3_mesh_facet_criteria_3<Tr>`

\sa `CGAL::make_periodic_3_mesh_3()`

*/
template< typename Tr >
class Periodic_3_mesh_cell_criteria_3 {
public:

/// \name Types
/// @{

/*!
Numerical type
*/
typedef Tr::FT FT;

/// @}

/// \name Creation
/// @{

/*!
Returns an object to serve as default criteria for cells.

\param periodic_domain the fundamental domain.
\param radius_edge_bound is the upper bound for the radius-edge ratio
                         of the tetrahedra.
\param radius_bound is a uniform upper bound for the circumradii of the tetrahedra in the mesh.
                    See section \ref introsecparam for further details.
Note that if one parameter is set to 0, then its corresponding criteria is ignored.

\tparam MD must be either a cuboid of type `Tr::Iso_cuboid` or a periodic mesh domain
        that is a model of `Periodic_3MeshDomain_3`.
*/
template<class MD>
Periodic_3_mesh_cell_criteria_3(const MD& periodic_domain,
                                FT radius_edge_bound,
                                FT radius_bound);

/*!
Returns an object to serve as default criteria for facets. The type `SizingField` must
be a model of the concept `MeshDomainField_3`. The behavior and semantic of the arguments are the same
as above, except that the radius bound parameter is a functional instead of a constant.
*/
template<class MD, class SizingField>
Periodic_3_mesh_cell_criteria_3(const MD& periodic_domain,
                                FT radius_edge_bound,
                                SizingField radius_bound);
/// @}

}; /* end Periodic_3_mesh_cell_criteria_3 */
} /* end namespace CGAL */
