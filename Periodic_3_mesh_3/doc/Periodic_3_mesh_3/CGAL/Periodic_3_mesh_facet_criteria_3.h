namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3MeshClasses

The class `Periodic_3_mesh_facet_criteria_3` is a model of `MeshFacetCriteria_3`.
It provides a uniform bound for the shape criterion, a uniform or variable sizing field
for the size criterion and/or a uniform or variable distance field
for the approximation error criterion.

This class is in essence identical to the non-periodic `CGAL::Mesh_cell_criteria_3`.
The difference is that the fundamental domain in which the construction
of the mesh is done (see \ref PkgPeriodic_3_mesh_3) must be known
by the cell criteria to be able to consider other copies of points or vertices.
For example, a facet of the periodic mesh might intersect the boundary of the
fundamental domain; to properly compute its size, some of the vertices
must be considered in a neighboring fundamental domain.

\tparam Tr must be identical to the nested type
`Triangulation` of the instance used as model of
`MeshComplex_3InTriangulation_3`.

\cgalModels `MeshFacetCriteria_3`

\sa `CGAL::Periodic_3_mesh_criteria_3<Tr>`
\sa `CGAL::Periodic_3_mesh_edge_criteria_3<Tr>`
\sa `CGAL::Periodic_3_mesh_cell_criteria_3<Tr>`
\sa `CGAL::make_periodic_3_mesh_3()`

\sa `CGAL::Mesh_facet_topology`
\sa `MeshDomainField_3`
*/
template< typename Tr >
class Periodic_3_mesh_facet_criteria_3 {
public:

/// \name Types
/// @{

/*!
Numerical type
*/
typedef Tr::Geom_traits::FT FT;

/// @}

/// \name Creation
/// @{

/*!
Returns an object to serve as criteria for facets.
\param periodic_domain the fundamental domain.
\param angle_bound is the lower bound for the angle in degrees of the
surface mesh facets.
\param radius_bound is a uniform upper bound
for the radius of the surface Delaunay balls.
\param distance_bound is an upper bound for the center-center distances
of the surface mesh facets.
\param topology is the set of topological constraints
which have to be verified by each surface facet. See
section \ref Mesh_3DelaunayRefinement for further details.
Note that if one parameter is set to 0, then its corresponding criteria is ignored.

\tparam MD must be either a cuboid of type `Tr::Iso_cuboid` or a periodic mesh domain
        that is a model of `Periodic_3MeshDomain_3`.
*/
template <class MD>
Periodic_3_mesh_facet_criteria_3(const MD& periodic_domain,
                                 FT angle_bound,
                                 FT radius_bound,
                                 FT distance_bound,
                                 Mesh_facet_topology topology = FACET_VERTICES_ON_SURFACE);

/*!
Returns an object to serve as criteria for facets. The types `SizingField` and
`DistanceField` must be models of the concept `MeshDomainField_3`.
The behavior and semantic of the arguments are the same
as above, except that the radius and distance bound parameters are functionals instead of constants.
*/
template <class MD, class SizingField, class DistanceField>
Periodic_3_mesh_facet_criteria_3(const MD& periodic_domain,
                                 FT angle_bound,
                                 SizingField radius_bound,
                                 DistanceField distance_bound,
                                 Mesh_facet_topology topology = FACET_VERTICES_ON_SURFACE);
/// @}

}; /* end Periodic_3_mesh_facet_criteria_3 */
} /* end namespace CGAL */
