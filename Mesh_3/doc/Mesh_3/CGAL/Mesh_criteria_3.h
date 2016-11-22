namespace CGAL {

/*!
\ingroup PkgMesh_3MeshClasses

The class `Mesh_criteria_3` is a model of both concepts `MeshCriteria_3` 
and `MeshCriteriaWithFeatures_3`. 
It gathers the refinement criteria for mesh tetrahedra and 
surface facets where 
surface facets are facets in the mesh approximating the domain surface patches. 
In addition, for domain with exposed 1-dimensional features, 
the class `Mesh_criteria_3` 
handles the definition of a sizing field to guide the discretization of 
1-dimensional features. 

\tparam Tr has to be instantiated with the type used for 
`C3T3::Triangulation`, 
where `C3T3` is the model of `MeshComplex_3InTriangulation_3` 
used in the mesh generation process, 
and `C3T3::Triangulation` its nested triangulation type. 

\cgalModels `MeshCriteria_3`

\cgalHeading{Example}

\code{.cpp}

// Create a Mesh_criteria_3<Tr> object with all cell and facet parameters set
Mesh_criteria_3<Tr> criteria (parameters::facet_angle=30,
                              parameters::facet_size=1,
                              parameters::facet_distance=0.1,
                              parameters::cell_radius_edge_ratio=2,
                              parameters::cell_size=1.5);

// Create a Mesh_criteria_3<Tr> object with size ignored (note that the order changed)
Mesh_criteria_3<Tr> criteria (parameters::cell_radius_edge_ratio=2,
                              parameters::facet_angle=30,
                              parameters::facet_distance=0.1);

\endcode

\sa `MeshCriteria_3` 
\sa `MeshCriteriaWithFeatures_3` 
\sa `MeshCellCriteria_3` 
\sa `MeshEdgeCriteria_3` 
\sa `MeshFacetCriteria_3` 
\sa `MeshDomainField_3` 
\sa `CGAL::Mesh_cell_criteria_3<Tr>` 
\sa `CGAL::Mesh_edge_criteria_3<Tr>` 
\sa `CGAL::Mesh_facet_criteria_3<Tr>` 
\sa `CGAL::Mesh_facet_topology` 

*/
template< typename Tr >
class Mesh_criteria_3 {
public:

/// \name Types 
/// @{

/*!
The criteria for edges. 
*/ 
typedef Mesh_edge_criteria_3<Tr> 
Edge_criteria; 

/*!
The criteria for facets. 
*/ 
typedef Mesh_facet_criteria_3<Tr> 
Facet_criteria; 

/*!
The 
criteria for cells. 
*/ 
typedef Mesh_cell_criteria_3<Tr> Cell_criteria; 

/// @} 

/// \name Creation 
/// @{

/*!
Construction from facet and cell criteria. The edge criteria are ignored 
in this case. 
*/ 
Mesh_criteria_3( 
Facet_criteria facet_criteria, 
Cell_criteria cell_criteria); 

/*!
Construction from edge, facet and cell criteria. 
*/ 
Mesh_criteria_3( 
Edge_criteria edge_criteria, 
Facet_criteria facet_criteria, 
Cell_criteria cell_criteria); 

/*!
\brief Construction from criteria parameters. This constructor uses named 
parameters (from <I>Boost.Parameter</I>) for convenient criteria 
construction.

\tparam FT must be a model of `Field`
\tparam Fieldi (`i`=1,..,4) should be either a model 
of the concept `Field` or a model of the concept `MeshDomainField_3`

The parameters are named parameters and can be passed in any order
provided their name is given (see example below).  The name of each
parameter is the one that is written in the description of the
function (e.g. `parameters::facet_size`).

The description of each parameter is as follows: 

- `edge_size`: a scalar field (resp. a constant) providing a space varying 
(resp. a uniform) 
upper bound for the lengths of curve segment edges. This parameter has to be set to a positive 
value when 1-dimensional features protection is used.

- `facet_angle`: a lower bound for the angles (in degrees) of the 
surface mesh facets. 

- `facet_size`: a scalar field (resp. a constant) describing 
a space varying (resp. a uniform) upper-bound or for the radii of the surface Delaunay balls. 

- `facet_distance`: a scalar field (resp. a constant) describing a space varying (resp. a uniform) 
upper bound for the distance between the facet circumcenter and the center of its surface
Delaunay ball.

- `facet_topology`: the set of topological constraints 
which have to be verified by each surface facet. The default value is 
`CGAL::FACET_VERTICES_ON_SURFACE`. See `Mesh_facet_topology` manual page to 
get all possible values. 

- `cell_radius_edge_ratio`: an upper bound for the radius-edge ratio of the mesh tetrahedra. 

- `cell_size`: a scalar field (resp. a constant) describing 
a space varying (resp. a uniform) upper-bound for the circumradii of the mesh tetrahedra. 

Note that each size or distance parameter can be specified using two ways: either as 
a scalar field or as a numerical value when the field is uniform. 

Each parameter has a special default value `ignored` which means that the 
corresponding criterion will be ignored. 
Numerical sizing or distance values, as well as scalar fields 
should be given in the unit used for coordinates of points in the mesh domain class 
of the mesh generation process. 

*/ 
template<typename FT, 
typename ...Fieldi> 
Mesh_criteria_3( 
Field1 parameters::edge_size = ignored, 
FT parameters::facet_angle = ignored, 
Field2 parameters::facet_size = ignored, 
Field3 parameters::facet_distance = ignored, 
Mesh_facet_topology parameters::facet_topology = CGAL::FACET_VERTICES_ON_SURFACE, 
FT parameters::cell_radius_edge_ratio = ignored, 
Field4 parameters::cell_size = ignored); 

/// @}

}; /* end Mesh_criteria_3 */
} /* end namespace CGAL */
