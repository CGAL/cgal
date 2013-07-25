
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `SurfaceMeshCellBase_3` describes the cell base type 
of the three dimensional triangulation used 
to embed the surface mesh. 

More precisely, 
the first template parameter `SurfaceMeshC2T3` of the function template
`CGAL::make_surface_mesh()` 
is a model of the concept 
`SurfaceMeshComplex_2InTriangulation_3` 
which describes a data structure to store 
a pure two dimensional complex 
embedded in a three dimensional triangulation. 
In particular, the type `SurfaceMeshC2T3` is required to provide 
a three dimensional triangulation type 
`SurfaceMeshC2T3::Triangulation_3`. 
The concept `SurfaceMeshCellBase_3` describes the cell base type 
required in this triangulation type. 

\cgalRefines `TriangulationCellBase_3` 
The concept `SurfaceMeshCellBase_3` adds four markers to mark
the facets of the triangulation that belong to the two dimensional
complex, and four markers that are helpers used in some operations to
mark for instance the facets that have been visited.
<br>
This concept also provides storage for the center of a
Delaunay surface ball. Given a surface and a 3D Delaunay
triangulation, a Delaunay surface ball is a ball circumscribed to a
facet of the triangulation and centered on the surface and empty of
triangulation vertices.  Such a ball does exist when the facet is part
of the restriction to the surface of a three dimensional
triangulation.  In the following we call surface center of a facet,
the center of its biggest Delaunay surface ball.

\cgalHasModel `CGAL::Surface_mesh_cell_base_3<Gt,Vb>`
\cgalHasModel `CGAL::Surface_mesh_default_triangulation_3::Cell`

\sa `SurfaceMeshTriangulation_3` 
\sa `SurfaceMeshComplex_2InTriangulation_3` 
\sa `CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr>` 
\sa `CGAL::Surface_mesh_default_triangulation_3` 
\sa `CGAL::make_surface_mesh()` 

*/

class SurfaceMeshCellBase_3 {
public:

/// \name Types 
/// @{

/*!
The point type, required to match the point type 
of the three dimensional 
triangulation in which the surface mesh is embedded. 
*/ 
typedef unspecified_type Point; 

/// @} 

/// \name Operations 
/// @{

/*!
returns `true`, if `facet(i)` is in the 2D complex. 
*/ 
bool is_facet_on_surface(int i); 

/*!
Sets `facet(i)` as part of the 2D complex, if `b` is `true`, 
and `NOT_IN_COMPLEX`, otherwise. 
*/ 
void set_facet_on_surface(int i, bool b ); 

/*!
Returns `true`, if `facet(i)` has been visited, 
`false` otherwise. 
*/ 
bool is_facet_visited (int i); 

/*!
Marks `facet(i)` as visited, if `b` is `true`, 
and non visited otherwise. 
*/ 
void set_facet_visited (int i, bool b); 

/*!
Returns a const reference to the surface center of `facet(i)`. 
*/ 
const Point& get_facet_surface_center(int i); 

/*!
Sets point `p` as the surface center of `facet(i)`. 
*/ 
void set_facet_surface_center (int i, Point p); 

/// @}

}; /* end SurfaceMeshCellBase_3 */

