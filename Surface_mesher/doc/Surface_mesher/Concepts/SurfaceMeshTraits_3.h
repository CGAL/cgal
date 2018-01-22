
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `SurfaceMeshTraits_3` describes the knowledge that is required on the 
surface to be meshed. A model of this concept 
implements an oracle that is able to tell whether a segment 
(or a ray, or a line) intersects the surface or not, 
and to compute some intersection 
points if any exists. The concept `SurfaceMeshTraits_3` also includes a funcctor able to provide 
a small set of initial points on the surface. 

\cgalHasModel `CGAL::Surface_mesh_traits_generator_3::Type` 

\sa `CGAL::make_surface_mesh()` 

*/

class SurfaceMeshTraits_3 {
public:

/// \name Types 
/// @{

/*!
The type of points. 
This type is required to match 
the point type of the 
three dimensional embedding triangulation 
`C2T3::Triangulation_3`. 
*/ 
typedef unspecified_type Point_3; 

/*!
The type of segments. 
*/ 
typedef unspecified_type Segment_3; 

/*!
The type of rays. 
*/ 
typedef unspecified_type Ray_3; 

/*!
The type of lines. 
*/ 
typedef unspecified_type Line_3; 

/*!
The surface type. 
*/ 
typedef unspecified_type Surface_3; 

/*!
A model of this type provides the operator 

`CGAL::Object operator()(Surface_3 surface, Type1 type1)` 

to compute the intersection of the surface 
with an object of type `Type1` which may be 
`Segment_3`, `Ray_3` or `Line_3` . 
*/ 
typedef unspecified_type Intersect_3; 

/*!
A model of this type provides the following operators 
to construct initial points on the surface: 

`template <class OutputIteratorPoints>` 
<br>
`OutputIteratorPoints operator()(OutputIteratorPoints pts)` 

which outputs a set of points on the surface, 

`template <class OutputIteratorPoints>` 
<br>
`OutputIteratorPoints operator() (OutputIteratorPoints pts, int n)` 

which outputs a set of `n` points on the surface. 
*/ 
typedef unspecified_type Construct_initial_points; 

/// @} 

/// \name Operations 
/// The following functions give access to the construction objects:
/// @{

/*!

*/ 
Intersect_3 intersect_3_object(); 

/*!

*/ 
Construct_initial_points construct_initial_points_object(); 

/// @}

}; /* end SurfaceMeshTraits_3 */

