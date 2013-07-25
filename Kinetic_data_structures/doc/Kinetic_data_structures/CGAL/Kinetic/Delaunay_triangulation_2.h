
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsTri2

The class `Kinetic::Delaunay_triangulation_2` maintains a Delaunay triangulation on top of the 
points contained in a `Kinetic::ActiveObjectsTable`. It has one 
main method of interest, `triangulation()`, which returns the 
triangulation it is maintaining. 

Note that the Delaunay triangulation is fully dynamic as it tracks points added to and removed from the `Kinetic::ActiveObjectsTable`. 

The class `Kinetic::Qt_triangulation_2<KineticTriangulation_2, QtWidget_2, QtMovingPoints_2>` displays a kinetic Delaunay 
triangulation using the Qt widget. 

This class is a good example of a simple, but non-trivial, kinetic 
data structure. 

The `Triangulation` template parameter must be a model of 
`CGAL::Delaunay_triangulation_2<Traits, Tds>` which uses 
`Traits::Default_instantaneous_kernel` as its geometric traits and a 
`Tds` whose face inherits from 
`Kinetic::Delaunay_triangulation_face_base_2<Traits, Base>`. 

The optional `Visitor` parameter takes a model of 
`Kinetic::DelaunayTriangulationVisitor_2`. Methods on this object will be called 
whenever the triangulation changes. 

\cgalModels `Ref_counted<T>`

\sa `Kinetic::DelaunayTriangulationVisitor_2`
\sa `Kinetic::Delaunay_triangulation_default_visitor_2`
\sa `Kinetic::Delaunay_triangulation_recent_edges_visitor_2<Triangulation>`
\sa `Kinetic::Delaunay_triangulation_event_log_visitor_2`
\sa `Kinetic::Qt_triangulation_2`

*/
template< typename Traits, typename Visitor, typename Triangulation >
class Delaunay_triangulation_2 {
public:

/// \name Types 
/// @{

/*!
The template argument triangulation. 
*/ 
typedef unspecified_type Triangulation; 

/*!
The template argument for the visitor. 
*/ 
typedef unspecified_type Visitor; 

/// @} 

/// \name Creation 
/// @{

/*!
Maintain the 
Delaunay triangulation of the points in 
`tr.active_points_2_handle()`. 
*/ 
Delaunay_triangulation_2(Traits tr); 

/// @} 

/// \name Operations 
/// @{

/*!
Access the triangulation that is maintained. 
*/ 
const Triangulation& triangulation() const; 

/*!
Access the visitor. 
*/ 
Visitor& visitor(); 

/*!
Insert the point. 
*/ 
Vertex_handle insert(Point_key k); 

/*!
Erase the vertex. 
*/ 
void erase(Vertex_handle h); 

/// @}

}; /* end Kinetic::Delaunay_triangulation_2 */
} /* end namespace Kinetic */
} /* end namespace CGAL */
