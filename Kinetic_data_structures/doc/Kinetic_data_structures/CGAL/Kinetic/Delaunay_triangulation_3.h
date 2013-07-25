
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsTri3

The class `Kinetic::Delaunay_triangulation_3` maintains a Delaunay triangulation on top of the 
points contained in a `Kinetic::ActiveObjectsTable`. It has one main method 
of interest. `triangulation()` which returns the triangulation it 
is maintaining. In addition, as an optimisation, you can turn on and 
off whether it is currently maintaining its certificates. This allows 
a large number of changes to the underlying points to be made at one 
time without recomputing the certificates each time a single point 
changes. 

Note that the Delaunay triangulation is fully dynamic as it tracks points added to and removed from the `Kinetic::ActiveObjectsTable`. 

The class `Kinetic::Qt_triangulation_3<Traits>`, included as part 
of the demo code, displays a kinetic Delaunay triangulation in three 
dimensions using the Coin library. 

The optional `Visitor` template argument is a model of 
`Kinetic::DelaunayTriangulationVisitor_3` and can be used to monitor 
changes in the kinetic data structure. 

The optional `Triangulation` template argument must be a model of 
a `CGAL::DelaunayTriangulation_3` which uses 
`Traits::Default_instantaneous_kernel` as its geometric traits and has 
`Kinetic::Delaunay_triangulation_cell_base_3<Traits, Base>` a the cell 
type. 

\sa `Kinetic::Regular_triangulation_3<Traits, Triangulation, Visitor>`
\sa `Kinetic::Delaunay_triangulation_2<Traits, Triangulation, Visitor`
\sa `Kinetic::Delaunay_triangulation_visitor_base_3`
\sa `Kinetic::Delaunay_triangulation_event_log_visitor_3`

*/
template< typename Traits, typename Visitor, typename Triangulation >
class Delaunay_triangulation_3 {
public:

/// \name Types 
/// @{

/*!
The template argument. 
*/ 
typedef unspecified_type Triangulation; 

/*!
The template argument. 
*/ 
typedef unspecified_type Visitor; 

/// @} 

/// \name Creation 
/// @{

/*!
Maintain the Delaunay triangulation of the points in `tr.active_points_3_handle()`. 
*/ 
Delaunay_triangulation_3(Traits tr); 

/// @} 

/// \name Operations 
/// @{

/*!
Access the triangulation that is maintained. 
*/ 
const Triangulation* triangulation(); 

/*!
This method returns true if the `Kinetic::Delaunay_triangulation_3` is currently maintaining certificates for a Delaunay triangulation. 
*/ 
bool has_certificates(); 

/*!
This method allows you to control whether the triangulation is maintaining certificates. 
*/ 
void set_has_certificates(bool tf); 

/*!
Access the visitor. 
*/ 
Visitor& visitor(); 

/// @}

}; /* end Kinetic::Delaunay_triangulation_3 */
} /* end namespace Kinetic */
} /* end namespace CGAL */
