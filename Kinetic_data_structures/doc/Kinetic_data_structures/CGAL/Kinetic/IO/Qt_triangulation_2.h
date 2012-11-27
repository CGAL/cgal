
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsSupport

The class draws a triangulation into a `CGAL::Qt_widget_2`. This class is 
very simple and a good one to look at if you want to see how to draw 
your own two dimensional kinetic data structure. 

See Section \ref seckds_delaunay_2_example for an example using this class. 

\sa `Kinetic::Qt_widget_2<Simulator>` 

*/
template< typename KineticTriangulation_2, typename QtWidget_2, typename QtMovingPoints_2 >
class Qt_triangulation_2 {
public:

/// \name Creation 
/// @{

/*! 
Construct the object and make all the connections with the appropriate other objects. 
*/ 
Qt_widget_2(KineticTriangulation_2::Handle,QtWidget_2::Handle,QtMovingPoints_2::Handle); 

/// @}

}; /* end Kinetic::Qt_triangulation_2 */
} /* end namespace Kinetic */
} /* end namespace CGAL */
