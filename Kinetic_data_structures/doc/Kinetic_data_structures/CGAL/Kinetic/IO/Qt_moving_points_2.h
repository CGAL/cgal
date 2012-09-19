
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsSupport

The class `Kinetic::Qt_moving_points_2` displays a set of moving points in 2D. 

See Section \ref seckds_delaunay_2_example for an example using this class. 

\sa `Kinetic::Qt_widget_2<Simulator>`

*/
template< typename Traits, typename QtWidget_2 >
class Qt_moving_points_2 {
public:

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Qt_moving_points_2(QtGui::Handle,Traits::Active_points_2_table::Handle); 

/// @}

}; /* end Kinetic::Qt_moving_points_2 */
} /* end namespace Kinetic */
} /* end namespace CGAL */
