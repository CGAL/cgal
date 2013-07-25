
namespace CGAL {
namespace Kinetic {

/*!
\ingroup PkgKdsSupport

This event erases a point from the `ActiveObjectsTable` when it 
is processed. 

\cgalModels `Kinetic::Simulator::Event`

\sa `Kinetic::ActiveObjectsTable`
\sa `Kinetic::Active_objects_vector<MovingObject>`

\cgalHeading{Example}

\code{.cpp} 

typedef CGAL::Kinetic::Exact_simulation_traits Simulation_traits; 
typedef Simulation_traits::Kinetic_kernel::Point_2 Moving_point_2; 
typedef CGAL::Kinetic::Insert_event<Simulation_traits::Active_points_2_table> Insert_event; 
typedef CGAL::Kinetic::Delaunay_triangulation_2<Simulation_traits> KDel; 

Simulation_traits tr; 

KDel kdel(tr); 

Moving_point_2 mp(Moving_point_2::NT(0), 
                  Moving_point_2::NT(0)); 
tr.simulator_handle()->new_event(Simulation_traits::Simulator::Time(3), 
                                 Erase_event(*tr.active_objects_table_handle()->keys_begin(), 
                                             tr.active_points_2_table_handle())); 

\endcode 

*/
template< typename ActiveObjectsTable >
class Erase_event {
public:

/// \name Creation 
/// @{

/*!
Erase the object `k` from the table `t` when processed. 
*/ 
Insert_event(ActiveObjectsTable::Key k, 
ActiveObjectsTable::Handle t); 

/// @}

}; /* end Kinetic::Erase_event */
} /* end namespace Kinetic */
} /* end namespace CGAL */
