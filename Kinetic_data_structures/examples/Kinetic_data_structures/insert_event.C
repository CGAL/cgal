#include <CGAL/KDS/Delaunay_triangulation_2.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Exact_simulation_traits_2.h>

int main(int, char *[])
{

    typedef CGAL::KDS::Exact_simulation_traits_2 Simulation_traits;
    typedef Simulation_traits::Kinetic_kernel::Point_2 Moving_point_2;
    typedef CGAL::KDS::Insert_event<Simulation_traits::Active_objects_table> Insert_event;
    typedef CGAL::KDS::Delaunay_triangulation_2<Simulation_traits> KDel;

    Simulation_traits tr;
    Simulation_traits::Simulator::Pointer sp= tr.simulator_pointer();

    KDel kdel(tr);
  
    std::ifstream in("data/points_2");
    in >> *tr.active_objects_table_pointer();

    sp->new_event(Simulation_traits::Simulator::Time(3), 
		  Insert_event(Moving_point_2(Moving_point_2::NT(0),
					      Moving_point_2::NT(0)),
			       tr.active_objects_table_pointer()));

    while (sp->next_event_time()
    < sp->end_time()) {
        sp->set_current_event_number(sp->current_event_number()+10);
        std::cout << "At time " << sp->current_time() << ":\n";
        std::cout << kdel.triangulation_data_structure();
    }

    return EXIT_SUCCESS;
};
