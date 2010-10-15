#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>

int main()
{

    typedef CGAL::Kinetic::Exact_simulation_traits Simulation_traits;
    typedef Simulation_traits::Kinetic_kernel::Point_1 Moving_point_1;
    typedef CGAL::Kinetic::Insert_event<Simulation_traits::Active_points_1_table> Insert_event;
    typedef CGAL::Kinetic::Sort<Simulation_traits> KDS;

    Simulation_traits tr(0, 10000.0);
    Simulation_traits::Simulator::Handle sp= tr.simulator_handle();

    KDS kds(tr);


    //CGAL_SET_LOG_LEVEL(CGAL::Kinetic::Log::LOTS);
    std::ifstream in("data/points_1");
    in >> *tr.active_points_1_table_handle();
    std::cout << *tr.active_points_1_table_handle();
    sp->new_event(Simulation_traits::Simulator::Time(3),
		  Insert_event(Moving_point_1(Moving_point_1::Coordinate(0)),
			       tr.active_points_1_table_handle()));

    while (sp->next_event_time()
    < sp->end_time()) {
        sp->set_current_event_number(sp->current_event_number()+10);
        std::cout << "At time " << sp->current_time() << ":\n";
        std::cout << kds;
    }

    return EXIT_SUCCESS;
}
