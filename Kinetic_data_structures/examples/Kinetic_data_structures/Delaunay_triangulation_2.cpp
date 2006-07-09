#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Exact_simulation_traits_2.h>

int main(int, char *[])
{

    typedef CGAL::Kinetic::Exact_simulation_traits_2 Simulation_traits;
    typedef Simulation_traits::Kinetic_kernel::Point_2 Moving_point_2;

    typedef CGAL::Kinetic::Delaunay_triangulation_2<Simulation_traits> KDel;

    Simulation_traits tr;
    Simulation_traits::Simulator::Handle sp= tr.simulator_handle();

    KDel kdel(tr);
  
    std::ifstream in("data/points_2");
    in >> *tr.active_points_2_table_handle();

    while (sp->next_event_time()
    < sp->end_time()) {
        sp->set_current_event_number(sp->current_event_number()+10);
        std::cout << "At time " << sp->current_time() << ":\n";
        std::cout << kdel.triangulation_data_structure();
    }

    return EXIT_SUCCESS;
};
