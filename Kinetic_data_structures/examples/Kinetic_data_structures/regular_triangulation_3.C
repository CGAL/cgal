#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits_3.h>
#include <CGAL/Kinetic/Regular_triangulation_3.h>

int main(int, char *[]) {
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits_3 Traits;
    typedef CGAL::Kinetic::Regular_triangulation_3<Traits> KDel;

    Traits tr;
    KDel kdel(tr);

    Traits::Simulator::Handle sp= tr.simulator_handle();
  
    std::ifstream in("data/weighted_points_3");
    in >> *tr.active_points_3_table_handle();
    
    std::cout <<  *tr.active_points_3_table_handle() << std::endl;

    kdel.set_has_certificates(true);
   
    sp->set_current_event_number(10000);
    return EXIT_SUCCESS;
};
