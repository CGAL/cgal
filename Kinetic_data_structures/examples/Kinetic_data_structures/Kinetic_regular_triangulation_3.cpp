#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Kinetic/Regular_triangulation_exact_simulation_traits.h>
#include <CGAL/Kinetic/Regular_triangulation_3.h>

int main()
{
    typedef CGAL::Kinetic::Regular_triangulation_exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Regular_triangulation_3<Traits> KDel;

    Traits tr(0,100000.0);
    KDel kdel(tr);

    Traits::Simulator::Handle sp= tr.simulator_handle();

    std::ifstream in("data/weighted_points_3");
    CGAL_assertion(in.good());
    in >> *tr.active_points_3_table_handle();
    CGAL_assertion(!in.fail());

    std::cout << *tr.active_points_3_table_handle();

    std::cout <<  *tr.active_points_3_table_handle() << std::endl;

    kdel.set_has_certificates(true);

    sp->set_current_event_number(10000);
    return EXIT_SUCCESS;
}
