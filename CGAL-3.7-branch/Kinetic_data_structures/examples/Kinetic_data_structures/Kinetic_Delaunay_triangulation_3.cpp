#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Delaunay_triangulation_event_log_visitor_3.h>

int main()
{

    typedef CGAL::Kinetic::Exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Delaunay_triangulation_3<Traits,
    CGAL::Kinetic::Delaunay_triangulation_event_log_visitor_3> KDel;

    Traits tr(0,10000);
    KDel kdel(tr);

    Traits::Simulator::Handle sp= tr.simulator_handle();

    std::ifstream in("data/points_3");
    in >> *tr.active_points_3_table_handle();

    kdel.set_has_certificates(true);

    sp->set_current_time(sp->end_time());

    std::cout << "Processed " << sp->current_event_number() << " events.\n";

    /*std::copy(kdel.visitor().events_begin(), kdel.visitor().events_end(),
      std::ostream_iterator<std::string>(std::cout, "\n"));*/

    return EXIT_SUCCESS;
}
