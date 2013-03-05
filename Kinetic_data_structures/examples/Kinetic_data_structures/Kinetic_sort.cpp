#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Sort.h>

int main()
{
    typedef CGAL::Kinetic::Exact_simulation_traits Traits;

    typedef CGAL::Kinetic::Sort<Traits> Sort;

    Traits tr(0,100000);
    Sort sort(tr);
    Traits::Simulator::Handle sp= tr.simulator_handle();

    std::ifstream in("data/points_1");
    in  >> *tr.active_points_1_table_handle();

    while (sp->next_event_time() != sp->end_time()) {
        sp->set_current_event_number(sp->current_event_number()+1);
    }

    return EXIT_SUCCESS;
}
