#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Exact_simulation_traits_1.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Sort.h>

int main(int, char *[])
{

    typedef CGAL::Kinetic::Exact_simulation_traits_1 Traits;

    typedef CGAL::Kinetic::Insert_event<Traits::Active_points_1_table> Insert_event;
    typedef Traits::Active_points_1_table::Data Moving_point;
    typedef CGAL::Kinetic::Sort<Traits> Sort;
    typedef Traits::Simulator::Time Time;

    Traits tr;
    Sort sort(tr);
    Traits::Simulator::Handle sp= tr.simulator_handle();

    std::ifstream in("data/points_1");
    in  >> *tr.active_points_1_table_handle();

    while (sp->next_event_time() != sp->end_time()) {
        sp->set_current_event_number(sp->current_event_number()+1);
    }

    return EXIT_SUCCESS;
};
