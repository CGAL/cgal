#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/Exact_simulation_traits_1.h>
#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Sort.h>

int main(int, char *[])
{

    typedef CGAL::Kinetic::Exact_simulation_traits_1 Traits;

    typedef CGAL::Kinetic::Insert_event<Traits::Active_objects_table> Insert_event;
    typedef Traits::Active_objects_table::Data Moving_point;
    typedef CGAL::Kinetic::Sort<Traits> Sort;
    typedef Traits::Simulator::Time Time;

    Traits tr;
    Sort sort(tr);
    Traits::Simulator::Pointer sp= tr.simulator_pointer();

    std::ifstream in("data/points_1");
    in  >> *tr.active_objects_table_pointer();
   

    while (sp->next_event_time() != sp->end_time()) {
        sp->set_current_event_number(sp->current_event_number()+1);
    }

   

    return EXIT_SUCCESS;
};
