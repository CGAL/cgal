#include <CGAL/KDS/Regular_triangulation_exact_simulation_traits_3.h>
#include <CGAL/KDS/Regular_triangulation_3.h>

int main(int, char *[])
{

    typedef CGAL::KDS::Regular_triangulation_exact_simulation_traits_3 Traits;

    typedef CGAL::KDS::Regular_triangulation_3<Traits> KDel;

    Traits tr;
    KDel kdel(tr);

    Traits::Simulator::Pointer sp= tr.simulator_pointer();

    Traits::Function_kernel::Construct_function cf
        = tr.function_kernel_object().construct_function_object();
   
  
    std::ifstream in("data/weighted_points_3");
    in >> *tr.active_objects_table_pointer();
    
    std::cout <<  *tr.active_objects_table_pointer() << std::endl;

    kdel.set_has_certificates(true);
   
    sp->set_current_event_number(10000);
    return EXIT_SUCCESS;
};
