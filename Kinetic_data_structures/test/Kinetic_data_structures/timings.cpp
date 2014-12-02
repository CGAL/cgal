#include <CGAL/config.h>
#define NDEBUG
// NDEBUG is defined after <CGAL/config.h> is included, to workaround the
// check of the testsuite that NDEBUG is not defined.

//#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS

#include <CGAL/Polynomial/internal/Isolating_interval.h>
#include <CGAL/Polynomial/internal/Simple_interval_root.h>
//#include <CGAL/Polynomial/Lazy_upper_bound_root_stack.h>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>
#include <CGAL/Kinetic/Sort.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

//template <class AO>
/*struct Slest_types
{
typedef CGAL::Simple_cartesian < CGAL::Kinetic::Default_field_nt > Static_kernel;
    typedef Static_kernel::FT NT;
    typedef CGAL::POLYNOMIAL::Polynomial < NT > Function;
    typedef CGAL::POLYNOMIAL::Upper_bound_root_stack_Descartes_traits <
        Function > Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Lazy_upper_bound_root_stack < Root_stack_traits >
        Root_stack;
    typedef CGAL::POLYNOMIAL::Kernel < Function, Root_stack > Function_kernel;
    typedef CGAL::Kinetic::Handle_degeneracy_function_kernel < Function_kernel >
        Simulator_function_kernel;
    typedef CGAL::Kinetic::Cartesian_kinetic_kernel < Function_kernel >
        Kinetic_kernel;
    typedef Simulator_function_kernel::Root Time;
    typedef CGAL::Kinetic::Two_list_pointer_event_queue < Time, double >Queue_base;

    struct Event_queue:public Queue_base
    {
        Event_queue(const Time & start, const Time & end):Queue_base(start, end) {
        }
    };

    typedef CGAL::Kinetic::Simulator < Simulator_function_kernel,
        Event_queue > Simulator;


};*/

 /*struct Lazy_exact_traits_1:public CGAL::Kinetic::Simulation_traits <
Slest_types::Static_kernel,
CGAL::Kinetic::Cartesian_instantaneous_kernel <
CGAL::Kinetic::Active_objects_vector < Slest_types::Kinetic_kernel::Point_1 >,
Slest_types::Static_kernel >, Slest_types::Kinetic_kernel,
Slest_types::Simulator,
CGAL::Kinetic::Active_objects_vector < Slest_types::Kinetic_kernel::Point_1 >
>
{
    typedef CGAL::Kinetic::Simulation_traits < Slest_types::Static_kernel,
        CGAL::Kinetic::Cartesian_instantaneous_kernel <
        CGAL::Kinetic::Active_objects_vector <
        Slest_types::Kinetic_kernel::Point_1 >, Slest_types::Static_kernel >,
        Slest_types::Kinetic_kernel, Slest_types::Simulator,
        CGAL::Kinetic::Active_objects_vector <
        Slest_types::Kinetic_kernel::Point_1 > >P;
    Lazy_exact_traits_1(const P::NT & lb, const P::NT & ub):P(lb, ub) {
    }
};
struct Lazy_exact_traits_3:public CGAL::Kinetic::Simulation_traits <
Slest_types::Static_kernel,
CGAL::Kinetic::Cartesian_instantaneous_kernel <
CGAL::Kinetic::Active_objects_vector < Slest_types::Kinetic_kernel::Point_3 >,
Slest_types::Static_kernel >, Slest_types::Kinetic_kernel,
Slest_types::Simulator,
CGAL::Kinetic::Active_objects_vector < Slest_types::Kinetic_kernel::Point_3 >
>
{
    typedef CGAL::Kinetic::Simulation_traits < Slest_types::Static_kernel,
        CGAL::Kinetic::Cartesian_instantaneous_kernel <
        CGAL::Kinetic::Active_objects_vector <
        Slest_types::Kinetic_kernel::Point_3 >, Slest_types::Static_kernel >,
        Slest_types::Kinetic_kernel, Slest_types::Simulator,
        CGAL::Kinetic::Active_objects_vector <
        Slest_types::Kinetic_kernel::Point_3 > >P;
    Lazy_exact_traits_3(const P::NT & lb, const P::NT & ub):P(lb, ub) {
    }
};*/

template < class Traits > double test_sort(unsigned int degree, unsigned int n)
{
    typedef CGAL::Kinetic::Sort < Traits > Sort;
    Traits tr(0, 10000);
    Sort sort(tr);
    CGAL::Random r;
    for (unsigned int i = 0; i < n; ++i) {
        std::vector < double >cf;
        for (unsigned int j = 0; j < degree + 1; ++j) {
            cf.push_back(r.get_double());
        }
        typename Traits::Kinetic_kernel::Motion_function fn(cf.begin(),
            cf.end());
        typename Traits::Kinetic_kernel::Point_1 pt(fn);
        tr.active_points_1_table_handle()->insert(pt);
    }
    CGAL::Timer timer;
    timer.start();
    int ne = 0;
    while (tr.simulator_handle()->next_event_time() !=
    tr.simulator_handle()->end_time()) {
        tr.simulator_handle()->set_current_event_number(tr.
            simulator_handle()->
            current_event_number()
            + 1);
        ++ne;
        if (ne == 1000)
            break;
    }
    timer.stop();
    return timer.time() / static_cast < double >(ne);
}


template < class Traits > void test_sort(const char *nm)
{
  unsigned int beg= 4, end=5;
    std::cout << "Solver: " << nm << std::endl;
    for (unsigned int i = beg; i < end; ++i) {
      std::printf("%6f\t", test_sort < Traits > (i, static_cast<int>(std::ceil(500.0/i))));
        std::cout << std::flush;
        if (i > 4)
            ++i;
    }
    std::cout << std::endl;
}


template < class Traits > double test_del(unsigned int degree, unsigned int n)
{
    typedef CGAL::Kinetic::Delaunay_triangulation_3 < Traits > Del;
    Traits tr;
    Del del(tr);
    CGAL::Random r;
    for (unsigned int i = 0; i < n; ++i) {
        std::vector < double >cf[3];
        for (unsigned int j = 0; j < degree + 1; ++j) {
            for (int k = 0; k < 3; ++k) {
                cf[k].push_back(r.get_double());
            }
        }
        typename Traits::Kinetic_kernel::Motion_function fn[3];
        for (unsigned int k = 0; k < 3; ++k)
            fn[k] =
                typename Traits::Kinetic_kernel::Motion_function(cf[k].begin(),
                cf[k].end());
        typename Traits::Kinetic_kernel::Point_3 pt(fn[0], fn[1], fn[2]);
        tr.active_objects_table_pointer()->insert(pt);
    }
    del.set_has_certificates(true);
    CGAL::Timer timer;
    timer.start();
    int ne = 0;
    while (tr.simulator_pointer()->next_event_time() !=
    tr.simulator_pointer()->end_time()) {
        tr.simulator_pointer()->set_current_event_number(tr.
            simulator_pointer()->
            current_event_number()
            + 1);
        ++ne;
        if (ne == 1000)
            break;
    }
    timer.stop();
    return timer.time() / static_cast < double >(ne);
}


template < class Traits > void test_del(const char *nm)
{
    std::cout << "Solver: " << nm << std::endl;
    for (unsigned int i = 1; i < 9; ++i) {
      printf("%6f\t", test_del < Traits > (i,static_cast<int>( std::ceil(20.0/std::sqrt(static_cast<double>(i))))));
        std::cout << std::flush;
        if (i > 4)
            ++i;
    }
    std::cout << std::endl;
}


int main(/* int argc, char *argv[] */)
{
  //std::cout << "Delaunay\n";
    //test_del < CGAL::Kinetic::Inexact_simulation_traits_3 > ("Numeric");
    std::cout << "Sort\n";
    test_sort < CGAL::Kinetic::Exact_simulation_traits > ("Upper bound");
   
    if (CGAL::Kinetic::internal::get_static_audit_failures() != 0) return EXIT_FAILURE;
    else return EXIT_SUCCESS;

    //test_sort < CGAL::Kinetic::Exact_simulation_traits_1 > ("Upper bound");
    //test_sort < CGAL::Kinetic::Inexact_simulation_traits_1 > ("Numeric");
    //test_sort<Lazy_exact_traits_1>("Lazy upper bound");
    /*std::cout << CGAL::POLYNOMIAL::internal::lazy_stats.created_ << " "
        << CGAL::POLYNOMIAL::internal::lazy_stats.isolated_ << " "
        << CGAL::POLYNOMIAL::internal::lazy_stats.refine_attempted_ << " "
        << CGAL::POLYNOMIAL::internal::lazy_stats.
        refine_succeeded_ << std::endl;*/
    
/*test_del<Lazy_exact_traits_3>("Lazy upper bound");
   std::cout << CGAL::POLYNOMIAL::lazy_stats.created_ << " "
   << CGAL::POLYNOMIAL::lazy_stats.isolated_ << " "
   << CGAL::POLYNOMIAL::lazy_stats.refine_attempted_ << " "
   << CGAL::POLYNOMIAL::lazy_stats.refine_succeeded_ << std::endl; */
}
