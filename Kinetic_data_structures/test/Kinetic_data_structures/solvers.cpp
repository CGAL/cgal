//#define NDEBUG
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Default_instantaneous_kernel.h>
#include <CGAL/Kinetic/Cartesian.h>
#include <CGAL/Kinetic/Active_objects_vector.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Simple_cartesian.h>

/*template <bool Skip>
struct Sest_types
{
  typedef CGAL::Simple_cartesian<CGAL::Kinetic::Default_field_nt> Static_kernel;
  typedef Static_kernel::FT NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack_Descartes_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  typedef CGAL::Kinetic::Handle_degeneracy_function_kernel<Function_kernel, false> Simulator_function_kernel;
  typedef CGAL::Kinetic::Cartesian_kinetic_kernel<Simulator_function_kernel> Kinetic_kernel;
  typedef typename Simulator_function_kernel::Root Time;
  typedef CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel> Queue_base;

  struct Event_queue: public Queue_base
  {
    Event_queue(const Time &start, const Time &end, Function_kernel fk, int num): Queue_base(start, end, fk, num){}
  };

  typedef CGAL::Kinetic::Default_simulator<Simulator_function_kernel, Event_queue > Simulator;
  typedef CGAL::Kinetic::Active_objects_vector<typename Kinetic_kernel::Point_1> Active_objects_table;
  typedef CGAL::Kinetic::Cartesian_instantaneous_kernel<Active_objects_table,
                                                        Static_kernel> Instantaneous_kernel;

							};*/

 /*template <bool Skip>
struct Exact_simulation_traits:
  public CGAL::Kinetic::Simulation_traits<typename Sest_types<Skip>::Static_kernel,
					  typename Sest_types<Skip>::Kinetic_kernel,
					  typename Sest_types<Skip>::Simulator >
{
};*/

template <class Traits, class Fn, class Rt>
void check_one(const Traits &tr, const Fn &fn, const Rt &lb, const Rt* rt)
{
  typename Traits::Kinetic_kernel::Function_kernel::Root_stack rs(fn, lb, std::numeric_limits<Rt>::infinity(),
								  tr.kinetic_kernel_object().function_kernel_object());

  while (*rt != std::numeric_limits<Rt>::infinity()) {
    if (rs.top() != *rt) {
      std::cerr << "ERROR For function " << fn << " expected " << *rt << " got " << rs.top() << std::endl;
    }
    ++rt;
    rs.pop();
  }
  if (!rs.empty()) {
    std::cerr << "ERROR For function " << fn << " expected " << *rt << " got " << rs.top() << std::endl;
  }
}


//check the KDS solvers

int main(int, char *[])
{

  /*{
    typedef Exact_simulation_traits<true> Traits;
    Traits tr;
    typedef Traits::Simulator::Root_stack::Root Root;
    Root inf= std::numeric_limits<Root>::infinity();
    Root zero= Root(0);
    Traits::Simulator::Function_kernel::Construct_function
    cf= tr.function_kernel_object().construct_function_object();

    {
    check_one(tr, cf(1,-2,1), zero, &inf);
    }

    {
    Root rts[]={Root(1), Root(1), Root(2), Root(2), Root(4), inf};
    Traits::Simulator::Function_kernel::Function f=  -cf(-1,1)*cf(-1,1)*cf(-2,1)*cf(-2,1)*cf(-4,1);
    check_one(tr,f , zero, rts);
    }
    {
    Root rts[]={Root(1), Root(1), Root(2),Root(3), Root(3), Root(4), inf};
    Traits::Simulator::Function_kernel::Function f=  cf(-1,1)*cf(-1,1)*cf(-2,1)*cf(-3,1)
    *cf(-3,1)*cf(-4,1);
    check_one(tr,f , zero, rts);
    }

    {
    Root rts[]={Root(0), Root(1), inf};
    Traits::Simulator::Function_kernel::Function f= -cf(0,1)*cf(1,-1);
    check_one(tr,f , zero, rts);
    }
    }*/

  /* {
    typedef Exact_simulation_traits<false> Traits;
    Traits tr;
    typedef Traits::Simulator::Time Root;
    Root inf= std::numeric_limits<Root>::infinity();
    Traits::Kinetic_kernel::Function_kernel::Construct_function
      cf= tr.kinetic_kernel_object().function_kernel_object().construct_function_object();
    Root zero(0);
    {
      Root rts[3]={Root(1), Root(1), inf};
      Traits::Simulator::Function_kernel::Function f=cf(1,-2,1);
      check_one(tr, f,zero, rts);
    }
    {
      Root rts[]={Root(1),Root(1),Root(2), Root(2), Root(4), inf};
      Traits::Simulator::Function_kernel::Function f=  -cf(-1,1)*cf(-1,1)*cf(-2,1)*cf(-2,1)*cf(-4,1);
      check_one(tr,f , zero, rts);
    }
    {
      Root rts[]={Root(1),Root(1),Root(2), Root(3), Root(3),Root(4), inf};
      Traits::Simulator::Function_kernel::Function f=  cf(-1,1)*cf(-1,1)*cf(-2,1)*cf(-3,1)
	*cf(-3,1)*cf(-4,1);
      check_one(tr,f , zero, rts);
    }
    {
      Root rts[]={Root(0), Root(1), inf};
      Traits::Simulator::Function_kernel::Function f= -cf(0,1)*cf(1,-1);
      check_one(tr,f , zero, rts);
    }
    }*/

  if (CGAL::Kinetic::internal::get_static_audit_failures() != 0) return EXIT_FAILURE;
  else return EXIT_SUCCESS;
}
