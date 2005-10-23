#define NDEBUG
#include <CGAL/KDS/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <CGAL/KDS/Cartesian_kinetic_kernel.h>
#include <CGAL/KDS/Notifying_table.h>
#include <CGAL/KDS/Simulation_traits.h>
#include <CGAL/KDS/Simulator.h>
#include <CGAL/KDS/Two_list_pointer_event_queue.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Simple_cartesian.h>


template <bool Skip>
struct Exact_simulation_types {
  typedef CGAL::Simple_cartesian<CGAL::Gmpq> Static_kernel;
  typedef Static_kernel::FT NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack_Descartes_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Upper_bound_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  typedef CGAL::KDS::Skip_even_roots_function_kernel<Function_kernel, Skip> Simulator_function_kernel;
  typedef CGAL::KDS::Cartesian_kinetic_kernel<Function_kernel> Kinetic_kernel;

  
  typedef  typename Simulator_function_kernel::Root Time;

  typedef CGAL::KDS::Two_list_pointer_event_queue<Time, double> Queue_base;
  struct Event_queue: public Queue_base{
    Event_queue(const Time &start): Queue_base(start){}
  };
  
  typedef CGAL::KDS::Simulator<Simulator_function_kernel, Event_queue > Simulator;
  typedef CGAL::KDS::Notifying_table<Kinetic_kernel::Point_1> Moving_point_table;
  typedef CGAL::KDS::Cartesian_instantaneous_kernel<Moving_point_table,
						    Static_kernel> Instantaneous_kernel;



};

template <bool Skip>
struct Exact_simulation_traits: 
  public CGAL::KDS::Simulation_traits<typename Exact_simulation_types<Skip>::Static_kernel,
			   typename Exact_simulation_types<Skip>::Kinetic_kernel,
			   typename Exact_simulation_types<Skip>::Simulator> {
private:
  typedef  CGAL::KDS::Simulation_traits<typename Exact_simulation_types<Skip>::Static_kernel,
			     typename Exact_simulation_types<Skip>::Kinetic_kernel,
			     typename Exact_simulation_types<Skip>::Simulator> P;
public:
 
  typedef typename Exact_simulation_types<Skip>::Instantaneous_kernel Instantaneous_kernel;
  typedef typename Exact_simulation_types<Skip>::Moving_point_table Moving_point_table;
 

  Exact_simulation_traits(): mp2_(new Moving_point_table()){

  }
  Moving_point_table* moving_point_table_pointer(){ return mp2_.get();}
  const Moving_point_table* moving_point_table_pointer() const { return mp2_.get();}
 
  Instantaneous_kernel instantaneous_kernel_object() const {
    return Instantaneous_kernel(mp2_,P::static_kernel_object());
  }
protected:
  typename Moving_point_table::Pointer mp2_;
};



template <class Traits, class Fn, class Rt>
void check_one(const Traits &tr, const Fn &fn, const Rt &lb, const Rt* rt) {
  typename Traits::Simulator::Root_stack rs(fn, lb, std::numeric_limits<Rt>::infinity(), 
					    tr.simulator_pointer()->function_kernel_object());
  
  while (*rt != std::numeric_limits<Rt>::infinity()) {
    if (rs.top() != *rt) {
      std::cerr << "ERROR For function " << fn << " expected " << *rt << " got " << rs.top() << std::endl;
    }
    ++rt;
    rs.pop();
  }
  if (!rs.empty()){
    std::cerr << "ERROR For function " << fn << " expected " << *rt << " got " << rs.top() << std::endl;
  }
}


//check the KDS solvers

int main(int, char *[]){


  {
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
      Root rts[]={Root(4), inf};
      Traits::Simulator::Function_kernel::Function f=  -cf(-1,1)*cf(-1,1)*cf(-2,1)*cf(-2,1)*cf(-4,1);
      check_one(tr,f , zero, rts);
    }
    {
      Root rts[]={Root(2), Root(4), inf};
      Traits::Simulator::Function_kernel::Function f=  cf(-1,1)*cf(-1,1)*cf(-2,1)*cf(-3,1)
	*cf(-3,1)*cf(-4,1);
      check_one(tr,f , zero, rts);
    }

    {
      Root rts[]={Root(0), Root(1), inf};
      Traits::Simulator::Function_kernel::Function f= -cf(0,1)*cf(1,-1);
      check_one(tr,f , zero, rts);
    }
  }

  {
    typedef Exact_simulation_traits<false> Traits;
    Traits tr;
    typedef Traits::Simulator::Root_stack::Root Root;
    Root inf= std::numeric_limits<Root>::infinity();
    Traits::Simulator::Function_kernel::Construct_function 
      cf= tr.function_kernel_object().construct_function_object();
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
  }

  {
    
  }

  return EXIT_SUCCESS;
}
