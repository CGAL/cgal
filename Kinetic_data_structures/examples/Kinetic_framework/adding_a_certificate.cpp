#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Polynomial/Sturm_root_stack_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
#include <CGAL/Kinetic/Active_objects_vector.h>
#include <CGAL/Kinetic/Default_instantaneous_kernel.h>
#include <CGAL/Kinetic/Cartesian.h>
#include <CGAL/Kinetic/Handle_degeneracy_function_kernel.h>
#include <CGAL/Kinetic/Default_simulator.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Certificate_generator.h>


template <class KineticKernel>
struct Positive_x_f_2 {
  typedef typename KineticKernel::Certificate_function result_type;
  typedef typename KineticKernel::Point_2 argument_type;
  result_type operator()(const argument_type &p){
    return result_type(p.x());
  }
};

template <class FunctionKernel>
class My_kinetic_kernel:
  public CGAL::Kinetic::Cartesian<FunctionKernel> {
  typedef CGAL::Kinetic::Cartesian<FunctionKernel> P;
  typedef My_kinetic_kernel<FunctionKernel> This;
public:
  typedef CGAL::Kinetic::Certificate_generator<This, Positive_x_f_2<This> > Positive_x_2;
  Positive_x_2 positive_x_2_object() const
  {
    return Positive_x_2(P::function_kernel_object());
  }
};


//#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

struct My_simulation_traits {
  typedef My_simulation_traits This;

  typedef CGAL::Exact_predicates_exact_constructions_kernel Static_kernel;
  //typedef CGAL::Regular_triangulation_euclidean_traits_3<Static_kernel_base> Static_kernel;
  typedef CGAL::POLYNOMIAL::Polynomial<Static_kernel::FT> Function;
  typedef CGAL::POLYNOMIAL::Sturm_root_stack_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Sturm_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;

  typedef CGAL::Kinetic::Handle_degeneracy_function_kernel<Function_kernel, false>  Simulator_function_kernel_base;
  struct Simulator_function_kernel: public Simulator_function_kernel_base{};

  typedef My_kinetic_kernel<Simulator_function_kernel> Kinetic_kernel;
  typedef CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel> Event_queue;
  typedef CGAL::Kinetic::Default_simulator<Simulator_function_kernel, Event_queue > Simulator;

  typedef CGAL::Kinetic::Active_objects_vector<Kinetic_kernel::Point_1> Active_points_1_table;
  typedef CGAL::Kinetic::Active_objects_vector<Kinetic_kernel::Point_2> Active_points_2_table;
  typedef CGAL::Kinetic::Active_objects_vector<Kinetic_kernel::Point_3> Active_points_3_table;
  // typedef Active_objects_vector<Kinetic_kernel::Weighted_point_3> Active_weighted_points_3_table;

  typedef CGAL::Kinetic::Default_instantaneous_kernel<This> Instantaneous_kernel;

  Active_points_1_table* active_points_1_table_handle() const { return ap1_.get();}
  Active_points_2_table* active_points_2_table_handle() const {return ap2_.get();}
  Active_points_3_table* active_points_3_table_handle() const {return ap3_.get();}
  //Active_weighted_points_3_table* active_weighted_points_3_table_handle() const {return awp3_.get();}

  Simulator* simulator_handle() const { return sim_.get();}
  const Static_kernel& static_kernel_object() const {return k_;}
  const Kinetic_kernel& kinetic_kernel_object() const {return kk_;}

  Instantaneous_kernel instantaneous_kernel_object() const {
    return Instantaneous_kernel(*this);
  }

  My_simulation_traits(const Simulator::Time &lb,
		       const Simulator::Time &ub): sim_(new Simulator(lb, ub)),
						   ap1_(new Active_points_1_table()),
						   ap2_(new Active_points_2_table()),
						   ap3_(new Active_points_3_table())
  {}


  bool is_exact() const {
    return true;
  }
protected:
  Simulator::Handle sim_;
  Active_points_1_table::Handle ap1_;
  Active_points_2_table::Handle ap2_;
  Active_points_3_table::Handle ap3_;
  //Active_weighted_points_3_table::Handle awp3_;
  Static_kernel k_;
  Kinetic_kernel kk_;
  Function_kernel fk_;
};




int main(int, char *[])
{
  typedef My_simulation_traits Traits;
  typedef CGAL::Kinetic::Delaunay_triangulation_2<Traits> KDel;

  Traits tr(0,100000.0);
  KDel kdel(tr);

  std::ifstream in("data/points_2");
  in >> *tr.active_points_2_table_handle();
  std::cout << "Read " << tr.active_points_2_table_handle()->size()
	    << " points" << std::endl;

  tr.simulator_handle()->set_current_time(tr.simulator_handle()->end_time());

  return EXIT_SUCCESS;
}
