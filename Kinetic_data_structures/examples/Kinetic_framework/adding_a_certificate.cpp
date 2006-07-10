#include <CGAL/Kinetic/Cartesian_kinetic_kernel.h>
#include <CGAL/Kinetic/Simulation_traits.h>
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
  public CGAL::Kinetic::Cartesian_kinetic_kernel<FunctionKernel> {
  typedef CGAL::Kinetic::Cartesian_kinetic_kernel<FunctionKernel> P;
  typedef My_kinetic_kernel<FunctionKernel> This;
public:
  typedef CGAL::Kinetic::Certificate_generator<This, Positive_x_f_2<This> > Positive_x_2;
  Positive_x_2 positive_x_2_object() const
  {
    return Positive_x_2(P::function_kernel_object());
  }
};


struct My_st_types: public CGAL::Kinetic::Suggested_exact_simulation_traits_types {
  typedef CGAL::Kinetic::Suggested_exact_simulation_traits_types P;
  typedef My_kinetic_kernel<P::Function_kernel>::Point_2 Active_object;
  typedef CGAL::Kinetic::Active_objects_vector<Active_object> Active_objects_table;
  typedef CGAL::Kinetic::Cartesian_instantaneous_kernel< Active_objects_table,
							 Static_kernel> Instantaneous_kernel;
};

struct My_simulation_traits:
  public  CGAL::Kinetic::Simulation_traits<My_st_types::Static_kernel,
			    My_st_types::Kinetic_kernel,
			    My_st_types::Simulator>
{
  typedef  CGAL::Kinetic::Simulation_traits<My_st_types::Static_kernel,
					    My_st_types::Kinetic_kernel,
					    My_st_types::Simulator> P;
  My_simulation_traits(const P::Time &lb= P::Time(0),
		       const P::Time &ub=std::numeric_limits<P::Time>::infinity()): 
    P(lb,ub), 
    ap_(new Active_points_2_table()) {}

  typedef My_st_types::Active_objects_table Active_points_2_table;
  Active_points_2_table* active_points_2_table_handle() {
    return ap_.get();
  }
  const Active_points_2_table* active_points_2_table_handle() const {
    return ap_.get();
  }
  typedef My_st_types::Instantaneous_kernel Instantaneous_kernel;
  Instantaneous_kernel instantaneous_kernel_object() const
  {
    return Instantaneous_kernel(ap_, static_kernel_object());
    }
protected:
  Active_points_2_table::Handle ap_;
};


int main(int, char *[])
{
  typedef My_simulation_traits Traits;
  typedef CGAL::Kinetic::Delaunay_triangulation_2<Traits> KDel;

  Traits tr;
  KDel kdel(tr);

  std::ifstream in("data/points_2");
  in >> *tr.active_points_2_table_handle();
  std::cout << "Read " << tr.active_points_2_table_handle()->size() 
	    << " points" << std::endl;
  
  tr.simulator_handle()->set_current_time(tr.simulator_handle()->end_time());

  return EXIT_SUCCESS;
}
