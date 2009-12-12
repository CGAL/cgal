#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/Kinetic/Default_simulator.h>

#include <CGAL/Polynomial/Sturm_root_stack_traits.h>
#include <CGAL/Polynomial/Sturm_root_stack.h>
#include <CGAL/Kinetic/Active_objects_vector.h>
#include <CGAL/Kinetic/Default_instantaneous_kernel.h>
#include <CGAL/Kinetic/Cartesian.h>
#include <CGAL/Kinetic/Handle_degeneracy_function_kernel.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Kinetic/Active_objects_vector.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template <class P>
struct Point_with_color: public P {
  typedef P Base_point;
  typedef typename P::Coordinate C;
  Point_with_color(const P &p, int c): P(p), color(c){}
  Point_with_color(const C &x, const C &y): P(x,y){}
  Point_with_color(){}
  int color;
};
template <class P>
std::istream &operator<<(std::istream &in, Point_with_color<P> &p) {
  typename P::Base_point bp;
  in >> bp;
  int c;
  in >> c;
  p= Point_with_color<P>(bp, c);
  return in;
}

template <class P>
std::ostream &operator>>(std::ostream &out, const Point_with_color<P> &p) {
  out << static_cast<typename P::Base_point>(p) << " " << p.color;
  return out;
}


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

  typedef CGAL::Kinetic::Cartesian<Simulator_function_kernel> Kinetic_kernel;
  typedef CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel> Event_queue;
  typedef CGAL::Kinetic::Default_simulator<Simulator_function_kernel, Event_queue > Simulator;

  typedef Point_with_color<Kinetic_kernel::Point_2> Point_2;


  typedef CGAL::Kinetic::Active_objects_vector<Kinetic_kernel::Point_1> Active_points_1_table;
  typedef CGAL::Kinetic::Active_objects_vector<Point_2> Active_points_2_table;
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


int main()
{
  typedef My_simulation_traits::Kinetic_kernel::Point_2 Moving_point_2;
  
  typedef CGAL::Kinetic::Delaunay_triangulation_2<My_simulation_traits> KDel;
  
  My_simulation_traits tr(0, 10000);
  My_simulation_traits::Simulator::Handle sp= tr.simulator_handle();
  
  KDel kdel(tr);
  
  kdel.set_has_certificates(false);
  std::ifstream in("data/points_with_color_2");
  in >> *tr.active_points_2_table_handle();
  
  

  kdel.set_has_certificates(true);
  
  
  std::cout << "Starting to run" << std::endl;
  while (sp->next_event_time()
	 < sp->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+10);
    std::cout << "At time " << sp->current_time() << ":\n";
    std::cout << kdel.triangulation_data_structure();
  }
  
  return EXIT_SUCCESS;
}
