#ifndef CGAL_KDS_EXACT_SIMULATION_1_H
#define CGAL_KDS_EXACT_SIMULATION_1_H
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

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
struct Exact_simulation_types_1: Suggested_exact_simulation_types {
  typedef CGAL::KDS::Notifying_table<Kinetic_kernel::Point_1> Moving_point_table;
  typedef CGAL::KDS::Cartesian_instantaneous_kernel<Moving_point_table,
						    Static_kernel> Instantaneous_kernel;
  
};

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

struct Exact_simulation_traits_1: 
  public Simulation_traits<internal::Exact_simulation_types_1::Static_kernel,
			   internal::Exact_simulation_types_1::Kinetic_kernel,
			   internal::Exact_simulation_types_1::Simulator> {
private:
  typedef  Simulation_traits<internal::Exact_simulation_types_1::Static_kernel,
			     internal::Exact_simulation_types_1::Kinetic_kernel,
			     internal::Exact_simulation_types_1::Simulator> P;
public:
 
  typedef internal::Exact_simulation_types_1::Instantaneous_kernel Instantaneous_kernel;
  typedef internal::Exact_simulation_types_1::Moving_point_table Moving_point_table;
 

  Exact_simulation_traits_1(): mp2_(new Moving_point_table()){

  }
  Moving_point_table* moving_point_table_pointer(){ return mp2_.get();}
  const Moving_point_table* moving_point_table_pointer() const { return mp2_.get();}
 
  Instantaneous_kernel instantaneous_kernel_object() const {
    return Instantaneous_kernel(mp2_,P::static_kernel_object());
  }
protected:
  Moving_point_table::Pointer mp2_;
};
CGAL_KDS_END_NAMESPACE

#endif
