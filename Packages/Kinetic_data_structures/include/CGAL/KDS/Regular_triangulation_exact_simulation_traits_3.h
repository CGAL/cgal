#ifndef CGAL_KDS_WEIGHTED_EXACT_SIMULATION_3_H
#define CGAL_KDS_WEIGHTED_EXACT_SIMULATION_3_H
#include <CGAL/KDS/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/KDS/Cartesian_kinetic_kernel.h>
#include <CGAL/KDS/Notifying_table.h>
#include <CGAL/KDS/Regular_triangulation_instantaneous_traits_3.h>
#include <CGAL/KDS/Simulation_traits.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/KDS/Exact_simulation_traits_3.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
struct Regular_triangulation_exact_simulation_types_3: public Exact_simulation_types_3  {
 
  typedef Exact_simulation_types_3 P;
  typedef CGAL::Regular_triangulation_euclidean_traits_3<P::Static_kernel> Static_kernel;

  typedef CGAL::KDS::Notifying_table<P::Kinetic_kernel::Weighted_point_3> Moving_point_table;
  typedef CGAL::KDS::Regular_triangulation_instantaneous_traits_3<Moving_point_table,
								  Static_kernel> Instantaneous_kernel;
  
};

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

struct Regular_triangulation_exact_simulation_traits_3: 
  public Simulation_traits<internal::Regular_triangulation_exact_simulation_types_3::Static_kernel,
			   internal::Regular_triangulation_exact_simulation_types_3::Kinetic_kernel,
			   internal::Regular_triangulation_exact_simulation_types_3::Simulator> {
private:
  typedef  Simulation_traits<internal::Regular_triangulation_exact_simulation_types_3::Static_kernel,
			     internal::Regular_triangulation_exact_simulation_types_3::Kinetic_kernel,
			     internal::Regular_triangulation_exact_simulation_types_3::Simulator> P;
public:
 
  typedef internal::Regular_triangulation_exact_simulation_types_3::Instantaneous_kernel Instantaneous_kernel;
  typedef internal::Regular_triangulation_exact_simulation_types_3::Moving_point_table Moving_point_table;
 

  Regular_triangulation_exact_simulation_traits_3(): mp2_(new Moving_point_table()){
  }
  Regular_triangulation_exact_simulation_traits_3(P::Simulator::Pointer p): P(p), mp2_(new Moving_point_table()){
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
