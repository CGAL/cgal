#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_DEFAULT_TRAITS_2_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_DEFAULT_TRAITS_2_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>

CGAL_KINETIC_BEGIN_NAMESPACE

template <class Simulation_traits_t, class Triangulation_t>
class Delaunay_triangulation_default_traits_2 {
protected:
  typedef Simulation_traits_t ST;

public:
  typedef typename ST::Kinetic_kernel::Positive_side_of_oriented_circle_2 SOC;
  typedef typename ST::Kinetic_kernel::Positive_orientation_2 O2;
  typedef typename Triangulation_t::Edge Edge;



  typedef Triangulation_t Triangulation;
  typedef internal::Triangulation_data_structure_helper_2<typename Triangulation::Triangulation_data_structure> TDS_helper;
  typedef typename ST::Simulator Simulator;
  typedef typename SOC::result_type Certificate_data;
  typedef typename Simulator::Time Time;


  typedef typename Simulator::Event_key Event_key;
  typedef typename ST::Active_points_2_table::Data Point_2;
  typedef typename ST::Active_points_2_table::Key Point_key;
  typedef typename ST::Instantaneous_kernel Instantaneous_kernel;
  typedef typename Simulator::NT NT;

 
  
  Delaunay_triangulation_default_traits_2(ST st): st_(st){
    soc_= st_.kinetic_kernel_object().positive_side_of_oriented_circle_2_object();
    o2_= st_.kinetic_kernel_object().positive_orientation_2_object();
  }
  
  typename Simulator::Handle simulator_handle() {
    return st_.simulator_handle();
  }

  typename Simulator::Const_handle simulator_handle() const {
    return st_.simulator_handle();
  }
  
  typedef typename ST::Active_points_2_table Active_points_2_table;
  typename Active_points_2_table::Handle active_points_2_table_handle() {
    return st_.active_points_2_table_handle();
  }


  NT rational_current_time() const {
    return st_.simulator_handle()->rational_current_time();
  }

  const Point_2& point(Point_key k) const {
    return st_.active_points_2_table_handle()->at(k);
  }

  Instantaneous_kernel instantaneous_kernel_object() const {
    return st_.instantaneous_kernel_object();
  }

  CGAL::Comparison_result compare_concurrent(typename Simulator::Event_key a,
					     Edge ea,
					     typename Simulator::Event_key b,
					     Edge eb) const {
    return CGAL::compare(a,b);
  }
 

  void point_changed(Point_key){}

  typedef std::pair<Time, Certificate_data> Certificate_pair;

  Certificate_pair certificate_failure_time(Edge e) {
    Certificate_data s;
    Point_key ks[4]; // must be 4, not 3 for both
     if (is_hull_edge(e)) {
       hull_points(e, ks);
       s= o2_(point(ks[0]), point(ks[1]), point(ks[2]),
	      st_.simulator_handle()->current_time(),
	      st_.simulator_handle()->end_time());
     } else {
       edge_points(e, ks);
       s= soc_(point(ks[0]), point(ks[1]),
	       point(ks[2]), point(ks[3]),
	       st_.simulator_handle()->current_time(),
	       st_.simulator_handle()->end_time());
     }
     return certificate_failure_time(e, s);
  }

  Certificate_pair certificate_failure_time(Edge e,
					    Certificate_data s) {
    Time t= s.failure_time();
    s.pop_failure_time();
    return std::make_pair(t, s);
  }



  bool is_hull_edge(const Edge &e) const {
    return ! TDS_helper::mirror_vertex(e)->point().is_valid()
      || ! TDS_helper::third_vertex(e)->point().is_valid()
      || ! TDS_helper::origin(e)->point().is_valid()
      || ! TDS_helper::destination(e)->point().is_valid();
  }


 
  void edge_points(const Edge &e, Point_key ks[4]) const {
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= TDS_helper::mirror_vertex(e)->point();
  }

  // very dangerous

  void hull_points(const Edge &e,
		   Point_key ks[4]) const {
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= TDS_helper::mirror_vertex(e)->point();

    bool odd_parity=false;
    bool infinity=false;
    for (unsigned int i=0; i<4; ++i) {
      if (infinity) {
	ks[i-1]=ks[i];
      }
      else {
	if (!ks[i].is_valid()) {
	  infinity=true;
	  odd_parity= ((i%2)==1);
	} 
      }
    }
    if (odd_parity) {
      std::swap(ks[0], ks[1]);
    }
  }

  SOC positive_side_of_oriented_circle_2_object() const {
    return soc_;
  }
protected:


  ST st_;
  SOC soc_;
  O2 o2_;
};

CGAL_KINETIC_END_NAMESPACE

#endif
