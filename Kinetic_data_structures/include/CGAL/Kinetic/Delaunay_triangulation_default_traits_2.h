// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_DELAUNAY_TRIANGULATION_DEFAULT_TRAITS_2_H
#define CGAL_KINETIC_DELAUNAY_TRIANGULATION_DEFAULT_TRAITS_2_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/tds_2_helpers.h>

namespace CGAL { namespace Kinetic {

template <class Simulation_traits_t, class Triangulation_t>
class Delaunay_triangulation_default_traits_2 {
protected:
  typedef Simulation_traits_t ST;

public:
  typedef typename ST::Kinetic_kernel::Side_of_oriented_circle_2 SOC;
  typedef typename ST::Kinetic_kernel::Orientation_2 O2;
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
    soc_= st_.kinetic_kernel_object().side_of_oriented_circle_2_object();
    o2_= st_.kinetic_kernel_object().orientation_2_object();
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


  NT rational_current_time() {
    NT nt= st_.simulator_handle()->next_time_representable_as_nt();
    /*std::cout << "Next time is " << nt << " " << CGAL::to_double(nt) << std::endl;
    std::cout << "Current time is " << st_.simulator_handle()->current_time() << " "
    << CGAL::to_double(st_.simulator_handle()->current_time()) << std::endl;*/
    CGAL_assertion(st_.simulator_handle()->current_time()== nt);
    //st_.simulator_handle()->set_current_time(nt);
    return nt;
  }

  const Point_2& point(Point_key k) const {
    return st_.active_points_2_table_handle()->at(k);
  }

  Instantaneous_kernel instantaneous_kernel_object() const {
    return st_.instantaneous_kernel_object();
  }

  CGAL::Comparison_result compare_concurrent(typename Simulator::Event_key a,
					     Edge ,
					     typename Simulator::Event_key b,
					     Edge ) const {
    if (a < b) return CGAL::SMALLER;
    else if (b < a) return CGAL::LARGER;
    else return CGAL::EQUAL;
    //return CGAL::compare(a,b);
  }
 

  void point_changed(Point_key){}

  typedef std::pair<Time, Certificate_data> Certificate_pair;

  bool internal_certificate_failure_time(Edge e, Point_key ks[4],  Time &t, Certificate_data &s) {
    s= soc_(point(ks[0]), point(ks[1]),
	    point(ks[2]), point(ks[3]),
	    st_.simulator_handle()->current_time(),
	    st_.simulator_handle()->end_time());
    return return_certificate_failure_time(e, t, s);
  }

  bool certificate_failure_time(Edge e, Certificate_data cd, Time &t, Certificate_data &s ) {
    s=cd;
    return return_certificate_failure_time(e, t, s);
  }


  bool hull_certificate_failure_time(Edge e,Point_key ks[3], Time &t, Certificate_data &s) {
    s= o2_(point(ks[0]), point(ks[1]), point(ks[2]),
	   st_.simulator_handle()->current_time(),
	   st_.simulator_handle()->end_time());
   
    return return_certificate_failure_time(e, t, s);
  }

  bool return_certificate_failure_time(Edge , Time &t, Certificate_data &s ) {
    if (s.will_fail()) {
      t= s.failure_time();
      s.pop_failure_time();
      return true;
    } else {
      return false;
    }
  }

  SOC positive_side_of_oriented_circle_2_object() const {
    return soc_;
  }
  
  bool is_exact() const {
    return st_.is_exact();
  }

protected:


  ST st_;
  SOC soc_;
  O2 o2_;
};

} } //namespace CGAL::Kinetic

#endif
