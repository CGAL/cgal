// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_WEIGHTED_EXACT_SIMULATION_3_H
#define CGAL_KDS_WEIGHTED_EXACT_SIMULATION_3_H
#include <CGAL/KDS/basic.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/KDS/Simulator.h>
#include <CGAL/KDS/Cartesian_kinetic_kernel.h>
#include <CGAL/Polynomial/Upper_bound_root_enumerator.h>
#include <CGAL/Polynomial/Upper_bound_enumerator_Descartes_traits.h>
#include <CGAL/KDS/Skip_even_roots_filter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/KDS/Notifying_table.h>
#include <CGAL/KDS/Regular_triangulation_instantaneous_traits_3.h>
#include <CGAL/KDS/Simulation_traits.h>
#include <CGAL/KDS/Two_list_pointer_queue.h>
#include <CGAL/KDS/Exact_simulation_traits_3.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
struct Exact_weighted_simulation_types_3: public Exact_simulation_types_3 {  
};

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

struct Exact_weighted_simulation_traits_3: 
  public Simulation_traits<internal::Exact_weighted_simulation_types_3::Static_kernel,
			   internal::Exact_weighted_simulation_types_3::Kinetic_kernel,
			   internal::Exact_weighted_simulation_types_3::Simulator> {
private:
  typedef  Simulation_traits<internal::Exact_weighted_simulation_types_3::Static_kernel,
			     internal::Exact_weighted_simulation_types_3::Kinetic_kernel,
			     internal::Exact_weighted_simulation_types_3::Simulator> P;
public:

  typedef CGAL::KDS::Notifying_table<internal::Exact_weighted_simulation_types_3::Kinetic_kernel::Weighted_point_3> Moving_point_table;
  typedef CGAL::KDS::Regular_triangulation_instantaneous_traits_3<Moving_point_table,
								  P::Static_kernel> Instantaneous_kernel;
  
 

  Exact_weighted_simulation_traits_3(): mp2_(new Moving_point_table()){
  }
  Exact_weighted_simulation_traits_3(P::Simulator::Pointer p): P(p), mp2_(new Moving_point_table()){
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
