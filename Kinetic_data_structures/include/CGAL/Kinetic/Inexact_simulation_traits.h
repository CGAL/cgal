// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kinetic_data_structures/include/CGAL/Kinetic/Exact_simulation_traits_2.h $
// $Id: Exact_simulation_traits_2.h 35766 2007-01-20 21:39:01Z drussel $
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_INEXACT_SIMULATION_H
#define CGAL_KINETIC_INEXACT_SIMULATION_H
#include <CGAL/Kinetic/basic.h>

#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/Kernel.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kinetic/Active_objects_vector.h>
#include <CGAL/Kinetic/Default_instantaneous_kernel.h>
#include <CGAL/Kinetic/Cartesian.h>
#include <CGAL/Kinetic/Derivitive_filter_function_kernel.h>
#include <CGAL/Kinetic/Default_simulator.h>
#include <CGAL/Kinetic/Heap_pointer_event_queue.h>

CGAL_KINETIC_BEGIN_NAMESPACE

struct Inexact_simulation_traits {
  typedef Inexact_simulation_traits This;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Static_kernel;
  typedef CGAL::POLYNOMIAL::Polynomial<Static_kernel::FT> Function;
  typedef CGAL::POLYNOMIAL::Root_stack_default_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Numeric_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;
  typedef CGAL::Kinetic::Derivitive_filter_function_kernel<Function_kernel> Simulator_function_kernel_base;
  struct Simulator_function_kernel: public Simulator_function_kernel_base{};

  typedef Cartesian<Simulator_function_kernel> Kinetic_kernel;
  typedef Heap_pointer_event_queue<Function_kernel> Event_queue;
  typedef Default_simulator<Simulator_function_kernel, Event_queue > Simulator;

  typedef Active_objects_vector<Kinetic_kernel::Point_1> Active_points_1_table;
  typedef Active_objects_vector<Kinetic_kernel::Point_2> Active_points_2_table;
  typedef Active_objects_vector<Kinetic_kernel::Point_3> Active_points_3_table;
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

  Inexact_simulation_traits(const Simulator::Time &lb,
			    const Simulator::Time &ub): sim_(new Simulator(lb, ub)),
							ap1_(new Active_points_1_table()),
							ap2_(new Active_points_2_table()),
							ap3_(new Active_points_3_table())
			   //awp3_(new Active_weighted_points_3_table())
{}
 
  
  bool is_exact() const {
    return false;
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
CGAL_KINETIC_END_NAMESPACE
#endif
