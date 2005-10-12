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

#ifndef CGAL_KDS_INEXACT_SIMULATION_1_H
#define CGAL_KDS_INEXACT_SIMULATION_1_H
#include <CGAL/KDS/basic.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/KDS/Simulator.h>
#include <CGAL/KDS/Cartesian_kinetic_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/KDS/Notifying_table.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <CGAL/KDS/Simulation_traits.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/KDS/Heap_pointer_event_queue.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
struct Inexact_simulation_types_1: public Suggested_inexact_simulation_types {
  typedef Suggested_inexact_simulation_types P;

  typedef CGAL::KDS::Notifying_table<Kinetic_kernel::Point_1> Moving_point_table;
  typedef CGAL::KDS::Cartesian_instantaneous_kernel<Moving_point_table,
						    Static_kernel> Instantaneous_kernel;
  
};

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

struct Inexact_simulation_traits_1: 
  public Simulation_traits<internal::Inexact_simulation_types_1::Static_kernel,
			   internal::Inexact_simulation_types_1::Kinetic_kernel,
			   internal::Inexact_simulation_types_1::Simulator> {
private:
  typedef  Simulation_traits<internal::Inexact_simulation_types_1::Static_kernel,
			     internal::Inexact_simulation_types_1::Kinetic_kernel,
			     internal::Inexact_simulation_types_1::Simulator> P;
public:
 
  typedef internal::Inexact_simulation_types_1::Instantaneous_kernel Instantaneous_kernel;
  typedef CGAL::KDS::Notifying_table<internal::Inexact_simulation_types_1::Kinetic_kernel::Point_1> Moving_point_table;
 

  Inexact_simulation_traits_1(): mp2_(new Moving_point_table()){

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
