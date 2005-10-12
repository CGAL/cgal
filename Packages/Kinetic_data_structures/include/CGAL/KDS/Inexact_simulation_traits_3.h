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

#ifndef CGAL_KDS_INEXACT_SIMULATION_3_H
#define CGAL_KDS_INEXACT_SIMULATION_3_H
#include <CGAL/KDS/basic.h>

#include <CGAL/KDS/Simulator.h>
#include <CGAL/KDS/Notifying_table.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <CGAL/KDS/Simulation_traits.h>
#include <CGAL/KDS/Two_list_pointer_event_queue.h>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE
struct Inexact_simulation_types_3: public Suggested_inexact_simulation_types {
  typedef Suggested_inexact_simulation_types P;
  typedef CGAL::KDS::Notifying_table<P::Kinetic_kernel::Point_3> Moving_point_table;
  typedef CGAL::KDS::Cartesian_instantaneous_kernel<Moving_point_table,
						    P::Static_kernel> Instantaneous_kernel;
};

CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

struct Inexact_simulation_traits_3: 
  public Simulation_traits<internal::Inexact_simulation_types_3::Static_kernel,
			   internal::Inexact_simulation_types_3::Kinetic_kernel,
			   internal::Inexact_simulation_types_3::Simulator> {
private:
  typedef  Simulation_traits<internal::Inexact_simulation_types_3::Static_kernel,
			     internal::Inexact_simulation_types_3::Kinetic_kernel,
			     internal::Inexact_simulation_types_3::Simulator> P;

 
  
public:
  /*struct Moving_point_table: public internal::Inexact_simulation_types_3::Moving_point_table {
    typedef boost::intrusive_ptr<Moving_point_table> Pointer;
    typedef boost::intrusive_ptr<Moving_point_table> Const_pointer;
    };*/
  typedef internal::Inexact_simulation_types_3::Moving_point_table Moving_point_table;

  struct Instantaneous_kernel: public internal::Inexact_simulation_types_3::Instantaneous_kernel{
    Instantaneous_kernel(Moving_point_table::Pointer mp, P::Static_kernel k): 
      internal::Inexact_simulation_types_3::Instantaneous_kernel(mp, k){}
  };
 
 

  Inexact_simulation_traits_3(): mp2_(new Moving_point_table()){
  }
  Inexact_simulation_traits_3(P::Simulator::Pointer p): P(p), mp2_(new Moving_point_table()){
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
