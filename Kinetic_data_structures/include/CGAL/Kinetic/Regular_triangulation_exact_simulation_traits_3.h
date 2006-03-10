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
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_WEIGHTED_EXACT_SIMULATION_3_H
#define CGAL_KINETIC_WEIGHTED_EXACT_SIMULATION_3_H
#include <CGAL/Kinetic/Regular_triangulation_instantaneous_traits_3.h>
#include <CGAL/Kinetic/Simulation_traits.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

CGAL_KINETIC_BEGIN_INTERNAL_NAMESPACE
struct Rest3_types: public Suggested_exact_simulation_traits_types
{
  typedef CGAL::Regular_triangulation_euclidean_traits_3< Suggested_exact_simulation_traits_types::Static_kernel> Static_kernel;
};

CGAL_KINETIC_END_INTERNAL_NAMESPACE

CGAL_KINETIC_BEGIN_NAMESPACE

struct Regular_triangulation_exact_simulation_traits_3:
  public  Simulation_traits<internal::Rest3_types::Static_kernel,
			    internal::Rest3_types::Kinetic_kernel,
			    internal::Rest3_types::Simulator> 
{
  typedef Simulation_traits<internal::Rest3_types::Static_kernel,
			    internal::Rest3_types::Kinetic_kernel,
			    internal::Rest3_types::Simulator> P;
  
  typedef Active_objects_vector<P::Kinetic_kernel::Weighted_point_3> Active_points_3_table;
  Active_points_3_table* active_points_3_table_handle() {
    return ap_.get();
  }
  const Active_points_3_table* active_points_3_table_handle() const {
    return ap_.get();
  }

  typedef Regular_triangulation_instantaneous_traits_3<Active_points_3_table, P::Static_kernel> Instantaneous_kernel;
  Instantaneous_kernel instantaneous_kernel_object() const {
    return Instantaneous_kernel(ap_, P::static_kernel_object());
  }
  Regular_triangulation_exact_simulation_traits_3(const P::Time &lb = 0,
						   const P::Time &ub = std::numeric_limits<P::Time>::infinity()): P(lb,ub), 
														  ap_(new Active_points_3_table()){}
 

protected:
  Active_points_3_table::Handle ap_;
};
CGAL_KINETIC_END_NAMESPACE
#endif
