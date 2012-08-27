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

#ifndef CGAL_CARTESIAN_REGULAR_INSTANTANEOUS_KERNEL_H
#define CGAL_CARTESIAN_REGULAR_INSTANTANEOUS_KERNEL_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Default_instantaneous_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <map>
#include <iostream>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/internal/To_static.h>
//#include <CGAL/Kinetic/Cartesian_static_converter.h>

/*#define CGAL_MSA(Pred, pred, Arg, d) typedef Instantaneous_adaptor<typename Static_kernel::Pred##_##d, Current_coordinates, Arg> Pred##_##d; \
  Pred##_##d pred##_##d##_object() const				\
  {									\
    typename Static_kernel::Pred##_##d sp= rep_->static_kernel().pred##_##d##_object();	\
    return Pred##_##d(current_coordinates_object(), sp);		\
    }*/

/*#define CGAL_TSO(name) typedef typename Static_kernel::name name*/

#define CGAL_MSA(Pred, pred, Arg, d) typedef Instantaneous_adaptor<typename P::Static_kernel::Pred##_##d, typename P::Kinetic_kernel::Pred##_##d, typename P::Rep, Arg> Pred##_##d; \
  Pred##_##d pred##_##d##_object() const				\
  {									\
    typename P::Static_kernel::Pred##_##d sp= P::rep_->static_kernel().pred##_##d##_object(); \
    typename P::Kinetic_kernel::Pred##_##d kp= P::rep_->kinetic_kernel().pred##_##d##_object(); \
    return Pred##_##d(P::rep_, sp, kp);					\
  }


namespace CGAL { namespace Kinetic {


template <class Traitst >
class Regular_triangulation_instantaneous_kernel: public Default_instantaneous_kernel<Traitst>
{
  typedef Regular_triangulation_instantaneous_kernel<Traitst> This;
public:
  typedef Traitst Traits;
  typedef Default_instantaneous_kernel<Traitst> P;

  //using P::Time;
  //using P::NT;

  Regular_triangulation_instantaneous_kernel(const Traits &tr): P(tr) {
  }

  /*using typename P::Point_1;
  using typename P::Point_2;
  using typename P::Point_3;*/
  typedef typename P::Point_3 Bare_point;
  typedef typename P::Point_3 Weighted_point_3;
  //typedef P::Static_kernel Static_kernel;

  CGAL_MSA(Power_test,power_test, Weighted_point_3, 3);
  CGAL_MSA(Equal, equal, Weighted_point_3, 3);
};
#undef CGAL_MSA

} } //namespace CGAL::Kinetic
#endif
