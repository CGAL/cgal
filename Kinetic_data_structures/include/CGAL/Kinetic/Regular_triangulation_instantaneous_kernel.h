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
#include <CGAL/Kinetic/Default_instantaneous_mapped_kernel.h>

namespace CGAL { namespace Kinetic {

template <class Traitst >
class Regular_triangulation_instantaneous_kernel
  : public Default_instantaneous_mapped_kernel<Traitst>
{
  typedef Regular_triangulation_instantaneous_kernel<Traitst> This;
public:
  typedef Traitst Traits;
  typedef Default_instantaneous_mapped_kernel<Traitst> P;

  Regular_triangulation_instantaneous_kernel(const Traits &tr): P(tr) { }
};

} } // namespace CGAL::Kinetic
#endif // CGAL_CARTESIAN_REGULAR_INSTANTANEOUS_KERNEL_H
