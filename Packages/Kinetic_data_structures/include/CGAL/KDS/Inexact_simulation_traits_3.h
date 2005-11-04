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

#include <CGAL/KDS/Simulation_traits.h>




CGAL_KDS_BEGIN_NAMESPACE

struct Inexact_simulation_traits_3: 
  public internal::Suggested_inexact_simulation_traits<internal::Sist_types::Kinetic_kernel::Point_3> {
};

CGAL_KDS_END_NAMESPACE

#endif
