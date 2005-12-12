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

#ifndef CGAL_KDS_EXACT_SIMULATION_1_H
#define CGAL_KDS_EXACT_SIMULATION_1_H
#include <CGAL/KDS/Simulation_traits.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Exact_simulation_traits_1:
public internal::Suggested_exact_simulation_traits<internal::Sest_types::Kinetic_kernel::Point_1>
{
    typedef internal::Suggested_exact_simulation_traits<internal::Sest_types::Kinetic_kernel::Point_1> P;
    Exact_simulation_traits_1(const P::Time &lb=0,
        const P::Time &ub=std::numeric_limits<P::Time>::infinity()): P(lb,ub){}
};
CGAL_KDS_END_NAMESPACE
#endif
