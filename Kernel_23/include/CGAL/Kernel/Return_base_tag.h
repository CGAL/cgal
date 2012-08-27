// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_RETURN_BASE_TAG_H
#define CGAL_KERNEL_RETURN_BASE_TAG_H

#include <CGAL/basic.h>

// This is a simple tag which is used as additional (first) argument in
// some kernel functors, to tell them to return the base (rep) class,
// instead of the main type (e.g. Kernel_base::Point_2 instead
// of Kernel::Point_2).  This is a minor optimization which prevents
// useless copies of the "reps".

// Those functors are only those used in the constructors of the kernel
// types like Point_2, so it's limited.

// The real solution will be to use "forwarding constructors", when they
// will be available in C++.
// In the mean time, this should be a mostly/hopefully internal hack.

namespace CGAL {

struct Return_base_tag {};

} //namespace CGAL

#endif
