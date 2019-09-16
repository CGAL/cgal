// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri <andreas.fabri@geometryfactory.com>

#ifndef CGAL_NEF_DEFAULT_ITEMS_H
#define CGAL_NEF_DEFAULT_ITEMS_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_3/SNC_items.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Extended_cartesian.h>

namespace CGAL {

template<class Kernel>
struct Default_items {
  typedef CGAL::SNC_indexed_items Items;
};

template<typename NT>
struct Default_items<CGAL::Extended_homogeneous<NT> > {
  typedef CGAL::SNC_items Items;
};

template<typename NT>
struct Default_items<CGAL::Extended_cartesian<NT> > {
  typedef CGAL::SNC_items Items;
};

} //namespace CGAL
#endif
