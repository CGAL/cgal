// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>
//                 Ophir Setter <ophirset@post.tau.ac.il>

#ifndef CGAL_GRAPH_TRAITS_DUAL_ARRANGEMENT_2_H
#define CGAL_GRAPH_TRAITS_DUAL_ARRANGEMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * Definition of the specialized Dual<Arrangement_2> class,
 * and the specialized boost::graph_traits<Dual<Arrangement_2> >class.
 */

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_on_surface_with_history_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>

namespace CGAL {

// Forward declaration.
template <class Type> class Dual;

}





#define CGAL_ARGT_CLASS CGAL::Arrangement_on_surface_2
#include <CGAL/Arrangement_2/graph_traits_Dual.h>
#define CGAL_ARGT_CLASS CGAL::Arrangement_2
#include <CGAL/Arrangement_2/graph_traits_Dual.h>
#define CGAL_ARGT_CLASS CGAL::Arrangement_on_surface_with_history_2
#include <CGAL/Arrangement_2/graph_traits_Dual.h>
#define CGAL_ARGT_CLASS CGAL::Arrangement_with_history_2
#include <CGAL/Arrangement_2/graph_traits_Dual.h>

#include <CGAL/enable_warnings.h>

#endif
