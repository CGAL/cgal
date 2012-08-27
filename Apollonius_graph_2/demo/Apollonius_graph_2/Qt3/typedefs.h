// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef _MK___TYPEDEFS_H
#define _MK___TYPEDEFS_H

#include <CGAL/basic.h>

#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float  exact_type;
//#include <CGAL/leda_real.h>
//typedef leda_real  exact_type;

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double>        Rep;


#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>


typedef CGAL::Apollonius_graph_filtered_traits_2<Rep,CGAL::Integral_domain_without_division_tag>  Gt;

typedef Gt::Point_2                           Point_2;
typedef Rep::Circle_2                         Circle_2;
typedef Gt::Site_2                            Apollonius_site_2;
typedef Gt::Site_2::Weight                    Weight;

//typedef CGAL::Apollonius_graph_2<Gt> AG_2;
typedef CGAL::Apollonius_graph_hierarchy_2<Gt> AG_2;


typedef AG_2::Vertex_handle                 Vertex_handle;



#endif  // _MK___TYPEDEFS_H
