// Copyright (c) 1999-2003  INRIA Sophia-Antipolis (France).
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
// $Revision$ 
// $Date$
// $Name$
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#define Simple_cartesian Sc

#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <CGAL/_test_weighted_alpha_shape_3.C>

// Filtered_kernel fails with Regular until weighted points are in the kernel.
typedef CGAL::Filtered_exact<double, CGAL::MP_Float> NT;
//typedef CGAL::Simple_cartesian<NT>  SC;

// Try to shorten symbol names (for VC++)
struct K : public CGAL::Simple_cartesian<NT> {};

typedef CGAL::Weighted_alpha_shape_euclidean_traits_3<K> Gt;
typedef Gt::Weighted_point   Weighted_point;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;

typedef CGAL::Regular_triangulation_3<Gt,Tds>        Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;

int main()
{
  _test_weighted_alpha_shape_3( Alpha_shape_3() );
  return 0;
}
