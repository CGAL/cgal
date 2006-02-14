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


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_euclidean_traits_3.h>
#include <CGAL/_test_cls_alpha_shape_3.C>


typedef double coord_type;
typedef CGAL::Simple_cartesian<coord_type>  SC;
struct K : public  CGAL::Filtered_kernel<SC> {};

typedef CGAL::Alpha_shape_euclidean_traits_3<K> Gt;
typedef CGAL::Alpha_shape_vertex_base_3<Gt> Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>   Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds> Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

int main()
{
  _test_cls_alpha_shape_3<Alpha_shape_3>();
  return 0;
}

// MipsPro prefers this after the other instantiations...
// Explicit instantiation of the whole class :
template class CGAL::Alpha_shape_3<Triangulation_3>;

