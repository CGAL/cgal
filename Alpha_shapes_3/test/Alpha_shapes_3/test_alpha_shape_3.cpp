// Copyright (c) 1999-2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/_test_cls_alpha_shape_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel           Gt;
typedef CGAL::Alpha_shape_vertex_base_3<Gt>                           Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>                             Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>                   Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>                        Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>                          Alpha_shape_3;

typedef CGAL::Alpha_shape_vertex_base_3<Gt,CGAL::Default,CGAL::Tag_true> EF_Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt,CGAL::Default,CGAL::Tag_true>   EF_Fb;
typedef CGAL::Triangulation_data_structure_3<EF_Vb,EF_Fb>                EF_Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,EF_Tds>                        EF_Triangulation_3;
typedef CGAL::Alpha_shape_3<EF_Triangulation_3,CGAL::Tag_true>           EF_Alpha_shape_3;

int main()
{
  _test_cls_alpha_shape_3<Alpha_shape_3>();
  _test_cls_alpha_shape_3<EF_Alpha_shape_3>();
  return 0;
}

// MipsPro prefers this after the other instantiations...
// Explicit instantiation of the whole class :
template class CGAL::Alpha_shape_3<Triangulation_3>;
template class CGAL::Alpha_shape_3<EF_Triangulation_3,CGAL::Tag_true>;

