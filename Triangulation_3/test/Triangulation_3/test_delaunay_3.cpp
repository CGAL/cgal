// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Francois Rebufat

#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include "test_dependencies.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

bool del=true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>
#include <CGAL/_test_cls_parallel_triangulation_3.h>

// Explicit instantiation of the whole class :
template class CGAL::Delaunay_triangulation_3<K>;

template <typename K, typename ConcurrencyTag = CGAL::Sequential_tag>
using Tds = CGAL::Triangulation_data_structure_3<CGAL::Triangulation_vertex_base_3<K>,
                                                 CGAL::Delaunay_triangulation_cell_base_3<K>,
                                                 ConcurrencyTag>;

template <typename K, typename ConcurrencyTag = CGAL::Sequential_tag>
using Tds_index = CGAL::Triangulation_data_structure_3<CGAL::VertexWithPoint<K>,
                                                       CGAL::Cell4Delaunay<K>,
                                                       ConcurrencyTag, CGAL::Index_tag>;

int main()
{
  using Cls = CGAL::Delaunay_triangulation_3<EPIC>;
  using Cls_with_epec = CGAL::Delaunay_triangulation_3<EPEC>;

  _test_cls_delaunay_3( Cls() );
  _test_cls_delaunay_3( Cls_with_epec() );

  using Cls_index = CGAL::Delaunay_triangulation_3<EPIC, Tds_index<EPIC>>;
  using Cls_with_epec_index = CGAL::Delaunay_triangulation_3<EPEC, Tds_index<EPEC>>;

  _test_cls_delaunay_3(Cls_index());
  _test_cls_delaunay_3(Cls_with_epec_index());

  using Tds_Delaunay_Cb = Tds<K>;
  using Cls_with_Delaunay_Cb = CGAL::Delaunay_triangulation_3<EPIC, Tds_Delaunay_Cb>;

  _test_cls_delaunay_3(Cls_with_Delaunay_Cb());

#ifdef CGAL_LINKED_WITH_TBB
  using Lock_ds = CGAL::Spatial_lock_grid_3<CGAL::Tag_priority_blocking>;
  using Tds_parallel = Tds<EPIC, CGAL::Parallel_tag>;
  using Cls_parallel = CGAL::Delaunay_triangulation_3<EPIC, Tds_parallel, CGAL::Default, Lock_ds>;
  // The following test won't do things in parallel since it doesn't provide
  // a lock data structure
  _test_cls_delaunay_3(Cls_parallel());
  // This test performs parallel operations
  _test_cls_parallel_triangulation_3(Cls_parallel());

  using Tds_parallel_index = Tds_index<EPIC, CGAL::Parallel_tag>;
  using Cls_parallel_index = CGAL::Delaunay_triangulation_3<EPIC, Tds_parallel_index, CGAL::Default, Lock_ds>;
  // The following test won't do things in parallel since it doesn't provide
  // a lock data structure
  _test_cls_delaunay_3(Cls_parallel_index());
  // This test performs parallel operations
  _test_cls_parallel_triangulation_3(Cls_parallel_index());
#endif

  // Second version for the circumcenter storing cell base class.
  using Vb = CGAL::Triangulation_vertex_base_3<K>;
  using Cb = CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K>;
  using TDS = CGAL::Triangulation_data_structure_3<Vb, Cb>;
  using Cls_circumcenter = CGAL::Delaunay_triangulation_3<K, TDS>;

  _test_cls_delaunay_3(Cls_circumcenter());

  using Vb4 = CGAL::VertexWithPoint<EPIC>;
  using Cb4 = CGAL::CellWithCircumcenter<EPIC>;
  using TDS4 = CGAL::Triangulation_data_structure_3<Vb4, Cb4, CGAL::Sequential_tag, CGAL::Index_tag>;
  using Cls4 = CGAL::Delaunay_triangulation_3<EPIC, TDS4>;

  _test_cls_delaunay_3(Cls4());

  return 0;
}
