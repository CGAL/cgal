// Copyright (c) 2023 Inria
// All rights reserved.
//
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Daniel Funke

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include <CGAL/CGAL_Ipelet_base.h>

namespace CGAL_nng {

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>              Triangulation;

const std::string Slab[] = {
  "NNG", "Help"
};

const std::string Hmsg[] = {
  "Draw the nearest-neighbor graph of a set of points"
};

struct nngIpelet
  : CGAL::Ipelet_base<Kernel, 2> {
  nngIpelet() : CGAL::Ipelet_base<Kernel, 2>("Nearest-neighbor graph", Slab, Hmsg){}
  void protected_run(int);
};

void nngIpelet::protected_run(int /*fn*/)
{
  std::list<Point_2> pt_list;

  read_active_objects(
    CGAL::dispatch_or_drop_output<Point_2>(
      std::back_inserter(pt_list)
    )
  );

  if (pt_list.empty()) {
    print_error_message("No mark selected");
    return;
  }

  Triangulation t(pt_list.begin(), pt_list.end());

  for(auto v = t.finite_vertices_begin();
       v != t.finite_vertices_end();
       ++v){

    auto nn = CGAL::nearest_neighbor(t, v);
    draw_in_ipe(Kernel::Segment_2(v->point(), nn->point()));
  }
}
}

CGAL_IPELET(CGAL_nng::nngIpelet)
