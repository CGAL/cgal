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

#include <boost/format.hpp>

namespace CGAL_nng {

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>              Triangulation;

const std::string Slab[] = {
    "k-nearest-neighbor graph", "Help"
};

const std::string Hmsg[] = {
  "Draw the k-nearest-neighbor graph of a set of points"
};

struct nngIpelet
  : CGAL::Ipelet_base<Kernel, 2> {
  nngIpelet() : CGAL::Ipelet_base<Kernel, 2>("k-nearest-neighbor graph", Slab, Hmsg){}
  void protected_run(int);
};

void nngIpelet::protected_run(int fn)
{

  if(fn == 1){
    show_help();
    return;
  }

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

  int ret_val;
  int kNeighbors=1;

  boost::tie(ret_val,kNeighbors)=request_value_from_user<int>((boost::format("Number of nearest neighbors (default : k=%1%)") % kNeighbors).str() );
  if (ret_val == -1) return;
  if (ret_val == 0) kNeighbors=1;

  Triangulation t(pt_list.begin(), pt_list.end());

  bool edgesDrawn = false;
  for(auto v = t.finite_vertices_begin();
       v != t.finite_vertices_end();
       ++v){

    std::vector<Triangulation::Vertex_handle> kNN;
    CGAL::nearest_neighbors(t, v, kNeighbors+1, std::back_inserter(kNN)); // +1 as v itself counts as its nearest neigbhor for CGAL::nearest_neighbors

    for(const auto & nn : kNN) {
      if(v->point() != nn->point()) {
        draw_in_ipe(Kernel::Segment_2(v->point(), nn->point()));
        edgesDrawn = true;
      }
    }
  }

  if(edgesDrawn) {
    group_selected_objects_();
    // don't create an empty group if no edges are drawn
  }
}
}

CGAL_IPELET(CGAL_nng::nngIpelet)
