// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_ADAPTATION_TRAITS_BASE_2_H
#define CGAL_VORONOI_DIAGRAM_2_ADAPTATION_TRAITS_BASE_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>
#include <CGAL/Voronoi_diagram_2/Cached_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class DG, class AS, class CVP, class NS>
class Adaptation_traits_base_2
{
private:
  typedef Adaptation_traits_base_2<DG,AS,CVP,NS>  Self;

public:
  typedef DG   Delaunay_graph;
  typedef AS   Access_site_2;
  typedef CVP  Construct_Voronoi_point_2;
  typedef NS   Nearest_site_2;

  typedef typename Delaunay_graph::Vertex_handle    Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Face_handle      Delaunay_face_handle;
  typedef typename Delaunay_graph::Edge             Delaunay_edge;

  typedef typename Functor_exists<Nearest_site_2>::Value  Has_nearest_site_2;

  Access_site_2 access_site_2_object() const {
    return Access_site_2();
  }

  Construct_Voronoi_point_2 construct_Voronoi_point_2_object() const {
    return Construct_Voronoi_point_2();
  }

  Nearest_site_2 nearest_site_2_object() const {
    return Nearest_site_2();
  }
};

//=========================================================================
//=========================================================================


} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_ADAPTATION_TRAITS_BASE_2_H
