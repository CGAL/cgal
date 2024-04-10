// Copyright(c) 2022  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_TRIANGULATION_GRAPHICS_SCENE_OPTIONS_H
#define CGAL_POLYGON_TRIANGULATION_GRAPHICS_SCENE_OPTIONS_H

#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Random.h>

template<class PT>
struct Polygon_triangulation_gs_options :
  public CGAL::Graphics_scene_options<typename PT::Triangulation,
                               typename PT::Vertex_handle,
                               typename PT::Finite_edges_iterator,
                               typename PT::Finite_faces_iterator>
{
  using T2=typename PT::Triangulation;
  template<class IPM>
  Polygon_triangulation_gs_options(IPM ipm)
  {
    this->colored_face =
      [](const T2&, const typename PT::Finite_faces_iterator) -> bool
      { return true; };

    this->face_color =
      [](const T2&, const typename PT::Finite_faces_iterator fh) -> CGAL::IO::Color
      {
        CGAL::Random random((unsigned int)(std::size_t)(&*fh));
        return get_random_color(random);
      };

    this->draw_face=
      [ipm](const T2&, const typename PT::Finite_faces_iterator fh) -> bool
      { return get(ipm, fh); };

    this->draw_edge=
      [ipm](const T2& pt, const typename PT::Finite_edges_iterator eh) -> bool
      {
        typename PT::Face_handle fh1=eh->first;
        typename PT::Face_handle fh2=pt.mirror_edge(*eh).first;
        return get(ipm, fh1) || get(ipm, fh2);
      };
  }
};

#endif // CGAL_POLYGON_TRIANGULATION_GRAPHICS_SCENE_OPTIONS_H
