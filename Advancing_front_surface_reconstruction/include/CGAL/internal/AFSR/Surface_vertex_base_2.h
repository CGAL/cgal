// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_AFSR_SURFACE_VERTEX_BASE_2_H
#define CGAL_AFSR_SURFACE_VERTEX_BASE_2_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


#include <CGAL/basic.h>
#include <CGAL/Triangulation_ds_vertex_base_2.h>

namespace CGAL {
  namespace AFSR {

    template < typename GT,
               typename V3,
               typename Vb = CGAL::Triangulation_ds_vertex_base_2<> >
    class Surface_vertex_base_2
      : public Vb

    {
      typedef typename Vb::Triangulation_data_structure    Tds;
    public:
      typedef GT                                    Geom_traits;
      typedef typename GT::Point_3                  Point;
      typedef Tds                                   Triangulation_data_structure;
      typedef typename Tds::Face_handle             Face_handle;
      typedef typename Tds::Vertex_handle           Vertex_handle;

      template < typename TDS2 >
      struct Rebind_TDS {
        typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
        typedef Surface_vertex_base_2<GT,V3, Vb2>           Other;
      };

    private:
      V3 _vertex;
    public:
      Surface_vertex_base_2() : Vb() {}
      Surface_vertex_base_2(Face_handle f) : Vb(f) {}

      void set_vertex(const V3& v)
      {
        _vertex = v;
      }

      V3 vertex_3() const
      {
        return _vertex;
      }

      const Point&  point() const { return _vertex->point(); }


    };





  } // namespace AFSR
} // namespace CGAL

#endif //CGAL::AFSR_SURFACE_VERTEX_BASE_2_H
