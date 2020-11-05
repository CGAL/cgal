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

#ifndef CGAL_AFSR_SURFACE_FACE_BASE_2_H
#define CGAL_AFSR_SURFACE_FACE_BASE_2_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


// This face class stores a facet in the tetrahedrization
// When it gets reoriented by the TDS, it also changes the facet

namespace CGAL {
  namespace AFSR {

    template < typename GT,
               typename F3,
               typename Fb = CGAL::Triangulation_ds_face_base_2<> >
    class Surface_face_base_2
      : public Fb
    {
      typedef typename Fb::Triangulation_data_structure    Tds;

    public:
      typedef typename Tds::Face_handle             Face_handle;
      typedef typename Tds::Vertex_handle           Vertex_handle;

      template < typename TDS2 >
      struct Rebind_TDS {
        typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
        typedef Surface_face_base_2<GT, F3, Fb2>           Other;
      };

    private:
      F3 _facet;
      bool _is_on_surface;

    public:
      Surface_face_base_2()
        : Fb(), _is_on_surface(true)
      {}

      Surface_face_base_2(Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2)
        : Fb(v0, v1, v2), _is_on_surface(true)
      {}

      Surface_face_base_2(Vertex_handle v0,
                          Vertex_handle v1,
                          Vertex_handle v2,
                          Face_handle n0,
                          Face_handle n1,
                          Face_handle n2)
        : Fb(v0, v1, v2, n0, n1, n2), _is_on_surface(true)
      {}

      void set_facet(const F3& facet)
      {
        _facet = facet;
      }

      const F3& facet() const
      {
        return _facet;
      }

      void set_is_on_surface(bool is_on_surface)
      {
        _is_on_surface = is_on_surface;
      }

      bool is_on_surface() const
      {
        return _is_on_surface;
      }


      void reorient()
      {
        Fb::reorient();
        if( is_on_surface()){
          _facet = std::make_pair(_facet.first->neighbor(_facet.second),
                                  _facet.first->neighbor(_facet.second)->index(_facet.first));
        }
      }

    };


  } // namespace AFSR
} // namespace CGAL

#endif // CGAL_AFSR_SURFACE_FACE_BASE_2_H
