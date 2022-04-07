// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_IMPLICIT_SURFACE_MESHER_VISITOR_H
#define CGAL_MESH_3_IMPLICIT_SURFACE_MESHER_VISITOR_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_2/Output_stream.h>

namespace CGAL {

  namespace Mesh_3 {

    template <
      typename Tr,
      typename Previous_level
      >
    class Visitor_for_surface {
      Previous_level* previous;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;

      typedef Previous_level Previous_visitor;

      Visitor_for_surface(Previous_visitor* p)
        : previous(p) {}

      template <typename E, typename P>
      void before_conflicts(E, P) const {}


      template <typename E, typename P, typename Z>
      void before_insertion(E, P, Z) const {}

      void after_insertion(const Vertex_handle& v)
      {
//         if(v->point().surface_index() == 0)
//         {
//           CGAL_MESHES_OUTPUT_STREAM << "?";
//           v->point().set_surface_index(1);
//         }
//         CGAL_MESHES_OUTPUT_STREAM << v->point().surface_index();
      }

      template <typename E, typename P, typename Z>
      void after_no_insertion(E, P, Z) const {}

      Previous_visitor& previous_level()
      {
        return *previous;
      }

    }; // end class Visitor_for_surface

 }  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_IMPLICIT_SURFACE_MESHER_VISITOR_H
