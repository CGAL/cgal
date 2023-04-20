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

#ifndef CGAL_AFSR_CONSTRUCT_POLYHEDRON_2
#define CGAL_AFSR_CONSTRUCT_POLYHEDRON_2

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

namespace CGAL {

  template <class Triangulation, class Filter>
  class Advancing_front_polyhedron_reconstruction;

  namespace AFSR {

    template <class HDS, class Surface>
    class Construct_polyhedron: public CGAL::Modifier_base<HDS> {

      const Surface& s;

    public:
      Construct_polyhedron(Surface& s)
        : s(s)
      {}

      void operator()( HDS& hds)
      {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( s.number_of_vertices(), s.number_of_facets(), 6* s.number_of_facets());

        typedef typename Surface::TDS_2 TDS_2;
        typedef typename TDS_2::Face_iterator Face_iterator;
        typedef typename TDS_2::Vertex_iterator Vertex_iterator;

        const TDS_2& tds = s.triangulation_data_structure_2();

        std::size_t index = 0;
        Vertex_iterator end = tds.vertices_end();

        for(Vertex_iterator vit = tds.vertices_begin(); vit != end; ++vit){
          if(vit->vertex_3() != s.triangulation_3().infinite_vertex()){
            B.add_vertex(vit->point());
            vit->vertex_3()->id() = index++;
          }
        }

        for(Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit){

          if(fit->is_on_surface()){
            B.begin_facet();
            for(int i=0; i < 3; i++){
              B.add_vertex_to_facet(fit->vertex(i)->vertex_3()->id());
            }
            B.end_facet();
          }
        }
        B.end_surface();
      }

    };

    template <class Polyhedron, class Surface>
    void
    construct_polyhedron(Polyhedron& P, Surface& surface)
    {
      typedef typename Polyhedron::HalfedgeDS             HalfedgeDS;
      Construct_polyhedron<HalfedgeDS, Surface> builder(surface);
      P.delegate(builder);
    }

  } // namespace AFSR
} // namespace CGAL
#endif // CGAL_AFSR_CONSTRUCT_POLYHEDRON_2
