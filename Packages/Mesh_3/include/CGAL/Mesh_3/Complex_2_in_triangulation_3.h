// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_3_H

namespace CGAL {

  namespace Mesh_3 {

    template <typename Tr, typename E_info, typename F_info>
    class Complex_2_in_triangulation_3
    {
    public:
      typedef Complex_2_in_triangulation_3<Tr,
                                           E_info,
                                           F_info> Self;
      typedef Tr Triangulation;
      typedef E_info Edge_info;
      typedef F_info Facet_info;

      typedef typename Triangulation::Geom_traits Geom_traits;

      typedef typename Triangulation::Vertex_handle Vertex_handle;
      typedef typename Triangulation::Edge Edge;
      typedef typename Triangulation::Facet Facet;
      typedef typename Triangulation::Cell_handle Cell_handle;

    public:
      enum Face_type { NOT_IN_COMPLEX , ISOLATED , BOUNDARY ,
                       REGULAR , SINGULAR };

    public:
      void set_in_complex(Facet f)
      {
        set_in_complex(f.first, f.second);
      }

      void set_in_complex(Cell_handle c, int index, bool b = true)
      {
        const Cell_handle& n = c->neighbor(index);
        
        c->set_constrained(index, b);
        n->set_constrained(n->index(c), b);
      }
      
      void set_in_complex(Edge e)
      {
        set_in_complex(e.first, e.second, e.third);
      }
      
      void set_in_complex(Cell_handle c, int i, int j, bool b = true)
      {
        const Vertex_handle& va = c->vertex(i);
        const Vertex_handle& vb = c->vertex(j);

        va->set_is_adjacent_by_constraint(vb, b);
	vb->set_is_adjacent_by_constraint(va, b);
      }

      void remove_from_complex(Facet f)
      {
        remove_from_complex(f.first, f.second);
      }
      

      void remove_from_complex(Cell_handle c, int index)
      {
        set_in_complex(c, index, false);
      }
      
      void remove_from_complex(Edge e)
      {
        remove_from_complex(e.first, e.second, e.third);
      }
      
      void remove_from_complex(Cell_handle c, int i, int j)
      {
        set_in_complex(c, i, j, false);
      }

      bool is_in_complex(Facet f)
      {
        return is_in_complex(f.first, f.second);
      }
      
      bool is_in_complex(Cell_handle c, int index)
      {
        return c->set_constrained(index);
      } 

      bool is_in_complex(Edge e)
      {
        return is_in_complex(e.first, e.second, e.third);
      }
      
      bool is_in_complex(Cell_handle c, int i, int j)
      {
        const Vertex_handle& va = c->vertex(i);
        const Vertex_handle& vb = c->vertex(j);

        const bool result = va->is_adjacent_by_constraint(vb);
        CGAL_assertion( vb->is_adjacent_by_constraint(va) == result );
        
        return result;
      }

    }; // end class Complex_2_in_triangulation_3

  } // end namespace Mesh_3

} // end namespace CGAL

#endif // CGAL_COMPLEX_2_IN_TRIANGULATION_3_H
