// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://jtournoi@scm.gforge.inria.fr/svnroot/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_complex_3_in_triangulation_3.h $
// $Id: comparison_operators.h  $
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef COMPARISON_OPERATORS_H
#define COMPARISON_OPERATORS_H

#include <CGAL/array.h>
#include <CGAL/Triangulation_utils_3.h>

namespace CGAL {
  namespace Mesh_3 {
    
    // http://en.wikipedia.org/wiki/Sorting_network
    template <typename T, typename Comparator>
    void
    sort4(T& q0, T& q1, T& q2, T& q3, const Comparator& less)
    {
      if(less(q2,q0)) std::swap(q0,q2);
      if(less(q3,q1)) std::swap(q1,q3);
      if(less(q1,q0)) std::swap(q0,q1);
      if(less(q3,q2)) std::swap(q2,q3);
      if(less(q2,q1)) std::swap(q1,q2);
    }

    template <typename T, typename Comparator>
    void
    sort3(T& q0, T& q1, T& q2, const Comparator& less)
    {
      if(less(q2,q0)) std::swap(q0,q2);
      if(less(q1,q0)) std::swap(q0,q1);
      if(less(q2,q1)) std::swap(q1,q2);
    }

    template <typename T, typename Comparator>
    void
    sort2(T& q0, T& q1, const Comparator& less)
    {
      if(less(q1,q0)) std::swap(q0,q1);
    }

    template<typename Vertex_handle>
    struct Vertex_handle_comparator 
      : public std::binary_function<Vertex_handle, Vertex_handle, bool> 
    {
      bool operator()(const Vertex_handle& v1, const Vertex_handle& v2) const 
      { 
        if(v1 == v2)
          return false;
        return v1->point() < v2->point(); 
      }
    };

    template<typename Tr>
    struct Cell_handle_comparator
      : public std::binary_function<typename Tr::Cell_handle,
                                    typename Tr::Cell_handle,
                                    bool>
    {
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Vertex_handle Vertex_handle;

      bool operator()(const Cell_handle& c1, const Cell_handle& c2) const
      {
        if(c1 == c2)
          return false;
        CGAL::cpp11::array<Vertex_handle,4> v1;
        CGAL::cpp11::array<Vertex_handle,4> v2;
        for(int i = 0; i < 4; ++i)
        {
          v1[i] = c1->vertex(i);
          v2[i] = c2->vertex(i);
        }
        Vertex_handle_comparator<Vertex_handle> vcomp;
        sort4(v1[0], v1[1], v1[2], v1[3], vcomp);
        sort4(v2[0], v2[1], v2[2], v2[3], vcomp);
        for(std::size_t i = 0; i < 4; ++i)
        {
          if(v1[i] == v2[i])
            continue;
          else return vcomp(v1[i], v2[i]);
        }
        return false;
      }
    };

    template<typename Tr>
    struct Triangulation_finite_facets_comparator
      : public std::binary_function<typename Tr::Facet,
                                    typename Tr::Facet,
                                    bool>
    {
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Vertex_handle Vertex_handle;

      bool operator()(const Facet& f1, const Facet& f2) const
      {
        if(f1 == f2)
          return false;
        CGAL_PROFILER("Compare facets");
        Vertex_handle_comparator<Vertex_handle> vcomp;
        CGAL::cpp11::array<Vertex_handle,3> vf1;
        CGAL::cpp11::array<Vertex_handle,3> vf2;
        for(int i = 0; i < 3; ++i)
        {
          vf1[i] = f1.first->vertex(
            Triangulation_utils_3::vertex_triple_index(f1.second,i));
          vf2[i] = f2.first->vertex(
            Triangulation_utils_3::vertex_triple_index(f2.second,i));
        }
        sort3(vf1[0], vf1[1], vf1[2], vcomp);
        sort3(vf2[0], vf2[1], vf2[2], vcomp);
        for(std::size_t i = 0; i < 3; ++i)
        {
          if(vf1[i] == vf2[i])
            continue;
          else return vcomp(vf1[i], vf2[i]);
        }
        return false;
      }
    };
    
    template<typename Facet>
    struct Polyhedron_Facet_handle_comparator  
      : public std::binary_function<typename Facet::Facet_handle, 
                                    typename Facet::Facet_handle, 
                                    bool> 
    {
      typedef typename Facet::Facet_handle  Facet_handle;
      typedef typename Facet::Vertex_handle Vertex_handle;
      typedef typename Facet::Halfedge_around_facet_circulator Facet_he_circ;

      bool operator()(const Facet_handle& pf1, const Facet_handle& pf2) const
      { 
        if(pf1 == pf2)
          return false;
        //collect vertices of both facets
        CGAL::cpp11::array<Vertex_handle, 3> vf1;
        CGAL::cpp11::array<Vertex_handle, 3> vf2;
        Facet_he_circ begin =	pf1->facet_begin();
        Facet_he_circ end = begin;
        std::size_t i = 0;
        do
        {
          vf1[i++] = begin->vertex(); 
          ++begin;
        }while(begin != end);

        begin =	pf2->facet_begin();
        end = begin;
        i = 0;
        do
        {
          vf2[i++] = begin->vertex(); 
          ++begin;
        }while(begin != end);

        //compare vertices
        Vertex_handle_comparator<Vertex_handle> vcomp;
        sort3(vf1[0], vf1[1], vf1[2], vcomp);
        sort3(vf2[0], vf2[1], vf2[2], vcomp);
        for(std::size_t i = 0; i < 3; i++)
        {
          if(vf1[i] == vf2[i])  
            continue;
          else 
            return vcomp(vf1[i], vf2[i]);
        }
        return false; //it is the same facet
      }
    };

    template<typename Halfedge>
    struct Halfedge_handle_comparator 
      : public std::binary_function<typename Halfedge::Halfedge_handle,
                                    typename Halfedge::Halfedge_handle,
                                    bool> 
    {
      typedef typename Halfedge::Halfedge_handle Halfedge_handle;
      typedef typename Halfedge::Vertex_handle Vertex_handle;

      bool operator()(const Halfedge_handle& he1, const Halfedge_handle& he2) const
      {     
        if(he1 == he2)
          return false;

        //collect vertices of both facets
        CGAL::cpp11::array<Vertex_handle,2> vhe1;
        vhe1[0] = he1->vertex();
        vhe1[1] = he1->opposite()->vertex();
      
        CGAL::cpp11::array<Vertex_handle,2> vhe2;
        vhe2[0] = he2->vertex();
        vhe2[1] = he2->opposite()->vertex();

        //compare vertices
        Vertex_handle_comparator<Vertex_handle> vcomp;
        sort2(vhe1[0], vhe1[1], vcomp);
        sort2(vhe2[0], vhe2[1], vcomp);

        //we want he and he->opposite() to be both in the set, for flooding
        if(he1 == he2->opposite())
          return vcomp(he1->vertex(), he2->vertex());
      
        if(vhe1[0] == vhe2[0])
          return vcomp(vhe1[1], vhe2[1]);
        else //either < or >
          return vcomp(vhe1[0], vhe2[0]);
      }
    };


  } //  namespace Mesh_3 {
} // namespace CGAL

#endif // COMPARISON_OPERATORS_H
