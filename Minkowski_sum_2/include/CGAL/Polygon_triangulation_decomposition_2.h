// Copyright (c) 2013  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel   <efifogel@gmail.com>

#ifndef CGAL_POLYGON_TRIANGULATION_DECOMPOSITION_2_H
#define CGAL_POLYGON_TRIANGULATION_DECOMPOSITION_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/General_polygon_set_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <iostream>
#include <vector>
#include <list>

namespace CGAL {

/*! \class Polygon_triangulation_decomposition Polygon_triangulation_decomposition.h
 * Constrained triangulation decomposition strategy.
 */
template <typename Kernel_,
          typename Container_ = std::vector<typename Kernel_::Point_2> >
class Polygon_triangulation_decomposition_2 {
public:
  typedef Kernel_                                       Kernel;
  typedef Container_                                    Container;

  typedef CGAL::Polygon_2<Kernel, Container>            Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel, Container> Polygon_with_holes_2;

private:
  struct Face_info {
    Face_info() {}
    int nesting_level;
    bool in_domain() { return nesting_level % 2 == 1; }
  };

  // Triangulation types
  typedef CGAL::Triangulation_vertex_base_2<Kernel>                     VB;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>  FBI;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel, FBI>      FB;
  typedef CGAL::Triangulation_data_structure_2<VB, FB>                  TDS;
  typedef CGAL::Exact_predicates_tag                                    Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

public:
  /*! Default constructor. */
  Polygon_triangulation_decomposition_2() {}

  // Destructor
  ~Polygon_triangulation_decomposition_2() {}

  /*! Decompose a polygon-with-holes into convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename OutputIterator_>
  OutputIterator_ operator()(const Polygon_2& pgm, OutputIterator_ oi) const
  {
    CDT cdt;

    // Insert boundary:
    insert_polygon(cdt, pgm);

    // Mark facets that are inside the domain bounded by the polygon
    mark_domains(cdt);

    typename CDT::Finite_faces_iterator fit;
    for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      if (! fit->info().in_domain()) continue;
      typename CDT::Vertex_handle vh0 = fit->vertex(0);
      typename CDT::Vertex_handle vh1 = fit->vertex(1);
      typename CDT::Vertex_handle vh2 = fit->vertex(2);
      Polygon_2 pgn;
      pgn.push_back(vh0->point());
      pgn.push_back(vh1->point());
      pgn.push_back(vh2->point());
      *oi++ = pgn;
    }
    return oi;
  }

  /*! Decompose a polygon-with-holes into convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename OutputIterator_>
  OutputIterator_
  operator()(const Polygon_with_holes_2& pgm, OutputIterator_ oi) const
  {
    CDT cdt;

    // Insert outer boundary:
    const Polygon_2& outer_pgn = pgm.outer_boundary();
    insert_polygon(cdt, outer_pgn);

    // Insert holes:
    typename Polygon_with_holes_2::Hole_const_iterator ith;
    for (ith = pgm.holes_begin(); ith != pgm.holes_end(); ++ith)
      insert_polygon(cdt, *ith);

    // Mark facets that are inside the domain bounded by the polygon
    mark_domains(cdt);

    typename CDT::Finite_faces_iterator fit;
    for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
      if (! fit->info().in_domain()) continue;
      typename CDT::Vertex_handle vh0 = fit->vertex(0);
      typename CDT::Vertex_handle vh1 = fit->vertex(1);
      typename CDT::Vertex_handle vh2 = fit->vertex(2);
      Polygon_2 pgn;
      pgn.push_back(vh0->point());
      pgn.push_back(vh1->point());
      pgn.push_back(vh2->point());
      *oi++ = pgn;
    }
    return oi;
  }

private:
  void mark_domains(CDT& cdt, typename CDT::Face_handle start, int index,
                    std::list<typename CDT::Edge>& border) const
  {
    if (start->info().nesting_level != -1) return;
    std::list<typename CDT::Face_handle> queue;
    queue.push_back(start);
    while (! queue.empty()) {
      typename CDT::Face_handle fh = queue.front();
      queue.pop_front();
      if (fh->info().nesting_level == -1) {
        fh->info().nesting_level = index;
        for (int i = 0; i < 3; i++) {
          typename CDT::Edge e(fh,i);
          typename CDT::Face_handle n = fh->neighbor(i);
          if (n->info().nesting_level == -1) {
            if (cdt.is_constrained(e)) border.push_back(e);
            else queue.push_back(n);
          }
        }
      }
    }
  }

  // Explore set of facets connected with non constrained edges,
  // and attribute to each such set a nesting level.
  // We start from facets incident to the infinite vertex, with a nesting
  // level of 0. Then we recursively consider the non-explored facets incident
  // to constrained edges bounding the former set and increase the nesting
  // level by 1.
  // Facets in the domain are those with an odd nesting level.
  void mark_domains(CDT& cdt) const
  {
    typename CDT::All_faces_iterator it;
    for (it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it)
      it->info().nesting_level = -1;

    std::list<typename CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while (! border.empty()) {
      typename CDT::Edge e = border.front();
      border.pop_front();
      typename CDT::Face_handle n = e.first->neighbor(e.second);
      if (n->info().nesting_level == -1)
        mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
  }

  template <typename Polygon_>
  void insert_polygon(CDT& cdt, const Polygon_& pgn) const
  {
    typedef Polygon_                                    Polygon_2;
    if (pgn.is_empty()) return;
    typedef typename Polygon_2::Vertex_circulator       Vertex_circulator;

    Vertex_circulator v_start = pgn.vertices_circulator();
    typename CDT::Vertex_handle vt_prev = cdt.insert(*v_start);
    Vertex_circulator v_curr = v_start;
    do {
      ++v_curr;
      typename CDT::Vertex_handle vt_curr = cdt.insert(*v_curr);
      cdt.insert_constraint(vt_prev, vt_curr);
      vt_prev = vt_curr;
    } while (v_curr != v_start);
  }
};

} //namespace CGAL

#endif
