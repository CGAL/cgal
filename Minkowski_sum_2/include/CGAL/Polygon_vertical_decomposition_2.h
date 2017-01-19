// Copyright (c) 2013  Tel-Aviv University (Israel).
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
// Author(s) : Efi Fogel   <efifogel@gmail.com>

#ifndef CGAL_POLYGON_VERTICAL_DECOMPOSITION_2_H
#define CGAL_POLYGON_VERTICAL_DECOMPOSITION_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Gps_segment_traits_2.h>

#include <vector>
#include <list>

namespace CGAL {


/*!
 * \class
 * Vertical decomposition strategy.
 */
template <typename Kernel_,
          typename Container_ = std::vector<typename Kernel_::Point_2> >
class Polygon_vertical_decomposition_2 {
public:
  typedef Kernel_                                        Kernel;
  typedef Container_                                     Container;

  typedef CGAL::Polygon_2<Kernel, Container>             Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel, Container>  Polygon_with_holes_2;

  typedef CGAL::Arr_segment_traits_2<Kernel>             Arr_segment_traits;
  typedef CGAL::Gps_segment_traits_2<Kernel, Container, Arr_segment_traits>
                                                         Traits_2;
  typedef CGAL::General_polygon_set_2<Traits_2>          General_polygon_set_2;
  typedef typename General_polygon_set_2::Arrangement_2  Arrangement_2;

  typedef typename Arrangement_2::Halfedge_const_iterator
    Halfedge_const_iterator;
  typedef typename Arrangement_2::Face_const_iterator    Face_const_iterator;

  typedef typename Arrangement_2::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle      Face_const_handle;

  typedef typename Arrangement_2::Vertex_handle          Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement_2::Face_handle            Face_handle;

  typedef std::pair<Vertex_const_handle, std::pair<CGAL::Object, CGAL::Object> >
                                                         Vert_decomp_entry;
  typedef std::list<Vert_decomp_entry>                   Vert_decomp_list;

  typedef typename Kernel::Point_2                       Point_2;

  typedef typename Arrangement_2::Outer_ccb_const_iterator
    Outer_ccb_const_iterator;

private:
  typedef typename Arrangement_2::X_monotone_curve_2     Segment_2;
  typedef typename Kernel::Line_2                        Line_2;

  typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;

  // An arrangement observer, used to receive notifications of face splits and
  // face mergers.
  class My_observer : public CGAL::Arr_observer<Arrangement_2> {
  public:
    My_observer(Arrangement_2& arr) : Arr_observer<Arrangement_2>(arr) {}

    virtual void after_split_face(Face_handle f, Face_handle new_f,
                                  bool /* is_hole */)
    { if (f->contained()) new_f->set_contained(true); }
  };

  // Kernel functors:
  typedef typename Kernel::Compare_x_2                   Compare_x_2;
  typedef typename Kernel::Intersect_2                   Intersect_2;
  typedef typename Kernel::Equal_2                       Equal_2;

  // Data members:
  const Traits_2* m_traits;
  bool m_own_traits;    // inidicates whether the kernel should be freed up.

  Compare_x_2 f_cmp_x;
  Intersect_2 f_intersect;
  Equal_2     f_equal;

public:
  /*! Default constructor. */
  Polygon_vertical_decomposition_2() :
    m_traits(NULL),
    m_own_traits(false)
  { init(); }

  /*! Constructor */
  Polygon_vertical_decomposition_2(const Traits_2& traits) :
    m_traits(&traits),
    m_own_traits(false)
  { init(); }

  /*! Initialize */
  void init()
  {
    // Allocate the traits if not provided.
    if (m_traits == NULL) {
      m_traits = new Traits_2;
      m_own_traits = true;
    }

    // Obtain kernel functors.
    const Kernel* kernel = m_traits;
    f_cmp_x = kernel->compare_x_2_object();
    f_intersect = kernel->intersect_2_object();
    f_equal = kernel->equal_2_object();
  }

  // Destructor
  ~Polygon_vertical_decomposition_2()
  {
    if (m_own_traits) {
      if (m_traits) {
        delete m_traits;
        m_traits = NULL;
      }
      m_own_traits = false;
    }
  }

  /*! Obtain the traits
   * \return the traits
   */
  const Traits_2& traits() const { return *m_traits; }

  /*! Decompose a polygon into convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename OutputIterator_>
  OutputIterator_ operator()(const Polygon_2& pgn, OutputIterator_ oi) const
  { return decomp(pgn, oi); }

  /*! Decompose a polygon with holes into convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename OutputIterator_>
  OutputIterator_
  operator()(const Polygon_with_holes_2& pgn, OutputIterator_ oi) const
  { return decomp(pgn, oi); }

private:
  /*!
   * Decompose a polygon-with-holes into convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename Polygon_, typename OutputIterator_>
  OutputIterator_ decomp(const Polygon_& pgn, OutputIterator_ oi) const
  {
    const Traits_2& traits = *m_traits;
    General_polygon_set_2 gps(traits);
    gps.insert(pgn);
    Arrangement_2& arr = gps.arrangement();
    My_observer obs(arr);
    vertical_decomposition(arr);
    Face_const_iterator fi;
    for (fi = arr.faces_begin(); fi != arr.faces_end(); ++fi) {
      if (! fi->contained()) continue;
      CGAL_assertion(fi->number_of_outer_ccbs() == 1);
      Outer_ccb_const_iterator oci = fi->outer_ccbs_begin();
      Halfedge_const_iterator first = *oci;
      Halfedge_const_iterator curr = first;
      Polygon_2 pgn;
      do {
        pgn.push_back(curr->target()->point());
        curr = curr->next();
      } while (curr != first);
      *oi++ = pgn;
    }
    return oi;
  }

  // Add a vertical segment from the given vertex to some other arrangement
  // feature.
  Halfedge_const_handle
  add_vertical_segment(Arrangement_2& arr, Vertex_handle v, CGAL::Object obj)
    const
  {
    Segment_2 seg;
    Vertex_const_handle vh;
    Halfedge_const_handle hh;
    Face_const_handle fh;
    Vertex_handle v2;

    if (CGAL::assign(vh, obj)) {
      // The given feature is a vertex.
      seg = Segment_2(v->point(), vh->point());
      v2 = arr.non_const_handle(vh);
    }
    else if (CGAL::assign(hh, obj)) {
      // The given feature is a halfedge. We ignore fictitious halfedges.
      if (hh->is_fictitious())
        return Halfedge_const_handle();

      // Check whether v lies in the interior of the x-range of the edge (in
      // which case this edge should be split).
      if (f_cmp_x(v->point(), hh->target()->point()) == CGAL::EQUAL) {
        // In case the target of the edge already has the same x-coordinate as
        // the vertex v, just connect these two vertices.
        seg = Segment_2(v->point(), hh->target()->point());
        v2 = arr.non_const_handle(hh->target());
      }
      else {
        // Compute the vertical projection of v onto the segment associated
        // with the halfedge. Split the edge and connect v with the split point.
        Line_2 supp_line(hh->source()->point(), hh->target()->point());
        Line_2 vert_line(v->point(),
                         Point_2(v->point().x(), v->point().y() + 1));
        Point_2 point;
        CGAL::assign(point, f_intersect(supp_line, vert_line));
        seg = Segment_2(v->point(), point);
        arr.split_edge(arr.non_const_handle(hh),
                       Segment_2(hh->source()->point(), point),
                       Segment_2(point, hh->target()->point()));
        v2 = arr.non_const_handle(hh->target());
      }
    }
    // Ignore faces and empty objects.
    else return Halfedge_const_handle();

    // Add the vertical segment to the arrangement using its two end vertices.
    return arr.insert_at_vertices(seg, v, v2);
  }

  // Construct the vertical decomposition of the given arrangement.
  void vertical_decomposition(Arrangement_2& arr) const
  {
    // For each vertex in the arrangment, locate the feature that lies
    // directly below it and the feature that lies directly above it.
    Vert_decomp_list vd_list;
    CGAL::decompose(arr, std::back_inserter(vd_list));

    // Go over the vertices (given in ascending lexicographical xy-order),
    // and add segements to the feautres below and above it.
    typename Vert_decomp_list::iterator it, prev = vd_list.end();
    for (it = vd_list.begin(); it != vd_list.end(); ++it) {
      // If the feature above the previous vertex is not the current vertex,
      // add a vertical segment to the feature below the vertex.
      Vertex_const_handle v;
      if ((prev == vd_list.end()) ||
          !CGAL::assign(v, prev->second.second) ||
          !f_equal(v->point(), it->first->point()))
        add_vertical_segment(arr, arr.non_const_handle(it->first),
                             it->second.first);
      // Add a vertical segment to the feature above the vertex.
      add_vertical_segment(arr, arr.non_const_handle(it->first),
                           it->second.second);
      prev = it;
    }
  }
};

} //namespace CGAL

#endif
