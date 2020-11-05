// Copyright (c) 2007  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>
//                                                 Guy Zucker        <guyzucke@post.tau.ac.il>

#ifndef CGAL_CONNECT_HOLES_H
#define CGAL_CONNECT_HOLES_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_default_dcel.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <list>
#include <set>
#include <CGAL/Unique_hash_map.h>
#include <iostream>

namespace CGAL {

template <class HANDLE>
struct _Less_handle
{
  bool operator() (const HANDLE& vh1, const HANDLE& vh2) const
  {
    return (&(*vh1) < &(*vh2));
  }
};

/*bool is_adjacent_hole(const Halfedge_handle he, Arrangement_2 arr) {
    return ((!he->twin()->face()->contained())
              && (he->twin()->face() != arr.unbounded_face()))
}*/

/*!
 * Connect the given polygon with holes, turning it into a sequence of
 * points, where the holes are connceted to the outer boundary using
 * zero-width passages.
 * For example:
 *              Input                             Output
 *  +----------------------------+    +-----*---------------*------+
 *  |                            |    |     |               |      |
 *  |     +------+        +--+   |    |     +------+        +--+   |
 *  |     |      |        |  |   |    |     |      |        |  |   |
 *  |     +------+         \ |   |    |     +----*-+         \ |   |
 *  |                       \|   |    |          |            \|   |
 *  |          +----+            |    |          +----+            |
 *  |         /     |            |    |         /     |            |
 *  |        +------+            |    |        +------+            |
 *  |                            |    |                            |
 *  +----------------------------+    +----------------------------+
 *
 * \param pwh The polygon with holes.
 * \param oi Output: An output iterator for the points.
 * \pre The polygons has an outer boundary.
 * \return A past-the-end iterator of the points.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator connect_holes(const Polygon_with_holes_2<Kernel,
                             Container>& pwh,
                             OutputIterator oi)
{
  typedef Polygon_2<Kernel,Container>              Polygon_2;
  typedef Arr_segment_traits_2<Kernel>             Traits_2;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Traits_2::X_monotone_curve_2    Segment_2;
  typedef CGAL::Gps_default_dcel<Traits_2>                  dcel;
  typedef General_polygon_set_2<Gps_segment_traits_2<Kernel, Container, Traits_2> , dcel>
                                                                                                                                          General_polygon_set_2;
  typedef Arrangement_2<Gps_segment_traits_2<Kernel, Container, Traits_2>, dcel>   Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

  CGAL_precondition (! pwh.is_unbounded());

  // In case the polygon has not holes, just go over its outer boundary
  // and report the points along it.
  const Polygon_2&                       outer_pgn = pwh.outer_boundary();
  typename Polygon_2::Vertex_circulator  first_v, curr_v, next_v;

  if (! pwh.has_holes())
  {
    first_v = curr_v = outer_pgn.vertices_circulator();

    do
    {
      *oi = *curr_v;
      ++oi;
      ++curr_v;

    } while (curr_v != first_v);

    return (oi);
  }

  // Construct the arrangement of all segments.
  General_polygon_set_2 gps(pwh);
  Arrangement_2         arr = gps.arrangement();

  // The resulting arrangment contains a single hole in the unbounded face,
  // which comprises a face f, with several holes in its interior.
  // Go over these holes and pick the topmost vertex in each hole.
  const Face_handle                      uf = arr.unbounded_face();
  typename Arrangement_2::Hole_iterator  f_hole_it = uf->holes_begin();
  const Face_handle                      f = (*f_hole_it)->twin()->face();
  typename Arrangement_2::Ccb_halfedge_circulator
                                         first, circ;
  Kernel                                 ker;
  typename Kernel::Compare_y_2           comp_y = ker.compare_y_2_object();
  typename Kernel::Compare_x_2           comp_x = ker.compare_x_2_object();
  Vertex_handle                          v_top;
  Comparison_result                      res;
  //construct vertex hash map (default data is false) and insert top vertices
  //handles as keys (key vertex data value is true)
  typedef typename CGAL::Unique_hash_map<Vertex_const_handle, bool> V_map;
  V_map top_vertices(false ,arr.number_of_faces());

  /*traversal of arrangement face holes - a hole in the face
  in arranements is disjoint from the outer boundary (different from BOP
  where the hole can have vertices along the outer boundary). We look for
  holes only inside faces that are part of the point set. This
  guarantees that if the input PWH had holes with vertices on the outer
  boundary they will not be connected by additional vertical segments to
  the outer boundary or another hole's vertex/edge.
  */
  typename Arrangement_2::Face_const_handle fit;
  //std::cout << arr.number_of_faces() << " faces:" << std::endl;
  typename Arrangement_2::Hole_const_iterator f_hole_itc;
  for (fit = arr.faces_begin(); fit != arr.faces_end() ; ++fit)
    if (fit->contained()) {
      for (f_hole_itc = fit->holes_begin(); f_hole_itc != fit->holes_end(); ++f_hole_itc){
        // Locate the topmost vertex in the current hole. In case of two (or more)
        // vertices with maximal y-coordinate, select the leftmost one.
        first = circ = (arr.non_const_handle(*f_hole_itc));
        v_top = circ->target();
        do
        {
          ++circ;
          res = comp_y (circ->target()->point(), v_top->point());
          if (res == CGAL::LARGER ||
            (res == CGAL::EQUAL &&
             comp_x (circ->target()->point(), v_top->point()) == CGAL::SMALLER))
          {
            v_top = circ->target();
          }
        } while (circ != first);
       // std::cout << "inserted top vertex at " << v_top->point() <<std::endl;
                        top_vertices[Vertex_const_handle (v_top)]=true;
    }
  }



  // Perform "vertical ray shooting" from each arrangement vertex, locating
  // the features that lie below and above it.
  typedef std::list<std::pair<Vertex_const_handle,
                              std::pair<Object, Object> > >   Ray_shoot_list;

  Ray_shoot_list                     vrs_list;
  typename Ray_shoot_list::iterator  vrs_iter;

  decompose (arr, std::back_inserter (vrs_list));

  // Go over the results of the batched vertical ray-shooting query.
  Vertex_const_handle                v;
  Vertex_handle                      v_above;
  Halfedge_const_handle              he;
  Halfedge_handle                    he_above;
  typename Kernel::Direction_2          dir_up (0, 1);
  typename Kernel::Construct_ray_2      ray = ker.construct_ray_2_object();
  typename Kernel::Construct_segment_2  segment =
                                            ker.construct_segment_2_object();
  typename Kernel::Intersect_2          intersect = ker.intersect_2_object();
  Object                                obj;
  Point_2                               ip;
  bool                                  assign_success;

  for (vrs_iter = vrs_list.begin(); vrs_iter != vrs_list.end(); ++vrs_iter)
  {
    if (!top_vertices.is_defined(vrs_iter->first))
      continue;

    v_top = arr.non_const_handle (vrs_iter->first);

    // In case the current vertex is a top vertex in one of the holes,
    // add a vertical segment connecting it with the feature above it.
    if (CGAL::assign (v, vrs_iter->second.second))
    {
      // v_top lies below a vertex v_above: Connect these two vertices.
      v_above = arr.non_const_handle (v);

      arr.insert_at_vertices (Segment_2 (v_top->point(), v_above->point()),
                              v_top, v_above);
      //std::cout << "connected ((" << v_top->point() << "),( " << v_above->point() << "))" <<std::endl;
    }
    else if (CGAL::assign (he, vrs_iter->second.second))
    {
      // v_top lies below the interior of the hafledge he_above:
      // Find the intersection of this halfegde with a vertical ray
      // emanating from v_top.
      he_above = arr.non_const_handle (he);

      obj = intersect (ray (v_top->point(), dir_up),
                       segment (he_above->source()->point(),
                                he_above->target()->point()));

      assign_success = CGAL::assign (ip, obj);
      CGAL_assertion (assign_success);

      if (assign_success)
      {
        // Split he_above at the computed intersection point.
        arr.split_edge (he_above,
                        Segment_2 (he_above->source()->point(), ip),
                        Segment_2 (ip, he_above->target()->point()));

        // Now he_above is split such that it becomes the predecessor
        // halfedge for the insertion of the vertical segment connecting
        // v_top and ip.
        arr.insert_at_vertices (Segment_2 (v_top->point(), ip),
                                he_above, v_top);
      //added for debugging
      //std::cout << "connected ((" << v_top->point() << "),( " << ip << "))" <<std::endl;
      }
    }
    else
    {
      // We should never reach here.
      CGAL_error_msg("top vertex is located in an unbounded face.");
    }
  }

  // The holes of the face f are now all connected to its outer boundary.
  // Go over this boundary and report the vertices along it.
  // Note that we start with a vertex located on the original PWH outer boundary.
  typedef typename Arrangement_2::Halfedge_const_iterator
                                              Halfedge_const_iterator;

  /*For the traversal we switch between two states using key vertices.
  A key vertex is a vertex on the outer boundary of the PWH that has  a degree>2
  (note that this includes new vertices added by the vertical ray shooting).
  In other words, it is a vertice that leads to holes. One state is searching for a key vertice that
  leads to a hole that has not been traversed. The other state is traversing (marking) a hole completely
  (starting and finishing at the same key vertice). Once the hole traversal is completed
  return to search state*/

  //start state is search state
  bool marking_hole_state = false;

  /*flags used for printing the edges for debugging purposes
  bool skip_print = false;
  bool antenna_trav = false; */

  //marker of first vertex of hole that is being traversed
  Vertex_handle hole_start, empty_handle;

  /*create a hash map container for outer boundary key vertices (degree>2
  and on outer boundary. Hash map is constructed with a
  size parameter so we'll traverse once around the outer boundary to
  count the key vertices and once to insert them*/

  f_hole_it = uf->holes_begin();
  Halfedge_const_handle he_han = *f_hole_it;
  if (he_han->face() != uf)
    he_han = he_han->twin();
  CGAL_assertion(he_han->face() == uf);
  //std::cout << "outer boundary:" <<std::endl;
  //std::cout << "(" << he_han->target()->point() << ")" <<std::endl;
  //calculate num of vertices on outer boundary for hash map creation.
  std::size_t size = 1;
  Halfedge_const_handle begin = he_han;
  he_han = he_han->next();
  CGAL_assertion(he_han->face() == uf);
  while (he_han != begin) {
    //std::cout << "(" << he_han->target()->point() << ")" <<std::endl;
    if (he_han->target()->degree()>2)
            size++;
    he_han = he_han->next();
    CGAL_assertion(he_han->face() == uf);
  }

   //construct vertex hash map (default data is false) and insert vertex
  //handles as keys (key vertex data value is true)
  V_map ver_map(false ,size);
  if (he_han->target()->degree()>2)
           ver_map[he_han->target()]=true;
  //std::cout << "(" << he_han->target()->point() << ")" <<std::endl;
  he_han = he_han->next();
  while (he_han != begin) {
    CGAL_assertion(he_han->face() == uf);
    if (he_han->target()->degree()>2)
      ver_map[he_han->target()]=true;
    //std::cout << "(" << he_han->target()->point() << ")" <<std::endl;
    he_han = he_han->next();
  }
  //std::cout << "outer boundary finished" <<std::endl;

  //get iterator to edge on outer boundary
  first = f->outer_ccb();
  //std::cout << "first edge is ((" << first->source()->point() << "),(" << first->target()->point() << "))" <<std::endl;
  Halfedge_const_iterator  start, curr, next;
  start = first;
  while (start->twin()->face() != uf)
    start = start->next();
  //std::cout << "start edge is ((" << start->source()->point() << "),(" << start->target()->point() << "))" <<std::endl;
  curr = start;
  do
  {
    /*traverse_hole (next for any vertice besides vertices on
    the outer boundary with deg>2. check hole on every step and
    mark new holes (if two holes are connected by a vertex)
    This is true also for traversing starting with an antenna halfedge*/
    if (marking_hole_state) {
       /*Add current hole to search structure.
       we traverse the hole fully an if we "incidently" traverse
       another hole as well (joined by a vertice) both will be traversed
       completely*/

       //case we are starting to traverse an antenna
      if (curr->face() == curr->twin()->face()) {
        //antenna_trav = true ;
        //std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
        *oi = curr->target()->point();
        ++oi;
        curr = curr->next();
      }
      Face_handle curr_face;
      //Traversal of the hole
      while (hole_start != arr.non_const_handle(curr->target())) {
        curr_face = arr.non_const_handle(curr->twin()->face());
        /*mark the hole as a hole that has been traversed to save
         multiple traversals*/
        curr_face->set_visited(true);
        //std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
        *oi = curr->target()->point();
        ++oi;
        //"turn inside" instead of next if target is located on outer boundary
        //(a key vertex encountered while traversing a hole)
        if (ver_map.is_defined(curr->target())) {
          curr = curr->twin()->prev()->twin();
        } else {//regular advance
          curr = curr->next();
        }
      } //exited loop target is the hole marking start vertex
      hole_start=empty_handle;
      marking_hole_state=false;
      /* next iteration willcheck all of the edges which have (curr->target())
      as their source to determine next move -
      does curr->target() have more holes to traverse or do we
      return to searching the outer boundary*/
    } else
    {//search for next hole to traverse

      /*Treatment of 4 possible cases (can be narrowed to 3) - the
       first two cannot co-exist here as vertices with degree of 2 that are
       not on the unbounded face will be encountered when traversing holes
       (different state) */
       next = curr->next();
      //unified 2 simple cases of traversal along outer boundary
      /* case target() is a simple vertice with a single possible edge on path.
      add it to output and continue searching for a hole. This is when traversing
      along the outer polygon boundary*/
      /*the case the next half edge is on the outer boundary meaning we finished
       handling holes connected to the target vertex.
       Add the target to the output, and continue searching for next hole*/

      /*if (curr->target()->degree()==2) {
        //std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
        //insert target point to result output iterator
        *oi = curr->target()->point();
        ++oi;
        curr = curr->next();
        //maintain the same state
        continue;
      }*/
      if ((curr->target()->degree()==2) ||(next->twin()->face() == uf)) {
        /*if (!skip_print)
          std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
        else {
           skip_print=false;
           antenna_trav=false;
        }*/
        *oi = curr->target()->point();
        ++oi;
        curr = curr->next();
        //maintain the same state
        continue;
      }

      /*if next is a boundary of a hole that has not been traversed,
       or an antenna (which should also be completely traversed until
       returning to the source vertex). an antenna to the outer boundary
       occurs only for polygons that are detached from the outer boundary
       and therefore these holes are connected to the boundary only by a
       single vertex and will be travrsed only once.
       Target must be added to output and state changed.
       */
      if  ((!next->twin()->face()->visited())
            //case of antenna
                      || (next->face()==next->twin()->face())) {
        /*if (!skip_print)
           std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
        else {
           skip_print=false;
           //can be reached after an antenna
           antenna_trav=false;
        }*/
        //insert target point to result output iterator
        *oi = curr->target()->point();
        ++oi;
        hole_start = arr.non_const_handle(curr->target());
        marking_hole_state=true;
        curr = curr->next();
        continue;
      } else {
        /*the case where target() is a part of several holes, and next is a
        boundary of a hole that has been traversed. This requires
        to continue "traversal alongside all edges whose target is also
        curr->target() to keep looking for a hole that hasn't been traversed.
        we do not insert the target to the output set to avoid duplication
        with cases 2 and 3*/

        /*bypass traversed hole from inside*/

        /*print debugging comment:
        the current edge needs to be printed now before moving on, but the
        current (target) vertex will be inserted to output iterator in case 2*/
        /*if (traversed_holes.find(arr.non_const_handle(curr->twin()->face())) == traversed_holes.end())
          std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
        */
        //case this is the last half edge of an antenna
        /*if (antenna_trav == true) {
          std::cout << "curr edge is ((" << curr->source()->point() << "),( " << curr->target()->point() << "))" <<std::endl;
          antenna_trav = false;
        }
        skip_print = true; */
        curr = next->twin()->next()->twin();
        //maintain the same state
        continue;
        /*now curr->target() has remained the same but curr->next() will
        be adjacent to a different hole or curr->next->twin->face()
        will be the unbounded face*/
      }
    }
  } while (curr != start);
  return (oi);
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
