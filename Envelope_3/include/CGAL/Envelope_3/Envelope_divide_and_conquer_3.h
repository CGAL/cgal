// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>
//                 Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_3_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_3_H

#define CGAL_ENVELOPE_SAVE_COMPARISONS
#define CGAL_ENVELOPE_USE_BFS_FACE_ORDER

#include <iostream>
#include <list>
#include <set>
#include <vector>
#include <map>
#include <time.h>

#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_3/Envelope_overlay_2.h>
#include <CGAL/Envelope_3/Envelope_element_visitor_3.h>
#include <CGAL/Envelope_3/set_dividors.h>

#ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER
#include <CGAL/Arr_face_index_map.h>
#include <CGAL/graph_traits_Dual_Arrangement_2.h>
#include <climits> // Needed because BGL forgot it.
#include <boost/graph/dijkstra_shortest_paths.hpp>
#endif

// this base divide & conquer algorithm splits the input into 2 groups,
// calculates the result over the 2 groups, and then merges the results like
// this:
// - overlays the maps, and set at each vertex, edge and face the data from
//   the 2 maps
// - foreach vertex, decide from the 2 envelopes data, which is the envelope
//   of all surfaces
// - foreach edge, do the same, using the Resolver class. an edge can split to
//   a constant number of parts, each with its own envelope data.
// - foreach face, do the same as for edges, using the Resolver class.
//   a face can split to a number of sub-faces that are linear with the face's
//   size, each with its own envelope data.
// - remove edges between faces with the same envelope data, which do not
//   contribute to the shape of the envelope (i.e. have the same envelope data
//   as their adjacent faces)
// - remove unneccessary vertices of two kinds:
//   a. vertices which have degree 2, the 2 incident edges can be geometrically
//      merged, and has the same envelope data as both these edges
//   b. isolated vertices which have the same envelope data as their incident
//      face

// the algorithm deals with some degenerate input including:
// 1. more than one surface on faces, edges, vertices
//    (the minimization diagram should also support this)
// 2. all degenerate cases in 2d (the minimization diagram model is
//    responsible for)

// some degenerate cases in the responsibility of the geometric traits
// 1. overlapping surfaces
// 2. a vertical surface that contributes only edge (or edges) to the envelope

namespace CGAL {

// The algorithm has 5 template parameters:
// 1. EnvelopeTraits_3        - the geometric traits class
// 2. MinimizationDiagram_2   - the type of the output, which is an arrangement
//                              with additional information (list of surfaces)
//                              in vertices, edges & faces
// 3. EnvelopeResolver_3      - part of the algorithm that solves the shape of
//                              the envelope between 2 surfaces over a feature
//                              of the arrangement
// 4. Overlay_2               - overlay of 2 MinimizationDiagram_2
template <class EnvelopeTraits_3, 
          class MinimizationDiagram_2,
          class EnvelopeResolver_3 =
            Envelope_element_visitor_3<EnvelopeTraits_3, MinimizationDiagram_2>,
          class Overlay_2 = Envelope_overlay_2<MinimizationDiagram_2> >
class Envelope_divide_and_conquer_3
{
public:
  typedef EnvelopeTraits_3                              Traits;
  typedef typename Traits::Surface_3                    Surface_3;
  typedef typename Traits::Xy_monotone_surface_3        Xy_monotone_surface_3;

  typedef MinimizationDiagram_2                         Minimization_diagram_2;

  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;

  typedef EnvelopeResolver_3                            Envelope_resolver;

  typedef Envelope_divide_and_conquer_3<Traits, Minimization_diagram_2,
                                        EnvelopeResolver_3,
                                        Overlay_2>      Self;

protected:

  typedef typename Minimization_diagram_2::Halfedge_const_iterator
    Halfedge_const_iterator;
  typedef typename Minimization_diagram_2::Halfedge_handle
    Halfedge_handle;
  typedef typename Minimization_diagram_2::Halfedge_iterator
    Halfedge_iterator;
  typedef typename Minimization_diagram_2::Face_handle
    Face_handle;
  typedef typename Minimization_diagram_2::Edge_iterator
    Edge_iterator;
  typedef typename Minimization_diagram_2::Face_iterator
    Face_iterator;
  typedef typename Minimization_diagram_2::Vertex_handle
    Vertex_handle;
  typedef typename Minimization_diagram_2::Vertex_iterator
    Vertex_iterator;
  typedef typename Minimization_diagram_2::Ccb_halfedge_circulator
    Ccb_halfedge_circulator;
  typedef typename Minimization_diagram_2::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;
  typedef typename Minimization_diagram_2::Outer_ccb_iterator 
    Outer_ccb_iterator;
  typedef typename Minimization_diagram_2::Inner_ccb_iterator 
    Inner_ccb_iterator;

  typedef Arr_observer<Minimization_diagram_2>          Md_observer;
  typedef typename Minimization_diagram_2::Dcel::Dcel_data_iterator
    Envelope_data_iterator;

#ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER
  typedef CGAL::Dual<Minimization_diagram_2>
    Dual_Minimization_diagram_2;
#endif

public:
  // c'tor
  Envelope_divide_and_conquer_3(Envelope_type type = ENVELOPE_LOWER)
  {    
    // Allocate the traits.
    m_geom_traits = new Traits;
    m_own_traits = true;

    // Allocate the Envelope resolver with our traits
    resolver = new Envelope_resolver(m_geom_traits, type);

    m_is_lower = ((type == ENVELOPE_LOWER) ? true : false);
  }

  Envelope_divide_and_conquer_3(const Traits* geom_traits,
                                Envelope_type type = ENVELOPE_LOWER)
  {
    // Set the traits.
    m_geom_traits = geom_traits;
    m_own_traits = false;

    // Allocate the Envelope resolver with our traits
    resolver = new Envelope_resolver(m_geom_traits, type);

    m_is_lower = ((type == ENVELOPE_LOWER) ? true : false);
  }
  
  // virtual destructor.
  virtual ~Envelope_divide_and_conquer_3()
  {
    // Free the traits object, if necessary.
    if (m_own_traits)
      delete m_geom_traits;

    // Free the resolver
    delete resolver;
  }

  // compute the envelope of surfaces in 3D, using the default arbitrary
  // dividor
  template <class SurfaceIterator>
  void construct_lu_envelope(SurfaceIterator begin, SurfaceIterator end,
                             Minimization_diagram_2& result)
  {
    Envelope_3::Arbitrary_dividor dividor;
    construct_lu_envelope(begin, end, result, dividor);
  }
    

  // compute the envelope of surfaces in 3D using the given set dividor
  template <class SurfaceIterator, class SetDividor>
  void construct_lu_envelope(SurfaceIterator begin, SurfaceIterator end, 
                             Minimization_diagram_2& result,
                             SetDividor& dividor)
  {
    if (begin == end)
    {
      return; // result is empty
    }
    
    // make the general surfaces xy-monotone
    std::list<Xy_monotone_surface_3> xy_monotones;
    for (; begin != end; ++begin)
      m_geom_traits->
        make_xy_monotone_3_object()(*begin, m_is_lower,
                                    std::back_inserter(xy_monotones));

    // recursively construct the envelope of the xy-monotone parts
    construct_lu_envelope_xy_monotones(xy_monotones.begin(), 
                                       xy_monotones.end(), result, dividor);

    CGAL_assertion(is_envelope_valid(result));     
  }

  // compute the envelope of xy-monotone surfaces in 3D,
  // using the default arbitrary dividor
  template <class SurfaceIterator>
  void construct_envelope_xy_monotone(SurfaceIterator begin,
                                      SurfaceIterator end,
                                      Minimization_diagram_2& result)
  {
    Envelope_3::Arbitrary_dividor dividor;
    construct_envelope_xy_monotone(begin, end, result, dividor);
  }

  // compute the envelope of xy-monotone surfaces in 3D using the given
  // set dividor
  template <class SurfaceIterator, class SetDividor>
  void construct_envelope_xy_monotone(SurfaceIterator begin,
                                      SurfaceIterator end,
                                      Minimization_diagram_2& result,
                                      SetDividor& dividor)
  {
    if (begin == end)
      return; // result is empty

    // recursively construct the envelope of the xy-monotone parts
    construct_lu_envelope_xy_monotones(begin, end, result, dividor);
    CGAL_assertion(is_envelope_valid(result));
  }

  /*! Access the traits object. */
  const Traits* get_traits() const
  {
    return m_geom_traits;
  }

  void reset()
  {
    resolver->reset();
  }

protected:

  // compute the envelope of xy-monotone surfaces in 3D 
  template <class SurfaceIterator, class SetDividor>
  void construct_lu_envelope_xy_monotones(SurfaceIterator begin,
                                          SurfaceIterator end,
                                          Minimization_diagram_2& result,
                                          SetDividor& dividor)
  {
    if (begin == end)
    {
      return; // result is empty
    }

    SurfaceIterator first = begin++;

    if (begin == end)
    {
      // only one surface is in the collection. insert it the result
      Xy_monotone_surface_3& surf = *first;
      
      deal_with_one_surface(surf, result);
      return;      
    }
   
    // divide the surfaces into 2 groups (insert surface to each group
    // alternately)
    // Efi: this copy is redundant. It is sufficient to determine the range
    std::list<Xy_monotone_surface_3> group1, group2;
    dividor(first, end,
            std::back_inserter(group1), std::back_inserter(group2));
    
    // recursively calculate the LU_envelope of the 2 groups
    Minimization_diagram_2 result1(m_geom_traits), result2(m_geom_traits);
    construct_lu_envelope_xy_monotones(group1.begin(), group1.end(),
                                       result1, dividor);
    construct_lu_envelope_xy_monotones(group2.begin(), group2.end(),
                                       result2, dividor);
        
    // merge the results:
    merge_envelopes(result1, result2, result);

    result1.clear();
    result2.clear();

    CGAL_assertion(is_envelope_valid(result));
  }

  void deal_with_one_surface(Xy_monotone_surface_3& surf,
                             Minimization_diagram_2& result)
  {
    typedef std::list<Object>                            Boundary_list;
    typedef std::pair<X_monotone_curve_2, Oriented_side> Boundary_xcurve;
    typedef Boundary_list::iterator                      Boundary_iterator;

    Boundary_list     boundary;
    m_geom_traits->
      construct_projected_boundary_2_object()(surf,
                                              std::back_inserter(boundary));
    
    if (boundary.empty())
    {
      //one infinite surface
      CGAL_assertion_msg(result.number_of_faces() == 1,
                         "In the beginning there should be only one face");
      result.faces_begin()->set_data(surf);
      return;
    }

    for (Boundary_iterator boundary_it = boundary.begin();
         boundary_it != boundary.end();
         ++boundary_it)
    {
      const Object& obj = *boundary_it;
      Boundary_xcurve boundary_cv;
      if (assign(boundary_cv, obj))
      {
        Oriented_side side = boundary_cv.second;
        Halfedge_handle he =
          insert_non_intersecting_curve(result, boundary_cv.first);
        
        if (side == ON_ORIENTED_BOUNDARY)
        {
          // vertical xy-surface
          he->face()->set_no_data();
          he->twin()->face()->set_no_data();
          
          continue;
        }

        if (he->face() != he->twin()->face())
        {
          // new face created.
          // 'he' is directed from left to right, so the face to the left
          // of 'he'  is above 'cv.
          Face_handle f;
          if (side == ON_NEGATIVE_SIDE) // the surface is below cv.
          {
            f = he->twin()->face();
            f->set_data(surf);
            he->face()->set_no_data();
          }
          else
          {
            CGAL_assertion(side == ON_POSITIVE_SIDE);
            f = he->face();
            f->set_data(surf);
            he->twin()->face()->set_no_data();
          }

          // init auxiliary data for f and its boundarys. 
          for(Outer_ccb_iterator ocit = f->outer_ccbs_begin();
              ocit != f->outer_ccbs_end(); ocit++){
            Ccb_halfedge_circulator face_hec = *ocit;
            Ccb_halfedge_circulator face_hec_begin = face_hec;
            do 
              {
                face_hec->set_is_equal_data_in_face(true);
                face_hec->set_has_equal_data_in_face(true);
                face_hec->set_has_equal_data_in_target_and_face(true);
                
                face_hec->twin()->set_is_equal_data_in_face(false);
                face_hec->twin()->set_has_equal_data_in_face(false);
                face_hec->twin()->set_has_equal_data_in_target_and_face(false);
                
                ++face_hec;
              } 
            while(face_hec != face_hec_begin);
          }
          for(Outer_ccb_iterator icit = f->inner_ccbs_begin();
              icit != f->inner_ccbs_end(); icit++){
            Ccb_halfedge_circulator face_hec = *icit;
            Ccb_halfedge_circulator face_hec_begin = face_hec;
            do 
              {
                face_hec->set_is_equal_data_in_face(true);
                face_hec->set_has_equal_data_in_face(true);
                face_hec->set_has_equal_data_in_target_and_face(true);
                
                face_hec->twin()->set_is_equal_data_in_face(false);
                face_hec->twin()->set_has_equal_data_in_face(false);
                face_hec->twin()->set_has_equal_data_in_target_and_face(false);
                
                ++face_hec;
              } 
            while(face_hec != face_hec_begin);
          }
        }
      }
      else
      {
        // the xy-surface is an isolated point
        Point_2 p;
        CGAL_assertion(assign(p, obj));
        assign(p, obj);
        insert_point(result, p);
      }
    }

    // update information in all the edges & vertices to indicate that
    // this surface is the envelope
    Halfedge_iterator hi = result.halfedges_begin();
    for (; hi != result.halfedges_end(); ++hi)
    {
      hi->set_data(surf);
      // since all the edges & vertices have their envelope data equal to the
      // current surface, we can set is/has equal_data_in_target of all
      // halfedges to true
      hi->set_is_equal_data_in_target(true);
      hi->set_has_equal_data_in_target(true);
    }

    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
    {
      vi->set_data(surf);
      if (vi->is_isolated())
      {
        // update the is/has equal_data_in_face flags according to the face data
        bool equal_data = !vi->face()->has_no_data();
        vi->set_is_equal_data_in_face(equal_data);
        vi->set_has_equal_data_in_face(equal_data);
      }
    }
  }
  

public:
  void merge_envelopes(Minimization_diagram_2& result1,
                       Minimization_diagram_2& result2,
                       Minimization_diagram_2& result)
  {
    // overlay the 2 arrangements
   
    Overlay_2 overlay;
    overlay(result1, result2, result);
  
    CGAL_expensive_assertion_msg(is_valid(result),
                                 "after overlay result is not valid");
       
    // make sure the aux flags are correctly set by the overlay
    //CGAL_assertion(verify_aux_flags(result));

    // for each face, edge and vertex in the result, should calculate
    // which surfaces are on the envelope
    // a face can be cut, or faces can be merged.
   
    // now the minimization diagram might change - we need to keep data in the
    // edges, when they're split
    Keep_edge_data_observer edge_observer(result, this);

    // compute the surface on the envelope for each edge
    // edge can be split as surfaces can intersect (or touch) over it
    std::list<Halfedge_handle> edges_to_resolve;
    Edge_iterator ei = result.edges_begin();
    for (; ei != result.edges_end(); ++ei)
    {
      Halfedge_handle hh = ei;
      // there must be data from at least one map, because all the surfaces
      // are continous
      if (!get_aux_is_set(hh, 0) || !get_aux_is_set(hh, 1))
        continue;
      CGAL_assertion(get_aux_is_set(hh, 0));
      CGAL_assertion(get_aux_is_set(hh, 1));
      CGAL_assertion(!aux_has_no_data(hh, 1) || !aux_has_no_data(hh, 0));
      if (aux_has_no_data(hh, 0) && !aux_has_no_data(hh, 1))
      {
        hh->set_decision(DAC_DECISION_SECOND);
        hh->twin()->set_decision(DAC_DECISION_SECOND);
        continue;
      }
      else if (!aux_has_no_data(hh, 0) && aux_has_no_data(hh, 1))
      {
        hh->set_decision(DAC_DECISION_FIRST);
        hh->twin()->set_decision(DAC_DECISION_FIRST);
        continue;
      }

      bool should_resolve = true;
#ifdef CGAL_ENVELOPE_SAVE_COMPARISONS
        if (hh->get_has_equal_aux_data_in_face(0) &&
            hh->get_has_equal_aux_data_in_face(1))
          should_resolve = false;

        if (hh->twin()->get_has_equal_aux_data_in_face(0) &&
            hh->twin()->get_has_equal_aux_data_in_face(1))
          should_resolve = false;
#endif

      // we collect the edges in a list to deal afterwards, because the resolve
      // can split edges, and destroy the iterator
      if (should_resolve)
        edges_to_resolve.push_back(hh);
    }
    // now deal with the edges
    typename std::list<Halfedge_handle>::iterator li;
    for (li = edges_to_resolve.begin(); li != edges_to_resolve.end(); ++li)
    {
      resolver->resolve(*li, result);
     
    }
    edges_to_resolve.clear();
    
    // decompose the result, to have faces without holes
   /* decompose(result);
    CGAL_expensive_assertion_msg(result.is_valid(), 
                       "after decomposition result is not valid");*/

    // compute the surface on the envelope for each face,
    // splitting faces if needed

    std::list<Face_handle> faces_to_split;

#ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER
    // we traverse the faces of result in BFS order to maximize the
    // efficiency gain by the conclusion mechanism of
    // compare_distance_to_envelope results
    // Create a mapping of result faces to indices.
    CGAL::Arr_face_index_map<Minimization_diagram_2> index_map(result);

    // Perform breadth-first search from the unbounded face, and use the BFS
    // visitor to associate each arrangement face with its discover time.
    Faces_order_bfs_visitor<CGAL::Arr_face_index_map<Minimization_diagram_2> >
      bfs_visitor(index_map, faces_to_split, this);
    Face_handle first_face = result.faces_begin();
    /* if (result.number_of_faces() > 1)
     *   first_face = ++(result.faces_begin());
     */

    boost::breadth_first_search(Dual_Minimization_diagram_2(result),
                                first_face,
                                boost::vertex_index_map(index_map).
                                visitor(bfs_visitor));
    index_map.detach();
#else
    // traverse the faces in arbitrary order  
    Face_iterator fi = result.faces_begin();
    for (; fi != result.faces_end(); ++fi)
    {
      Face_handle fh = fi;
        // if a surface of one map doesn't exist, then we set the second surface
        if (aux_has_no_data(fh, 0) && !aux_has_no_data(fh, 1))
        {
          fh->set_decision(DAC_DECISION_SECOND);
          continue;
        }
        else if (aux_has_no_data(fh, 0) && aux_has_no_data(fh, 1))
        {
          fh->set_decision(EQUAL);
          fh->set_no_data();

          continue;
        }
        else if (!aux_has_no_data(fh, 0) && aux_has_no_data(fh, 1))
        {
          fh->set_decision(DAC_DECISION_FIRST);
          continue;
        }

        // here, we have both surfaces.
        // we save the face in a list for a later treatment, because the
        // face can change and destroy the iterator
        faces_to_split.push_back(fh);
      }
#endif

    deal_with_faces_to_split(faces_to_split, result);

//    #ifndef CGAL_ENVELOPE_SAVE_COMPARISONS
//      hi = result.halfedges_begin();
//      for (; hi != result.halfedges_end(); ++hi, ++hi)
//      {
//        if (!hi->is_decision_set())
//          resolver->resolve(hi, result);
//      }
//    #endif
    
    // detach the edge_observer from result, since no need for it anymore
    edge_observer.detach();

    // compute the surface on the envelope for each vertex
    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      if (vh->is_decision_set())
        continue;
      // there must be data from at least one map, because all the surfaces
      // are continous
      CGAL_assertion(get_aux_is_set(vh, 0));
      CGAL_assertion(get_aux_is_set(vh, 1));
      CGAL_assertion(!aux_has_no_data(vh, 1) || !aux_has_no_data(vh, 0));
      if (aux_has_no_data(vh, 0) && !aux_has_no_data(vh, 1))
      {
        vh->set_decision(DAC_DECISION_SECOND);
        continue;
      }
      else if (!aux_has_no_data(vh, 0) && aux_has_no_data(vh, 1))
      {
        vh->set_decision(DAC_DECISION_FIRST);
        continue;
      }
      resolver->resolve(vh);
     
    }

    CGAL_expensive_assertion_msg(result.is_valid(),
                                 "after resolve result is not valid");

    // make sure that aux_source and decision are set at all features
    // after all resolvings
    CGAL_assertion(check_resolve_was_ok(result));
 
    // make sure the aux flags are correctly after all resolvings
    //CGAL_assertion(verify_aux_flags(result));

    // finally, remove unneccessary edges, between faces  with the same surface
    // (and which are not degenerate)
  
    remove_unneccessary_edges(result);
    CGAL_expensive_assertion_msg(result.is_valid(), 
                       "after remove edges result is not valid");

    // also remove unneccessary vertices (that were created in the process of
    // vertical decomposition but the vertical edge was removed)
    remove_unneccessary_vertices(result);
    CGAL_expensive_assertion_msg(result.is_valid(), 
                       "after remove vertices result is not valid");

    // update is_equal_data and has_equal_data of halfedge->face and
    // vertex->face relations, according to the decision, and the aux
    // similar flags
    update_flags(result);

    // update the envelope surfaces according to the decision and the aux
    // surfaces in aux source
    update_envelope_surfaces_by_decision(result);
    
    // make sure that all the flags are correctly set on the envelope result
    //CGAL_assertion(verify_flags(result));
    CGAL_expensive_assertion_msg(is_valid(result),
                                 "after merge result is not valid");
  }


protected:
 
  void deal_with_faces_to_split(std::list<Face_handle>& faces_to_split,
                                Minimization_diagram_2& result)
  {
    // for each face in faces_to_split, find the intersection over the face,
    // and split the face
    typename std::list<Face_handle>::iterator li;
    for (li = faces_to_split.begin(); li != faces_to_split.end(); ++li)
      resolver->resolve(*li, result);
    faces_to_split.clear();
  }

  template <class InputIterator>
  bool is_equal_data(const InputIterator& begin1,
                     const InputIterator& end1,
                     const InputIterator& begin2,
                     const InputIterator& end2)
  {
    // insert the input data objects into a set
    std::set<Xy_monotone_surface_3> first(begin1, end1);
    std::set<Xy_monotone_surface_3> second(begin2, end2);

    if (first.size() != second.size())
      return false;

    return (first == second);
  }

  // todo: should remove the uses of this method from this class
  template <class InputIterator>
  bool has_equal_data(const InputIterator& begin1,
                      const InputIterator& end1,
                      const InputIterator& begin2,
                      const InputIterator& end2)
  {
    // insert the input data objects into a set
    std::set<Xy_monotone_surface_3> first(begin1, end1);
    std::set<Xy_monotone_surface_3> second(begin2, end2);
    std::list<Xy_monotone_surface_3> intersection;
    std::set_intersection(first.begin(), first.end(),
                          second.begin(), second.end(),
                          std::back_inserter(intersection));
    return (intersection.size() > 0);
    return true;
  }
  // todo: should remove the uses of this method from this class
  template <class FeatureHandle1, class FeatureHandle2>
  bool has_equal_aux_data(unsigned int id, FeatureHandle1 fh1,
                          FeatureHandle2 fh2)
  {
    Envelope_data_iterator begin1, end1, begin2, end2;
    get_aux_data_iterators(id, fh1, begin1, end1);
    get_aux_data_iterators(id, fh2, begin2, end2);
    bool has_eq = has_equal_data(begin1, end1, begin2, end2);
    return has_eq;
  }

  // Remove unneccessary edges, between faces with the same surface
  // (and which are not degenerate)
  void remove_unneccessary_edges(Minimization_diagram_2& result)
  {
    // collect all those edges in this list, and remove them all at the end
    // (thus, not destroying the iterator)
    std::list<Halfedge_handle> edges;
    Halfedge_iterator hi = result.halfedges_begin();
    for (; hi != result.halfedges_end(); ++hi, ++hi)
    {
      Halfedge_handle hh = hi;
      if (can_remove_edge(hh))
        edges.push_back(hh);
    }
    
    for (typename std::list<Halfedge_handle>::iterator ci = edges.begin();
         ci != edges.end(); ++ci)
    {
      // if the endpoints become isolated after the removal we need to remove
      // them if they have the same data as the edge
      Halfedge_handle h = *ci;
      Vertex_handle src = h->source();
      Vertex_handle trg = h->target();
      CGAL_assertion_code(Face_handle h_face = h->face());
      bool remove_src = can_remove_edge_target(h->twin());
      bool remove_trg = can_remove_edge_target(h);
      bool src_is_equal_0 = (h->twin()->get_is_equal_aux_data_in_face(0) &&
                             h->twin()->get_is_equal_aux_data_in_target(0));
      bool src_is_equal_1 = (h->twin()->get_is_equal_aux_data_in_face(1) &&
                             h->twin()->get_is_equal_aux_data_in_target(1));
      bool trg_is_equal_0 = (h->get_is_equal_aux_data_in_face(0) &&
                             h->get_is_equal_aux_data_in_target(0));
      bool trg_is_equal_1 = (h->get_is_equal_aux_data_in_face(1) &&
                             h->get_is_equal_aux_data_in_target(1));
      bool src_has_equal_0 =
        h->twin()->get_has_equal_aux_data_in_target_and_face(0);
      bool src_has_equal_1 =
        h->twin()->get_has_equal_aux_data_in_target_and_face(1);
      bool trg_has_equal_0 = h->get_has_equal_aux_data_in_target_and_face(0);
      bool trg_has_equal_1 = h->get_has_equal_aux_data_in_target_and_face(1);

      /* A vertex at an open boundary is removed once it becomes redundant
       * regardless of the boolean values passed as the 2nd and 3rd argument
       * to the remove_edge() member function.
       */
      if (src->is_at_open_boundary()) remove_src = true;
      if (trg->is_at_open_boundary()) remove_trg = true;
      result.remove_edge(*ci, remove_src, remove_trg);
      // otherwise, we should make sure, they will not be removed
      // the first check is needed since if the vertex was removed, then the
      // handle is invalid
      if (!remove_src && src->is_isolated()) {
        // to be precise we copy from the halfedge-face and halfedge-target
        // relations
        src->set_is_equal_aux_data_in_face(0, src_is_equal_0);
        src->set_is_equal_aux_data_in_face(1, src_is_equal_1);
        // todo: the has_equal flags should be updated also
        // make sure h_face is also src face
        CGAL_assertion(h_face == src->face());
        // CGAL_assertion(src_has_equal_0 ==
        // has_equal_aux_data(0, src, h_face));
        // CGAL_assertion(src_has_equal_1 ==
        // has_equal_aux_data(1, src, h_face));
        src->set_has_equal_aux_data_in_face(0, src_has_equal_0);
        src->set_has_equal_aux_data_in_face(1, src_has_equal_1);
      }
      if (!remove_trg && trg->is_isolated())
      {
        trg->set_is_equal_aux_data_in_face(0, trg_is_equal_0);
        trg->set_is_equal_aux_data_in_face(1, trg_is_equal_1);
        // make sure h_face is also trg face
        CGAL_assertion(h_face == trg->face());
        // CGAL_assertion(trg_has_equal_0 ==
        // has_equal_aux_data(0, trg, h_face));
        // CGAL_assertion(trg_has_equal_1 ==
        // has_equal_aux_data(1, trg, h_face));
        trg->set_has_equal_aux_data_in_face(0, trg_has_equal_0);
        trg->set_has_equal_aux_data_in_face(1, trg_has_equal_1);
      }
    }
  }

  template <class FeatureHandle>
  void get_aux_data_iterators(unsigned int id, FeatureHandle fh,
                              Envelope_data_iterator& begin,
                              Envelope_data_iterator& end)
  {
    Halfedge_handle h;
    Vertex_handle v;
    Face_handle f;

    const Object& o = fh->get_aux_source(id);
    CGAL_assertion(!o.is_empty());

    // aux source of a face must be a face!
    // aux source of a halfedge can be face or halfedge
    // aux source of a vertex can be face, halfedge or vertex

    // this is why we start with a check for a face, then halfedge
    // and last vertex

    if (assign(f, o))
    {
      begin = f->begin_data();
      end = f->end_data();
    }
    else if (assign(h, o))
    {
      begin = h->begin_data();
      end = h->end_data();
    }
    else
    {
      CGAL_assertion_code(bool b = )
      assign(v, o);
      CGAL_assertion(b);
      begin = v->begin_data();
      end = v->end_data();
    }
  }
  
  // check if we can remove the edge from the envelope
  // this can be done if the envelope surfaces on the edge are the same as
  // the envelope surfaces on both sides of the edge
  // (or if the edge is fake, i.e. created in the vd process)
  bool can_remove_edge(Halfedge_handle hh)
  {
    Face_handle f1 = hh->face(), f2 = hh->twin()->face();
      
    // we check if the decision done on the edge is equal to the decision
    // done on the faces. if not, then the envelope surfaces must differ
    CGAL_assertion(hh->is_decision_set() && f1->is_decision_set() &&
                   f2->is_decision_set());
    if (hh->get_decision() != f1->get_decision() ||
        hh->get_decision() != f2->get_decision())
    {
      return false;
    }

    // now, check the equality of the surfaces list according to the decision
    CGAL::Dac_decision decision = hh->get_decision();
    bool equal_first = (hh->get_is_equal_aux_data_in_face(0) &&
                        hh->twin()->get_is_equal_aux_data_in_face(0));
    bool equal_second = (hh->get_is_equal_aux_data_in_face(1) &&
                         hh->twin()->get_is_equal_aux_data_in_face(1));

    if (decision == DAC_DECISION_FIRST)
      return equal_first;
   
    if (decision == DAC_DECISION_SECOND)
      return equal_second;
    
    return (equal_first && equal_second);
  }


  // check if can remove the edge's target if we remove the edge
  // this means that the target has the same envelope information as the edge
  bool can_remove_edge_target(Halfedge_handle h)
  {
    // \todo Use new design
    /* The code below uses the is_at_open_boundary.
       The comment from the calling function says:
       "if the endpoints become isolated after the removal we need to remove
        them if they have the same data as the edge."

       When I tried to use the following code I got a Segmentation Fault when
       trying to compute power diagram:
       
       if ((v->parameter_space_in_x() != ARR_INTERIOR) ||
          (v->parameter_space_in_y() != ARR_INTERIOR))
          return false;

     */
    Vertex_handle v = h->target();
    if (v->is_at_open_boundary())
      return false;


    /*if (v->get_is_fake() && !v->is_decision_set())
      return true;
   

    if (h->get_is_fake() && !h->is_decision_set())
    {
      h->set_decision(h->face()->get_decision());
      h->twin()->set_decision(h->get_decision());
    }*/
    CGAL_assertion(v->is_decision_set());
    CGAL_assertion(h->is_decision_set());

    // if the decision done on the vertex and edge are different,
    // the envelope differs too.
    if (h->get_decision() != v->get_decision())
      return false;

    // now, check the equality of the surfaces list according to the decision

    CGAL::Dac_decision decision = h->get_decision();
    bool equal_first = (h->get_is_equal_aux_data_in_target(0));
    bool equal_second = (h->get_is_equal_aux_data_in_target(1));

    if (decision == DAC_DECISION_FIRST)
      return equal_first;
   
    if (decision == DAC_DECISION_SECOND)
      return equal_second;
    
    return (equal_first && equal_second);
  }
  
  // check if we can remove an isolated vertex from the envelope
  // this can be done if the envelope surfaces on the vertex are the same as
  // the envelope surfaces on its incident face
  bool can_remove_isolated_vertex(Vertex_handle vh)
  {
    Face_handle f = vh->face();
    CGAL_assertion(vh->is_decision_set() && f->is_decision_set());
    // if the decision done on the vertex and face are different,
    // the envelope differs too.
    if (vh->get_decision() != f->get_decision())
      return false;

    // now, check the equality of the surfaces list according to the decision
    CGAL::Dac_decision decision = vh->get_decision();

    bool equal_first = (vh->get_is_equal_aux_data_in_face(0));
    bool equal_second = (vh->get_is_equal_aux_data_in_face(1));

    if (decision == DAC_DECISION_FIRST)
      return equal_first;
    
    if (decision == DAC_DECISION_SECOND)
      return equal_second;
    
    return (equal_first && equal_second);
  }

  // check if we can remove a (non-isolated) vertex from the envelope
  // we check only from the combinatoric aspect of the envelope, not from
  // the geometric point of view (i.e. if the curves can merge)
  // this can be done if the envelope surfaces on the vertex are the same as
  // the envelope surfaces on its 2 incident halfedges
  bool combinatorically_can_remove_vertex(Vertex_handle vh)
  {
    Halfedge_around_vertex_circulator hec1 = vh->incident_halfedges();
    Halfedge_around_vertex_circulator hec2 = hec1++;
    Halfedge_handle he1 = hec1, he2 = hec2;
    CGAL_assertion(he1 != he2);
    CGAL_assertion(he1->is_decision_set() && he2->is_decision_set());
    
    /* if (vh->get_is_fake()) {
     *   CGAL_assertion(he1->get_decision() == he2->get_decision());
     *   return true;
     * }
     */
    
    CGAL_assertion(vh->is_decision_set());
    // if the decision done on the vertex and its incident halfedges are
    // different, the envelope differs too.
    CGAL_assertion(vh == he1->target() && vh == he2->target());
    if (vh->get_decision() != he1->get_decision() ||
        vh->get_decision() != he2->get_decision())
      return false;

    // now, check the equality of the surfaces list according to the decision
    CGAL::Dac_decision decision = vh->get_decision();
    bool equal_first = (he1->get_is_equal_aux_data_in_target(0) &&
                        he2->get_is_equal_aux_data_in_target(0));
    bool equal_second = (he1->get_is_equal_aux_data_in_target(1) &&
                         he2->get_is_equal_aux_data_in_target(1));

    if (decision == DAC_DECISION_FIRST)
      return equal_first;
    
    if (decision == DAC_DECISION_SECOND)
      return equal_second;

    return (equal_first && equal_second);
  }
  
  // Remove unneccessary vertices, which have degree 2, and the 2 curves 
  // can be merged
  // (and which are not degenerate)
  void remove_unneccessary_vertices(Minimization_diagram_2& result)
  {
    // we have 2 types of unneccessary vertices: those with degree 2 (that
    // satisfy all the conditions below), and isolated vertices that have the
    // same envelope information as the face they're contained in.
    
    // all the vertices that don't have their data set, are those vertices
    // on vertical edges, created in the decomposition process,
    // and are not neccessary

    // also those vertices with degree 2, that can merge their 2 edges and 
    // with same data as both these edges, can be removed
    
    // collect all vertices candidate to remove in this list,
    // and remove the correct ones at the end
    // (thus, not destroying the iterator)
    std::list<Vertex_handle> candidates_to_remove;
    std::list<Vertex_handle> isolated_to_remove;

    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      if (!vh->is_decision_set() || vh->degree() == 2)
        candidates_to_remove.push_back(vh);
      else if (vh->is_isolated() && can_remove_isolated_vertex(vh))
        isolated_to_remove.push_back(vh);
    }

    typename Traits::Merge_2 curves_merge = m_geom_traits->merge_2_object();
    typename Traits::Are_mergeable_2 curves_can_merge = 
      m_geom_traits->are_mergeable_2_object();

    // check the candidates and remove if necessary
    typename std::list<Vertex_handle>::iterator ci;
    for (ci = candidates_to_remove.begin(); 
        ci != candidates_to_remove.end(); ++ci)
    {
      Vertex_handle vh = *ci;
      CGAL_assertion(vh->degree() == 2);
      
      // we can remove this vertex only if the data on its halfedges is the 
      // same
      if (!combinatorically_can_remove_vertex(vh))
        continue;
                     
      // merge the edges, if geometrically possible (if data on vertex is not
      // set, then it must be geometrically possible)
      Halfedge_around_vertex_circulator hec1 = vh->incident_halfedges();
      Halfedge_around_vertex_circulator hec2 = hec1++;
      Halfedge_handle he1 = hec1, he2 = hec2;
      

      const X_monotone_curve_2& a = he1->curve(), b = he2->curve();
      CGAL_assertion(vh->is_decision_set() || curves_can_merge(a,b));
      if (vh->is_decision_set() && !curves_can_merge(a,b))
        continue;

      X_monotone_curve_2 c;
      curves_merge(a,b,c);

      // the decisions on he1 and he2 were the same, so the decision on 
      // the edge that will be left after the merge will be ok
      // but we need to take care of the bool flags of the target relation
      // of the edge that will be left
      // he1 and he2 both point to the removed vertex, so each should copy
      // the other's twin flags. the flags in the twins don't change
      //    he1            he2
      // x----------->x<------------x
      //  <----------- ------------>
      //   he1->twin()   he2->twin()
      he1->set_is_equal_aux_data_in_target
        (0, he2->twin()->get_is_equal_aux_data_in_target(0));
      he1->set_is_equal_aux_data_in_target
        (1, he2->twin()->get_is_equal_aux_data_in_target(1));
      he1->set_has_equal_aux_data_in_target
        (0, he2->twin()->get_has_equal_aux_data_in_target(0));
      he1->set_has_equal_aux_data_in_target
        (1, he2->twin()->get_has_equal_aux_data_in_target(1));
      he1->set_has_equal_aux_data_in_target_and_face
        (0, he2->twin()->get_has_equal_aux_data_in_target_and_face(0));
      he1->set_has_equal_aux_data_in_target_and_face
        (1, he2->twin()->get_has_equal_aux_data_in_target_and_face(1));

      he2->set_is_equal_aux_data_in_target
        (0, he1->twin()->get_is_equal_aux_data_in_target(0));
      he2->set_is_equal_aux_data_in_target
        (1, he1->twin()->get_is_equal_aux_data_in_target(1));
      he2->set_has_equal_aux_data_in_target
        (0, he1->twin()->get_has_equal_aux_data_in_target(0));
      he2->set_has_equal_aux_data_in_target
        (1, he1->twin()->get_has_equal_aux_data_in_target(1));
      he2->set_has_equal_aux_data_in_target_and_face
        (0, he1->twin()->get_has_equal_aux_data_in_target_and_face(0));
      he2->set_has_equal_aux_data_in_target_and_face
        (1, he1->twin()->get_has_equal_aux_data_in_target_and_face(1));

      // order of halfedges for merge doesn't matter
#if !defined(CGAL_NO_ASSERTIONS)
      Halfedge_handle new_edge =
#endif
      result.merge_edge(he1, he2 ,c);
      CGAL_assertion(new_edge->is_decision_set());

      CGAL_expensive_assertion_msg(result.is_valid(), 
                                   "after remove vertex result is not valid");
    }

    // remove isolated vertices

    typename std::list<Vertex_handle>::iterator li;
    for (li = isolated_to_remove.begin(); li != isolated_to_remove.end(); ++li)
    {
      Vertex_handle vh = *li;
      CGAL_assertion(vh->is_isolated());
      CGAL_assertion(can_remove_isolated_vertex(vh));

      result.remove_isolated_vertex(vh);
    }
  }

  template <class FeatureHandle>
  void update_envelope_surfaces_by_decision(FeatureHandle fh)
  {
    CGAL::Dac_decision decision = fh->get_decision();

    Halfedge_handle h;
    Vertex_handle v;
    Face_handle f;

    if (decision == DAC_DECISION_FIRST || decision == DAC_DECISION_BOTH)
    {
      Envelope_data_iterator begin, end;
      get_aux_data_iterators(0, fh, begin, end);
      fh->set_data(begin, end);
    }
    if (decision == DAC_DECISION_SECOND || decision == DAC_DECISION_BOTH)
    {
      // copy data from second envelope
      Envelope_data_iterator begin, end;
      get_aux_data_iterators(1, fh, begin, end);

      if (decision == DAC_DECISION_SECOND)
        fh->set_data(begin, end);
      else
        fh->add_data(begin, end);
    }    
  }

  // foreach feature of result, update the envelope surfaces, according
  // to the decision done over this feature, and its aux sources
  void update_envelope_surfaces_by_decision(Minimization_diagram_2& result)
  {
    // vertices
    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
      update_envelope_surfaces_by_decision(vi);
    
    // edges
    Halfedge_iterator hi = result.halfedges_begin();
    for (; hi != result.halfedges_end(); ++hi)
      update_envelope_surfaces_by_decision(hi);
    
    // faces
    Face_iterator fi = result.faces_begin();
    for (; fi != result.faces_end(); ++fi)
      update_envelope_surfaces_by_decision(fi);
  }

  // update the is_equal/has_equal flags of the result envelope
  void update_edge_face_flags(Halfedge_handle h)
  {
    bool is_equal, has_equal;
    is_equal = (h->get_decision() == h->face()->get_decision());
    // has equal can be true even if the decision is not the same,
    // but has same surfaces, i.e. one of the features got DAC_DECISION_BOTH
    // decision, and the other didn't
    has_equal = (h->get_decision() == h->face()->get_decision() ||
                 h->get_decision() == DAC_DECISION_BOTH ||
                 h->face()->get_decision() == DAC_DECISION_BOTH);

    CGAL::Dac_decision decision = h->face()->get_decision();
    bool is_equal_first = (h->get_is_equal_aux_data_in_face(0));
    bool has_equal_first = (h->get_has_equal_aux_data_in_face(0));
    bool is_equal_second = (h->get_is_equal_aux_data_in_face(1));
    bool has_equal_second = (h->get_has_equal_aux_data_in_face(1));
    if (decision == DAC_DECISION_FIRST)
    {
      is_equal &= is_equal_first;
      has_equal &= has_equal_first;
    }
    else if (decision == DAC_DECISION_SECOND)
    {
      is_equal &= is_equal_second;
      has_equal &= has_equal_second;
    }
    else
    {
      is_equal &= (is_equal_first & is_equal_second);
      // we check if the halfedge has a different decision, and if so,

      // we update the flag according to the halfedge decision
      decision = h->get_decision();
      if (decision == DAC_DECISION_FIRST)
      has_equal &= has_equal_first;
      else if (decision == DAC_DECISION_SECOND)
      has_equal &= has_equal_second;
      else      
      has_equal &= (has_equal_first & has_equal_second);
    }

    h->set_is_equal_data_in_face(is_equal);
    h->set_has_equal_data_in_face(has_equal);
  }
  void update_edge_target_flags(Halfedge_handle h)
  {
    bool is_equal, has_equal;

    is_equal = (h->get_decision() == h->target()->get_decision());
    // has equal can be true even if the decision is not the same,
    // but has same surfaces, i.e. one of the features got DAC_DECISION_BOTH
    // decision, and the other didn't
    has_equal = (h->get_decision() == h->target()->get_decision() ||
                 h->get_decision() == DAC_DECISION_BOTH ||
                 h->target()->get_decision() == DAC_DECISION_BOTH);

    CGAL::Dac_decision decision = h->get_decision();
    bool is_equal_first = (h->get_is_equal_aux_data_in_target(0));
    bool has_equal_first = (h->get_has_equal_aux_data_in_target(0));
    bool is_equal_second = (h->get_is_equal_aux_data_in_target(1));
    bool has_equal_second = (h->get_has_equal_aux_data_in_target(1));
    if (decision == DAC_DECISION_FIRST)
    {
      is_equal &= is_equal_first;
      has_equal &= has_equal_first;
    }

    else if (decision == DAC_DECISION_SECOND)
    {
      is_equal &= is_equal_second;
      has_equal &= has_equal_second;
    }
    else
    {
      is_equal &= (is_equal_first & is_equal_second);
      // we check if the vertex has a different decision, and if so,
      // we update the flag according to the vertex decision
      decision = h->target()->get_decision();
      if (decision == DAC_DECISION_FIRST)
      has_equal &= has_equal_first;
      else if (decision == DAC_DECISION_SECOND)
      has_equal &= has_equal_second;

      else      
      has_equal &= (has_equal_first & has_equal_second);
    }
    h->set_is_equal_data_in_target(is_equal);
    h->set_has_equal_data_in_target(has_equal);
  }
  void update_target_face_flags(Halfedge_handle h)
  {
    bool has_equal;
    // has equal can be true even if the decision is not the same,
    // but has same surfaces, i.e. one of the features got DAC_DECISION_BOTH
    // decision, and the other didn't
    has_equal = (h->face()->get_decision() == h->target()->get_decision() ||
                 h->face()->get_decision() == DAC_DECISION_BOTH ||
                 h->target()->get_decision() == DAC_DECISION_BOTH);

    CGAL::Dac_decision decision = h->face()->get_decision();
    bool has_equal_first = (h->get_has_equal_aux_data_in_target_and_face(0));
    bool has_equal_second = (h->get_has_equal_aux_data_in_target_and_face(1));
    if (decision == DAC_DECISION_FIRST)
      has_equal &= has_equal_first;
    else if (decision == DAC_DECISION_SECOND)
      has_equal &= has_equal_second;
    else
    {
      // we check if the vertex has a different decision, and if so,
      // we update the flag according to the vertex decision
      decision = h->target()->get_decision();
      if (decision == DAC_DECISION_FIRST)
        has_equal &= has_equal_first;
      else if (decision == DAC_DECISION_SECOND)
        has_equal &= has_equal_second;

      else      
        has_equal &= (has_equal_first & has_equal_second);
    }
    h->set_has_equal_data_in_target_and_face(has_equal);
  }

  void update_vertex_face_flags(Vertex_handle v, Face_handle f)
  {
    bool is_equal, has_equal;
    is_equal = (v->get_decision() == f->get_decision());
    // has equal can be true even if the decision is not the same,
    // but has same surfaces, i.e. one of the features got DAC_DECISION_BOTH
    // decision, and the other didn't

    has_equal = (v->get_decision() == f->get_decision() ||
                 v->get_decision() == DAC_DECISION_BOTH ||
                 f->get_decision() == DAC_DECISION_BOTH);

    CGAL::Dac_decision decision = v->get_decision();
    bool is_equal_first = (v->get_is_equal_aux_data_in_face(0));
    bool has_equal_first = (v->get_has_equal_aux_data_in_face(0));
    bool is_equal_second = (v->get_is_equal_aux_data_in_face(1));
    bool has_equal_second = (v->get_has_equal_aux_data_in_face(1));
    if (decision == DAC_DECISION_FIRST)
    {
      is_equal &= is_equal_first;
      has_equal &= has_equal_first;
    }
    else if (decision == DAC_DECISION_SECOND)
    {
      is_equal &= is_equal_second;
      has_equal &= has_equal_second;
    }
    else
    {
      is_equal &= (is_equal_first & is_equal_second);
      // we check if the face has a different decision, and if so,
      // we update the flag according to the face decision
      decision = f->get_decision();
      if (decision == DAC_DECISION_FIRST)
      has_equal &= has_equal_first;
      else if (decision == DAC_DECISION_SECOND)
        has_equal &= has_equal_second;
      else      
      has_equal &= (has_equal_first & has_equal_second);
    }
    v->set_is_equal_data_in_face(is_equal);
    v->set_has_equal_data_in_face(has_equal);
  }

  void update_flags(Minimization_diagram_2& result)
  {    
    // edges
    Halfedge_iterator hi = result.halfedges_begin();
    for (; hi != result.halfedges_end(); ++hi)
    {
      update_edge_face_flags(hi);
      update_edge_target_flags(hi);
      update_target_face_flags(hi);
    }

    // vertices
    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
      if (vi->is_isolated())
    update_vertex_face_flags(vi, vi->face());

  }

  template <class FeatureHandle>
  bool get_aux_is_set(FeatureHandle fh, unsigned int id)
  {
    return fh->get_aux_is_set(id);
  }
 
  template <class FeatureHandle>
  bool aux_has_no_data(FeatureHandle fh, unsigned int id)
  {
    const Object& o = fh->get_aux_source(id);
    Halfedge_handle h;
    Vertex_handle v;
    Face_handle f;
   
    // aux source of a face must be a face!
    // aux source of a halfedge can be face or halfedge
    // aux source of a vertex can be face, halfedge or vertex
    // this is why we start with a check for a face, then halfedge
    // and last vertex
    if (assign(f, o))
      return f->has_no_data();
    else if (assign(h, o))
      return h->has_no_data();
    else
    {
      CGAL_assertion_code(bool b =)
      assign(v, o);
      CGAL_assertion(b);
      return v->has_no_data();
    }
  }

protected:

  //***************************************************************************
  // methods for assertion and checking
  //***************************************************************************

  void get_data_iterators(Object aux_src,
                          Envelope_data_iterator& begin,
                          Envelope_data_iterator& end)
  {
    CGAL_assertion(!aux_src.is_empty());
    Vertex_handle v;
    Halfedge_handle h;
    Face_handle f;

    if (assign(v, aux_src))
    {
      begin = v->begin_data();
      end = v->end_data();
    }
    else if (assign(h, aux_src))
    {
      begin = h->begin_data();
      end = h->end_data();
    }
    else
    {
      CGAL_assertion(assign(f, aux_src));
      assign(f, aux_src);
      begin = f->begin_data();
      end = f->end_data();
    }
  }
  bool is_equal_data(Object o,
                     Envelope_data_iterator begin,
                     Envelope_data_iterator end)
  {
    CGAL_assertion(!o.is_empty());
    Vertex_handle v;
    Halfedge_handle h;
    Face_handle f;

    if (assign(v, o))
      return v->is_equal_data(begin, end);
    else if (assign(h, o))
      return h->is_equal_data(begin, end);
    else
    {
      CGAL_assertion(assign(f, o));
      assign(f, o);
      return f->is_equal_data(begin, end);
    }
  }
  bool has_equal_data(Object o,
                      Envelope_data_iterator begin,
                      Envelope_data_iterator end)
  {
    CGAL_assertion(!o.is_empty());
    Vertex_handle v;
    Halfedge_handle h;
    Face_handle f;

    if (assign(v, o))
      return v->has_equal_data(begin, end);
    else if (assign(h, o))
      return h->has_equal_data(begin, end);
    else
    {
      CGAL_assertion(assign(f, o));
      assign(f, o);
      return f->has_equal_data(begin, end);
    }
  }

  void print_decisions(Minimization_diagram_2& result)
  {
    Face_iterator fi = result.faces_begin();
    for (; fi != result.faces_end(); ++fi)
    {
      if (fi->is_unbounded())
        continue;
      Face_handle fh = fi;
      Ccb_halfedge_circulator face_hec = fh->outer_ccb();
      Ccb_halfedge_circulator face_hec_begin = face_hec;
      do {
        ++face_hec;
      } while(face_hec != face_hec_begin);

      Inner_ccb_iterator inner_iter = fh->holes_begin();
      for (; inner_iter != fh->holes_end(); ++inner_iter)

      {
        face_hec = face_hec_begin = (*inner_iter);
        do {
          ++face_hec;
        } while(face_hec != face_hec_begin);
      }
    }
  }
   

  // confirm that aux source and decision are set over all minimization 
  // diagram features
  bool check_resolve_was_ok(Minimization_diagram_2& result)
  {
    bool all_ok = true;
    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      /* if (vh->get_is_fake())
       *   continue;
       */

      all_ok &= (vh->get_aux_is_set(0));
      CGAL_assertion_msg(all_ok, "aux source (0) not set over vertex");
      all_ok &= (vh->get_aux_is_set(1));
      CGAL_assertion_msg(all_ok, "aux source (1) not set over vertex");
      all_ok &= (vh->is_decision_set());
      CGAL_assertion_msg(all_ok, "decision was not set over vertex");

    }
    Halfedge_iterator hi = result.halfedges_begin();
    for (; hi != result.halfedges_end(); ++hi)
    {
      Halfedge_handle hh = hi;
     /* if (hh->get_is_fake())
      *   continue;
      */

      all_ok &= (hh->get_aux_is_set(0));
      CGAL_assertion_msg(all_ok, "aux source (0) not set over edge");
      all_ok &= (hh->get_aux_is_set(1));
      CGAL_assertion_msg(all_ok, "aux source (1) not set over edge");
      all_ok &= (hh->is_decision_set());
      CGAL_assertion_msg(all_ok, "decision was not set over edge");
    }
    Face_iterator fi = result.faces_begin();
    for (; fi != result.faces_end(); ++fi)
    {
      Face_handle fh = fi;

      all_ok &= (fh->get_aux_is_set(0));
      CGAL_assertion_msg(all_ok, "aux source (0) not set over face");
      all_ok &= (fh->get_aux_is_set(1));
      CGAL_assertion_msg(all_ok, "aux source (1) not set over face");
      all_ok &= (fh->is_decision_set());
      CGAL_assertion_msg(all_ok, "decision was not set over face");
    }
    return all_ok;
  }

  // confirm that envelope data is set over all minimization diagram features
  // and that no fake feature are left
  bool is_envelope_valid(Minimization_diagram_2& result)
  {
    bool all_ok = true;
    Vertex_iterator vi = result.vertices_begin();
    for (; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      
      all_ok &= (vh->get_is_set());
      CGAL_assertion_msg(all_ok, "data not set over vertex");
      all_ok &= (!vh->has_no_data());

      CGAL_assertion_msg(all_ok, "data empty over vertex");      
     /* all_ok &= (!vh->get_is_fake());*/
      CGAL_assertion_msg(all_ok, "fake vertex in envelope");
    }
    Halfedge_iterator hi = result.halfedges_begin();
    for (; hi != result.halfedges_end(); ++hi)
    {
      Halfedge_handle hh = hi;
      
      all_ok &= (hh->get_is_set());
      if (!all_ok)
        std::cerr << "edge: " << hh->curve() << std::endl;
      CGAL_assertion_msg(all_ok, "data not set over edge");
      all_ok &= (!hh->has_no_data());
      if (!all_ok)
        std::cerr << "edge: " << hh->curve() << std::endl;
      CGAL_assertion_msg(all_ok, "data empty over edge");

      /*all_ok &= (!hh->get_is_fake());*/
      CGAL_assertion_msg(all_ok, "fake edge in envelope");
    }
    Face_iterator fi = result.faces_begin();
    for (; fi != result.faces_end(); ++fi)
    {
      Face_handle fh = fi;
      all_ok &= (fh->get_is_set());
      CGAL_assertion_msg(all_ok, "data not set over face");
    }
    return all_ok;
  }
    
  // observer for the minimization diagram
  // keeps the relevant data in the new faces
  class Keep_face_data_observer : public Md_observer
  {
  public:
    typedef typename Minimization_diagram_2::Face_handle Face_handle;
    
    Keep_face_data_observer(Minimization_diagram_2& arr) :
      Md_observer(arr)
    {}

    virtual void after_split_face(Face_handle org_f,
                                  Face_handle new_f,
                                  bool /* is_hole*/)
    {
      // update data in the new face from the original face
      if (org_f->get_aux_is_set(0))
        new_f->set_aux_source(0, org_f->get_aux_source(0));
      if (org_f->get_aux_is_set(1))
        new_f->set_aux_source(1, org_f->get_aux_source(1));
      if (org_f->is_decision_set())
        new_f->set_decision(org_f->get_decision());
    }
  };
  

  // observer for the minimization diagram
  // keeps the relevant data in the new edges & vertices
  class Keep_edge_data_observer : public Md_observer
  {
  public:
    typedef typename Minimization_diagram_2::Halfedge_handle   Halfedge_handle;
    typedef typename Minimization_diagram_2::Vertex_handle     Vertex_handle;
    typedef typename Minimization_diagram_2::X_monotone_curve_2 
                                                            X_monotone_curve_2;

    typedef typename Envelope_divide_and_conquer_3<Traits, 
                                                   Minimization_diagram_2,
                                                   EnvelopeResolver_3, 
                                                   Overlay_2>::Self Self;
    Keep_edge_data_observer(Minimization_diagram_2& arr,
                            Self* b) :
      Md_observer(arr), base(b)
    {
      CGAL_assertion(base != NULL);
    }

    /* virtual void before_split_edge (Halfedge_handle e,
     *                                 Vertex_handle v,
     *                                 const X_monotone_curve_2& c1,
     *                                 const X_monotone_curve_2& c2)
     * {}
     */

    virtual void after_split_edge(Halfedge_handle he1, Halfedge_handle he2)
    {
      // update data of the new vertex, which is the common vertex of he1 and 
      // he2, and of the new edge according to the data in the original edge
      CGAL_assertion(he2->source() == he1->target());

      Vertex_handle new_vertex;
//      if (he2->source() == he1->target() ||
//          he2->source() == he1->source())
      new_vertex = he2->source();
//      else
//        new_vertex = he2->target();

      CGAL_assertion(!new_vertex->is_decision_set());
      CGAL_assertion(!new_vertex->get_aux_is_set(0));
      CGAL_assertion(!new_vertex->get_aux_is_set(1));

      // find the halfedge with the additional information, to be copied into
      // the second halfedge
      Halfedge_handle org_he = he1, new_he = he2;
      
      if (org_he->is_decision_set())
      {
        new_he->set_decision(org_he->get_decision());
        new_he->twin()->set_decision(org_he->get_decision());
        new_vertex->set_decision(org_he->get_decision());
      }        
      if (org_he->get_aux_is_set(0))
      {
        new_vertex->set_aux_source(0, org_he->get_aux_source(0));
        new_he->set_aux_source(0, org_he->get_aux_source(0));
        new_he->twin()->set_aux_source(0, org_he->twin()->get_aux_source(0));
      }      
      if (org_he->get_aux_is_set(1))
      {
        new_vertex->set_aux_source(1, org_he->get_aux_source(1));
        new_he->set_aux_source(1, org_he->get_aux_source(1));
        new_he->twin()->set_aux_source(1, org_he->twin()->get_aux_source(1));
      }

      /*new_he->set_is_fake(org_he->get_is_fake());*/
      /*new_he->twin()->set_is_fake(org_he->get_is_fake());*/
      /*new_vertex->set_is_fake(org_he->get_is_fake());*/

      // update all new bools
      new_he->set_is_equal_aux_data_in_face
        (0, org_he->get_is_equal_aux_data_in_face(0));
      new_he->twin()->set_is_equal_aux_data_in_face
        (0, org_he->twin()->get_is_equal_aux_data_in_face(0));
      new_he->set_is_equal_aux_data_in_face
        (1, org_he->get_is_equal_aux_data_in_face(1));
      new_he->twin()->set_is_equal_aux_data_in_face
        (1, org_he->twin()->get_is_equal_aux_data_in_face(1));

      new_he->set_has_equal_aux_data_in_face
        (0, org_he->get_has_equal_aux_data_in_face(0));
      new_he->twin()->set_has_equal_aux_data_in_face
        (0, org_he->twin()->get_has_equal_aux_data_in_face(0));
      new_he->set_has_equal_aux_data_in_face
        (1, org_he->get_has_equal_aux_data_in_face(1));
      new_he->twin()->set_has_equal_aux_data_in_face
        (1, org_he->twin()->get_has_equal_aux_data_in_face(1));

      // new_he->target is the original edge's target, and org_he->target
      // is the new vertex
      new_he->set_is_equal_aux_data_in_target
        (0, org_he->get_is_equal_aux_data_in_target(0));
      new_he->set_is_equal_aux_data_in_target
        (1, org_he->get_is_equal_aux_data_in_target(1));
      org_he->set_is_equal_aux_data_in_target(0, true);
      org_he->set_is_equal_aux_data_in_target(1, true);
      new_he->set_has_equal_aux_data_in_target
        (0, org_he->get_has_equal_aux_data_in_target(0));
      new_he->set_has_equal_aux_data_in_target
        (1, org_he->get_has_equal_aux_data_in_target(1));
      org_he->set_has_equal_aux_data_in_target
        (0, !base->aux_has_no_data(org_he, 0));
      org_he->set_has_equal_aux_data_in_target
        (1, !base->aux_has_no_data(org_he, 1));
      new_he->set_has_equal_aux_data_in_target_and_face
        (0, org_he->get_has_equal_aux_data_in_target_and_face(0));
      new_he->set_has_equal_aux_data_in_target_and_face
        (1, org_he->get_has_equal_aux_data_in_target_and_face(1));
      org_he->set_has_equal_aux_data_in_target_and_face
        (0, org_he->get_has_equal_aux_data_in_face(0));
      org_he->set_has_equal_aux_data_in_target_and_face
        (1, org_he->get_has_equal_aux_data_in_face(1));

      // new_he->source is the new vertex, and org_he->source is the
      // original vertex
      new_he->twin()->set_is_equal_aux_data_in_target(0, true);
      new_he->twin()->set_is_equal_aux_data_in_target(1, true);
      new_he->twin()->set_has_equal_aux_data_in_target
        (0,!base->aux_has_no_data(org_he, 0));
      new_he->twin()->set_has_equal_aux_data_in_target
        (1,!base->aux_has_no_data(org_he, 1));

      new_he->twin()->set_has_equal_aux_data_in_target_and_face
        (0, org_he->twin()->get_has_equal_aux_data_in_face(0));
      new_he->twin()->set_has_equal_aux_data_in_target_and_face
        (1, org_he->twin()->get_has_equal_aux_data_in_face(1));
    }

  protected:
    Self* base;
  };

#ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER

  // A BFS visitor class which collects the faces that need resolving
  // in a list according to the discover time
  // In our case graph vertices represent minimization diagram faces
  template <class IndexMap>
  class Faces_order_bfs_visitor : public boost::default_bfs_visitor
  {
  public:
    typedef typename Envelope_divide_and_conquer_3<Traits,
                                                   Minimization_diagram_2,
                                                   EnvelopeResolver_3,
                                                   Overlay_2>::Self Self;

  protected:
    const IndexMap* index_map;          // Mapping vertices to indices
    std::list<Face_handle>& faces;      // The ordered faces list
    Self* base;

  public:

    // Constructor.
    Faces_order_bfs_visitor(const IndexMap& imap, std::list<Face_handle>& f,
                            Self* b) :
      index_map(&imap),
      faces(f),
      base(b)
    {
      CGAL_assertion(base != NULL);
    }

    // Write the discover time for a given vertex.
    template <typename Vertex, typename Graph>

    void discover_vertex(Vertex fh, const Graph&)
    {
      // first we check if we can set the decision immediately
      // if a surface of one map doesn't exist, then we set the second surface
      if (base->aux_has_no_data(fh, 0) && !base->aux_has_no_data(fh, 1))
        fh->set_decision(DAC_DECISION_SECOND);
      else if (base->aux_has_no_data(fh, 0) && base->aux_has_no_data(fh, 1))
      {

        fh->set_decision(EQUAL);
        fh->set_no_data();
      }
      else if (!base->aux_has_no_data(fh, 0) && base->aux_has_no_data(fh, 1))
        fh->set_decision(DAC_DECISION_FIRST);
      else
        // here, we have both surfaces.
        // we save the face in a list for a later treatment, because the
        // face can change and destroy the iterator
        faces.push_back(fh);
    }
  };
#endif
  
protected:
  Envelope_resolver* resolver;
  const Traits* m_geom_traits;
  bool m_own_traits; 
  bool m_is_lower;
};

} //namespace CGAL

#endif //CGAL_ENVELOPE_DIVIDE_AND_CONQUER_3_H
