// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source: /CVSROOT/CGAL/Packages/Envelope_3/include/CGAL/Envelope_divide_and_conquer_3.h,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_3_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_3_H

#define CGAL_ENVELOPE_SAVE_COMPARISONS
#define CGAL_ENVELOPE_USE_BFS_FACE_ORDER

#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Envelope_3/Env_overlay_2.h>
#include <CGAL/Envelope_3/Envelope_element_visitor_3.h>
#include <CGAL/Envelope_3/No_vertical_decomposition_2.h>
#include <CGAL/Envelope_3/set_dividors.h>

#ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <CGAL/graph_traits_Dual_Arrangement_2.h>
#include <CGAL/Arr_face_map.h>
#endif

#include <iostream>
#include <cassert>
#include <list>
#include <set>
#include <vector>
#include <map>
#include <time.h>

// this base divide & conquer algorithm splits the input into 2 groups,
// calculates the result over the 2 groups, and then merges the results like
// this:
// - overlays the maps, and set at each vertex, edge and face the data from
//   the 2 maps
// - foreach vertex, decide from the 2 envelopes data, which is the envelope
//   of all surfaces
// - foreach edge, do the same, using the Resolver class. an edge can split to
//   a constant number of parts, each with its own envelope data.
// - possibly do a decomposition (vertical / partial or other) to make faces
//   with no holes
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

CGAL_BEGIN_NAMESPACE

// The algorithm has 5 template parameters:
// 1. EnvelopeTraits_3        - the geometric traits class
// 2. MinimizationDiagram_2   - the type of the output, which is an arrangement
//                              with additional information (list of surfaces)
//                              in vertices, edges & faces
// 3. VerticalDecomposition_2 - a vertical decomposition for 2D
//                              MinimizationDiagram_2 so that we can switch
//                              from partial/full/other/none decomposition.
// 4. EnvelopeResolver_3      - part of the algorithm that solves the shape of
//                              the envelope between 2 surfaces over a feature
//                              of the arrangement
// 5. Overlay_2               - overlay of 2 MinimizationDiagram_2
template <class EnvelopeTraits_3, class MinimizationDiagram_2,
          class VerticalDecomposition_2 = No_vertical_decomposition_2<MinimizationDiagram_2>,
          class EnvelopeResolver_3 =
                Envelope_element_visitor_3<EnvelopeTraits_3, MinimizationDiagram_2>,
          class Overlay_2 = Envelope_overlay_2<MinimizationDiagram_2> >
class Envelope_divide_and_conquer_3
{
public:
  typedef EnvelopeTraits_3                                          Traits;
  typedef typename Traits::Surface_3                                Surface_3;
  typedef typename Traits::Xy_monotone_surface_3                    Xy_monotone_surface_3;

  typedef MinimizationDiagram_2                                     Minimization_diagram_2;

  typedef typename Traits::Point_2                                  Point_2;
  typedef typename Traits::X_monotone_curve_2                       X_monotone_curve_2;
  typedef typename Traits::Curve_2                                  Curve_2;

  typedef EnvelopeResolver_3                                        Envelope_resolver;
  typedef VerticalDecomposition_2                                   Vertical_decomposition_2;

  typedef Envelope_divide_and_conquer_3<Traits, Minimization_diagram_2,
                                        Vertical_decomposition_2,
                                        EnvelopeResolver_3,
                                        Overlay_2>                  Self;

protected:

  typedef typename Minimization_diagram_2::Halfedge_const_iterator  Halfedge_const_iterator;
  typedef typename Minimization_diagram_2::Halfedge_handle          Halfedge_handle;
  typedef typename Minimization_diagram_2::Halfedge_iterator        Halfedge_iterator;
  typedef typename Minimization_diagram_2::Face_handle              Face_handle;
  typedef typename Minimization_diagram_2::Face_iterator            Face_iterator;
  typedef typename Minimization_diagram_2::Vertex_handle            Vertex_handle;
  typedef typename Minimization_diagram_2::Vertex_iterator          Vertex_iterator;
  typedef typename Minimization_diagram_2::Hole_iterator            Hole_iterator;
  typedef typename Minimization_diagram_2::Ccb_halfedge_circulator  Ccb_halfedge_circulator;
	typedef typename Minimization_diagram_2::Halfedge_around_vertex_circulator
                                                                    Halfedge_around_vertex_circulator;

  typedef Arr_observer<Minimization_diagram_2>                      Md_observer;
  typedef typename Minimization_diagram_2::Dcel::Dcel_data_iterator Envelope_data_iterator;

#ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER
  typedef CGAL::Dual<Minimization_diagram_2>                        Dual_Minimization_diagram_2;
#endif

public:
  // c'tor
  Envelope_divide_and_conquer_3()
  {    
    // Allocate the traits.
    traits = new Traits;                                                     

    own_traits = true;

    // Allocate the Envelope resolver with our traits
    resolver = new Envelope_resolver(traits);
  }

  Envelope_divide_and_conquer_3(Traits* tr)
  {
    // Set the traits.
    traits = tr;
    own_traits = false;

    // Allocate the Envelope resolver with our traits
    resolver = new Envelope_resolver(traits);
  }
  
  // virtual destructor.
  virtual ~Envelope_divide_and_conquer_3()
  {
    // Free the traits object, if necessary.
    if (own_traits)
      delete traits;

    // Free the resolver
    delete resolver;
  }

  // compute the envelope of surfaces in 3D, using the default arbitrary
  // dividor
  template <class SurfaceIterator>
  void construct_lu_envelope(SurfaceIterator begin, SurfaceIterator end,
                             Minimization_diagram_2 &result)
  {
    Arbitrary_dividor dividor;
    construct_lu_envelope(begin, end, result, dividor);
  }
    

  // compute the envelope of surfaces in 3D using the given set dividor
  template <class SurfaceIterator, class SetDividor>
  void construct_lu_envelope(SurfaceIterator begin, SurfaceIterator end, 
                             Minimization_diagram_2 &result,
                             SetDividor& dividor)
  {
    if (begin == end)
    {
      return; // result is empty
    }
    
    // make the general surfaces xy-monotone
    std::list<Xy_monotone_surface_3> xy_monotones;
    for(; begin != end; ++begin)
      traits->construct_envelope_xy_monotone_parts_3_object()
                        (*begin, std::back_inserter(xy_monotones));

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
                                      Minimization_diagram_2 &result)
  {
    Arbitrary_dividor dividor;
    construct_envelope_xy_monotone(begin, end, result, dividor);
  }

  // compute the envelope of xy-monotone surfaces in 3D using the given
  // set dividor
  template <class SurfaceIterator, class SetDividor>
  void construct_envelope_xy_monotone(SurfaceIterator begin,
                                      SurfaceIterator end,
                                      Minimization_diagram_2 &result,
                                      SetDividor& dividor)
  {
    if (begin == end)
      return; // result is empty

    init_stats();

    // recursively construct the envelope of the xy-monotone parts
    construct_lu_envelope_xy_monotones(begin, end, result, dividor);
    CGAL_assertion(is_envelope_valid(result));
  }

  /*! Access the traits object (const version). */
  const Traits* get_traits () const
  {
    return (traits);
  }

  /*! Access the traits object (non-const version). */
  Traits* get_traits () 
  {
    return (traits);
  }

  void reset()
  {
    resolver->reset();
    // reset statistical measures
    init_stats();
  }

  
  
protected:

  // compute the envelope of xy-monotone surfaces in 3D 
  template <class SurfaceIterator, class SetDividor>
  void construct_lu_envelope_xy_monotones(SurfaceIterator begin,
                                          SurfaceIterator end,
                                          Minimization_diagram_2 &result,
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
   
    // divide the surfaces into 2 groups (insert surface to each group alternately)
    std::list<Xy_monotone_surface_3> group1, group2;
    dividor(first, end,
            std::back_inserter(group1), std::back_inserter(group2));
    
    // recursively calculate the LU_envelope of the 2 groups
    Minimization_diagram_2 result1(traits), result2(traits);
    construct_lu_envelope_xy_monotones(group1.begin(), group1.end(), result1, dividor);
    construct_lu_envelope_xy_monotones(group2.begin(), group2.end(), result2, dividor);
        
    // merge the results:
    merge_envelopes(result1, result2, result);

    result1.clear();
    result2.clear();
    
    CGAL_assertion(is_envelope_valid(result));   
  }

  // insert one surface into an empty minimization diagram
  // to create the envelope on this surface only
  void deal_with_one_surface(Xy_monotone_surface_3& surf, Minimization_diagram_2& result)
  {
   
    std::list<Curve_2> boundary_curves;
    traits->construct_projected_boundary_curves_2_object()(surf, std::back_inserter(boundary_curves));

    typename std::list<Curve_2>::iterator boundary_it = boundary_curves.begin();
    for(; boundary_it != boundary_curves.end(); ++boundary_it)
    {
      insert_curve(result, *boundary_it);
    }
    // todo: should we do incremental/aggregate insert?
    //  insert(result, boundary_curves.begin(), boundary_curves.end());
    
    // update the information in the faces

    // update the holes of the unbounded face to indicate that this surface is the envelope
    // over them
    Face_handle uface = result.unbounded_face();
    Hole_iterator hole = uface->holes_begin();
    for (; hole != uface->holes_end(); ++hole)
    {
      Ccb_halfedge_circulator hec = (*hole);
      Halfedge_handle ccb_he = hec;

      // unbounded face has no data
      ccb_he->set_is_equal_data_in_face(false);
      ccb_he->set_has_equal_data_in_face(false);
      ccb_he->set_has_equal_data_in_target_and_face(false);
      
      Face_handle face = ccb_he->twin()->face();
      if (face->is_unbounded())
        continue;
      face->set_data(surf);
      // update the is/has equal_data_in_face in all the halfedges of the face's
      // boundary to true
  	  // and also target-face has_equal flag to true
      Ccb_halfedge_circulator face_hec = face->outer_ccb();
      Ccb_halfedge_circulator face_hec_begin = face_hec;
      do {
        face_hec->set_is_equal_data_in_face(true);
        face_hec->set_has_equal_data_in_face(true);
        face_hec->set_has_equal_data_in_target_and_face(true);

        ++face_hec;
      } while(face_hec != face_hec_begin);

      Hole_iterator inner_iter = face->holes_begin();
      for (; inner_iter != face->holes_end(); ++inner_iter)
      {
        face_hec = face_hec_begin = (*inner_iter);
        do {
          face_hec->set_is_equal_data_in_face(true);
          face_hec->set_has_equal_data_in_face(true);
          face_hec->set_has_equal_data_in_target_and_face(true);
          ++face_hec;
        } while(face_hec != face_hec_begin);
      }
    }
    
    // update other faces to indicate that no surface is over them
    Face_iterator fi = result.faces_begin();
    for(; fi != result.faces_end(); ++fi)
    {
      Face_handle fh = fi;
      if (!fh->get_is_set())
      {
        fh->set_no_data();
        // update the is/has equal_data_in_face in all the halfedges of the face's
        // boundary to false
        Ccb_halfedge_circulator face_hec, face_hec_begin;
        if (!fh->is_unbounded())
        {
          face_hec = fh->outer_ccb();
          face_hec_begin = face_hec;
          do {
            face_hec->set_is_equal_data_in_face(false);
            face_hec->set_has_equal_data_in_face(false);
            face_hec->set_has_equal_data_in_target_and_face(false);
            ++face_hec;
          } while(face_hec != face_hec_begin);
        }
        Hole_iterator inner_iter = fh->holes_begin();
        for (; inner_iter != fh->holes_end(); ++inner_iter)
        {
          face_hec = face_hec_begin = (*inner_iter);
          do {
            face_hec->set_is_equal_data_in_face(false);
            face_hec->set_has_equal_data_in_face(false);
            face_hec->set_has_equal_data_in_target_and_face(false);
            ++face_hec;
          } while(face_hec != face_hec_begin);
        }
      }
    }
    uface->set_no_data();

    // update information in all the edges & vertices to indicate that
    // this surface is the envelope
    Halfedge_iterator hi = result.halfedges_begin();
    for(; hi != result.halfedges_end(); ++hi)
    {
      Halfedge_handle hh = hi;
      hi->set_data(surf);
      // since all the edges & vertices have their envelope data equal to the
      // current surface, we can set is/has equal_data_in_target of all halfedges
      // to true
      hi->set_is_equal_data_in_target(true);
      hi->set_has_equal_data_in_target(true);
    }

    Vertex_iterator vi = result.vertices_begin();
    for(; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      vh->set_data(surf);
      if (vh->is_isolated())
      {
        // update the is/has equal_data_in_face flags according to the face data
        bool equal_data = !vh->face()->has_no_data();
        vh->set_is_equal_data_in_face(equal_data);
        vh->set_has_equal_data_in_face(equal_data);
      }
    }

	  CGAL_assertion(verify_flags(result));
  }

public:
  
  void merge_envelopes(Minimization_diagram_2& result1,
                       Minimization_diagram_2& result2,
                       Minimization_diagram_2& result)
  {
    // overlay the 2 arrangements
   
    overlay(result1, result2, result);
  
    CGAL_expensive_assertion_msg(is_valid(result), "after overlay result is not valid");
       
    // make sure the aux flags are correctly set by the overlay
	  CGAL_assertion(verify_aux_flags(result));

    // for each face, edge and vertex in the result, should calculate
    // which surfaces are on the envelope
    // a face can be cut, or faces can be merged.
   
    // now the minimization diagram might change - we need to keep data in the
    // edges, when they're split
    Keep_edge_data_observer edge_observer(result, this);

    // compute the surface on the envelope for each edge
    // edge can be split as surfaces can intersect (or touch) over it
    std::list<Halfedge_handle> edges_to_resolve;
    Halfedge_iterator hi = result.halfedges_begin();
    for(; hi != result.halfedges_end(); ++hi, ++hi)
    {
      Halfedge_handle hh = hi;
      // there must be data from at least one map, because all the surfaces are continous
      if (!get_aux_is_set(hh, 0) || !get_aux_is_set(hh, 1))
        continue;
      CGAL_assertion(get_aux_is_set(hh, 0));
      CGAL_assertion(get_aux_is_set(hh, 1));
      CGAL_assertion(!aux_has_no_data(hh, 1) || !aux_has_no_data(hh, 0));
      if (aux_has_no_data(hh, 0) && !aux_has_no_data(hh, 1))
      {
    		hh->set_decision(SECOND);
    		hh->twin()->set_decision(SECOND);
        continue;
      }
      else if (!aux_has_no_data(hh, 0) && aux_has_no_data(hh, 1))
      {
    		hh->set_decision(FIRST);
    		hh->twin()->set_decision(FIRST);
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
    decompose(result);
    CGAL_expensive_assertion_msg(result.is_valid(), 
                       "after decomposition result is not valid");

    // compute the surface on the envelope for each face,
    // splitting faces if needed

    std::list<Face_handle> faces_to_split;

    #ifdef CGAL_ENVELOPE_USE_BFS_FACE_ORDER
      // we traverse the faces of result in BFS order to maximize the
      // efficiency gain by the conclusion mechanism of
      // compare_distance_to_envelope results
      // Create a mapping of result faces to indices.
      CGAL::Arr_face_index_map<Minimization_diagram_2> index_map (result);

      // Perform breadth-first search from the unbounded face, and use the BFS
      // visitor to associate each arrangement face with its discover time.
      Faces_order_bfs_visitor<CGAL::Arr_face_index_map<Minimization_diagram_2> >
                                      bfs_visitor (index_map, faces_to_split, this);
      Face_handle first_face = result.faces_begin();
      if (result.number_of_faces() > 1)
        first_face = ++(result.faces_begin());

      boost::breadth_first_search (Dual_Minimization_diagram_2(result),
                                   first_face,
                                   boost::vertex_index_map(index_map).
                                   visitor (bfs_visitor));
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
          fh->set_decision(SECOND);
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
          fh->set_decision(FIRST);
          continue;
        }

        // here, we have both surfaces.
        // we save the face in a list for a later treatment, because the face can change
        // and destroy the iterator
        faces_to_split.push_back(fh);
      }
    #endif
        
    deal_with_faces_to_split(faces_to_split, result);

//    #ifndef CGAL_ENVELOPE_SAVE_COMPARISONS
//      hi = result.halfedges_begin();
//      for(; hi != result.halfedges_end(); ++hi, ++hi)
//      {
//        if (!hi->is_decision_set())
//          resolver->resolve(hi, result);
//      }
//    #endif
    
    // detach the edge_observer from result, since no need for it anymore
    edge_observer.detach();

    // compute the surface on the envelope for each vertex
    Vertex_iterator vi = result.vertices_begin();
    for(; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      if (vh->is_decision_set() || vh->get_is_fake())
        continue;
      // there must be data from at least one map, because all the surfaces are continous
      CGAL_assertion(get_aux_is_set(vh, 0));
      CGAL_assertion(get_aux_is_set(vh, 1));
      CGAL_assertion(!aux_has_no_data(vh, 1) || !aux_has_no_data(vh, 0));
      if (aux_has_no_data(vh, 0) && !aux_has_no_data(vh, 1))
      {
        vh->set_decision(SECOND);
        continue;
      }
      else if (!aux_has_no_data(vh, 0) && aux_has_no_data(vh, 1))
      {
        vh->set_decision(FIRST);
        continue;
      }
      resolver->resolve(vh);
     
    }

    CGAL_expensive_assertion_msg(result.is_valid(), "after resolve result is not valid");

    // make sure that aux_source and decision are set at all features
	  // after all resolvings
    CGAL_assertion(check_resolve_was_ok(result));
 
    // make sure the aux flags are correctly after all resolvings
	  CGAL_assertion(verify_aux_flags(result));

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

    // update is_equal_data and has_equal_data of halfedge->face and vertex->face
    // relations, according to the decision, and the aux similar flags
	  update_flags(result);

    // update the envelope surfaces according to the decision and the aux
	  // surfaces in aux source
    update_envelope_surfaces_by_decision(result);
    
	  // make sure that all the flags are correctly set on the envelope result
	  CGAL_assertion(verify_flags(result));
    CGAL_expensive_assertion_msg(is_valid(result), "after merge result is not valid");
  }


protected:

  // do a vertical decomposition
  // TODO: want to do a partial decomposition, or any other decomposition that
  //       afterwards we have simple faces (with no holes)
  void decompose(Minimization_diagram_2& result)
  {
    // before decomposing, we attach an observer that will copy data
    // into newly created faces (from the original faces they are split from)
    Keep_face_data_observer obs(result);
    Decomposition_observer dec_obs(result, this);
    
    obs.detach();
    dec_obs.detach();
  }
  
  void deal_with_faces_to_split(std::list<Face_handle>& faces_to_split, Minimization_diagram_2& result)
  {
    // for each face in faces_to_split, find the intersection over the face, and split the face
    typename std::list<Face_handle>::iterator li;
    for (li = faces_to_split.begin(); li != faces_to_split.end(); ++li)
    {
     
      resolver->resolve(*li, result);
    
    }
    faces_to_split.clear();
  }

  template <class InputIterator>
  bool is_equal_data(const InputIterator & begin1,
                     const InputIterator & end1,
                     const InputIterator & begin2,
                     const InputIterator & end2)
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
  bool has_equal_data(const InputIterator & begin1,
                      const InputIterator & end1,
                      const InputIterator & begin2,
                      const InputIterator & end2)
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
  bool has_equal_aux_data(unsigned int id, FeatureHandle1 fh1, FeatureHandle2 fh2)
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
      {
        edges.push_back(hh);
      }
    }

    // this is for counting the number of extra faces
    Stats_observer count_faces(result);
    
    for(typename std::list<Halfedge_handle>::iterator ci = edges.begin(); ci != edges.end(); ++ci)
    {
      // if the endpoints become isolated after the removal we need to remove
      // them if they have the same data as the edge
      Halfedge_handle h = *ci;
      Vertex_handle src = h->source(), trg = h->target();
      Face_handle h_face = h->face();
      bool remove_src = can_remove_edge_target(h->twin()),
           remove_trg = can_remove_edge_target(h);
      bool src_is_equal_0 = (h->twin()->get_is_equal_aux_data_in_face(0) &&
                             h->twin()->get_is_equal_aux_data_in_target(0));
      bool src_is_equal_1 = (h->twin()->get_is_equal_aux_data_in_face(1) &&
                             h->twin()->get_is_equal_aux_data_in_target(1));
      bool trg_is_equal_0 = (h->get_is_equal_aux_data_in_face(0) &&
                             h->get_is_equal_aux_data_in_target(0));
      bool trg_is_equal_1 = (h->get_is_equal_aux_data_in_face(1) &&
                             h->get_is_equal_aux_data_in_target(1));
      bool src_has_equal_0 = h->twin()->get_has_equal_aux_data_in_target_and_face(0);
      bool src_has_equal_1 = h->twin()->get_has_equal_aux_data_in_target_and_face(1);
      bool trg_has_equal_0 = h->get_has_equal_aux_data_in_target_and_face(0);
      bool trg_has_equal_1 = h->get_has_equal_aux_data_in_target_and_face(1);

      result.remove_edge(*ci, remove_src, remove_trg);
      // otherwise, we should make sure, they will not be removed
      // the first check is needed since if the vertex was removed, then the
      // handle is invalid
      if (!remove_src && src->is_isolated())
      {
        // to be precise we copy from the halfedge-face and halfedge-target relations
        src->set_is_equal_aux_data_in_face(0, src_is_equal_0);
        src->set_is_equal_aux_data_in_face(1, src_is_equal_1);
        // todo: the has_equal flags should be updated also
        // make sure h_face is also src face
        CGAL_assertion(h_face == src->face());
        // CGAL_assertion(src_has_equal_0 == has_equal_aux_data(0, src, h_face));
        // CGAL_assertion(src_has_equal_1 == has_equal_aux_data(1, src, h_face));
        src->set_has_equal_aux_data_in_face(0, src_has_equal_0);
        src->set_has_equal_aux_data_in_face(1, src_has_equal_1);
      }
      if (!remove_trg && trg->is_isolated())
      {
        trg->set_is_equal_aux_data_in_face(0, trg_is_equal_0);
        trg->set_is_equal_aux_data_in_face(1, trg_is_equal_1);
        // make sure h_face is also trg face
        CGAL_assertion(h_face == trg->face());
        // CGAL_assertion(trg_has_equal_0 == has_equal_aux_data(0, trg, h_face));
        // CGAL_assertion(trg_has_equal_1 == has_equal_aux_data(1, trg, h_face));
        trg->set_has_equal_aux_data_in_face(0, trg_has_equal_0);
        trg->set_has_equal_aux_data_in_face(1, trg_has_equal_1);
      }
    }

    count_faces.detach();    
  }

  template <class FeatureHandle>
  void get_aux_data_iterators(unsigned int id, FeatureHandle fh,
                              Envelope_data_iterator& begin,
                              Envelope_data_iterator& end)
  {
    Halfedge_handle h;
    Vertex_handle v;
  	Face_handle f;

    Object o = fh->get_aux_source(id);
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
    if (hh->get_is_fake())
    {
      CGAL_assertion(f2->get_decision() == f1->get_decision());
      return true;
    }
      
    // we check if the decision done on the edge is equal to the decision
    // done on the faces. if not, then the envelope surfaces must differ
    CGAL_assertion(hh->is_decision_set() && f1->is_decision_set() && f2->is_decision_set());
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

    // we assert that the flags' values are correct using comparison of data
    // as in the old way (using set operations)
    CGAL_assertion_code (
      Envelope_data_iterator begin1;
      Envelope_data_iterator end1;
      Envelope_data_iterator begin2;
      Envelope_data_iterator end2;

      get_aux_data_iterators(0, hh, begin1, end1);
      get_aux_data_iterators(0, f1, begin2, end2);
      bool b1 = is_equal_data(begin1, end1, begin2, end2);


      get_aux_data_iterators(0, f2, begin2, end2);
      bool b2 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(1, hh, begin1, end1);
      get_aux_data_iterators(1, f1, begin2, end2);
      bool b3 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(1, f2, begin2, end2);
      bool b4 = is_equal_data(begin1, end1, begin2, end2);
    );
    // after removes, the aux_source might be wrong for source that 

	  // has no connection to the decision, so cannot use assertions here
    // todo (after return) - is it correct to put it in comment
    // CGAL_assertion(equal_first == (b1 && b2));
    // CGAL_assertion(equal_second == (b3 && b4));
    if (decision == FIRST)
    {
      CGAL_assertion(equal_first == (b1 && b2));
      return equal_first;
    }
    else if (decision == SECOND)
	  {
      CGAL_assertion(equal_second == (b3 && b4));
      return equal_second;
	  }
    else
	  {
      CGAL_assertion(equal_first == (b1 && b2));
      CGAL_assertion(equal_second == (b3 && b4));
      return (equal_first && equal_second);
	  }
  }


  // check if can remove the edge's target if we remove the edge
  // this means that the target has the same envelope information as the edge
  bool can_remove_edge_target(Halfedge_handle h)
  {
    Vertex_handle v = h->target();
    if (v->get_is_fake() && !v->is_decision_set())
      return true;
    CGAL_assertion(v->is_decision_set());

    if (h->get_is_fake() && !h->is_decision_set())
    {
      h->set_decision(h->face()->get_decision());
      h->twin()->set_decision(h->get_decision());
    }
    CGAL_assertion(h->is_decision_set());

    // if the decision done on the vertex and edge are different,
    // the envelope differs too.
    if (h->get_decision() != v->get_decision())
      return false;

    // now, check the equality of the surfaces list according to the decision

    CGAL::Dac_decision decision = h->get_decision();
    bool equal_first = (h->get_is_equal_aux_data_in_target(0));
    bool equal_second = (h->get_is_equal_aux_data_in_target(1));

    // we assert that the flags' values are correct using comparison of data
    // as in the old way (using set operations)
    CGAL_assertion_code (
      Envelope_data_iterator begin1;
      Envelope_data_iterator end1;
      Envelope_data_iterator begin2;
      Envelope_data_iterator end2;

      get_aux_data_iterators(0, h, begin1, end1);
      get_aux_data_iterators(0, v, begin2, end2);
      bool b1 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(1, h, begin1, end1);
      get_aux_data_iterators(1, v, begin2, end2);
      bool b2 = is_equal_data(begin1, end1, begin2, end2);
    );
    
    // we don't assert here because the flags may only be true for
    // aux data that is related to the decision.
    //CGAL_assertion(equal_first == b1);
    //CGAL_assertion(equal_second == b2);

    if (decision == FIRST)
    {
      CGAL_assertion(equal_first == b1) ;
      return equal_first;
    }
    else if (decision == SECOND)
    {
      CGAL_assertion(equal_second == b2);
      return equal_second;
    }
    else
    {
      CGAL_assertion(equal_first == b1);
      CGAL_assertion(equal_second == b2);
      return (equal_first && equal_second);
    }
  }
  
  // check if we can remove an isolated vertex from the envelope
  // this can be done if the envelope surfaces on the vertex are the same as
  // the envelope surfaces on its incident face
  bool can_remove_isolated_vertex(Minimization_diagram_2& result, Vertex_handle vh)
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

    // we assert that the flags' values are correct using comparison of data
    // as in the old way (using set operations)
    CGAL_assertion_code (
      Envelope_data_iterator begin1;
      Envelope_data_iterator end1;
      Envelope_data_iterator begin2;
      Envelope_data_iterator end2;

      get_aux_data_iterators(0, vh, begin1, end1);
      get_aux_data_iterators(0, f, begin2, end2);
      bool b1 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(1, vh, begin1, end1);
      get_aux_data_iterators(1, f, begin2, end2);
      bool b2 = is_equal_data(begin1, end1, begin2, end2);
    );
    
    if (decision == FIRST)
  	{
      CGAL_assertion(equal_first == b1);
      return equal_first;
	  }
    else if (decision == SECOND)
	  {
      CGAL_assertion(equal_second == b2);
      return equal_second;
  	}
    else
	  {
      CGAL_assertion(equal_first == b1);
      CGAL_assertion(equal_second == b2);
      return (equal_first && equal_second);      
	  }
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
    
    if (vh->get_is_fake())
    {
      CGAL_assertion(he1->get_decision() == he2->get_decision());
      return true;
    }
    
    CGAL_assertion(vh->is_decision_set());
    // if the decision done on the vertex and its incident halfedges are different,
    // the envelope differs too.
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

    // we assert that the flags' values are correct using comparison of data
    // as in the old way (using set operations)
    CGAL_assertion_code (
      Envelope_data_iterator begin1;
      Envelope_data_iterator end1;
      Envelope_data_iterator begin2;
      Envelope_data_iterator end2;

      get_aux_data_iterators(0, vh, begin1, end1);
      get_aux_data_iterators(0, he1, begin2, end2);
      bool b1 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(0, he2, begin2, end2);
      bool b2 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(1, vh, begin1, end1);
      get_aux_data_iterators(1, he1, begin2, end2);
      bool b3 = is_equal_data(begin1, end1, begin2, end2);

      get_aux_data_iterators(1, he2, begin2, end2);
      bool b4 = is_equal_data(begin1, end1, begin2, end2);
    );
    // we cannot have the assertions here, since after merge,
	  // aux_source might not be correct for a source which is not related
    // to the decision
    //CGAL_assertion(equal_first == (b1 && b2));
    //CGAL_assertion(equal_second == (b3 && b4));

    if (decision == FIRST)
  	{
      CGAL_assertion(equal_first == (b1 && b2));
      return equal_first;
	  }
    else if (decision == SECOND)
    {
      CGAL_assertion(equal_second == (b3 && b4));
      return equal_second;
	  }
    else
	  {
      CGAL_assertion(equal_first == (b1 && b2));
      CGAL_assertion(equal_second == (b3 && b4));
      return (equal_first && equal_second);
  	}
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
    for(; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      if (!vh->is_decision_set() || vh->degree() == 2)
        candidates_to_remove.push_back(vh);
      else if (vh->is_isolated() &&
               can_remove_isolated_vertex(result, vh))
        isolated_to_remove.push_back(vh);
    }

    typename Traits::Merge_2 curves_merge = traits->merge_2_object();
    typename Traits::Are_mergeable_2 curves_can_merge = 
                                              traits->are_mergeable_2_object();

    // check the candidates and remove if necessary
    typename std::list<Vertex_handle>::iterator ci;
    for(ci = candidates_to_remove.begin(); 
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
      he1->set_is_equal_aux_data_in_target(0, he2->twin()->get_is_equal_aux_data_in_target(0));
      he1->set_is_equal_aux_data_in_target(1, he2->twin()->get_is_equal_aux_data_in_target(1));
      he1->set_has_equal_aux_data_in_target(0, he2->twin()->get_has_equal_aux_data_in_target(0));
      he1->set_has_equal_aux_data_in_target(1, he2->twin()->get_has_equal_aux_data_in_target(1));
      he1->set_has_equal_aux_data_in_target_and_face
		       (0, he2->twin()->get_has_equal_aux_data_in_target_and_face(0));
      he1->set_has_equal_aux_data_in_target_and_face
		       (1, he2->twin()->get_has_equal_aux_data_in_target_and_face(1));

      he2->set_is_equal_aux_data_in_target(0, he1->twin()->get_is_equal_aux_data_in_target(0));
      he2->set_is_equal_aux_data_in_target(1, he1->twin()->get_is_equal_aux_data_in_target(1));
      he2->set_has_equal_aux_data_in_target(0, he1->twin()->get_has_equal_aux_data_in_target(0));
      he2->set_has_equal_aux_data_in_target(1, he1->twin()->get_has_equal_aux_data_in_target(1));
      he2->set_has_equal_aux_data_in_target_and_face
		       (0, he1->twin()->get_has_equal_aux_data_in_target_and_face(0));
      he2->set_has_equal_aux_data_in_target_and_face
		       (1, he1->twin()->get_has_equal_aux_data_in_target_and_face(1));

      // order of halfedges for merge doesn't matter
      Halfedge_handle new_edge = result.merge_edge(he1, he2 ,c);

      CGAL_assertion(new_edge->is_decision_set());

      CGAL_expensive_assertion_msg(result.is_valid(), 
                         "after remove vertex result is not valid");
    }

    // remove isolated vertices

    typename std::list<Vertex_handle>::iterator li;
    for(li = isolated_to_remove.begin(); li != isolated_to_remove.end(); ++li)
    {
      Vertex_handle vh = *li;
      CGAL_assertion(vh->is_isolated());
      CGAL_assertion(can_remove_isolated_vertex(result, vh));

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

    if (decision == FIRST || decision == BOTH)
    {
      Envelope_data_iterator begin, end;
      get_aux_data_iterators(0, fh, begin, end);
      fh->set_data(begin, end);
    }
    if (decision == SECOND || decision == BOTH)
    {
      // copy data from second envelope
      Envelope_data_iterator begin, end;
      get_aux_data_iterators(1, fh, begin, end);

      if (decision == SECOND)
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
    for(; vi != result.vertices_end(); ++vi)
    {
      update_envelope_surfaces_by_decision(vi);
    }
    
    // edges
    Halfedge_iterator hi = result.halfedges_begin();
    for(; hi != result.halfedges_end(); ++hi)
    {
      update_envelope_surfaces_by_decision(hi);
    }
    
    // faces
    Face_iterator fi = result.faces_begin();
    for(; fi != result.faces_end(); ++fi)
      update_envelope_surfaces_by_decision(fi);
  }

  // update the is_equal/has_equal flags of the result envelope
  void update_edge_face_flags(Halfedge_handle h)
  {
    bool is_equal, has_equal;
	  is_equal = (h->get_decision() == h->face()->get_decision());
    // has equal can be true even if the decision is not the same,
	  // but has same surfaces, i.e. one of the features got BOTH
	  // decision, and the other didn't
	  has_equal = (h->get_decision() == h->face()->get_decision() ||
		             h->get_decision() == BOTH ||
				         h->face()->get_decision() == BOTH);

    CGAL::Dac_decision decision = h->face()->get_decision();
    bool is_equal_first = (h->get_is_equal_aux_data_in_face(0));
    bool has_equal_first = (h->get_has_equal_aux_data_in_face(0));
    bool is_equal_second = (h->get_is_equal_aux_data_in_face(1));
    bool has_equal_second = (h->get_has_equal_aux_data_in_face(1));
    if (decision == FIRST)
  	{
      is_equal &= is_equal_first;
	    has_equal &= has_equal_first;
	  }
    else if (decision == SECOND)
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
      if (decision == FIRST)
	      has_equal &= has_equal_first;
      else if (decision == SECOND)
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
	  // but has same surfaces, i.e. one of the features got BOTH
	  // decision, and the other didn't
	  has_equal = (h->get_decision() == h->target()->get_decision() ||
		             h->get_decision() == BOTH ||
				         h->target()->get_decision() == BOTH);

    CGAL::Dac_decision decision = h->get_decision();
    bool is_equal_first = (h->get_is_equal_aux_data_in_target(0));
    bool has_equal_first = (h->get_has_equal_aux_data_in_target(0));
    bool is_equal_second = (h->get_is_equal_aux_data_in_target(1));
    bool has_equal_second = (h->get_has_equal_aux_data_in_target(1));
    if (decision == FIRST)
  	{
      is_equal &= is_equal_first;
	    has_equal &= has_equal_first;
  	}

    else if (decision == SECOND)
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
      if (decision == FIRST)
	      has_equal &= has_equal_first;
      else if (decision == SECOND)
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
	  // but has same surfaces, i.e. one of the features got BOTH
	  // decision, and the other didn't
	  has_equal = (h->face()->get_decision() == h->target()->get_decision() ||
		           h->face()->get_decision() == BOTH ||
				   h->target()->get_decision() == BOTH);

    CGAL::Dac_decision decision = h->face()->get_decision();
    bool has_equal_first = (h->get_has_equal_aux_data_in_target_and_face(0));
    bool has_equal_second = (h->get_has_equal_aux_data_in_target_and_face(1));
    if (decision == FIRST)
	    has_equal &= has_equal_first;
    else if (decision == SECOND)
	    has_equal &= has_equal_second;
    else
	  {
      // we check if the vertex has a different decision, and if so,
	    // we update the flag according to the vertex decision
	    decision = h->target()->get_decision();	
      if (decision == FIRST)
	      has_equal &= has_equal_first;
      else if (decision == SECOND)
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
	  // but has same surfaces, i.e. one of the features got BOTH
	  // decision, and the other didn't

	  has_equal = (v->get_decision() == f->get_decision() ||
		             v->get_decision() == BOTH ||
				         f->get_decision() == BOTH);

    CGAL::Dac_decision decision = v->get_decision();
    bool is_equal_first = (v->get_is_equal_aux_data_in_face(0));
    bool has_equal_first = (v->get_has_equal_aux_data_in_face(0));
    bool is_equal_second = (v->get_is_equal_aux_data_in_face(1));
    bool has_equal_second = (v->get_has_equal_aux_data_in_face(1));
    if (decision == FIRST)
  	{
      is_equal &= is_equal_first;
	    has_equal &= has_equal_first;
	  }
    else if (decision == SECOND)
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
      if (decision == FIRST)
	      has_equal &= has_equal_first;
      else if (decision == SECOND)
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
    for(; hi != result.halfedges_end(); ++hi)
    {
      update_edge_face_flags(hi);
	    update_edge_target_flags(hi);
		update_target_face_flags(hi);
	  }

    // vertices
    Vertex_iterator vi = result.vertices_begin();
    for(; vi != result.vertices_end(); ++vi)
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
	  Object o = fh->get_aux_source(id);
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
    for(; fi != result.faces_end(); ++fi)
    {
      if (fi->is_unbounded())
        continue;
      Face_handle fh = fi;
      Ccb_halfedge_circulator face_hec = fh->outer_ccb();
      Ccb_halfedge_circulator face_hec_begin = face_hec;
      do {
        ++face_hec;
      } while(face_hec != face_hec_begin);

      Hole_iterator inner_iter = fh->holes_begin();
      for (; inner_iter != fh->holes_end(); ++inner_iter)

      {
        face_hec = face_hec_begin = (*inner_iter);
        do {
          ++face_hec;
        } while(face_hec != face_hec_begin);
      }
    }
  }
   
  // check that all aux flags are correctly set, by checking the sets of 
  // surfaces
  bool verify_aux_flags(Minimization_diagram_2& result)
  {
    Envelope_data_iterator begin, end;
    bool all_ok = true;
    Halfedge_iterator hi = result.halfedges_begin();
    for(; hi != result.halfedges_end(); ++hi)
    {
      Halfedge_handle h = hi;
      if (h->get_is_fake())
	      continue;

      Object h_src1 = h->get_aux_source(0);
      Object h_src2 = h->get_aux_source(1);

      // check halfedge-face flags
      Face_handle f = h->face();
      Object f_src1 = f->get_aux_source(0);
      Object f_src2 = f->get_aux_source(1);

      get_data_iterators(h_src1, begin, end);
      all_ok = (is_equal_data(f_src1, begin, end) ==
  	           h->get_is_equal_aux_data_in_face(0));
      if (!all_ok)
      {
        std::cout << "real is_equal = " << is_equal_data(f_src1, begin, end) << std::endl;
  	    std::cout << "flags is_equal = " << h->get_is_equal_aux_data_in_face(0) << std::endl;
      }
      CGAL_assertion(all_ok);

      all_ok = (has_equal_data(f_src1, begin, end) ==
     	            h->get_has_equal_aux_data_in_face(0));
      CGAL_assertion(all_ok);

      get_data_iterators(h_src2, begin, end);
      all_ok = (is_equal_data(f_src2, begin, end) ==
  	        h->get_is_equal_aux_data_in_face(1));
      CGAL_assertion(all_ok);

      if (!all_ok)
      {
        std::cout << "real has equal = " << has_equal_data(f_src2, begin, end) << std::endl;
      	std::cout << "flags has equal = " << h->get_has_equal_aux_data_in_face(1) << std::endl;
      }
      all_ok = (has_equal_data(f_src2, begin, end) ==
  	              h->get_has_equal_aux_data_in_face(1));
      CGAL_assertion(all_ok);

      // check halfedge-target flags
      Vertex_handle v = h->target();
      Object v_src1 = v->get_aux_source(0);
      Object v_src2 = v->get_aux_source(1);

      get_data_iterators(h_src1, begin, end);

      all_ok = (is_equal_data(v_src1, begin, end) ==
  	              h->get_is_equal_aux_data_in_target(0));
      CGAL_assertion(all_ok);

      all_ok = (has_equal_data(v_src1, begin, end) ==
                  h->get_has_equal_aux_data_in_target(0));
      CGAL_assertion(all_ok);

      get_data_iterators(h_src2, begin, end);
      all_ok = (is_equal_data(v_src2, begin, end) ==
  	              h->get_is_equal_aux_data_in_target(1));
      CGAL_assertion(all_ok);

      all_ok = (has_equal_data(v_src2, begin, end) ==
  	              h->get_has_equal_aux_data_in_target(1));
      CGAL_assertion(all_ok);

      // check target-face flags
      get_data_iterators(v_src1, begin, end);
      all_ok = (has_equal_data(f_src1, begin, end) ==
  	            h->get_has_equal_aux_data_in_target_and_face(0));
      CGAL_assertion(all_ok);

      get_data_iterators(v_src2, begin, end);
      all_ok = (has_equal_data(f_src2, begin, end) ==
  	            h->get_has_equal_aux_data_in_target_and_face(1));
      CGAL_assertion(all_ok);
    }

    Vertex_iterator vi = result.vertices_begin();
    for(; vi != result.vertices_end(); ++vi)
      if (vi->is_isolated() && !vi->get_is_fake())
      {
    	// check face flags
        Vertex_handle v = vi;
        Object v_src1 = v->get_aux_source(0);
        Object v_src2 = v->get_aux_source(1);

        Face_handle f = v->face();
        Object f_src1 = f->get_aux_source(0);
        Object f_src2 = f->get_aux_source(1);

        get_data_iterators(v_src1, begin, end);
        all_ok = (is_equal_data(f_src1, begin, end) ==
  	              v->get_is_equal_aux_data_in_face(0));
        if (!all_ok)
      	{
      	  std::cout << "problem with is_equal flags of isolated point "
      		    << vi->point() << std::endl;
      	  std::cout << "real is equal = "
                    << is_equal_data(f_src1, begin, end) << std::endl;
      	  std::cout << "flags is equal = "
                    << v->get_has_equal_aux_data_in_face(0) << std::endl;
      	}
      	CGAL_assertion(all_ok);

        all_ok = (has_equal_data(f_src1, begin, end) ==
  	                v->get_has_equal_aux_data_in_face(0));
        CGAL_assertion(all_ok);

        get_data_iterators(v_src2, begin, end);
        all_ok = (is_equal_data(f_src2, begin, end) ==
  	                v->get_is_equal_aux_data_in_face(1));
        CGAL_assertion(all_ok);

        all_ok = (has_equal_data(f_src2, begin, end) ==
  	                v->get_has_equal_aux_data_in_face(1));
        CGAL_assertion(all_ok);
    }

    return all_ok;
  }

  // check that all flags are correctly set, by checking the sets of surfaces
  bool verify_flags(Minimization_diagram_2& result)
  {
    Envelope_data_iterator begin, end;
	  bool all_ok = true;

    Halfedge_iterator hi = result.halfedges_begin();
	  for(; hi != result.halfedges_end(); ++hi)
  	{
	    Halfedge_handle h = hi;
      Face_handle f = h->face();
	    // check halfedge-face flags
	    all_ok = (h->is_equal_data(f->begin_data(), f->end_data()) ==
		            h->get_is_equal_data_in_face());
	    CGAL_assertion(all_ok);

	    all_ok = (h->has_equal_data(f->begin_data(), f->end_data()) ==
		            h->get_has_equal_data_in_face());

	    CGAL_assertion(all_ok);

	    // check halfedge-target flags
      Vertex_handle v = h->target();

	    all_ok = (h->is_equal_data(v->begin_data(), v->end_data()) ==
		            h->get_is_equal_data_in_target());
	    if (!all_ok)
		    std::cout << "flag value: " << h->get_is_equal_data_in_target()
		              << " real value " << h->is_equal_data(v->begin_data(), v->end_data())
					        << std::endl;
	    CGAL_assertion(all_ok);

	    all_ok = (h->has_equal_data(v->begin_data(), v->end_data()) ==
		            h->get_has_equal_data_in_target());
	    if (!all_ok)
		    std::cout << "flag value: " << h->get_has_equal_data_in_target()
		              << " real value " << h->has_equal_data(v->begin_data(), v->end_data())
					        << std::endl;
	    CGAL_assertion(all_ok);

		// check target-face flags
	    all_ok = (h->face()->has_equal_data(v->begin_data(), v->end_data()) ==
		          h->get_has_equal_data_in_target_and_face());
	    if (!all_ok)
		    std::cout << "flag value: " << h->get_has_equal_data_in_target_and_face()
		              << " real value " << h->face()->has_equal_data(v->begin_data(), v->end_data())
					        << std::endl;
	    CGAL_assertion(all_ok);

	  }

    Vertex_iterator vi = result.vertices_begin();
	  for(; vi != result.vertices_end(); ++vi)
      if (vi->is_isolated())
	    {
		    // check face flags
        Vertex_handle v = vi;

        Face_handle f = v->face();

  	    all_ok = (v->is_equal_data(f->begin_data(), f->end_data()) ==
  		            v->get_is_equal_data_in_face());
  	    CGAL_assertion(all_ok);

  	    all_ok = (v->has_equal_data(f->begin_data(), f->end_data()) ==
  		            v->get_has_equal_data_in_face());
        if (!all_ok)
		      std::cout << "flag value: " << v->get_has_equal_data_in_face()
		                << " real value " << v->has_equal_data(f->begin_data(), f->end_data())
					          << std::endl;
	      CGAL_assertion(all_ok);
	    }

    return all_ok;
  }

  // confirm that aux source and decision are set over all minimization 
  // diagram features
  bool check_resolve_was_ok(Minimization_diagram_2& result)
  {
    bool all_ok = true;
    Vertex_iterator vi = result.vertices_begin();
    for(; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
  	  if (vh->get_is_fake())
	    	continue;

      
      all_ok &= (vh->get_aux_is_set(0));
      CGAL_assertion_msg(all_ok, "aux source (0) not set over vertex");
      all_ok &= (vh->get_aux_is_set(1));
      CGAL_assertion_msg(all_ok, "aux source (1) not set over vertex");
      all_ok &= (vh->is_decision_set());
      CGAL_assertion_msg(all_ok, "decision was not set over vertex");

    }
    Halfedge_iterator hi = result.halfedges_begin();
    for(; hi != result.halfedges_end(); ++hi)
    {
      Halfedge_handle hh = hi;
      if (hh->get_is_fake())
      	continue;
      
      all_ok &= (hh->get_aux_is_set(0));
      CGAL_assertion_msg(all_ok, "aux source (0) not set over edge");
      all_ok &= (hh->get_aux_is_set(1));
      CGAL_assertion_msg(all_ok, "aux source (1) not set over edge");
      all_ok &= (hh->is_decision_set());
      CGAL_assertion_msg(all_ok, "decision was not set over edge");
    }
    Face_iterator fi = result.faces_begin();
    for(; fi != result.faces_end(); ++fi)
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
    for(; vi != result.vertices_end(); ++vi)
    {
      Vertex_handle vh = vi;
      
      all_ok &= (vh->get_is_set());
      CGAL_assertion_msg(all_ok, "data not set over vertex");
      all_ok &= (!vh->has_no_data());

      CGAL_assertion_msg(all_ok, "data empty over vertex");      
      all_ok &= (!vh->get_is_fake());
      CGAL_assertion_msg(all_ok, "fake vertex in envelope");
    }
    Halfedge_iterator hi = result.halfedges_begin();
    for(; hi != result.halfedges_end(); ++hi)
    {
      Halfedge_handle hh = hi;
      
      all_ok &= (hh->get_is_set());
      if (!all_ok)
        std::cout << "edge: " << hh->curve() << std::endl;
      CGAL_assertion_msg(all_ok, "data not set over edge");
      all_ok &= (!hh->has_no_data());
      if (!all_ok)

        std::cout << "edge: " << hh->curve() << std::endl;
      CGAL_assertion_msg(all_ok, "data empty over edge");

      all_ok &= (!hh->get_is_fake());
      CGAL_assertion_msg(all_ok, "fake edge in envelope");
    }
    Face_iterator fi = result.faces_begin();
    for(; fi != result.faces_end(); ++fi)
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
                                  bool is_hole)
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
  // sets the relevant data in the new edges of the decomposition
  class Decomposition_observer : public Md_observer
  {
  public:
    typedef typename Minimization_diagram_2::Halfedge_handle Halfedge_handle;
    typedef typename Envelope_divide_and_conquer_3<Traits, 
                                                   Minimization_diagram_2,
                                                   Vertical_decomposition_2,
                                                   EnvelopeResolver_3, 
                                                   Overlay_2>::Self Self;

    Decomposition_observer(Minimization_diagram_2& arr, Self* b = NULL) :
      Md_observer(arr), base(b)
    { CGAL_assertion(base); }

    virtual void after_create_edge(Halfedge_handle e)
    {
      e->set_is_fake(true);
      e->twin()->set_is_fake(true);

      e->set_aux_source(0, e->face()->get_aux_source(0));
      e->set_aux_source(1, e->face()->get_aux_source(1));

      e->twin()->set_aux_source(0, e->twin()->face()->get_aux_source(0));
      e->twin()->set_aux_source(1, e->twin()->face()->get_aux_source(1));

      //  set decision, if possible
      if (base->aux_has_no_data(e, 0) && !base->aux_has_no_data(e, 1))
      {
        e->set_decision(SECOND);
        e->twin()->set_decision(SECOND);
      }
      else if (!base->aux_has_no_data(e, 0) && base->aux_has_no_data(e, 1))
      {
  	    e->set_decision(FIRST);
        e->twin()->set_decision(FIRST);
      }
      else if(base->aux_has_no_data(e, 0) && base->aux_has_no_data(e, 1))
      {
  	    e->set_decision(EQUAL);
        e->twin()->set_decision(EQUAL);
      } 

      // set halfedge-face flags
      e->set_is_equal_aux_data_in_face(0, true);
      e->twin()->set_is_equal_aux_data_in_face(0, true);
      e->set_is_equal_aux_data_in_face(1, true);
      e->twin()->set_is_equal_aux_data_in_face(1, true);

      e->set_has_equal_aux_data_in_face(0, !base->aux_has_no_data(e->face(), 0));
      e->twin()->set_has_equal_aux_data_in_face(0,!base->aux_has_no_data(e->face(), 0));
      e->set_has_equal_aux_data_in_face(1, !base->aux_has_no_data(e->face(), 1));
      e->twin()->set_has_equal_aux_data_in_face(1, !base->aux_has_no_data(e->face(), 1));

      // we need to set halfedge-target is_equal flags
      // we need the other halfedge that points to the face and to the vertex 
      // if exists (this is the twin's prev halfedge), or the isolated vertex info
	  // we also set has equal flags
      if (e->twin()->prev() != e)
      {
      	CGAL_assertion(e->twin()->face() == e->twin()->prev()->face());
        Halfedge_handle prev = e->twin()->prev();
        e->set_is_equal_aux_data_in_target(0, prev->get_is_equal_aux_data_in_face(0) &&
                         				      prev->get_is_equal_aux_data_in_target(0));
        e->set_is_equal_aux_data_in_target(1, prev->get_is_equal_aux_data_in_face(1) &&
                                              prev->get_is_equal_aux_data_in_target(1));

		e->set_has_equal_aux_data_in_target(0, 
			    prev->get_has_equal_aux_data_in_target_and_face(0));
		e->set_has_equal_aux_data_in_target(1, 
			    prev->get_has_equal_aux_data_in_target_and_face(1));

		e->set_has_equal_aux_data_in_target_and_face(0, 
			    prev->get_has_equal_aux_data_in_target_and_face(0));
		e->set_has_equal_aux_data_in_target_and_face(1, 
			    prev->get_has_equal_aux_data_in_target_and_face(1));

      }
      else
      {
       	// the target of e was isolated, before we added e
       	e->set_is_equal_aux_data_in_target(0, 
			e->target()->get_is_equal_aux_data_in_face(0));
       	e->set_is_equal_aux_data_in_target(1, 
			e->target()->get_is_equal_aux_data_in_face(1));

       	e->set_has_equal_aux_data_in_target(0, 
			e->target()->get_has_equal_aux_data_in_face(0));
       	e->set_has_equal_aux_data_in_target(1, 
			e->target()->get_has_equal_aux_data_in_face(1));

       	e->set_has_equal_aux_data_in_target_and_face(0, 
			e->target()->get_has_equal_aux_data_in_face(0));
       	e->set_has_equal_aux_data_in_target_and_face(1, 
			e->target()->get_has_equal_aux_data_in_face(1));
      }	

      if (e->prev() != e->twin())
      {
       	CGAL_assertion(e->face() == e->prev()->face());

       	Halfedge_handle prev = e->prev();
       	e->twin()->set_is_equal_aux_data_in_target(0, prev->get_is_equal_aux_data_in_face(0) &&
			 		      prev->get_is_equal_aux_data_in_target(0));
      	e->twin()->set_is_equal_aux_data_in_target(1, prev->get_is_equal_aux_data_in_face(1) &&
					      prev->get_is_equal_aux_data_in_target(1));
	
		e->twin()->set_has_equal_aux_data_in_target(0, 
			    prev->get_has_equal_aux_data_in_target_and_face(0));
		e->twin()->set_has_equal_aux_data_in_target(1, 
			    prev->get_has_equal_aux_data_in_target_and_face(1));

		e->twin()->set_has_equal_aux_data_in_target_and_face(0, 
			    prev->get_has_equal_aux_data_in_target_and_face(0));
		e->twin()->set_has_equal_aux_data_in_target_and_face(1, 
			    prev->get_has_equal_aux_data_in_target_and_face(1));
      }
      else
      {
      	// the source of e was isolated, before we added e
      	e->twin()->set_is_equal_aux_data_in_target(0,
                               e->source()->get_is_equal_aux_data_in_face(0));
      	e->twin()->set_is_equal_aux_data_in_target(1,
                               e->source()->get_is_equal_aux_data_in_face(1));

       	e->twin()->set_has_equal_aux_data_in_target(0, 
			                   e->source()->get_has_equal_aux_data_in_face(0));
       	e->twin()->set_has_equal_aux_data_in_target(1, 
			                   e->source()->get_has_equal_aux_data_in_face(1));

       	e->twin()->set_has_equal_aux_data_in_target_and_face(0, 
			                   e->source()->get_has_equal_aux_data_in_face(0));
       	e->twin()->set_has_equal_aux_data_in_target_and_face(1, 
			                   e->source()->get_has_equal_aux_data_in_face(1));
      }	

      // we don't set the halfedge-target has_equal flags, because the setting will not 
      // improve performance (it will not save geometric operations)

      // TODO: if we set correctly all has_equal falgs, maybe we can get rid
      // of "fake" flags
      // TODO: if a fake edge overlaps a projected intersection, and thus
      // becomes not fake (and we will not want to remove it at the end) -
      // what happens in the code? check and fix!
    }
  protected:
    Self *base;
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
                                                   Vertical_decomposition_2,
                                                   EnvelopeResolver_3, 
                                                   Overlay_2>::Self Self;
    Keep_edge_data_observer(Minimization_diagram_2& arr,
			    Self* b = NULL) :
      Md_observer(arr), base(b)
    {
      CGAL_assertion(base);
    }

    virtual void before_split_edge (Halfedge_handle e,
                                    Vertex_handle v,
                                    const X_monotone_curve_2& c1,
                                    const X_monotone_curve_2& c2)
    {
    }

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

      new_he->set_is_fake(org_he->get_is_fake());
      new_he->twin()->set_is_fake(org_he->get_is_fake());
      new_vertex->set_is_fake(org_he->get_is_fake());

      // update all new bools
      new_he->set_is_equal_aux_data_in_face(0, org_he->get_is_equal_aux_data_in_face(0));
      new_he->twin()->set_is_equal_aux_data_in_face(0, org_he->twin()->get_is_equal_aux_data_in_face(0));
      new_he->set_is_equal_aux_data_in_face(1, org_he->get_is_equal_aux_data_in_face(1));
      new_he->twin()->set_is_equal_aux_data_in_face(1, org_he->twin()->get_is_equal_aux_data_in_face(1));

      new_he->set_has_equal_aux_data_in_face(0, org_he->get_has_equal_aux_data_in_face(0));
      new_he->twin()->set_has_equal_aux_data_in_face(0, org_he->twin()->get_has_equal_aux_data_in_face(0));
      new_he->set_has_equal_aux_data_in_face(1, org_he->get_has_equal_aux_data_in_face(1));
      new_he->twin()->set_has_equal_aux_data_in_face(1, org_he->twin()->get_has_equal_aux_data_in_face(1));

      // new_he->target is the original edge's target, and org_he->target is the new vertex
      new_he->set_is_equal_aux_data_in_target(0, org_he->get_is_equal_aux_data_in_target(0));
      new_he->set_is_equal_aux_data_in_target(1, org_he->get_is_equal_aux_data_in_target(1));
      org_he->set_is_equal_aux_data_in_target(0, true);
      org_he->set_is_equal_aux_data_in_target(1, true);
      new_he->set_has_equal_aux_data_in_target(0, org_he->get_has_equal_aux_data_in_target(0));
      new_he->set_has_equal_aux_data_in_target(1, org_he->get_has_equal_aux_data_in_target(1));
      org_he->set_has_equal_aux_data_in_target(0, !base->aux_has_no_data(org_he, 0));
      org_he->set_has_equal_aux_data_in_target(1, !base->aux_has_no_data(org_he, 1));
      new_he->set_has_equal_aux_data_in_target_and_face
		          (0, org_he->get_has_equal_aux_data_in_target_and_face(0));
      new_he->set_has_equal_aux_data_in_target_and_face
		          (1, org_he->get_has_equal_aux_data_in_target_and_face(1));
      org_he->set_has_equal_aux_data_in_target_and_face
		          (0, org_he->get_has_equal_aux_data_in_face(0));
      org_he->set_has_equal_aux_data_in_target_and_face
		          (1, org_he->get_has_equal_aux_data_in_face(1));

      // new_he->source is the new vertex, and org_he->source is the original vertex
      new_he->twin()->set_is_equal_aux_data_in_target(0, true);
      new_he->twin()->set_is_equal_aux_data_in_target(1, true);
      new_he->twin()->set_has_equal_aux_data_in_target(0,!base->aux_has_no_data(org_he, 0));
      new_he->twin()->set_has_equal_aux_data_in_target(1,!base->aux_has_no_data(org_he, 1));

      new_he->twin()->set_has_equal_aux_data_in_target_and_face
		          (0, org_he->twin()->get_has_equal_aux_data_in_face(0));
      new_he->twin()->set_has_equal_aux_data_in_target_and_face
		          (1, org_he->twin()->get_has_equal_aux_data_in_face(1));
    }
  protected:
    Self *base;
  };

  // this observer counts the number of merge_face events
  // in order to count the number of unnecessary faces that were computed
  class Stats_observer : public Md_observer
  {
  public:
    typedef typename Minimization_diagram_2::Face_handle Face_handle;

    Stats_observer(Minimization_diagram_2& arr) :
      Md_observer(arr), counter(0)
    {}

    virtual void after_merge_face (Face_handle /* f */)
    {
      ++counter;      
    }

    unsigned int get_counter()
    {
      return counter;
    }
    
  protected:
    unsigned int counter;
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
                                                 Vertical_decomposition_2,
                                                 EnvelopeResolver_3,
                                                 Overlay_2>::Self Self;

  protected:
    const IndexMap          *index_map; // Mapping vertices to indices
    std::list<Face_handle>& faces;      // The ordered faces list
    Self                    *base;

  public:

    // Constructor.
    Faces_order_bfs_visitor(const IndexMap& imap,
                            std::list<Face_handle>& f,
                            Self* b = NULL) :
      index_map (&imap),
      faces(f),
      base(b)
    {
      CGAL_assertion(base);
    }

    // Write the discover time for a given vertex.
    template <typename Vertex, typename Graph>

    void discover_vertex (Vertex fh, const Graph& g)
    {
      // first we check if we can set the decision immediately
      // if a surface of one map doesn't exist, then we set the second surface
      if (base->aux_has_no_data(fh, 0) && !base->aux_has_no_data(fh, 1))
        fh->set_decision(SECOND);
      else if (base->aux_has_no_data(fh, 0) && base->aux_has_no_data(fh, 1))
      {

        fh->set_decision(EQUAL);
        fh->set_no_data();
      }
      else if (!base->aux_has_no_data(fh, 0) && base->aux_has_no_data(fh, 1))
        fh->set_decision(FIRST);
      else
        // here, we have both surfaces.
        // we save the face in a list for a later treatment, because the face can change
        // and destroy the iterator
        faces.push_back(fh);
    }
  };
#endif
  
protected:
  Overlay_2                 overlay;
  Vertical_decomposition_2  vertical_decomposition;
  Envelope_resolver         *resolver;
  Traits                    *traits;
  // Should we evetually free the traits object
  bool                      own_traits; 
};

CGAL_END_NAMESPACE

#endif //CGAL_ENVELOPE_DIVIDE_AND_CONQUER_3_H
