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
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_PARTIAL_VERTICAL_DECOMPOSITION_2_H
#define CGAL_PARTIAL_VERTICAL_DECOMPOSITION_2_H

//#define CGAL_DEBUG_PARTIAL_VD

#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Partial_vd_visitor.h>
#include <CGAL/Partial_vd_meta_traits.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Timer.h>
#include <vector>
#include <iostream>

CGAL_BEGIN_NAMESPACE

// To use this partial vertical decomposition, the arrangement's traits
// should supply 2 additional methods:
// 1. construct_vertical_2 - to construct a vertical X_monotone_curve_2
//    from 2 points with same x coordinate
// 2. vertical_ray_shoot_2 - to get a point on a X_monotone_curve_2 with
//    a given x coordinate (assuming it is in the curve's x-range)


template <class Arrangement_>
class Partial_vertical_decomposition_2
{
public:
  typedef Arrangement_                                 Arrangement;
  // Arrangement types:
  typedef typename Arrangement::Traits_2               Traits_2;
  typedef typename Traits_2::X_monotone_curve_2        Base_X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Base_Point_2;

  typedef typename Arrangement::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement::Vertex_const_iterator  Vertex_const_iterator;
  typedef typename Arrangement::Edge_const_iterator    Edge_const_iterator;
  typedef typename Arrangement::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Arrangement::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement::Vertex_handle          Vertex_handle;
  typedef typename Arrangement::Size                   Size;


  // Define meta-traits class for the batched point location:
  typedef Partial_vd_meta_traits<Traits_2, Arrangement>
                                                         Meta_traits_2;

  typedef typename Meta_traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Meta_traits_2::Point_2                Point_2;

  typedef std::pair<Object, Object>                      Vd_pair;
  typedef std::list<Vd_pair>                             Vd_pairs_container;
  typedef std::list<Vd_pair>::iterator                   Vd_pairs_iter;
  typedef std::back_insert_iterator<Vd_pairs_container>  Vd_pairs_oiter;
  // Define the sweep-line visitor:
  typedef Partial_vd_visitor<Meta_traits_2,
                             Arrangement,
                             Vd_pairs_oiter>             Visitor;

  typedef Basic_sweep_line_2<Meta_traits_2, Visitor>     Sweep_line;

  typedef Unique_hash_map<Halfedge_handle, Halfedge_handle> Halfedges_map;

  // Do a partial vertical decomposition on existing arrangement "arr"
  void operator()(Arrangement& arr)
  {
    #ifdef CGAL_DEBUG_PARTIAL_VD  
      cout << "before partial vd: print edges of arr" << endl;
      for(Halfedge_const_iterator hi = arr.halfedges_begin();
          hi != arr.halfedges_end(); ++hi, ++hi)
          cout << hi->curve() << endl;
      cout << "before partial vd: print isolated vertices of arr" << endl;
      for(Vertex_const_iterator vi = arr.vertices_begin();
          vi != arr.vertices_end(); ++vi)
          if (vi->is_isolated())
            cout << vi->point() << endl;
    #endif
    
    // Go over all arrangement edges.
    std::vector<X_monotone_curve_2>  xcurves_vec;
    xcurves_vec.resize(arr.number_of_edges());
    typename Traits_2::Compare_xy_2    comp_xy =
      arr.get_traits()->compare_xy_2_object();
    Edge_const_iterator       eit;

    Size i = 0;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit, ++i)
    {
      // Associate each x-monotone curve with the halfedge that represent it
      // that is directed from left to right.
      if(comp_xy(eit->source()->point(),
  	             eit->target()->point()) == SMALLER)
        xcurves_vec[i] = X_monotone_curve_2(eit->curve(),eit);
      else
        xcurves_vec[i] = X_monotone_curve_2(eit->curve(),eit->twin());
    }

    // Associate each isolated point with the vertex that represents it
    std::vector<Point_2>    iso_points;
    for(Vertex_const_iterator v_itr = arr.vertices_begin();
        v_itr != arr.vertices_end();
        ++v_itr)
    {
      if(v_itr->is_isolated())
        iso_points.push_back(Point_2(v_itr->point(), v_itr));
    }

    // Perform the sweep
    Vd_pairs_container vd_pairs;
    Visitor     visitor (arr, std::back_inserter(vd_pairs));
    Sweep_line  sweep_line (&visitor);

    sweep_timer.start();
    sweep_line.sweep(xcurves_vec.begin(),
  		               xcurves_vec.end(),
                     iso_points.begin(),
                     iso_points.end());
    sweep_timer.stop();

    add_timer.start();
    _add_vertical_edges(arr, vd_pairs);
    add_timer.stop();
    
    #ifdef CGAL_DEBUG_PARTIAL_VD
      std::cout << "finish partial vd" << std::endl;
      print_times();
    #endif
  }
protected:
  // add the vertical edges (that were determined by the sweep) to the
  // arrangement
  void _add_vertical_edges(Arrangement& arr, Vd_pairs_container& vd_pairs)
  {
    Vertex_const_handle invalid_v(NULL);
    Halfedge_const_handle invalid_he(NULL);

    Halfedge_handle prev_split_he(NULL);
    Base_Point_2 prev_split_pt;
    Vertex_handle prev_split_v(NULL);

    // map original halfedge in the arrangement to its current rightmost part
    // which should be split when more than one split of this halfedge is needed
    Halfedges_map map_orig_to_rightmost;
        
    Vd_pairs_iter it = vd_pairs.begin();
    for(; it != vd_pairs.end(); ++it)
    {
      Vd_pair cur_pair = *it;
      Vertex_const_handle v1, v2;
      Halfedge_const_handle h1, h2;
      if (!CGAL::assign(v1, cur_pair.first))
      {
        CGAL_assertion(CGAL::assign(h1, cur_pair.first));
        CGAL::assign(h1, cur_pair.first);
      }
      if (!CGAL::assign(v2, cur_pair.second))
      {
        CGAL_assertion(CGAL::assign(h2, cur_pair.second));
        CGAL::assign(h2, cur_pair.second);
      }
      // we have 2 vertices, no split is needed
      if (v1 != invalid_v && v2 != invalid_v)
      {
        if (should_add_vertical_edge(arr.non_const_handle(v1), arr.non_const_handle(v2)))
        {
          #ifdef CGAL_DEBUG_PARTIAL_VD
            std::cout << "got vertex-vertex pair: " << v1->point()
                      << " and " << v2->point() << std::endl;
          #endif
          insert_timer.start();
          #ifdef CGAL_DEBUG_PARTIAL_VD
            std::cout << "before insert at vertices: " << v1->point() << " , "
                      << v2->point() << std::endl;
          #endif
          arr.insert_at_vertices(
               arr.get_traits()->construct_vertical_2_object()(v1->point(), v2->point()),
               arr.non_const_handle(v1),
               arr.non_const_handle(v2));
          insert_timer.stop();
        }
        continue;
      }
      Vertex_handle v;
      Halfedge_handle orig_split_he, split_he;
      
      if (v1 != invalid_v)
      {
        // we must have h2 valid
        CGAL_assertion(h2 != invalid_he);
        v = arr.non_const_handle(v1);
        orig_split_he = arr.non_const_handle(h2);
      }
      else
      {
        // we must have v2 and h1 valid
        CGAL_assertion(v2 != invalid_v);
        CGAL_assertion(h1 != invalid_he);
        v = arr.non_const_handle(v2);
        orig_split_he = arr.non_const_handle(h1);
      }

      if (!should_add_vertical_edge(v, orig_split_he))
        continue;
        
      #ifdef CGAL_DEBUG_PARTIAL_VD
        std::cout << "got vertex-halfedge pair: " << v->point()
                  << " and " << orig_split_he->curve() << std::endl;
      #endif      
      // split split_he and connect the split point with v
      if (map_orig_to_rightmost.is_defined(orig_split_he))
        // we should split the rightmost halfedge instead
        split_he = map_orig_to_rightmost[orig_split_he];
      else
        split_he = orig_split_he;

      Base_Point_2 split_p;
      Vertex_handle split_v;
      if (prev_split_he == orig_split_he &&
          arr.get_traits()->compare_x_2_object()(prev_split_pt, v->point()) == EQUAL)
      {
        split_p = prev_split_pt;
        split_v = prev_split_v;
      }
      else
      {
        split_p = arr.get_traits()->vertical_ray_shoot_2
                                         (v->point(), split_he->curve());
        Base_X_monotone_curve_2 a,b;
        #ifdef CGAL_DEBUG_PARTIAL_VD
          std::cout << "before edge_split, curve=" << split_he->curve()
                    << std::endl << " point= " << split_p << std::endl;
        #endif
        arr.get_traits()->split_2_object()(split_he->curve(), split_p, a, b);
        split_he = arr.split_edge(split_he, a, b);
        // split always returns the halfedge with source =  original source
        // so the current rightmost part is split_he->next()
        map_orig_to_rightmost[orig_split_he] = split_he->next();
        split_v = split_he->target();

      }

      prev_split_he = orig_split_he;
      prev_split_pt = split_p;
      prev_split_v = split_v;                                   

      // insert the vertical edge
      insert_timer.start();
      #ifdef CGAL_DEBUG_PARTIAL_VD
        std::cout << "before insert at vertices: " << v->point() << " , "
                  << split_p << std::endl;
      #endif
      arr.insert_at_vertices(
            arr.get_traits()->construct_vertical_2_object()(v->point(), split_p),
            v,
            split_v);
      insert_timer.stop();
//      insert_x_monotone(arr, Base_X_monotone_curve_2(v->point(), split_p));   
    }
    vd_pairs.clear();
  }

  bool should_add_vertical_edge(Vertex_handle v1, Vertex_handle v2)
  {
    return true;
//    return (!v1->get_is_intersection() || !v2->get_is_intersection());
  }
  bool should_add_vertical_edge(Vertex_handle v, Halfedge_handle he)
  {
    return true;
//    return (!v->get_is_intersection());
  }
  void print_times()
  {
      std::cout << "Partial vd times: " << std::endl;
      std::cout << "sweep: " << sweep_timer.time() << " seconds" << std::endl;
      std::cout << "add vertical edges: " << add_timer.time() << " seconds" << std::endl;
      std::cout << "insert edges: " << insert_timer.time() << " seconds" << std::endl;
  }

protected:
  mutable Timer sweep_timer;
  mutable Timer add_timer;
  mutable Timer insert_timer;

};

CGAL_END_NAMESPACE

#endif
