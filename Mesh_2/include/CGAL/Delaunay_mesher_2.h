// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_MESHER_2_H
#define CGAL_DELAUNAY_MESHER_2_H

#include <CGAL/license/Mesh_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_2/Refine_edges_with_clusters.h>
#include <CGAL/Mesh_2/Refine_edges_visitor.h>
#include <CGAL/Mesh_2/Refine_faces.h>

namespace CGAL {

template <typename Tr, typename Crit>
class Delaunay_mesher_2 
{

  /** \name \c Tr types */
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Edge Edge;

  typedef typename Tr::Point Point;

  /** \name Types needed for private member datas */
  typedef Mesh_2::Refine_edges_with_clusters<Tr,
    Mesh_2::Is_locally_conforming_Gabriel<Tr> > Edges_level;

  typedef Mesh_2::Refine_faces<Tr, Crit, Edges_level> Faces_level;

public:

  typedef Tr Triangulation;
  typedef Crit Criteria;
  typedef Mesh_2::Clusters<Tr> Clusters;

  /** \name Types needed to mark the domain for Delaunay_mesh_2<Tr> */

  typedef std::list<Point> Seeds;
  typedef typename Seeds::const_iterator Seeds_iterator;
  typedef Seeds_iterator Seeds_const_iterator;

private:
  // --- PRIVATE MEMBER DATAS ---
  Tr& tr;
  Criteria criteria;
  Null_mesher_level null_level;
  Null_mesh_visitor null_visitor;
  Clusters clusters_;
  Edges_level edges_level;
  Faces_level faces_level;
  Mesh_2::Refine_edges_visitor_from_faces<Faces_level> visitor;

  Seeds seeds;
  bool seeds_mark;
public:
  /** \name CONSTRUCTORS */
  Delaunay_mesher_2(Tr& tr_, const Criteria& criteria_ = Criteria())
    : tr(tr_),
      criteria(criteria_), 
      null_level(),
      null_visitor(),
      clusters_(tr),
      edges_level(tr, clusters_, null_level),
      faces_level(tr, criteria, edges_level),
      visitor(faces_level, edges_level, null_visitor),
      initialized(false)
  {
  }

  Delaunay_mesher_2(Tr& tr_, Edges_level& edges_level_,
                    const Criteria& criteria_ = Criteria())
    : tr(tr_),
      criteria(criteria_), 
      null_level(),
      null_visitor(),
      clusters_(tr),
      edges_level(edges_level_),
      faces_level(tr, criteria, edges_level),
      visitor(faces_level, null_visitor),
      initialized(false)
  {
  }

  /** SEEDS HANDLING FUNCTIONS */

  Seeds_const_iterator seeds_begin() const
  {
    return seeds.begin();
  }
  
  Seeds_const_iterator seeds_end() const
  {
    return seeds.end();
  }

private:
  /** \name INITIALIZED */

  bool initialized;

public:
  /** \name MARKING FUNCTIONS */

  /** The value type of \a InputIterator should be \c Point, and represents
      seeds. Connected components of seeds are marked with the value of
      \a mark. Other components are marked with \c !mark. The connected
      component of infinite faces is always marked with \c false.
  */
  template <class InputIterator>
  void set_seeds(InputIterator b, InputIterator e,
                 const bool mark = false,
                 const bool do_it_now = false)
  {
    seeds.clear();
    std::copy(b, e, std::back_inserter(seeds));
    seeds_mark=mark;
    if(do_it_now) mark_facets();
  }

  void clear_seeds()
  {
    seeds.clear();
    seeds_mark = false;
  }

  void mark_facets(bool domain_specified = false)
  {
    if(!domain_specified) {
      mark_facets(tr, seeds.begin(), seeds.end(), seeds_mark);
    }
    else {
      propagate_marks(tr.infinite_face(), false);
    }
  }

  /** Procedure that marks facets according to a list of seeds. */
  template <typename Seeds_it>
  static void mark_facets(Tr& tr, 
                          Seeds_it begin,
                          Seeds_it end,
                          bool mark = false)
  {
    if (tr.dimension()<2) return;
    if( begin != end )
      {
        for(typename Tr::All_faces_iterator it=tr.all_faces_begin();
            it!=tr.all_faces_end();
            ++it)
          it->set_in_domain(!mark);

        for(Seeds_it sit=begin; sit!=end; ++sit)
          {
            Face_handle fh=tr.locate(*sit);
            if(fh!=NULL)
              propagate_marks(fh, mark);
          }
	propagate_marks(tr.infinite_face(), false);
      }
    else
      mark_convex_hull(tr);
  }

  /**
   * Marks all faces of the convex hull but those connected to the
   * infinite faces.
   */
  static void mark_convex_hull(Tr& tr)
  {
    for(typename Tr::All_faces_iterator fit=tr.all_faces_begin();
        fit!=tr.all_faces_end();
        ++fit)
      fit->set_in_domain(true);
    propagate_marks(tr.infinite_face(), false);
  }

  /** Propagates the mark \c mark recursivly. */
  static void propagate_marks(const Face_handle fh, bool mark)
  {
    // std::queue only works with std::list on VC++6, and not with
    // std::deque, which is the default
    // But it should be fixed by VC++7 know. [Laurent Rineau 2003/03/24]
    std::queue<Face_handle/*, std::list<Face_handle>*/> face_queue;
    fh->set_in_domain(mark);
    face_queue.push(fh);
    while( !face_queue.empty() )
      {
        Face_handle fh = face_queue.front();
        face_queue.pop();
        for(int i=0;i<3;i++)
          {
            const Face_handle& nb = fh->neighbor(i);
            if( !fh->is_constrained(i) && (mark != nb->is_in_domain()) )
              {
                nb->set_in_domain(mark);
                face_queue.push(nb);
              }
          }
      }
  }
  /** \name MESHING FUNCTIONS */

  void refine_mesh()
  {
    if(initialized != true) init();
    faces_level.refine(visitor);
  }

  /** \name REMESHING FUNCTIONS */

  void set_criteria(const Criteria& criteria_,
                    bool recalculate_bad_faces = true)
  {
    criteria = criteria_;
    if (recalculate_bad_faces) faces_level.scan_triangulation();  
  }

  const Criteria& get_criteria() const 
  {
    return criteria;  
  }

  template <class Fh_it>
  void set_bad_faces(Fh_it begin, Fh_it end)
  {
    faces_level.set_bad_faces(begin, end);
  }

  /** \name STEP BY STEP FUNCTIONS */

  /**
     Initialize the data structures
     (The call of this function is REQUIRED before any step by step
     operation).
  */
  void init(bool domain_specified = false)
  {
    mark_facets(domain_specified);
    clusters_.create_clusters();
    edges_level.scan_triangulation();
    faces_level.scan_triangulation();
    initialized = true;
  }

  bool is_refinement_done ()
  {
    return faces_level.is_algorithm_done();  
  }

  bool
  step_by_step_refine_mesh()
  {
    return faces_level.try_to_insert_one_point(visitor);
  }

  bool try_one_step_refine_mesh()
  {
    return faces_level.one_step(visitor);
  }

  /** \name ACCESS FUNCTIONS */

  const Mesh_2::Clusters<Tr>& clusters() const 
  {
    return clusters_;
  }

  const Triangulation& triangulation() const
  {
    return tr;
  }

  /** \name DEBUGGING FUNCTIONS */

  typedef typename Edges_level::Constrained_edge Constrained_edge;

  bool is_edges_refinement_done()
  {
    return edges_level.is_algorithm_done();
  }

  Edge next_encroached_edge() 
  {
    return edges_level.get_next_element();
  }

  const Face_handle next_bad_face() 
  {
    return faces_level.get_next_element();
  }

  const Point next_refinement_point() 
  {
    if( !edges_level.is_algorithm_done() )
      return edges_level.refinement_point(next_encroached_edge());
    else
      return faces_level.refinement_point(next_bad_face());
  }

  typedef typename Edges_level::Edges_const_iterator
    Encroached_edges_const_iterator;

  typedef typename Faces_level::Bad_faces_const_iterator
    Bad_faces_const_iterator;
  
  Encroached_edges_const_iterator encroached_edges_begin() const
  {
    return edges_level.begin();
  }

  Encroached_edges_const_iterator encroached_edges_end() const
  {
    return edges_level.end();
  }

  Bad_faces_const_iterator bad_faces_begin() const
  {
    return faces_level.begin();
  }
  
  Bad_faces_const_iterator bad_faces_end() const
  {
    return faces_level.end();
  }
}; // end class Delaunay_mesher_2

// --- GLOBAL FUNCTIONS ---

template <typename Tr, typename Criteria>
void
refine_Delaunay_mesh_2(Tr& t,
                       const Criteria& criteria = Criteria(), bool domain_specified=false)
{
  typedef Delaunay_mesher_2<Tr, Criteria> Mesher;

  Mesher mesher(t, criteria);
  mesher.init(domain_specified);
  mesher.refine_mesh();
}


template <typename Tr, typename Criteria, typename InputIterator>
void
refine_Delaunay_mesh_2(Tr& t,
                       InputIterator b, InputIterator e,
                       const Criteria& criteria = Criteria(),
                       bool mark = false)
{
  typedef Delaunay_mesher_2<Tr, Criteria> Mesher;

  Mesher mesher(t, criteria);
  mesher.set_seeds(b, e, mark);
  mesher.refine_mesh();
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_DELAUNAY_MESHER_2_H
