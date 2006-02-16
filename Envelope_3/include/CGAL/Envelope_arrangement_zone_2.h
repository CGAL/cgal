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
// $Id$
// 
//
// Author(s)     : Michal Meyerovitch <gorgymic@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_ARRANGEMENT_ZONE_2_H
#define CGAL_ENVELOPE_ARRANGEMENT_ZONE_2_H

/*! \file
 * Defintion of the Envelope_arrangement_zone_2 class.
 */

#include <CGAL/Arrangement_zone_2.h>
#include <CGAL/Arr_observer.h>
#include <list>
#include <map>
#include <set>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class for computing the zone of a given $x$-monotone curve in a given
 * arrangement.
 * The arrangement parameter corresponds to the underlying arrangement, and
 * the zone-visitor parameter corresponds to a visitor class which is capable
 * of receiving notifications on the arrangment features the query curve
 * traverses. The visitor has to support the following functions:
 * - init(), for initializing the visitor with a given arrangement.
 * - found_subcurve(), called when a non-intersecting x-monotone curve is
 *                     computed and located in the arrangement.
 * - found_overlap(), called when an x-monotone curve overlaps an existing
 *                    halfedge in the arrangement.
 * Both the second and the third functions return pair<Halfedge_handle, bool>,
 * where the halfedge handle corresponds to the halfedge created or modified
 * by the visitor (if valid), and the Boolean value indicates whether we
 * should halt the zone-computation process.
 *
 * This class improves the complexity of the zone algorithm of
 * Arrangement_zone_2 for the cases where we enter a face many times
 */
template <class Arrangement_, class ZoneVisitor_>
class Envelope_arrangement_zone_2 :
        public Arrangement_zone_2< Arrangement_, ZoneVisitor_ >,
        public Arr_observer< Arrangement_ >
{
public:

  typedef Arrangement_                            Arrangement_2;
  typedef typename Arrangement_2::Traits_2        Traits_2;

  typedef ZoneVisitor_                            Visitor;

  typedef typename Arrangement_2::Vertex_handle   Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement_2::Face_handle     Face_handle;

  typedef std::pair<Halfedge_handle, bool>        Visitor_result;

  typedef typename Traits_2::Point_2              Point_2;
  typedef typename Traits_2::X_monotone_curve_2   X_monotone_curve_2;

  typedef Arrangement_zone_2< Arrangement_2, Visitor >
                                                  Base_zone_2;
  typedef Arr_observer< Arrangement_2 >           Base_observer;                                                                           
protected:

  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                  Ccb_halfedge_circulator;
  typedef Arr_traits_adaptor_2<Traits_2>          Traits_adaptor_2;

  // Types used for caching intersection points:
  typedef std::pair<Point_2, unsigned int>        Intersect_point_2;
  typedef std::list<CGAL::Object>                 Intersect_list;
  typedef std::map<const X_monotone_curve_2*,
                   Intersect_list>                Intersect_map;
  typedef typename Intersect_map::iterator        Intersect_map_iterator;

  typedef std::set<const X_monotone_curve_2*>     Curves_set;
  typedef typename Curves_set::iterator           Curves_set_iterator;

  // Data member:

  // Set of faces that we dealt with already
  Unique_hash_map<Face_handle, bool>              discovered_faces;
  
  // map halfedges that were split to the right part of it
  // (which is the only part that might be encountered later again.
  class Less_halfedge
  {
  public:
    bool operator()(Halfedge_handle h1, Halfedge_handle h2)
    {
      return (&*h1 < &*h2);
    }
  };
  
  typedef std::map<Halfedge_handle, Halfedge_handle,Less_halfedge> Halfedges_map;
  typedef typename Halfedges_map::iterator           Halfedges_map_iter;
  Halfedges_map                                      split_map;

  // Set of features (halfedges and vertices) that intersect cv
  // (the current portion of the query curve), ordered according to the
  // intersection order on cv from left to right

  class Less_inter_vertex
  {
  public:
    Less_inter_vertex(Traits_adaptor_2* t) : traits(t)
    {}
    
    bool operator() (const Vertex_handle& v1,
                     const Vertex_handle& v2) const
    {
      return (traits->compare_xy_2_object()(v1->point(), v2->point())
              == SMALLER);
    }
  protected:
    Traits_adaptor_2       *traits;
  };
  
  class Less_inter_edge
  {
  public:
    Less_inter_edge(Traits_adaptor_2 *t, Intersect_map& imap) :
      traits(t), zone_inter_map(imap)
    {}

    bool operator() (const Halfedge_handle& h1,
                     const Halfedge_handle& h2) const
    {
      // find the leftmost intersection with each curve
      // (the curves must be in the inter_map)
      const Intersect_point_2  *ip;
      const X_monotone_curve_2 *icv;
      Point_2                   ip1, ip2;
      bool                      no_inter1, no_inter2;

      Intersect_map_iterator iter = zone_inter_map.find (&(h1->curve()));
      CGAL_assertion(iter != zone_inter_map.end());
      if (iter != zone_inter_map.end())
      {
        // Retrieve the intersections list from the map.
        Intersect_list& inter_list = iter->second;

        if (inter_list.empty())
          no_inter1 = true;
        else
        {
          no_inter1 = false;
          // Locate the first intersection 
          // Compare that current object with left_pt.
          ip = object_cast<Intersect_point_2> (&(inter_list.front()));

          if (ip != NULL)
            ip1 = ip->first;
          else
          {
    	      icv = object_cast<X_monotone_curve_2> (&(inter_list.front()));
    	      CGAL_assertion (icv != NULL);

            ip1 = traits->construct_min_vertex_2_object()(*icv);
          }
        }
      }

      iter = zone_inter_map.find (&(h2->curve()));
      CGAL_assertion(iter != zone_inter_map.end());
      if (iter != zone_inter_map.end())
      {
        // Retrieve the intersections list from the map.
        Intersect_list& inter_list = iter->second;

        if (inter_list.empty())
          no_inter2 = true;
        else
        {
          no_inter2 = false;
          // Locate the first intersection
          // Compare that current object with left_pt.
          ip = object_cast<Intersect_point_2> (&(inter_list.front()));

          if (ip != NULL)
            ip2 = ip->first;
          else
          {
    	      icv = object_cast<X_monotone_curve_2> (&(inter_list.front()));
    	      CGAL_assertion (icv != NULL);

            ip2 = traits->construct_min_vertex_2_object()(*icv);
          }
        }
      }

//      if (no_inter1 &&  no_inter2)
//        return (h1 < h2);
//      else if (no_inter1)
//        return false;
//      else if (no_inter2)
//        return true;
//      else
      CGAL_assertion(!no_inter1 && !no_inter2);
      return (traits->compare_xy_2_object()(ip1, ip2) == SMALLER);
    }
    
  protected:
    Traits_adaptor_2       *traits;
    Intersect_map &zone_inter_map;
  };
  
  typedef std::set<Halfedge_handle, Less_inter_edge>   Halfedge_sorted_set;
  typedef typename Halfedge_sorted_set::iterator       He_sorted_set_iter;
  Halfedge_sorted_set                                  intersect_he_sorted_set;
  
  typedef std::set<Vertex_handle, Less_inter_vertex>   Vertices_sorted_set;
  typedef typename Vertices_sorted_set::iterator       V_sorted_set_iter;
  Vertices_sorted_set                                  intersect_v_sorted_set;
  
public:

  /*!
   * Constructor.                                              typename Arrangement_2::Ccb_halfedge_circulator
   * \param _arr The arrangement for which we compute the zone.
   * \param _visitor A pointer to a zone-visitor object.
   */
  Envelope_arrangement_zone_2 (Arrangement_2& _arr,
		      Visitor *_visitor) :
    Base_zone_2(_arr, _visitor),
    Base_observer(_arr),
    intersect_he_sorted_set(Less_inter_edge(traits, inter_map)),
    intersect_v_sorted_set(Less_inter_vertex(traits))   
  {
  }

  virtual ~Envelope_arrangement_zone_2(){}

  /*!
   * Compute the zone of the given curve and issue the apporpriate
   * notifications for the visitor.
   */
  virtual void compute_zone ()
  {
    Base_zone_2::compute_zone();
    intersect_v_sorted_set.clear();
    intersect_he_sorted_set.clear();
    discovered_faces.clear();
    split_map.clear();
  }

  /*!
   * Notification after an edge was split.
   * \param e1 A handle to one of the twin halfedges forming the first edge.
   * \param e2 A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_edge (Halfedge_handle e1,
                                 Halfedge_handle e2)
  {
    // we assume that e1 is the original edge that was split
    // update split_map with the pair e1 and the rightmost
    // halfedge part
    if (e1->direction() == LARGER)
    {
      split_map[e1] = e2;
      split_map[e1->twin()] = e2;
    }
    // in the other case, don't need the map
  }

  /*!
   * Notification after a face was split.
   * \param f A handle to the face we have just split.
   * \param new_f A handle to the new face that has been created.
   * \param is_hole Whether the new face forms a hole inside f.
   */
  virtual void after_split_face (Face_handle f,
                                 Face_handle new_f,
                                 bool is_hole)
  {
    // update the set of discovered faces
    if (discovered_faces.is_defined(f))
      discovered_faces[new_f] = discovered_faces.default_value();
  }

protected:

  /*!
   * Compute the (lexicographically) leftmost intersection of the query
   * curve with the boundary of a given face in the arrangement.
   * The function computes sets intersect_p, intersect_he (or alternatively
   * overlap_cv and intersect_he) and set the flags found_intersect and
   * found_overlap accordingly.
   * \param face A handle to the face.
   * \param on_boundary Specifies whether the left endpoint of the curve lies
   *                    on the face boundary.
   */
  virtual void _leftmost_intersection_with_face_boundary (Face_handle face,
                                                          bool on_boundary)
  {
    found_intersect = false;
    found_overlap = false;
    found_iso_vert = false;

    // Go over the outer boundary of the face (if one exists), and try to
    // locate intersections of cv with the edges along the boundary.
    typename Traits_adaptor_2::Compare_xy_2            compare_xy =
                                        traits->compare_xy_2_object();
    typename Traits_adaptor_2::Is_in_x_range_2         is_in_x_range =
                                        traits->is_in_x_range_2_object();
    typename Traits_adaptor_2::Construct_min_vertex_2  min_vertex =
                                        traits->construct_min_vertex_2_object();
    typename Traits_adaptor_2::Construct_max_vertex_2  max_vertex =
                                        traits->construct_max_vertex_2_object();
    typename Traits_adaptor_2::Compare_y_at_x_2        compare_y_at_x =
                                        traits->compare_y_at_x_2_object();

    //Base_zone_2::_leftmost_intersection_with_face_boundary(face, on_boundary);
    if (!discovered_faces.is_defined(face))
    {
      // find all intersections with face boundary, and insert into the
      // intersection sets
      Ccb_halfedge_circulator  he_first;
      if (! face->is_unbounded())
      {
        // Get circulators for the outer boundary of the face.
        he_first = face->outer_ccb();
        _intersect_with_ccb(he_first, on_boundary);
      }

      typename Arrangement_2::Hole_iterator   holes_it;
      for (holes_it = face->holes_begin();
           holes_it != face->holes_end(); ++holes_it)
      {
        // Get circulators for the boundary of the current hole.
        he_first = *holes_it;
        _intersect_with_ccb(he_first, on_boundary);
      }

      typename Arrangement_2::Isolated_vertex_iterator   iso_verts_it;
      for (iso_verts_it = face->isolated_vertices_begin();
           iso_verts_it != face->isolated_vertices_end(); ++iso_verts_it)
      {
        // If the isolated vertex is not in the x-range of our curve, disregard it.
        if (! is_in_x_range (cv, iso_verts_it->point()))
          continue;

        // In case the isolated vertex lies on the curve, add it to the
        // intersection vertex set
        if (compare_y_at_x (iso_verts_it->point(), cv) == EQUAL)
        {
          intersect_v_sorted_set.insert(iso_verts_it);
        }

      } // End:: traversal of the isolated vertices inside the face.
      
      // mark that we discovered the current face
      discovered_faces[face] = discovered_faces.default_value(); 
    }
    // find the leftmost intersection (it should relate to the current face)

    // should find the halfedge with the leftmost intersection
    // and the isolated point with the leftmost intersection
    // and compare them
    bool he_exist, v_exist;
    // the leftmost intersections
    Point_2 vp, hep;
    V_sorted_set_iter  v_iter;
    He_sorted_set_iter he_iter;
    bool he_inter_is_point;
    Halfedge_handle inter_he;
    
    if (intersect_v_sorted_set.empty())
      v_exist = false;
    else
    {
      v_exist = true;
      v_iter = intersect_v_sorted_set.begin();
      vp = (*v_iter)->point();
    }

    CGAL::Object               obj;
    const Intersect_point_2   *int_p;
    const X_monotone_curve_2  *icv;
    Point_2                    ip;

    if (intersect_he_sorted_set.empty())
      he_exist = false;
    else
    {
      do
      {
        he_iter = intersect_he_sorted_set.begin();
        //  now find the next intersection with this halfedge
        // todo: is the false ok?
        Halfedges_map_iter split_iter = split_map.find(inter_he);
        if (split_iter == split_map.end())
          inter_he = *he_iter;
        else
        {
          Halfedge_handle old = inter_he;
          inter_he = split_iter->second;
          split_map.erase(old);
        }

        obj = _compute_next_intersection (inter_he, false);
        if (obj.is_empty())
          intersect_he_sorted_set.erase(he_iter);

      }
      while(obj.is_empty() && !intersect_he_sorted_set.empty());

      if (intersect_he_sorted_set.empty())
        he_exist = false;
      else
      {
        he_exist = true;
        CGAL_assertion (! obj.is_empty());

        // We have found an intersection (either a simple point or an
        // overlapping x-monotone curve).
	      int_p = object_cast<Intersect_point_2> (&obj);
        if (int_p != NULL)
        {
          hep = int_p->first;
          he_inter_is_point = true;
        }
        else
        {
          // We have located an overlapping curve. Assign ip as its left
          // endpoint.
	        icv = object_cast<X_monotone_curve_2> (&obj);
	        CGAL_assertion (icv != NULL);

          hep = min_vertex (*icv);
          he_inter_is_point = false;
        }
      }
    }
    
    if (!v_exist && !he_exist)
      return; // no intersection at all

    if (v_exist && he_exist)
    {
      // compare the intersections of the vertex and the halfedge
      // and return the leftmost between them
      Comparison_result res = compare_xy(hep, vp);
      CGAL_assertion(res != EQUAL); // imposible intersection between halfedge
      // and isolated vertex

      // change v_exist/he_exist in order for the below conditiona to work
      if (res == SMALLER)
        // halfedge wins
        v_exist = false;
      else
        he_exist = false;
    }
    
    if (v_exist && !he_exist)
    {
      // return the intersection with isolated vertex
      intersect_v = *v_iter;
      intersect_p = intersect_v->point();
      ip_mult = 0;
      found_intersect = true;
      found_iso_vert = true;

      // remove the vertex from the sorted set
      intersect_v_sorted_set.erase(v_iter);

      return;
    }

    if (!v_exist && he_exist)
    {
      // return the intersection with a halfedge 
      if (he_inter_is_point)
      {
        intersect_p = hep;
        ip_mult = int_p->second;
        intersect_he = inter_he;
        found_intersect = true;
      }
      else
      {
        // begin of overlapping curve
        intersect_p = hep;
        ip_mult = 0;
        overlap_cv = *icv;
        intersect_he = inter_he;
        found_overlap = true;
        found_intersect = true;
      }
      // remove the halfedge from the sorted intersection set
      intersect_he_sorted_set.erase(he_iter);
      // remove the found intersection from the list, and if the list is not
      // empty, reenter the halfedge into the sorted set
      _remove_next_intersection(intersect_he);
      Intersect_map_iterator iter = inter_map.find(&(intersect_he->curve()));
      Intersect_list& inter_list = iter->second;
      if (!inter_list.empty())
        intersect_he_sorted_set.insert(intersect_he);
        
      // todo: should remove all halfedges that intersect in the same point from the set   
    }
  }

  void _intersect_with_ccb(Ccb_halfedge_circulator he_first, bool on_boundary)
  {
    typename Traits_adaptor_2::Compare_xy_2            compare_xy =
                                        traits->compare_xy_2_object();
    typename Traits_adaptor_2::Is_in_x_range_2         is_in_x_range =
                                        traits->is_in_x_range_2_object();
    typename Traits_adaptor_2::Construct_min_vertex_2  min_vertex =
                                        traits->construct_min_vertex_2_object();
    typename Traits_adaptor_2::Construct_max_vertex_2  max_vertex =
                                        traits->construct_max_vertex_2_object();
    typename Traits_adaptor_2::Compare_y_at_x_2        compare_y_at_x =
                                        traits->compare_y_at_x_2_object();

    CGAL::Object               obj;
    Point_2                    ip;
    bool                       left_equals_curr_endpoint;

    Ccb_halfedge_circulator he_curr = he_first;
    do
    {
      // If we have already found an intersection with the twin halfedge,
      // we do not have to compute intersections with the current halfedge.
      // This happens if we already discovered the twin's face
      // todo: this doesn't work for antennas!
      if (discovered_faces.is_defined(he_curr->twin()->face()) ||
          inter_map.find(&(he_curr->curve())) != inter_map.end())
      {
        ++he_curr;
        continue;
      }

      left_equals_curr_endpoint = false;
      if (on_boundary)
      {
        // Check if the left endpoint of the inserted curve (which is located
        // on the boundary of our face) equals one of the endpoints of the
        // current halfedge. If it equals the right endpoint of the current
        // halfedge, we can skip this edge, as there is no true overlap in
        // the x-range. Otherwise, we keep track of the fact that left_v is
        // the left end-vertex of the current halfedge.
        if (he_curr->target() == left_v)
        {
          left_equals_curr_endpoint = true;

          if (he_curr->direction() == SMALLER)
          {
            ++he_curr;
            continue;
          }
        }
        else if (he_curr->source() == left_v)
        {
          left_equals_curr_endpoint = true;

          if (he_curr->direction() == LARGER)
          {
            ++he_curr;
            continue;
          }
        }
      }

      // Check whether the two curves overlap in their x-range (in order
      // to avoid unnecessary intersection computations).
      if (! left_equals_curr_endpoint &&
          (compare_xy (max_vertex (he_curr->curve()), left_pt) != LARGER ||
           ! is_in_x_range (cv, he_curr->curve())))
      {
        // In case there is no overlap, the two x-monotone curves obviously
        // do not intersect.
        ++he_curr;
        continue;
      }

      // The intersection of the halfedge with the curve have not been
      // computed yet, so we have to compute them now.
      // Note that the first curve we intersect is
      // always the subcurve associated with the given halfegde and the second
      // curve is the one we insert. Even though the order seems unimportant, we
      // exploit this fact in some of the traits classes in order to optimize
      // computations.
      Intersect_list           inter_list;

      traits->intersect_2_object() (he_curr->curve(), cv,
                                    std::back_inserter(inter_list));

      // if there is intersection with the halfedge endpoint, we remove it
      if (! inter_list.empty() && left_equals_curr_endpoint)
        inter_list.pop_front();

      // Insert the list of valid intersections into the map.
      inter_map[&(he_curr->curve())] = inter_list;

      if (! inter_list.empty())
      {
        // insert the intersection to the sorted est
        intersect_he_sorted_set.insert(he_curr);
      }
      
      // Move to the next edge along the ccb
      ++he_curr;
    } while (he_curr != he_first);


  }

};

CGAL_END_NAMESPACE

#endif
