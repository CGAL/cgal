// Copyright (c) 2000  Tel-Aviv University (Israel).
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
// Author(s)     : Eyal flato <flato@math.tau.ac.il>

#ifndef CGAL_PLANAR_MAP_WITH_INTERSECTIONS_H
#define CGAL_PLANAR_MAP_WITH_INTERSECTIONS_H

#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections_misc.h>
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif
#include <CGAL/Planar_map_2/Pm_change_notification.h>

#include <CGAL/Sweep_line_2/Pmwx_aggregate_insert.h>
#include <CGAL/Sweep_line_2/Pmwx_aggregate_insert_old.h>

#ifndef CGAL_PM_COUNT_OPERATIONS_TIMES
#define CGAL_PM_START_OP(x) 
#define CGAL_PM_END_OP(x)  
#define CGAL_DEFINE_COUNT_OP_TIMES_OBJECT
#endif

#if defined(CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG)
#define CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG(msg) std::cout << msg;
#define CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_CODE(code) code
#else
#define CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG(msg)
#define CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_CODE(code)
#endif

CGAL_BEGIN_NAMESPACE

template<class Planar_map_>
class Planar_map_with_intersections_2 : public Planar_map_
{
public:
  typedef Planar_map_                                       Planar_map;
  typedef Planar_map_with_intersections_2<Planar_map>       Self;
  typedef typename Planar_map::Traits                       Traits;
  typedef typename Planar_map::Traits_wrap                  Pm_traits_wrap;
  typedef Planar_map_with_intersections_traits_wrap<Traits> Pmwx_traits_wrap;
  typedef typename Planar_map::Halfedge_handle              Halfedge_handle;
  typedef typename Planar_map::Vertex_handle                Vertex_handle;
  typedef typename Planar_map::Face_handle                  Face_handle;
  typedef typename Traits::X_monotone_curve_2               X_monotone_curve_2;
  typedef typename Traits::Curve_2                          Curve_2;
  typedef typename Traits::Point_2                          Point_2;

  typedef typename Planar_map::Change_notification
    Change_notification;

  typedef typename Planar_map::Halfedge_iterator            Halfedge_iterator;
    
  // Obsolete, for backward compatability
  typedef Point_2                               Point;
  typedef X_monotone_curve_2                    X_curve;
  typedef Curve_2                               Curve;
  typedef Change_notification                   Pmwx_change_notification;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Planar_map::traits;
  using Planar_map::halfedges_end;
  using Planar_map::halfedges_begin;
  using Planar_map::unbounded_face;
  using Planar_map::pl;
#endif

  CGAL_DEFINE_COUNT_OP_TIMES_OBJECT

  // Constructors
  // ------------

  Planar_map_with_intersections_2() :
#ifndef CGAL_NO_PM_DEFAULT_POINT_LOCATION
    Planar_map(new Pmwx_traits_wrap,
               new Pm_walk_along_line_point_location<Planar_map>, NULL),
#else
    Planar_map(new Pmwx_traits_wrap, NULL, NULL),
#endif
    pmwx_use_delete_traits(true), 
    pmwx_use_delete_pl(true)
  { 
    pmwx_traits = (Pmwx_traits_wrap*)traits;
  }    
        
  Planar_map_with_intersections_2(Pm_point_location_base<Planar_map> *pl_ptr) :
    Planar_map(new Pmwx_traits_wrap, pl_ptr, NULL),
    pmwx_use_delete_traits(true), 
    pmwx_use_delete_pl(false)
  {
    pmwx_traits = (Pmwx_traits_wrap*)traits;
  }

  Planar_map_with_intersections_2(const Traits& tr_, 
                                  Pm_point_location_base<Planar_map> *pl_ptr) :
    Planar_map(new Pmwx_traits_wrap(tr_), pl_ptr, NULL),
    pmwx_use_delete_traits(true), 
    pmwx_use_delete_pl(false)
  {
    pmwx_traits = (Pmwx_traits_wrap*)traits;
  }

  Planar_map_with_intersections_2(Pmwx_traits_wrap *tr_ptr,
                                  Pm_point_location_base<Planar_map> *pl_ptr) :
    Planar_map(tr_ptr, pl_ptr, NULL),
    pmwx_use_delete_traits(false), 
    pmwx_use_delete_pl(false)
  {
    pmwx_traits = (Pmwx_traits_wrap*)traits;
  }

  /*! Copy Constructor */
  Planar_map_with_intersections_2(const Self & rhs);

  /*! Planar_map converter */
  Planar_map_with_intersections_2(const Planar_map & pm);

  /*! Destructor */
  ~Planar_map_with_intersections_2();
        
  // Operations
  // ----------

  /*! Obtain the supporting curve of the given halfedge
   * \param he the halfedge
   * \param en the notification object
   * \return the supporting curve of he if exists
   */
  const X_monotone_curve_2 & curve(Halfedge_handle he,
                                   Change_notification * en)
  { 
    if (en == NULL)
      return he->curve();
    if (!en->have_support_curve())
      return he->curve();
    return en->edge_support_curve(he); 
  }

  // finds the intersection of <cv> directed <direction_right> with the
  // halfedge <he>. The returned intersections is based on the intersection
  // of the supporting curves (if they exist).
  bool
  directed_nearest_intersection_with_halfedge(const X_monotone_curve_2 & cv,
                                              Halfedge_handle he,
                                              const Point_2 & ref_point,
                                              bool direction_right,
                                              Point_2 & xp1,
                                              Point_2 & xp2,
                                              X_monotone_curve_2 & overlap_cv,
                                              Change_notification * en)
  {
    Object res = (direction_right) ?
      pmwx_traits->nearest_intersection_to_right(cv, curve(he, en), ref_point)
      :
      pmwx_traits->nearest_intersection_to_left(cv, curve(he, en), ref_point);

    if (res.is_empty())
      return (false);
                
    /*
     * check for an intersection on the real curves. assume there is none.
     *
     * since we are checking on the parent, we should make sure that the 
     * intersection point is on the halfedge_cv and not only on the parent.
     * do not worry: we will get the same intersection point for the correct
     * halfedge_cv as well, and therefore we can throw it away if it's 
     * not on halfedge_cv
     * no need to check for cv because the checked side of it is not 
     * in the arrangement yet so there is no possibility for an 
     * intersection point not on cv.
     */

    const X_monotone_curve_2 & he_cv = he->curve();
    
    if (CGAL::assign(xp1, res)) {
      // The intersection is a point:
      xp2 = xp1;        //! \todo is this really needed?
      if (traits->point_in_x_range(he_cv, xp1) &&
 	  traits->curve_compare_y_at_x(xp1, he_cv) == EQUAL) {
        return true;
      }
      return false;
    }

    if (CGAL::assign(overlap_cv, res)) {
      // There is an overlap, the intersection is a subcurve:

      xp1 = pmwx_traits->curve_source(overlap_cv);
      xp2 = pmwx_traits->curve_target(overlap_cv);

      if (!is_left_low(xp1, xp2))
        std::swap(xp1, xp2);

      Point_2 left_point =
        traits->point_leftlow_most(pmwx_traits->curve_source(he_cv),
                                   pmwx_traits->curve_target(he_cv));
      Point_2 right_point =
        traits->point_righttop_most(pmwx_traits->curve_source(he_cv),
                                    pmwx_traits->curve_target(he_cv));
  
      if (is_left_low(xp1, left_point))
        xp1 = left_point;

      if (is_left_low(right_point, xp2))
        xp2 = right_point;

      if (is_left_low(xp1, xp2))
        return (true);

      return (false);
    }

    return (false);
  }


  // input: cv.source is on vh
  // is_overlap is true if cv overlaps prev_halfedge around vh.
  void find_face_of_curve_around_vertex(const X_monotone_curve_2 & cv, 
                                        const X_monotone_curve_2 & orig_cv, 
                                        const Vertex_handle & vh,
                                        Halfedge_handle & prev_halfedge,
                                        Point_2 & overlap_end_pt,
                                        bool & is_overlap,
                                        X_monotone_curve_2 & overlap_cv,
                                        Change_notification * en)
  {
    CGAL_PM_START_OP(1);
     
    typename Planar_map::Halfedge_around_vertex_circulator next, prev, start, 
      last_next_checked;
    X_monotone_curve_2 xcv;
    Point_2 xp1, xp2;
    const Point_2 & source_point = pmwx_traits->curve_source(cv);
    const Point_2 & target_point = pmwx_traits->curve_target(cv);
    bool direction_right = is_left_low(source_point, target_point);
    bool b;
                
    last_next_checked = halfedges_end();        
    start = vh->incident_halfedges();
    prev = start;
    next = prev;
    ++next;
    do {
      if ((last_next_checked != prev) && 
          (traits->curves_overlap(cv, prev->curve())))
      {
        // cv overlapps prev->curve()
        b = directed_nearest_intersection_with_halfedge(orig_cv, prev,
                                                        source_point,
                                                        direction_right,
                                                        xp1, xp2, overlap_cv,
                                                        en);
        // Verify that there is indeed an intersection
        CGAL_assertion(b);
        // Verify that there is indeed an overlap
        CGAL_assertion(!point_equal(xp1, xp2));

        // the overlapping part might not start from the source 
        // vertex (in polyline for example), so if this is the case, 
        // we ignore the overlap. 
        if ((point_equal(xp1, source_point)) ||
            (point_equal(xp2, source_point)))
        {
          if (point_equal(vh->point(), xp1)) {
            overlap_end_pt = xp2;
          }
          else {
            CGAL_assertion(point_equal(vh->point(), xp2));
            overlap_end_pt = xp1;
          }
          
          prev_halfedge = prev;
          is_overlap = true;
          CGAL_PM_END_OP(1);
          return;
        }
      }

      // Same check as above is done now to next->curve() to prevent
      // unexpected behavior of is_between_cw when one of the curves
      // overlaps at vh.
      // last_next_checked is set to next so if the test fails
      // (namely, there is no overlap at vh) then it will not 
      // be executed again at the next loop when prev will be this 
      // next. this is for efficiency only.
      if (traits->curves_overlap(cv, next->curve())) {
        last_next_checked = next;
        // cv overlapps next->curve()
        b = directed_nearest_intersection_with_halfedge(orig_cv, next,
                                                        source_point,
                                                        direction_right,
                                                        xp1, xp2, overlap_cv,
                                                        en);
        // Verify that there is indeed an intersection
        CGAL_assertion(b);
        // Verify that there is indeed an overlap
        CGAL_assertion(!point_equal(xp1, xp2));

        // the overlapping part might not start from the source 
        // vertex (in polyline for example), so if this is the case, 
        // we ignore the overlap. 
        if ((point_equal(xp1, source_point)) ||
            (point_equal(xp2, source_point)))
        {
          if (point_equal(vh->point(), xp1)) {
            overlap_end_pt = xp2;
          }
          else {
            CGAL_assertion(point_equal(vh->point(), xp2));
            overlap_end_pt = xp1;
          }
          
          prev_halfedge = next;
          is_overlap = true;
          CGAL_PM_END_OP(1);
          return;
        }
      }
  
      /////////***** Eyal end
      // The following remarked test is not fully correct
      // since there can be an overlap that does not start at vh
      // but elsewhere on cv (e.g., polylines).
      // this condition is redundant if curve_is_between_cw does
      // not return true when overlap occurs, but it is needed
      // since it is not a defined behaviour (in specs).
      //// if (!traits->curves_overlap(cv, next->curve())) ----
      bool b1,b2;
      if (next != prev) 
      {
        if ((pmwx_traits->curve_is_between_cw(cv, prev->curve(), next->curve(),
                                              vh->point(), b1, b2)))
        {
          prev_halfedge = prev;
          is_overlap = false;
          CGAL_PM_END_OP(1);
          return;
        }
      }
      ++next;
      ++prev;
    } while (prev != start);
                
    // assuming no overlapping occurs, the only situation in 
    // which we reached this line is when there is only one 
    // edge from vh
    CGAL_assertion(next == prev);
    prev_halfedge = prev;
    is_overlap = false;
    CGAL_PM_END_OP(1);
    return;
  }
        
  /*! Find the first intersection in the direction cv.source --> cv.target
   * in <face>.
   * The returned <intersection> is the intersection curve of 
   * <cv> and <face>'s boundary.
   * returned <halfedge> is the halfedge on which the intersection occurs,
   * in case of intersection-point halfedge->source will contain this point
   * \return false if no intersection is found and true otherwise.
   */
  bool find_first_intersection_in_face(const Face_handle & face,
                                       const X_monotone_curve_2 & cv,
                                       const X_monotone_curve_2 & orig_cv, 
                                       Halfedge_handle & halfedge,
                                       Point_2 & best_xpnt1,
                                       Point_2 & best_xpnt2,
                                       Change_notification * en)
  {
    CGAL_PM_START_OP(2);
    Halfedge_handle best_halfedge_x;
    typename Planar_map::Ccb_halfedge_circulator che, che_beg;
    Point_2 xpnt, xp1, xp2;
    bool b, better_xcv;
    bool intersection_exist = false;
    const Point_2 & source_point = pmwx_traits->curve_source(cv);
    const Point_2 & target_point = pmwx_traits->curve_target(cv);
    bool direction_right = is_left_low(source_point, target_point);
    
    Point_2 best_point = pmwx_traits->curve_target(cv);
    
    if (!face->is_unbounded()) {
      che = face->outer_ccb();
      che_beg = che;
      do {
        if (intersection_exist) {
          // optimization: comparing the
          // endpoints of the halfedge and the best intersection 
          // point found so far, Improve performance.
          if (direction_right) {
            if (!is_left(leftmost(che->source()->point(),
                                  che->target()->point()),
                         best_point))
            {
              ++che;
              continue;
            }
          } else {
            if (!is_right(rightmost(che->source()->point(),
                                    che->target()->point()),
                          best_point))
            {
              ++che;
              continue;
            }
          }
        }
          
        // optimization: checking the x-range, Highly improve performance.
        if (!in_x_range(orig_cv,che->curve())) {
          ++che;
          continue;
        }
        
        X_monotone_curve_2 overlap_cv;
        b = directed_nearest_intersection_with_halfedge(orig_cv, che,
                                                        source_point,
                                                        direction_right,
                                                        xp1, xp2, overlap_cv,
                                                        en);
        if (b) {
          // direct xp1-xp2 to be like cv
          if (direction_right) {
            if (!is_left_low(xp1, xp2))
              std::swap(xp1, xp2);
          } else {
            if (!is_right_top(xp1, xp2))
              std::swap(xp1, xp2);
          }
          
          xpnt = xp1;
          if (direction_right)
            better_xcv = is_left_low(xpnt, best_point);
          else
            better_xcv = is_right_top(xpnt, best_point);
          if (better_xcv || (!intersection_exist)) {
            // xcv is better
            best_halfedge_x = che;
            best_point = xpnt;
            best_xpnt1 = xp1;
            best_xpnt2 = xp2;
            intersection_exist = true;
          }
        }
        ++che;
      } while (che != che_beg);
    }

    //if (intersection_exist)
    //  {
    //    halfedge = best_halfedge_x;
    
    //    CGAL_PM_END_OP(2);
    //      return intersection_exist;
    //  }

    typename Planar_map::Holes_iterator holes;
    for (holes = face->holes_begin(); holes != face->holes_end(); holes++) {
      che = *holes;
      che_beg = che;
      do {
        if (intersection_exist) {
           if (direction_right) {
             if (!(is_left_low(leftmost(che->source()->point(),
                                        che->target()->point()),
                               best_point)))
             {
               ++che;
               continue;
             }
           } else {
             if (!(is_right_top(rightmost(che->source()->point(),
                                          che->target()->point()),
                                best_point)))
             {
               ++che;
               continue;
             }
           }
        }
        
        if (!in_x_range(orig_cv,che->curve())) {
          ++che;
          continue;
        }

        X_monotone_curve_2 overlap_cv;
        b = directed_nearest_intersection_with_halfedge(orig_cv, che,
                                                        source_point,
                                                        direction_right,
                                                        xp1, xp2, overlap_cv,
                                                        en);
        if (b) {
          // direct xp1-xp2 to be like cv
          if (direction_right) {
            if (!is_left_low(xp1, xp2))
              std::swap(xp1, xp2);
          } else {
            if (!is_right_top(xp1, xp2))
              std::swap(xp1, xp2);
          }
          
          xpnt = xp1;
          if (direction_right)
            better_xcv = is_left_low(xpnt, best_point);
          else
            better_xcv = is_right_top(xpnt, best_point);
          if (better_xcv || (!intersection_exist)) {
            // xcv is better
            best_halfedge_x = che;
            best_point = xpnt;
            best_xpnt1 = xp1;
            best_xpnt2 = xp2;
            intersection_exist = true;
          }
        }
        ++che;
      } while (che != che_beg);
    }
    
    
    halfedge = best_halfedge_x;
    
    CGAL_PM_END_OP(2);
    return intersection_exist;
  }

  /*!
   */
  void get_vertex_of_point_on_halfedge
  (
    const Point_2 &point,                     
    Halfedge_handle halfedge,
    Vertex_handle &vertex_of_point, 
    // in case of split it is easy to compute:
    Halfedge_handle &vertex_of_point_prev_halfedge,
    // true if vertex_of_point_prev_halfedge is set :
    bool &vertex_of_point_prev_halfedge_set,
    Change_notification *en)
  {
    CGAL_PM_START_OP(3);
    if (point_equal(point, halfedge->source()->point())) {
      vertex_of_point = halfedge->source();
      vertex_of_point_prev_halfedge_set = false;
    } else if (point_equal(point, halfedge->target()->point())) {
      vertex_of_point = halfedge->target();
      vertex_of_point_prev_halfedge_set = false;
    } else {
      // intersection in the interior - split
      X_monotone_curve_2 split1, split2;
      Halfedge_handle he_split;
      if (point_equal(halfedge->source()->point(), 
                      pmwx_traits->curve_source(halfedge->curve())))
      {
        pmwx_traits->directed_curve_split(halfedge->curve(), 
                                          halfedge->source()->point(), 
                                          point, split1, split2);
        he_split = Planar_map::split_edge(halfedge, split1, split2);
        if (en != NULL) 
          en->split_edge(he_split, he_split->next_halfedge(), split1, split2);
                            
        vertex_of_point = he_split->target();
        vertex_of_point_prev_halfedge = he_split->next_halfedge()->twin();
        vertex_of_point_prev_halfedge_set = true;
      } else {
        Halfedge_handle twin_halfedge = halfedge->twin();
        pmwx_traits->directed_curve_split(twin_halfedge->curve(),
                                          twin_halfedge->source()->point(), 
                                          point, split1, split2);
        he_split = Planar_map::split_edge(twin_halfedge, split1, split2);
        if (en != NULL)
          en->split_edge(he_split, he_split->next_halfedge(), split1, split2);
                            
        // change he_split to be on the current face becase we spliited 
        // the twin_halfedge and not halfedge - this is to be 
        // consistent with 
        // the arrangement code that is based on the fact the the 
        // splitted halfedge
        // is the one that is directed like the curve
        he_split = he_split->next_halfedge()->twin();
        vertex_of_point = he_split->target();
        vertex_of_point_prev_halfedge = he_split->next_halfedge()->twin();
        vertex_of_point_prev_halfedge_set = true;
      }
    }
                
    // in order to enable polyline overlap
    // (a little ugly solution. should be made nicer) : 
    vertex_of_point_prev_halfedge_set = false; 
    CGAL_PM_END_OP(3);
  }

 
  /*! Insert the first part of cv into the face source_face 
   * returning: 
   *   1. inserted edge
   *   2. remaining curve (the part that was not inserted)
   *   3. remaining curve source vertex
   */
  void insert_intersecting_xcurve_in_face_interior
  (const X_monotone_curve_2 & cv,                     // inserted curve
   const X_monotone_curve_2 & orig_cv, 
   Face_handle source_face,
   Halfedge_handle & inserted_halfedge,
   bool & remaining_curve_trivial,
   X_monotone_curve_2 & remaining_curve,
   Vertex_handle & remaining_curve_source_vertex, 
   // in case of split it is easy to compute :
   Halfedge_handle & remaining_curve_prev_halfedge,
   // true if remaining_curve_face is set :
   bool &remaining_curve_prev_halfedge_set,  
   Change_notification * en)
  {
    CGAL_PM_START_OP(4);
    remaining_curve_trivial = false;
    //std::cout << "iisifi " 
    // << "insert_intersecting_xcurve_in_face_interior: " << cv << std::endl;
    Halfedge_handle intersection_he;
    Point_2 xpnt1, xpnt2;
                
    if (find_first_intersection_in_face(source_face, cv, orig_cv,
                                        intersection_he, xpnt1, xpnt2, en))
    {
      Point_2 insert_point;
      X_monotone_curve_2 cv_first_part;
      bool direction_right = is_left_low(pmwx_traits->curve_source(cv),
                                         pmwx_traits->curve_target(cv));

      insert_point = (direction_right) ?
        pmwx_traits->point_leftlow_most(xpnt1, xpnt2) :
        pmwx_traits->point_righttop_most(xpnt1, xpnt2);
                        
      CGAL_assertion(!point_equal(pmwx_traits->curve_source(cv),
                                  insert_point));

      remaining_curve_trivial = point_equal(pmwx_traits->curve_target(cv),
                                            insert_point);
      if (remaining_curve_trivial)
        cv_first_part = cv;
      else
        pmwx_traits->directed_curve_split(cv, pmwx_traits->curve_source(cv),
                                          insert_point, cv_first_part,
                                          remaining_curve);

                        
      CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_CODE(
        std::cout << "cv " << cv << std::endl;
        std::cout << "xpnt1 " << xpnt1 << std::endl;
        std::cout << "xpnt2 " << xpnt2 << std::endl;
        std::cout << "insert_point " << insert_point << std::endl;
        std::cout << "cv_first_part " << cv_first_part << std::endl;
        std::cout << "remaining_curve " << remaining_curve << std::endl;
        std::cout << "intersection_he->curve() " 
                  << intersection_he->curve() << std::endl;
        );
        
      get_vertex_of_point_on_halfedge(insert_point,
                                      intersection_he,
                                      remaining_curve_source_vertex,
                                      remaining_curve_prev_halfedge,
                                      remaining_curve_prev_halfedge_set, 
                                      en);

      inserted_halfedge =
        Planar_map::insert_from_vertex(cv_first_part, 
                                       remaining_curve_source_vertex);
      if (en != NULL) en->add_edge(cv_first_part, inserted_halfedge,
                                   false, false);
      CGAL_PM_END_OP(4);
      return;
    }

    inserted_halfedge = Planar_map::insert_in_face_interior(cv, source_face);
    if (en != NULL) {
      en->add_edge(cv, inserted_halfedge, true, false);
      en->add_hole(source_face, inserted_halfedge);
    }
    //std::cout << "iisifi inserted_halfedge " 
    // << inserted_halfedge->source()->point() <<
    //  " ----> " << inserted_halfedge->target()->point() << std::endl;

    if (point_equal(inserted_halfedge->source()->point(),
                    pmwx_traits->curve_source(cv)))
      remaining_curve_source_vertex = inserted_halfedge->target();
    else
      remaining_curve_source_vertex = inserted_halfedge->source();
    remaining_curve_trivial = true;
    remaining_curve_prev_halfedge_set = false;
    CGAL_PM_END_OP(4);
  }

  // source_prev_halfedge->face() is the face we are working on
  // source_prev_halfedge->target() is the vertex on whice cv's source is
  // assuming cv does not overlap an edge from source_prev_halfedge->target()
  // insert the first part of cv into face where cv's source is on the 
  // boundary of face
  // returning: 
  //   1. inserted edge
  //   2. remaining curve (the part that was not inserted)
  //   3. remaining curve source vertex 
  void insert_intersecting_xcurve_from_boundary_of_face
  (const X_monotone_curve_2 & cv,                     // inserted curve
   const X_monotone_curve_2 & orig_cv, 
   Halfedge_handle source_prev_halfedge,
   Halfedge_handle & inserted_halfedge,
   bool & remaining_curve_trivial,
   X_monotone_curve_2 & remaining_curve,
   Vertex_handle & remaining_curve_source_vertex, 
   // in case of split it is easy to compute :
   Halfedge_handle & remaining_curve_prev_halfedge, 
   // true if remaining_curve_face is set :
   bool &remaining_curve_prev_halfedge_set,  
   Change_notification * en)
  {
    CGAL_PM_START_OP(5);
    remaining_curve_trivial = false;
    //std::cout << "iifbof " 
    // << "insert_intersecting_xcurve_from_boundary_of_face: " 
    // << cv << std::endl;
    //std::cout << "iifbof " 
    // << "         source_prev_halfedge: " 
    // << source_prev_halfedge->source()->point() <<
    //   " ----> " << source_prev_halfedge->target()->point() << std::endl;
    Halfedge_handle intersection_he;
    Point_2 xpnt1, xpnt2;
    Face_handle orig_face = source_prev_halfedge->face();
                
    if (find_first_intersection_in_face(source_prev_halfedge->face(),
                                        cv, orig_cv, intersection_he,
                                        xpnt1, xpnt2, en))
    {
      Point_2 insert_point;
      X_monotone_curve_2 cv_first_part;

      // keep the source vertex. We can't rely on 
      // source_prev_halfedge->target() because it might be changed if 
      // source_prev_halfedge is intersected by the new curve
      Vertex_handle source_vertex(source_prev_halfedge->target());

      bool direction_right = is_left_low(pmwx_traits->curve_source(cv),
                                         pmwx_traits->curve_target(cv));
      insert_point = (direction_right) ?
        pmwx_traits->point_leftlow_most(xpnt1, xpnt2) :
        pmwx_traits->point_righttop_most(xpnt1, xpnt2);

      CGAL_assertion(!point_equal(pmwx_traits->curve_source(cv),
                                  insert_point));

      remaining_curve_trivial = point_equal(pmwx_traits->curve_target(cv),
                                            insert_point);
      if (remaining_curve_trivial)
        cv_first_part = cv;
      else
        pmwx_traits->directed_curve_split(cv, pmwx_traits->curve_source(cv),
                                          insert_point, cv_first_part,
                                          remaining_curve);
                        
      get_vertex_of_point_on_halfedge(insert_point, intersection_he,
                                      remaining_curve_source_vertex,
                                      remaining_curve_prev_halfedge,
                                      remaining_curve_prev_halfedge_set, en);
                        
      //std::cout << "iifbob " << "insert_at_vertices: " << 
      //  cv_first_part << 
      //  "   vx1: "<< source_prev_halfedge->target()->point() <<
      //  "   prev_src: "<< source_prev_halfedge->source()->point() <<
      //  "   vx2: "<< remaining_curve_source_vertex->point() << std::endl;
      CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_CODE(
        std::cout << "cv " << cv << std::endl;
        std::cout << "xpnt1 " << xpnt1 << std::endl;
        std::cout << "xpnt2 " << xpnt2 << std::endl;
        std::cout << "insert_point " << insert_point << std::endl;
        std::cout << "cv_first_part " << cv_first_part << std::endl;
        std::cout << "remaining_curve " << remaining_curve << std::endl;
        std::cout << "intersection_he->curve() " 
                  << intersection_he->curve() << std::endl;
        );
      
      inserted_halfedge = 
        Planar_map::insert_at_vertices(cv_first_part, source_vertex, 
                                       remaining_curve_source_vertex);
      if (en != NULL) 
      {
        en->add_edge(cv_first_part, inserted_halfedge, true, false);
        if (inserted_halfedge->face() == orig_face)
          en->split_face(inserted_halfedge->face(), 
                         inserted_halfedge->twin()->face());
        else
          en->split_face(inserted_halfedge->twin()->face(), 
                         inserted_halfedge->face());
      }
      CGAL_PM_END_OP(5);
      return;
    }

    inserted_halfedge = 
      Planar_map::insert_from_vertex(cv, source_prev_halfedge->target()); 
    if (en != NULL) 
      en->add_edge(cv, inserted_halfedge,true, false);
    if (point_equal(inserted_halfedge->source()->point(),
                    pmwx_traits->curve_source(cv)))
      remaining_curve_source_vertex = inserted_halfedge->target();
    else
      remaining_curve_source_vertex = inserted_halfedge->source();
    remaining_curve_trivial = true;
    remaining_curve_prev_halfedge_set = false;
    CGAL_PM_END_OP(5);
  }

  /*!
   */
  Halfedge_handle
  insert_intersecting_xcurve(const X_monotone_curve_2 & cv_, // inserted curve
                             Vertex_handle & source_vertex,
                             // to be set by the function :  
                             Vertex_handle & target_vertex, 
                             bool source_valid,
                             Change_notification * en = NULL)
  {
    CGAL_PM_START_OP(6);
    //if a vertex on which an endpoint of cv_ is known then set cv to 
    //have this endpoint as it's source
    // at the end source_vertex and target_vertex will contain the
    // end-vertices of cv
    Vertex_handle curr_vertex = source_vertex;
    X_monotone_curve_2 cv = cv_; 
    X_monotone_curve_2 orig_cv = cv;
    X_monotone_curve_2 split1, split2;
    Point_2 overlap_end_pt;
    Halfedge_handle inserted_he, prev_halfedge, last_edge;
    bool is_overlap;
    bool next_face_valid = false;
    bool remaining_curve_trivial = false; 

    const Point_2 & source_point = pmwx_traits->curve_source(cv);
    
    if (!source_valid) {
      /* If the source is invalid, we need to locate the source if exists,
       * or create it by splitting an existing halfedge or insertion of the
       * new curve
       */
      CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("!");

      typename Planar_map::Locate_type lt;
      Halfedge_handle he;
      //-- locate pmwx_traits->curve_source(cv)
      CGAL_PM_START_OP(7);
      he = locate(source_point, lt);
      CGAL_PM_END_OP(7);
                        
      if (lt == Planar_map::VERTEX) {
        // The source point exists! We simply set the source appropriately
        source_vertex = he->target();
        curr_vertex = source_vertex;
      } else if (lt == Planar_map::EDGE) {
        /* The source point does not exist, but it lies on an existng
         * halfedge. In this case we split the halfedge and set the source
         * vertex appropriately
         */
        Halfedge_handle he_split;
        X_monotone_curve_2 split1, split2;
        if (point_equal(he->source()->point(),
                        pmwx_traits->curve_source(he->curve())))
        {
          pmwx_traits->directed_curve_split(he->curve(),
                                            he->source()->point(),
                                            source_point, split1, split2);
          he_split = Planar_map::split_edge(he, split1, split2);
        } else {
          Halfedge_handle twin_he = he->twin();
          pmwx_traits->directed_curve_split(twin_he->curve(),
                                            twin_he->source()->point(),
                                            source_point, split1, split2);
          he_split = Planar_map::split_edge(twin_he, split1, split2);
        }
        if (en != NULL) 
          en->split_edge(he_split, he_split->next_halfedge(), split1, split2);
        source_vertex = he_split->target();
        curr_vertex = source_vertex;
      } else {
        // The source point lies inside the unbounded or a general face.
        Face_handle face = (lt == Planar_map::UNBOUNDED_FACE) ?
          unbounded_face() : he->face();
        // insert_intersecting_xcurve_in_face_interior(curr_vertex)
        X_monotone_curve_2 remaining_curve; 
        insert_intersecting_xcurve_in_face_interior(cv, orig_cv, face,
                                                    inserted_he,
                                                    remaining_curve_trivial,
                                                    remaining_curve,    
                                                    curr_vertex,
                                                    prev_halfedge,
                                                    next_face_valid, en);
        if (!remaining_curve_trivial)
          cv = remaining_curve;
        last_edge = inserted_he;
        target_vertex = curr_vertex; // temporary - can be changed later
        source_vertex = (inserted_he->source() == curr_vertex) ?
          inserted_he->target() : inserted_he->source();
        //inserted_edges.push_back(inserted_he);
      }
                        
      // by now: curr_vertex and source_vertex are set
    }
                
    while (!remaining_curve_trivial) {
      CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("@");

      Halfedge_handle he_split;
      X_monotone_curve_2 split1, split2;

      if (!next_face_valid) { 
        CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("#");

        //-- locate cv around curr_vertex
        //std::cout << "iis " << " ---- curr_vertex: " 
        //<< curr_vertex->point() << std::endl;
        X_monotone_curve_2 overlap_cv;
        find_face_of_curve_around_vertex(cv, orig_cv, curr_vertex,
                                         prev_halfedge, overlap_end_pt,
                                         is_overlap, overlap_cv, en);
        //std::cout << "iis " << " ---- prev_halfedge: " 
        //<< prev_halfedge->source()->point() <<
        //   "  ---->  " << prev_halfedge->target()->point() << std::endl;
        if (is_overlap) {
          // if overlaps an edge from curr_vertex 
          CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("$");

          // rem: prev_halfedge->target == curr_vertex
          if (point_equal(prev_halfedge->source()->point(), overlap_end_pt)) {
            // whole edge is overlapped, proceed to its other end
            CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("%");

            if (en != NULL) 
              en->add_edge(overlap_cv, prev_halfedge, false, true);
            last_edge = prev_halfedge;
            // update cv
                        
            CGAL_assertion(point_equal(pmwx_traits->curve_source(cv),
                                       curr_vertex->point()));
            remaining_curve_trivial =
              point_equal(pmwx_traits->curve_target(cv),
                          prev_halfedge->source()->point());
            if (!remaining_curve_trivial) {
              pmwx_traits->directed_curve_split(cv, curr_vertex->point(), 
                 prev_halfedge->source()->point(), split1, split2);
              cv = split2;
            }
            curr_vertex = prev_halfedge->source();
          } else { 
            // otherwise - split prev_halfedge and let curr_vertex 
            //be the splitting point
            CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("^");

            if (point_equal(prev_halfedge->twin()->source()->point(), 
                            pmwx_traits->curve_source(prev_halfedge->twin()
                                                      ->curve())))
            {
              pmwx_traits->directed_curve_split(prev_halfedge->curve(),
                                                curr_vertex->point(), 
                 overlap_end_pt, split1, split2);
                                
              // split prev_halfedge->twin() so that the splitted edge 
              // source will be the same as split1's source
              he_split = 
                Planar_map::split_edge(prev_halfedge->twin(), split1, split2);
              if (en != NULL) {
                en->split_edge(he_split, he_split->next_halfedge(), 
                               split1, split2);
                en->add_edge(overlap_cv, he_split, true, true);
              }
              last_edge = he_split;
                            
              // update cv
              remaining_curve_trivial =
                point_equal(pmwx_traits->curve_target(cv),
                            he_split->target()->point());
              if (!remaining_curve_trivial) {
                pmwx_traits->directed_curve_split(cv,
                                                  pmwx_traits->
                                                  curve_source(cv),
                                                  he_split->target()->point(),
                                                  split1, split2);
                cv = split2;
              }
              curr_vertex = he_split->target();
            } else {
              pmwx_traits->
                directed_curve_split(prev_halfedge->twin()->curve(), 
                                     pmwx_traits->
                                     curve_source(prev_halfedge->twin()->
                                                  curve()), overlap_end_pt,
                                     split1, split2);
                        
              // split prev_halfedge->twin() so that the splitted
              //  edge 
              // source will be the same as split1's source
              he_split = Planar_map::split_edge(prev_halfedge, split1, split2);
              if (en != NULL) {
                en->split_edge(he_split, he_split->next_halfedge(), 
                               split1, split2);
                en->add_edge(overlap_cv, he_split->next_halfedge()->twin(),
                             true, true);
              }
              last_edge = he_split->next_halfedge()->twin();

              // update cv
              remaining_curve_trivial =
                point_equal(pmwx_traits->curve_target(cv),
                            he_split->
                            next_halfedge()->twin()->target()->point());
              if (!remaining_curve_trivial) {
                pmwx_traits->directed_curve_split(cv,
                                                  pmwx_traits->
                                                  curve_source(cv),
                                                  he_split->next_halfedge()->
                                                  twin()->target()->point(),
                                                  split1, split2);
                cv = split2;
              }
              curr_vertex = he_split->next_halfedge()->twin()->target();

            }
          }
                                        
          target_vertex = curr_vertex; // temporary -can be changed later
          next_face_valid = false; 
          continue;
        }
        // if not overlap then in face
        next_face_valid = true;
      }
                        
      // now - we are sure that cv is in the interior of
      // prev_halfedge->face (does not overlaps an halfedge around 
      // curr_vertex)
                        
      CGAL_PM_WITH_INTERSECTIONS_PRINT_DEBUG_MSG("&");

      X_monotone_curve_2 remaining_curve;
      insert_intersecting_xcurve_from_boundary_of_face(cv, orig_cv,
                                                       prev_halfedge,
                                                       inserted_he,
                                                       remaining_curve_trivial,
                                                       remaining_curve,
                                                       curr_vertex,
                                                       prev_halfedge,
                                                       next_face_valid, en);
      if (!remaining_curve_trivial)
        cv = remaining_curve;
      last_edge = inserted_he;

      target_vertex = curr_vertex; // temporary - can be changed later
    }

    if (!point_equal(last_edge->target()->point(),
                     pmwx_traits->curve_target(cv_)))
      last_edge = last_edge->twin();
    CGAL_PM_END_OP(6);
    return last_edge; 
  }

  /*! Insert a given curve into the map. One of its endpoints is conditionally
   * the mapping of a given vertex
   */
  Halfedge_handle insert_intersecting_curve(const Curve_2 & c,
                                            Vertex_handle & source_vertex,
                                            bool source_valid,
                                            Change_notification * en = NULL)
  {
    typedef typename Traits::X_monotone_curve_2         X_monotone_curve_2;
    typedef std::list<X_monotone_curve_2>               X_monotone_curve_list;
    typedef typename X_monotone_curve_list::const_iterator
      X_monotone_curve_iter;
    
    Vertex_handle src = source_vertex;
    Vertex_handle trg;
    Halfedge_handle last_he;
    X_monotone_curve_list x_list;
    traits->curve_make_x_monotone(c, std::back_inserter(x_list));

    X_monotone_curve_iter it = x_list.begin();
    last_he = insert_intersecting_xcurve(*it, src, trg, source_valid, en); 
    src = trg;
    for (it++; it != x_list.end(); it++) {
      last_he = insert_intersecting_xcurve(*it, src, trg, true, en); 
      src = trg;
    }

    return last_he;
  }

  /*! Insert a given curve that one of its endpoints is the mapping of a given
   * vertex into the map.
   * \param cv the curve to insert
   * \param src an existing vertex in the map and one the endpoint of cv
   * \param en the notification
   * \return the last inserted halfedge. Its target maps to the target point of
   * the last inserted X-monotone curve
   */
  Halfedge_handle insert_from_vertex(const Curve_2 & cv, Vertex_handle src, 
                                     Change_notification * en = NULL)
  {
    CGAL_precondition(!point_equal(pmwx_traits->curve_source(cv),
                                   pmwx_traits->curve_target(cv)));
    return insert_intersecting_curve(cv, src,  true, en);
  }

  // return the last inserted halfedge whose target points to the last 
  // point of the inserted xcurve
  Halfedge_handle insert(const Curve_2 & c, Change_notification * en = NULL)
  {
    // If curve is x-monotone then its source is different from its target.
    // (which is not true for non x-monotone curves, e.g, triangles.)

    Vertex_handle src;
    return insert_intersecting_curve(c, src, false, en);
  }


  // inserts a given curve container into the map using Eti's sweep
  template <class X_monotone_curve_2_iterator>
  Halfedge_iterator insert_old(const X_monotone_curve_2_iterator & begin,
			       const X_monotone_curve_2_iterator & end,
			       Change_notification * en = NULL)
  {
    typedef Pmwx_aggregate_insert_old<X_monotone_curve_2_iterator, Traits, 
      Self ,Change_notification> Pmwx_agg_insert;
    Pmwx_agg_insert p(traits);
    p.insert_curves(begin, end, *this, en);

    return halfedges_begin();
  }

  // inserts a given curve container into the map using Tali's sweep
  template <class X_monotone_curve_2_iterator>
  Halfedge_iterator insert(const X_monotone_curve_2_iterator & begin,
			   const X_monotone_curve_2_iterator & end,
			   Change_notification * en = NULL)
  {
    typedef Pmwx_aggregate_insert<X_monotone_curve_2_iterator, Traits, 
      Self ,Change_notification> Pmwx_agg_insert;
    Pmwx_agg_insert p(traits);
    p.insert_curves(begin, end, *this, en);

    return halfedges_begin();
  }

  // Non intersecting insert methods:

  //! inserts a given curve into the map.
  Halfedge_handle non_intersecting_insert(const X_monotone_curve_2 & cv,
                                          Change_notification * en = NULL)
  { return Planar_map::insert(cv, en); }

  //! iterates through a given range of curves, inserting the curves into the
  // map.
  template <class X_monotone_curve_2_iterator> Halfedge_iterator
  non_intersecting_insert(const X_monotone_curve_2_iterator & begin,
                          const X_monotone_curve_2_iterator & end,
                          Change_notification * en = NULL)
  { return Planar_map::insert(begin, end, en); }

  //! inserts a given curve as a new inner component of a given face.
  Halfedge_handle
  non_intersecting_insert_in_face_interior(const X_monotone_curve_2 & cv, 
                                           Face_handle f,
                                           Change_notification * en = NULL)
  { return Planar_map::insert_in_face_interior(cv, f, en); }

  //! inserts a given curve that one of its endpoints is held by the target
  // vertex of a given halfedge into the map.
  Halfedge_handle
  non_intersecting_insert_from_vertex(const X_monotone_curve_2 & cv, 
                                      Halfedge_handle h,
                                      Change_notification * en = NULL) 
  { return Planar_map::insert_from_vertex(cv, h, en); }

  //! inserts a given curve that both of its endpoints are held by the target
  // vertices of two given halfedges respectively into the map.
  Halfedge_handle
  non_intersecting_insert_at_vertices(const X_monotone_curve_2 & cv, 
                                      Halfedge_handle h1, 
                                      Halfedge_handle h2,
                                      Change_notification * en = NULL)
  { return Planar_map::insert_at_vertices(cv, h1, h2, en); } 

  //! inserts a given curve that one of its endpoints is held by a given vertex
  // into the map.
  Halfedge_handle
  non_intersecting_insert_from_vertex(const X_monotone_curve_2 & cv, 
                                      Vertex_handle v1,
                                      Change_notification * en = NULL) 
  { return Planar_map::insert_from_vertex(cv, v1, en); }

  //! inserts a given curve that both of its endpoints are held by two given
  // vertices respectively into the map.
  Halfedge_handle
  non_intersecting_insert_at_vertices(const X_monotone_curve_2 & cv, 
                                      Vertex_handle v1, Vertex_handle v2,
                                      Change_notification * en = NULL)
  { return Planar_map::insert_at_vertices(cv, v1, v2, en); } 
    
protected:
  bool in_x_range(const X_monotone_curve_2 & cv1,
                  const X_monotone_curve_2 & cv2)
  {
    return  (curve_in_x_range(cv1,cv2) || curve_in_x_range(cv2,cv1));   
  }
  
  bool curve_in_x_range(const X_monotone_curve_2 & cv1,
                        const X_monotone_curve_2 & cv2)
  {
    return ((traits->point_in_x_range(cv1, pmwx_traits->curve_source(cv2)) ||
             traits->point_in_x_range(cv1, pmwx_traits->curve_target(cv2))));
  
  }

  // Data fields:
  Pmwx_traits_wrap * pmwx_traits;
  bool pmwx_use_delete_traits;
  bool pmwx_use_delete_pl;

private:
  bool is_left(const Point_2 &p1, const Point_2 &p2) const
  { return (traits->compare_x(p1, p2) == SMALLER); }

  bool is_right(const Point_2 &p1, const Point_2 &p2) const 
  { return (traits->compare_x(p1, p2) == LARGER); }

  bool is_left_low(const Point_2 &p1, const Point_2 &p2) const
  { return (traits->compare_xy(p1, p2) == SMALLER); }

  bool is_right_top(const Point_2 & p1, const Point_2 & p2) const
  { return is_left_low(p2, p1); }
  
  const Point_2& leftmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_left(p1, p2) ? p1 : p2); }

  const Point_2& rightmost(const Point_2 &p1, const Point_2 &p2) const
  { return (is_right(p1, p2) ? p1 : p2); }

  /*! returns true iff the two points are the same */
  bool point_equal(const Point_2 & p1, const Point_2 & p2) const
  { return traits->point_equal(p1, p2); }
};

/*! Copy Constructor */
template<class Pm>
Planar_map_with_intersections_2<Pm>::
Planar_map_with_intersections_2(const Self & rhs) :
  Planar_map(rhs),
  pmwx_use_delete_traits(false),
  pmwx_use_delete_pl(false)
{
  pmwx_traits = (Pmwx_traits_wrap*)traits;
}

/*! Copy Constructor */
template<class Pm>
Planar_map_with_intersections_2<Pm>::
Planar_map_with_intersections_2(const Planar_map & pm) : 
  Planar_map(pm),
  pmwx_use_delete_traits(false),
  pmwx_use_delete_pl(false)
{
  pmwx_traits = (Pmwx_traits_wrap*)traits;
}

/*! Destructor */
template<class Pm>
Planar_map_with_intersections_2<Pm>::~Planar_map_with_intersections_2()
{
  if (pmwx_use_delete_traits) delete traits;
  if (pmwx_use_delete_pl) delete pl;
}

CGAL_END_NAMESPACE

#endif
