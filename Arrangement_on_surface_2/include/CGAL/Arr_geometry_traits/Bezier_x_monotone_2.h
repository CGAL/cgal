// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>
//                 Iddo Hanniel <iddoh@cs.technion.ac.il>

#ifndef CGAL_BEZIER_X_MONOTONE_2_H
#define CGAL_BEZIER_X_MONOTONE_2_H

/*! \file
 * Header file for the _Bezier_x_monotone_2 class.
 */

#include <CGAL/Arr_geometry_traits/Bezier_curve_2.h>
#include <CGAL/Arr_geometry_traits/Bezier_point_2.h>
#include <CGAL/Arr_geometry_traits/Bezier_cache.h>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of an x-monotone Bezier subcurve, specified by a Bezier curve
 * and two end points.
 */
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_,
          class Bounding_traits_>
class _Bezier_x_monotone_2
{
public:

  typedef Rat_kernel_                             Rat_kernel;
  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  typedef Bounding_traits_                        Bounding_traits;
  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>              Curve_2;
  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>              Point_2;
  typedef _Bezier_x_monotone_2<Rat_kernel,
                               Alg_kernel,
                               Nt_traits,
                               Bounding_traits>         Self;

  typedef _Bezier_cache<Nt_traits>                      Bezier_cache;

private:

  typedef typename Alg_kernel::Point_2            Alg_point_2;
  typedef typename Rat_kernel::Point_2            Rat_point_2;

  typedef typename Nt_traits::Integer             Integer;
  typedef typename Nt_traits::Rational            Rational;
  typedef typename Nt_traits::Algebraic           Algebraic;
  typedef typename Nt_traits::Polynomial          Polynomial;

  typedef typename Point_2::Originator               Originator;
  typedef typename Point_2::Originator_iterator      Originator_iterator;
  typedef typename Bounding_traits::Bez_point_bound  Bez_point_bound;
  typedef typename Bounding_traits::Bez_point_bbox   Bez_point_bbox;

  // Type definition for the vertical tangency-point mapping.
  typedef typename Bezier_cache::Curve_id                 Curve_id;
  typedef std::pair<Curve_id, Curve_id>                   Curve_pair;
  typedef typename Bezier_cache::Vertical_tangency_list   Vert_tang_list;
  typedef typename Bezier_cache::Vertical_tangency_iter   Vert_tang_iter;

  // Type definition for the intersection-point mapping.
  typedef typename Bezier_cache::Intersection_list        Intersect_list;
  typedef typename Bezier_cache::Intersection_iter        Intersect_iter;

  // Representation of an intersection point with its multiplicity:
  typedef std::pair<Point_2, unsigned int>                Intersection_point_2;

  /*! \class Less_intersection_point
   * Comparison functor for intersection points.
   */
  class Less_intersection_point
  {
  private:
    Bezier_cache    *p_cache;
    
  public:

    Less_intersection_point (Bezier_cache& cache) :
      p_cache (&cache)
    {}

    bool operator() (const Intersection_point_2& ip1,
                     const Intersection_point_2& ip2) const
    {
      return (ip1.first.compare_xy (ip2.first, *p_cache) == SMALLER);
    }
  };

  // Type definition for the bounded intersection-point mapping.
  typedef std::list<Point_2>                          Intersection_list;

  /*! \struct
   * An auxiliary functor for comparing pair of curve IDs.
   */
  struct Less_curve_pair
  {
    bool operator() (const Curve_pair& cp1, const Curve_pair& cp2) const
    {
      // Compare the pairs of IDs lexicographically.
      return (cp1.first < cp2.first ||
              (cp1.first == cp2.first && cp1.second < cp2.second));
    }
  };

  /*! \struct Subcurve
   * For the usage of the _exact_vertical_position() function.
   */
  struct Subcurve
  {
    std::list<Rat_point_2>    control_points;
    Rational                  t_min;
    Rational                  t_max;

    /*! Get the rational bounding box of the subcurve. */
    void bbox (Rational& x_min, Rational& y_min,
               Rational& x_max, Rational& y_max) const
    {
      typename std::list<Rat_point_2>::const_iterator  pit =
        control_points.begin();

      CGAL_assertion (pit != control_points.end());

      x_min = x_max = pit->x();
      y_min = y_max = pit->y();

      for (++pit; pit != control_points.end(); ++pit)
      {
        if (CGAL::compare (x_min, pit->x()) == LARGER)
          x_min = pit->x();
        else if (CGAL::compare (x_max, pit->x()) == SMALLER)
          x_max = pit->x();

        if (CGAL::compare (y_min, pit->y()) == LARGER)
          y_min = pit->y();
        else if (CGAL::compare (y_max, pit->y()) == SMALLER)
          y_max = pit->y();
      }

      return;
    }
  };

public:

  typedef std::map<Curve_pair,
                   Intersection_list,
                   Less_curve_pair>               Intersection_map;
  typedef typename Intersection_map::value_type   Intersection_map_entry;
  typedef typename Intersection_map::iterator     Intersection_map_iterator;

private:

  // Data members:
  Curve_2     _curve;         /*! The supporting Bezier curve. */
  Point_2     _ps;            /*! The source point. */
  Point_2     _pt;            /*! The target point. */
  bool        _dir_right;     /*! Is the subcurve directed right (or left). */
  bool        _inc_to_right;  /*! Does the parameter value increase when
                                  traversing the subcurve from left to
                                  right. */
  bool        _is_vert;       /*! Is the subcurve a vertical segment. */

public:

  /*! Default constructor. */
  _Bezier_x_monotone_2 () :
    _dir_right (false),
    _is_vert (false)
  {}

  /*!
   * Constructor given two endpoints.
   * \param B The supporting Bezier curve.
   * \param ps The source point.
   * \param pt The target point.
   * \param cache Caches the vertical tangency points and intersection points.
   * \pre B should be an originator of both ps and pt.
   */
  _Bezier_x_monotone_2 (const Curve_2& B,
                        const Point_2& ps, const Point_2& pt,
                        Bezier_cache& cache);

  /*!
   * Get the supporting Bezier curve.
   */
  const Curve_2& supporting_curve () const
  {
    return (_curve);
  }

  /*!
   * Get the source point.
   */
  const Point_2& source () const
  {
    return (_ps);
  }

  /*!
   * Get the target point.
   */
  const Point_2& target () const
  {
    return (_pt);
  }

  /*!
   * Get the left endpoint (the lexicographically smaller one).
   */
  const Point_2& left () const
  {
    return (_dir_right ? _ps : _pt);
  }

  /*!
   * Get the right endpoint (the lexicographically larger one).
   */
  const Point_2& right () const
  {
    return (_dir_right ? _pt : _ps);
  }

  /*!
   * Check if the subcurve is a vertical segment.
   */
  bool is_vertical () const
  {
    return (_is_vert);
  }

  /*!
   * Check if the subcurve is directed from left to right.
   */
  bool is_directed_right () const
  {
    return (_dir_right);
  }

  /*!
   * Get the approximate parameter range defining the curve.
   * \return A pair of t_src and t_trg, where B(t_src) is the source point
   *         and B(t_trg) is the target point.
   */
  std::pair<double, double> parameter_range () const;

  /*!
   * Get the relative position of the query point with respect to the subcurve.
   * \param p The query point.
   * \param cache Caches the vertical tangency points and intersection points.
   * \pre p is in the x-range of the arc.
   * \return SMALLER if the point is below the arc;
   *         LARGER if the point is above the arc;
   *         EQUAL if p lies on the arc.
   */
  Comparison_result point_position (const Point_2& p,
                                    Bezier_cache& cache) const;

  /*!
   * Compare the relative y-position of two x-monotone subcurve to the right
   * of their intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \param cache Caches the vertical tangency points and intersection points.
   * \pre p is the common left endpoint of both subcurves.
   * \return SMALLER if (*this) lies below cv to the right of p;
   *         EQUAL in case of an overlap (should not happen);
   *         LARGER if (*this) lies above cv to the right of p.
   */
  Comparison_result compare_to_right (const Self& cv,
                                      const Point_2& p,
                                      Bezier_cache& cache) const;

  /*!
   * Compare the relative y-position of two x-monotone subcurve to the left
   * of their intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \param cache Caches the vertical tangency points and intersection points.
   * \pre p is the common right endpoint of both subcurves.
   * \return SMALLER if (*this) lies below cv to the right of p;
   *         EQUAL in case of an overlap (should not happen);
   *         LARGER if (*this) lies above cv to the right of p.
   */
  Comparison_result compare_to_left (const Self& cv,
                                     const Point_2& p,
                                     Bezier_cache& cache) const;

  /*!
   * Check whether the two subcurves are equal (have the same graph).
   * \param cv The other subcurve.
   * \param cache Caches the vertical tangency points and intersection points.
   * \return (true) if the two subcurves have the same graph;
   *         (false) otherwise.
   */
  bool equals (const Self& cv,
               Bezier_cache& cache) const;

  /*!
   * Compute the intersections with the given subcurve.
   * \param cv The other subcurve.
   * \param inter_map Caches the bounded intersection points.
   * \param cache Caches the vertical tangency points and intersection points.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template<class OutputIterator>
  OutputIterator intersect (const Self& cv,
                            Intersection_map& inter_map,
                            Bezier_cache& cache,
                            OutputIterator oi) const
  {
    // In case we have two x-monotone subcurves of the same Bezier curve,
    // the only intersections that may occur are common endpoints.
    if (_curve.is_same (cv._curve))
    {
      if (left().is_same (cv.left()) || left().is_same (cv.right()))
      {
        *oi = CGAL::make_object (Intersection_point_2 (left(), 0));
        ++oi;
      }

      if (right().is_same (cv.left()) || right().is_same (cv.right()))
      {
        *oi = CGAL::make_object (Intersection_point_2 (right(), 0));
        ++oi;
      }

      return (oi);
    }
   
    // Compute the intersections of the two sucurves. Note that for caching
    // purposes we always apply the _intersect() function on the subcurve whose
    // curve ID is smaller.
    std::vector<Intersection_point_2>   ipts;
    Self                                ovlp_cv;
    bool                                do_ovlp;

    CGAL_assertion (_curve.id() != cv._curve.id());
    if (_curve.id() < cv._curve.id())
      do_ovlp = _intersect (cv, inter_map, cache, ipts, ovlp_cv);
    else
      do_ovlp = cv._intersect (*this, inter_map, cache, ipts, ovlp_cv);

    // In case of overlap, just report the overlapping subcurve.
    if (do_ovlp)
    {
      *oi = CGAL::make_object (ovlp_cv);
      ++oi;
      return (oi);
    }

    // If we have a set of intersection points, sort them in ascending
    // xy-lexicorgraphical order, and insert them to the output iterator.
    typename std::vector<Intersection_point_2>::const_iterator  ip_it;

    std::sort (ipts.begin(), ipts.end(), Less_intersection_point (cache));
    for (ip_it = ipts.begin(); ip_it != ipts.end(); ++ip_it)
    {
      *oi = CGAL::make_object (*ip_it);
      ++oi;
    }
    return (oi);
  }

  /*!
   * Split the subcurve into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the subcurve (not one of its endpoints).
   */
  void split (const Point_2& p,
              Self& c1, Self& c2) const;

  /*!
   * Check if the two subcurves are mergeable.
   * \param cv The other subcurve.
   * \return Whether the two subcurves can be merged.
   */
  bool can_merge_with (const Self& cv) const;

  /*!
   * Merge the current arc with the given arc.
   * \param cv The other subcurve.
   * \pre The two arcs are mergeable.
   * \return The merged arc.
   */
  Self merge (const Self& cv) const;

  /*!
   * Flip the subcurve (swap its source and target points).
   * \return The flipped subcurve.
   */
  Self flip () const
  {
    // TODO: Is this "legal"? Should we touch the Bezier curve instead
    // so that _trg > _src in all cases?
    Self  cv = *this;

    cv._ps = this->_pt;
    cv._pt = this->_ps;
    cv._dir_right = ! this->_dir_right;

    return (cv);
  }

private:

  /*!
   * Check if the given t-value is in the range of the subcurve.
   * \param t The parameter value.
   * \param cache Caches the vertical tangency points and intersection points.
   * \return If t in the parameter-range of the subcurve.
   */
  bool _is_in_range (const Algebraic& t,
                     Bezier_cache& cache) const;

  /*!
   * Check if the given point lies in the range of this x-monotone subcurve.
   * \param p The point, which lies on the supporting Bezier curve.
   * \param is_certain Output: Is the answer we provide is certain.
   * \return Whether p is on the x-monotone subcurve.
   */
  bool _is_in_range (const Point_2& p, bool& is_certain) const;

  /*!
   * Compute a y-coordinate of a point on the x-monotone subcurve with a
   * given x-coordinate.
   * \param x0 The given x-coodinate.
   * \param cache Caches the vertical tangency points and intersection points.
   * \return The y-coordinate.
   */
  Algebraic _get_y (const Rational& x0,
                    Bezier_cache& cache) const;

  /*!
   * Compare the slopes of the subcurve with another given Bezier subcurve at
   * their given intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \pre p lies of both subcurves.
   * \pre Neither of the subcurves is a vertical segment.
   * \return SMALLER if (*this) slope is less than cv's;
   *         EQUAL if the two slopes are equal;
   *         LARGER if (*this) slope is greater than cv's.
   */
  Comparison_result _compare_slopes (const Self& cv,
                                     const Point_2& p) const;

  /*!
   * Get the range of t-value over which the subcurve is defined.
   * \param cache Caches the vertical tangency points and intersection points.
   * \return A pair comprised of the t-value for the source point and the
   *         t-value for the target point.
   */
  std::pair<Algebraic, Algebraic> _t_range (Bezier_cache& cache) const;

  /*!
   * Compare the relative y-position of two x-monotone subcurve to the right
   * (or to the left) of their intersection point, whose multiplicity is
   * greater than 1.
   * \param cv The other subcurve.
   * \param p The query point.
   * \param to_right Should we compare to p's right or to p's left.
   * \param cache Caches the vertical tangency points and intersection points.
   * \pre The x-coordinate of the right endpoint of *this is smaller than
   *      (or equal to) the x-coordinate of the right endpoint of cv.
   * \pre p is the common left endpoint of both subcurves.
   * \pre Neither of the subcurves is a vertical segment.
   * \return SMALLER if (*this) lies below cv to the right of p;
   *         EQUAL in case of an overlap (should not happen);
   *         LARGER if (*this) lies above cv to the right of p.
   */
  Comparison_result _compare_to_side (const Self& cv,
                                      const Point_2& p,
                                      bool to_right,
                                      Bezier_cache& cache) const;

  /*!
   * Approximate the intersection points between the two given curves.
   * \param B1 The first Bezier curve.
   * \param B2 The second Bezier curve.
   * \param inter_pts Output: An output list of intersection points.
   * \return Whether all intersection points where successfully approximated.
   */
  bool _approximate_intersection_points (const Curve_2& B1,
                                         const Curve_2& B2,
                                         std::list<Point_2>& inter_pts) const;

  /*!
   * Compute the intersections with the given subcurve.
   * \param cv The other subcurve.
   * \param inter_map Caches the bounded intersection points.
   * \param cache Caches the vertical tangency points and intersection points.
   * \param ipts Output: A vector of intersection points + multiplicities.
   * \param ovlp_cv Output: An overlapping subcurve (if exists).
   * \return Whether an overlap has occured.
   */
  bool _intersect (const Self& cv,
                   Intersection_map& inter_map,
                   Bezier_cache& cache,
                   std::vector<Intersection_point_2>& ipts,
                   Self& ovlp_cv) const;

  /*!
   * Compute the exact vertical position of the given point with respect to
   * the x-monotone curve.
   * \param p The point.
   * \param force_exact Sould we force an exact result.
   * \return SMALLER if the point is below the arc;
   *         LARGER if the point is above the arc;
   *         EQUAL if p lies on the arc.
   */
  Comparison_result _exact_vertical_position (const Point_2& p,
                                              bool force_exact) const;

};

/*!
 * Exporter for Bezier curves.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits, 
          class Bounding_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Bezier_x_monotone_2<Rat_kernel, Alg_kernel, Nt_traits, 
                                       Bounding_traits>& cv)
{
  os << cv.supporting_curve()
     << " | " << cv.source() 
     << " --> " << cv.target();

  return (os);
}

// ---------------------------------------------------------------------------
// Constructor given two endpoints.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_Bezier_x_monotone_2
        (const Curve_2& B,
         const Point_2& ps, const Point_2& pt,
         Bezier_cache& cache) :
  _curve (B),
  _ps (ps),
  _pt (pt),
  _is_vert (false)
{
  // Get the originators of the point that correspond to the curve B.
  Originator_iterator   ps_org = ps.originator(B);
  CGAL_precondition (ps_org != ps.originators_end()); 
  
  Originator_iterator   pt_org = pt.originator(B);
  CGAL_precondition (pt_org != pt.originators_end());
  
  // Check if the subcurve is directed left or right.
  const Comparison_result    res = _ps.compare_x (_pt, cache);
  // Iddo: TODO - alternative check if the original curve is vertical,
  //       a check that can work on intervals.
  
  if (res == EQUAL)
  {
    // We have a vertical segment. Check if the source is below the target.
    _is_vert = true;
    // Iddo: TODO change the check to use compare_xy (or point bbox).
    _dir_right = (CGAL::compare (_ps.y(), _pt.y()) == SMALLER);
  }
  else
  {
    _dir_right = (res == SMALLER);
  }
  
  // Check if the value of the parameter t increases when we traverse the
  // curve from left to right: If the curve is directed to the right, we
  // check if t_src < t_trg, otherwise we check whether t_src > t_trg.
  Comparison_result      t_res;
  
  if (CGAL::compare (ps_org->point_bound().t_max,
                     pt_org->point_bound().t_min) == SMALLER ||
      CGAL::compare (ps_org->point_bound().t_min,
                     pt_org->point_bound().t_max) == LARGER)
  {
    // Perform the comparison assuming that the possible parameter
    // values do not overlap.
    t_res = CGAL::compare (ps_org->point_bound().t_min, 
                           pt_org->point_bound().t_min);
  }
  else
  {
    // In this case both exact parameter values must be known.
    // We use them to perform an exact comparison.
    CGAL_assertion (ps_org->has_parameter() && pt_org->has_parameter());
    
    t_res = CGAL::compare (ps_org->parameter(), pt_org->parameter());
  }
  
  CGAL_precondition (t_res != EQUAL);
  
  if (_dir_right)
    _inc_to_right = (t_res == SMALLER);
  else
    _inc_to_right = (t_res == LARGER);
}

// ---------------------------------------------------------------------------
// Get the approximate parameter range defining the curve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
  std::pair<double, double> 
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::parameter_range () const
{
  // First try to use the approximate representation of the endpoints.
  Originator_iterator  s_org = _ps.originator (_curve);
  CGAL_assertion (s_org != _ps.originators_end());

  Originator_iterator  t_org = _pt.originator (_curve);
  CGAL_assertion (t_org != _pt.originators_end());

  double  t_src = (CGAL::to_double (s_org->point_bound().t_min) +
                   CGAL::to_double (s_org->point_bound().t_max)) / 2;
  double  t_trg = (CGAL::to_double (t_org->point_bound().t_min) +
                   CGAL::to_double (t_org->point_bound().t_max)) / 2;

  return (std::make_pair (t_src, t_trg));
}

// ---------------------------------------------------------------------------
// Get the relative position of the query point with respect to the subcurve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::point_position
    (const Point_2& p,
     Bezier_cache& cache) const
{
  // First check whether p has the same x-coordinate as one of the endpoints.
  const Comparison_result  res1 = p.compare_x (_ps, cache);
  
  if (res1 == EQUAL)
  {
    // In this case both points must be exact.
    CGAL_assertion (p.is_exact() && _ps.is_exact());

    // If both point are rational, compare their rational y-coordinates.
    if (p.is_rational() && _ps.is_rational())
    {
      const Rat_point_2&  rat_p = (Rat_point_2) p;
      const Rat_point_2&  rat_ps = (Rat_point_2) _ps;

      return (CGAL::compare (rat_p.y(), rat_ps.y()));
    }

    // Compare the algebraic y-coordinates.
    return (CGAL::compare (p.y(), _ps.y()));
  }
  
  const Comparison_result  res2 = p.compare_x (_pt, cache);
  
  if (res2 == EQUAL)
  {
    // In this case both points must be exact.
    CGAL_assertion (p.is_exact() && _pt.is_exact());

    // If both point are rational, compare their rational y-coordinates.
    if (p.is_rational() && _pt.is_rational())
    {
      const Rat_point_2&  rat_p = (Rat_point_2) p;
      const Rat_point_2&  rat_pt = (Rat_point_2) _pt;

      return (CGAL::compare (rat_p.y(), rat_pt.y()));
    }

    // Compare the algebraic y-coordinates.
    return (CGAL::compare (p.y(), _pt.y()));
  }
  
  // Make sure that p is in the x-range of our subcurve.
  CGAL_precondition (res1 != res2);
  
  // Check for the case when curve is an originator of the point.
  Originator_iterator   p_org = p.originator(_curve);
 
  if (p_org != p.originators_end())
  {
    Originator_iterator      ps_org = _ps.originator(_curve);
    CGAL_assertion (ps_org != _ps.originators_end());
  
    Originator_iterator      pt_org = _pt.originator(_curve);
    CGAL_assertion (pt_org != _pt.originators_end());

    // Check if the point is in the parameter range of this subcurve.
    // First try an approximate check of the parameter bounds.
    bool  correct_res;
    bool  in_range = false;

    in_range = _is_in_range (p, correct_res);

    if (! correct_res)
    {
      // Perform the comparsion in an exact manner.
      if (! p.is_exact())
        p.make_exact (cache);

      CGAL_assertion (p_org->has_parameter());

      in_range = _is_in_range (p_org->parameter(), cache);
    }

    if (in_range)
      return (EQUAL);
  }
  
  // Call the vertical-position function that uses the bounding-boxes
  // to evaluate the comparsion result.
  typename Bounding_traits::Control_point_vec cp;

  std::copy (_curve.control_points_begin(), _curve.control_points_end(),
             std::back_inserter(cp));

  Originator_iterator           ps_org = _ps.originator(_curve);
  CGAL_assertion (ps_org != _ps.originators_end());
  
  Originator_iterator           pt_org = _pt.originator(_curve);
  CGAL_assertion (pt_org != _pt.originators_end());
  
  Comparison_result             res_bound = EQUAL;
  typename Bounding_traits::NT  x_min, y_min, x_max, y_max;
  bool                          can_refine;

  p.get_bbox (x_min, y_min, x_max, y_max);

  if (CGAL::compare (ps_org->point_bound().t_max,
                     pt_org->point_bound().t_min) == SMALLER)
  {
    // Examine the parameter range of the originator of the source point
    // with respect to the current subcurve B, and make sure that B(t_max)
    // lies to the left of p if the curve is directed from left to right
    // (or to the right of p, if the subcurve is directed from right to left).
    can_refine = ! _ps.is_exact();
    do
    {
      const Rat_point_2&  ps = _curve (ps_org->point_bound().t_max);

      if ((_dir_right && CGAL::compare (ps.x(), x_min) != LARGER) ||
          (! _dir_right && CGAL::compare (ps.x(), x_max) != SMALLER))
        break;

      if (can_refine)
        can_refine = _ps.refine();
    } while (can_refine);

    // Examine the parameter range of the originator of the target point
    // with respect to the current subcurve B, and make sure that B(t_min)
    // lies to the right of p if the curve is directed from left to right
    // (or to the left of p, if the subcurve is directed from right to left).
    can_refine = ! _pt.is_exact();
    do
    {
      const Rat_point_2&  pt = _curve (pt_org->point_bound().t_min);

      if ((_dir_right && CGAL::compare (pt.x(), x_max) != SMALLER) ||
          (! _dir_right && CGAL::compare (pt.x(), x_min) != LARGER))
        break;

      if (can_refine)
        can_refine = _pt.refine();
    } while (can_refine);

    // In this case the parameter value of the source is smaller than the
    // target's, so we compare the point with the subcurve of B defined over
    // the proper parameter range.
    res_bound = p.vertical_position (cp,
                                     ps_org->point_bound().t_max,
                                     pt_org->point_bound().t_min);
  }
  else if (CGAL::compare (pt_org->point_bound().t_max,
                          ps_org->point_bound().t_min) == SMALLER)
  {
    // Examine the parameter range of the originator of the source point
    // with respect to the current subcurve B, and make sure that B(t_min)
    // lies to the left of p if the curve is directed from left to right
    // (or to the right of p, if the subcurve is directed from right to left).
    can_refine = ! _ps.is_exact();
    do
    {
      const Rat_point_2&  ps = _curve (ps_org->point_bound().t_min);

      if ((_dir_right && CGAL::compare (ps.x(), x_min) != LARGER) ||
          (! _dir_right && CGAL::compare (ps.x(), x_max) != SMALLER))
        break;

      if (can_refine)
        can_refine = _ps.refine();
    } while (can_refine);

    // Examine the parameter range of the originator of the target point
    // with respect to the current subcurve B, and make sure that B(t_max)
    // lies to the right of p if the curve is directed from left to right
    // (or to the left of p, if the subcurve is directed from right to left).
    can_refine = ! _pt.is_exact();
    do
    {
      const Rat_point_2&  pt = _curve (pt_org->point_bound().t_max);

      if ((_dir_right && CGAL::compare (pt.x(), x_max) != SMALLER) ||
          (! _dir_right && CGAL::compare (pt.x(), x_min) != LARGER))
        break;

      if (can_refine)
        can_refine = _pt.refine();
    } while (can_refine);

    // In this case the parameter value of the source is large than the
    // target's, so we compare the point with the subcurve of B defined over
    // the proper parameter range.
    res_bound = p.vertical_position (cp,
                                     pt_org->point_bound().t_max,
                                     ps_org->point_bound().t_min);
  }

  if (res_bound != EQUAL)
    return (res_bound);
 
  // In this case we have to switch to exact computations and check whether
  // p lies of the given subcurve. We take one of p's originating curves and
  // compute its intersections with our x-monotone curve.
  if (! p.is_exact())
    p.make_exact (cache);

  CGAL_assertion (p.originators_begin() != p.originators_end());
    
  // Iddo: If the point is a rational point (e.g., ray shooting)
  //       use comparison between Y(root_of(X0-X(t))) and Y0.
  Originator   org = *(p.originators_begin());
  bool         do_ovlp;
  bool         swap_order = (_curve.id() > org.curve().id());
  const Intersect_list&  inter_list = (! swap_order ?
    (cache.intersections (_curve.id(),
                              _curve.x_polynomial(), _curve.x_norm(),
                              _curve.y_polynomial(), _curve.y_norm(),
                              org.curve().id(),
                              org.curve().x_polynomial(), org.curve().x_norm(),
                              org.curve().y_polynomial(), org.curve().y_norm(),
                              do_ovlp)) :
    (cache.intersections (org.curve().id(),
                              org.curve().x_polynomial(), org.curve().x_norm(),
                              org.curve().y_polynomial(), org.curve().y_norm(),
                              _curve.id(),
                              _curve.x_polynomial(), _curve.x_norm(),
                              _curve.y_polynomial(), _curve.y_norm(),
                              do_ovlp)));

  if (do_ovlp)
    return (EQUAL);
    
  // Go over the intersection points and look for p there.
  Intersect_iter       iit;
    
  for (iit = inter_list.begin(); iit != inter_list.end(); ++iit)
  {
    // Get the parameter of the originator and compare it to p's parameter.
    const Algebraic&  s = swap_order ? iit->s : iit->t;
      
    if (CGAL::compare (s, org.parameter()) == EQUAL)
    {
      // Add this curve as an originator for p.
      const Algebraic&  t = swap_order ? iit->t : iit->s;
    
      CGAL_assertion (_is_in_range (t, cache));

      Point_2&  pt = const_cast<Point_2&> (p);
      pt.add_originator (Originator (_curve, t));

      // The point p lies on the subcurve.
      return (EQUAL);
    }
  }
  
  // We now that p is not on the subcurve. We therefore subdivide the curve
  // using exact rational arithmetic, until we reach a separation between the
  // curve and the point. (This case should be very rare.)
  // Note that we first try to work with inexact endpoint representation, and
  // only if we fail we make the endpoints of the x-monotone curves exact.
  if (! p.is_exact())
    p.make_exact (cache);

  Comparison_result  exact_res = _exact_vertical_position (p, false);

  if (exact_res != EQUAL)
    return (exact_res);

  if (! _ps.is_exact())
    _ps.make_exact (cache);
  
  if (! _pt.is_exact())
    _pt.make_exact (cache);

  return (_exact_vertical_position (p, true));
}

// ---------------------------------------------------------------------------
// Compare the relative y-position of two x-monotone subcurves to the right
// of their intersection point.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::compare_to_right
        (const Self& cv,
         const Point_2& p,
         Bezier_cache& cache) const
{
  CGAL_precondition (p.compare_xy (right(), cache) != LARGER);
  CGAL_precondition (p.compare_xy (cv.right(), cache) != LARGER);
  
  if (this == &cv)
    return (EQUAL);
  
  // Make sure that p equals the left endpoint of both curves. Note that this
  // is important to carry out this test, as it assures us the eventually both
  // curves are originators of p.
  if (! p.equals (left(), cache) || ! p.equals (cv.left(), cache))
  {
    CGAL_assertion (p.originator(_curve) != p.originators_end() &&
                    p.originator(cv._curve) != p.originators_end());
  }

  // Check for vertical subcurves. A vertical segment is above any other
  // x-monotone subcurve to the right of their common endpoint.
  if (is_vertical())
  {
    if (cv.is_vertical())
      // Both are vertical segments with a common endpoint, so they overlap:
      return (EQUAL);
    
    return (LARGER);
  }
  else if (cv.is_vertical())
  {
    return (SMALLER);
  }
    
  // Check if both subcurves originate from the same Bezier curve.
  Nt_traits       nt_traits;

  if (_curve.is_same (cv._curve))
  {
    // Get the originator, and make sure that p is a vertical tangency
    // point of this originator.
    Originator_iterator  org = p.originator(_curve);
    
    CGAL_assertion (org != p.originators_end());    
    CGAL_assertion (_inc_to_right != cv._inc_to_right);
    
    if (! p.is_exact())
    {
      CGAL_assertion (org->point_bound().point_type ==
                      Bez_point_bound::VERTICAL_TANGENCY_PT);

      // Comparison based on the control polygon of the bounded vertical
      // tangency point, using the fact this polygon is y-monotone.
      const typename Bounding_traits::Control_point_vec& cp = 
        org->point_bound().bounding_polyline;
      
      if (_inc_to_right)
      {
        return (CGAL::compare (cp.back().y(), cp.front().y()));
      }
      else
      {
        return (CGAL::compare (cp.front().y(), cp.back().y()));
      }
    }
    
    // Iddo: Handle rational points (using de Casteljau derivative)?
    
    // In this case we know that we have a vertical tangency at t0, so
    // X'(t0) = 0. We evaluate the sign of Y'(t0) in order to find the
    // vertical position of the two subcurves to the right of this point.
    CGAL_assertion (org->has_parameter());
    
    const Algebraic&  t0 = org->parameter();
    Polynomial        polyY_der = nt_traits.derive (_curve.y_polynomial());
    const CGAL::Sign  sign_der =
      CGAL::sign (nt_traits.evaluate_at (polyY_der, t0));
    
    CGAL_assertion (sign_der != CGAL::ZERO);
    
    if (_inc_to_right)
      return ((sign_der == CGAL::POSITIVE) ? LARGER : SMALLER);
    else
      return ((sign_der == CGAL::NEGATIVE) ? LARGER : SMALLER);
  }
  
  // Compare the slopes of the two supporting curves at p. In the general
  // case, the slopes are not equal and their comparison gives us the
  // vertical order to p's right. 
  Comparison_result   slope_res = _compare_slopes (cv, p);
  
  if (slope_res != EQUAL)
    return (slope_res);

  // Compare the two subcurves by choosing some point to the right of p
  // and comparing the vertical position there.
  Comparison_result   right_res;
  
  if (right().compare_x (cv.right(), cache) != LARGER)
  {
    right_res = _compare_to_side (cv, p,
                                  true,           // Compare to p's right.
                                  cache);
  }
  else
  {
    right_res = cv._compare_to_side (*this, p,
                                     true,        // Compare to p's right.
                                     cache);
    
    right_res = CGAL::opposite (right_res);
  }
  
  return (right_res);
}

// ---------------------------------------------------------------------------
// Compare the relative y-position of two x-monotone subcurve to the left
// of their intersection point.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::compare_to_left
        (const Self& cv,
         const Point_2& p,
         Bezier_cache& cache) const
{
  CGAL_precondition (p.compare_xy (left(), cache) != SMALLER);
  CGAL_precondition (p.compare_xy (cv.left(), cache) != SMALLER);
  
  if (this == &cv)
    return (EQUAL);

  // Make sure that p equals the right endpoint of both curves. Note that this
  // is important to carry out this test, as it assures us the eventually both
  // curves are originators of p.
  if (! p.equals (right(), cache) || ! p.equals (cv.right(), cache))
  {
    CGAL_assertion (p.originator(_curve) != p.originators_end() &&
                    p.originator(cv._curve) != p.originators_end());
  }
  
  // Check for vertical subcurves. A vertical segment is below any other
  // x-monotone subcurve to the left of their common endpoint.
  if (is_vertical())
  {
    if (cv.is_vertical())
      // Both are vertical segments with a common endpoint, so they overlap:
      return (EQUAL);
    
    return (SMALLER);
  }
  else if (cv.is_vertical())
  {
    return (LARGER);
  }
    
  // Check if both subcurves originate from the same Bezier curve.
  Nt_traits       nt_traits;

  if (_curve.is_same (cv._curve))
  {
    // Get the originator, and make sure that p is a vertical tangency
    // point of this originator.
    Originator_iterator  org = p.originator(_curve);
    
    CGAL_assertion (org != p.originators_end());    
    CGAL_assertion (_inc_to_right != cv._inc_to_right);
    
    if (! p.is_exact())
    {
      CGAL_assertion (org->point_bound().point_type ==
                      Bez_point_bound::VERTICAL_TANGENCY_PT);

      // Comparison based on the control polygon of the bounded vertical
      // tangency point, using the fact this polygon is y-monotone.
      const typename Bounding_traits::Control_point_vec& cp = 
        org->point_bound().bounding_polyline;
      
      if (_inc_to_right)
      {
        return (CGAL::compare(cp.front().y(), cp.back().y()));
      }
      else
      {
        return (CGAL::compare(cp.back().y(), cp.front().y()));
      }
    }
    
    // Iddo: Handle rational points (using de Casteljau derivative)?
    
    // In this case we know that we have a vertical tangency at t0, so
    // X'(t0) = 0. We evaluate the sign of Y'(t0) in order to find the
    // vertical position of the two subcurves to the right of this point.
    CGAL_assertion (org->has_parameter());
    
    const Algebraic&  t0 = org->parameter();
    Polynomial        polyY_der = nt_traits.derive (_curve.y_polynomial());
    const CGAL::Sign  sign_der =
      CGAL::sign (nt_traits.evaluate_at (polyY_der, t0));
    
    CGAL_assertion (sign_der != CGAL::ZERO);
    
    if (_inc_to_right)
      return ((sign_der == CGAL::NEGATIVE) ? LARGER : SMALLER);
    else
      return ((sign_der == CGAL::POSITIVE) ? LARGER : SMALLER);
  }
  
  // Compare the slopes of the two supporting curves at p. In the general
  // case, the slopes are not equal and their comparison gives us the
  // vertical order to p's right; note that we swap the order of the curves
  // to obtain their position to the left.
  Comparison_result   slope_res = cv._compare_slopes (*this, p);
  
  if (slope_res != EQUAL)
    return (slope_res);
  
  // Compare the two subcurves by choosing some point to the left of p
  // and compareing the vertical position there.
  Comparison_result   left_res;
  
  if (left().compare_x (cv.left(), cache) != SMALLER)
  {
    left_res = _compare_to_side (cv, p,
                                 false,          // Compare to p's left.
                                 cache);
  }
  else
  {
    left_res = cv._compare_to_side (*this, p,
                                    false,       // Compare to p's left.
                                    cache);
    left_res = CGAL::opposite (left_res);
  }
  
  return (left_res);
}

// ---------------------------------------------------------------------------
// Check whether the two subcurves are equal (have the same graph).
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::equals
        (const Self& cv,
         Bezier_cache& cache) const
{
  // Check if the two subcurve have overlapping supporting curves.
  if (! _curve.is_same (cv._curve))
  {
    // Ron: For the time being, we do not deal with this case ...
    /*
    if (! _curve.has_same_support (cv._curve))
      return (false);
    
    // Mark that the two curves overlap.
    const Curve_id               id1 = _curve.id();
    const Curve_id               id2 = cv._curve.id();
    Curve_pair                   curve_pair;
    
    if (id1 < id2)
      curve_pair = Curve_pair (id1, id2);
    else
      curve_pair = Curve_pair (id2, id1);
    
    if (inter_map.find (curve_pair) == inter_map.end())
    {
      Intersection_info&  info = inter_map[curve_pair];
      info.second = true;
    }
    */
  }

  // Check for equality of the endpoints.
  return (left().equals (cv.left(), cache) &&
          right().equals (cv.right(), cache));
}

// ---------------------------------------------------------------------------
// Split the subcurve into two at a given split point.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
void _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::split
        (const Point_2& p,
         Self& c1, Self& c2) const
{
  CGAL_precondition (p.originator(_curve) != p.originators_end());

  // Duplicate the curve.
  c1 = c2 = *this;
    
  // Perform the split.
  if (_dir_right)
  {
    c1._pt = p;      
    c2._ps = p;
  }
  else
  {
    c1._ps = p;
    c2._pt = p;
  }
  
  return;
}

// ---------------------------------------------------------------------------
// Check if the two subcurves are mergeable.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::can_merge_with
        (const Self& cv) const
{
  // Note that we only allow merging subcurves of the same originating
  // Bezier curve (overlapping curves will not do in this case).
  return (_curve.is_same (cv._curve) &&
          (right().is_same (cv.left()) || left().is_same (cv.right())));
  
  return (false);
}

// ---------------------------------------------------------------------------
// Merge the current arc with the given arc.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
typename _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::Self
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::merge
        (const Self& cv) const
{
  CGAL_precondition (_curve.is_same (cv._curve));
  
  Self    res = cv;
  
  if (right().is_same (cv.left()))
  {
    // Extend the subcurve to the right.
    if (_dir_right)
    {
      res._pt = cv.right();
    }
    else
    {
      res._ps = cv.right();
    }
  }
  else
  {
    CGAL_precondition (left().is_same (cv.right()));
    
    // Extend the subcurve to the left.
    if (_dir_right)
    {
      res._ps = cv.left();
    }
    else
    {
      res._pt = cv.left();
    }
  }

  return (res);
}

// ---------------------------------------------------------------------------
// Check if the given t-value is in the range of the subcurve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_is_in_range
        (const Algebraic& t,
         Bezier_cache& cache) const
{
  // First try to use the approximate representation of the endpoints.
  Originator_iterator  s_org = _ps.originator (_curve);
  CGAL_assertion (s_org != _ps.originators_end());

  Originator_iterator  t_org = _pt.originator (_curve);
  CGAL_assertion (t_org != _pt.originators_end());

  bool  p_lt_ps = (CGAL::compare (t,
                                  Algebraic (s_org->point_bound().t_min)) == 
                   SMALLER);
  bool  p_gt_ps = (CGAL::compare (t,
                                  Algebraic (s_org->point_bound().t_max)) ==
                   LARGER);
  bool  p_lt_pt = (CGAL::compare (t,
                                  Algebraic (t_org->point_bound().t_min)) ==
                   SMALLER);
  bool  p_gt_pt = (CGAL::compare (t,
                                  Algebraic (t_org->point_bound().t_max)) ==
                   LARGER);

  if ((p_gt_ps && p_lt_pt) || (p_lt_ps && p_gt_pt))
  {
    // The point p is definately in the x-range of the subcurve, as its
    // parameter is between the source and target parameters.
    return (true);
  }
  
  if ((p_lt_ps && p_lt_pt) || (p_gt_ps && p_gt_pt))
  {
    // The point p is definately not in the x-range of the subcurve,
    // as its parameter is smaller than both source and target parameter
    // (or greater than both of them).
    return (false);
  }
  
  // Obtain the exact t-range of the curve and peform an exact comparison.
  std::pair<Algebraic, Algebraic> range = _t_range (cache);
  const Algebraic&                t_src = range.first;
  const Algebraic&                t_trg = range.second;
  
  const Comparison_result  res1 = CGAL::compare (t, t_src);
  const Comparison_result  res2 = CGAL::compare (t, t_trg);
  
  return (res1 == EQUAL || res2 == EQUAL || res1 != res2);
}

// ---------------------------------------------------------------------------
// Check if the given point lies in the range of this x-monotone subcurve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_is_in_range
        (const Point_2& p,
         bool& is_certain) const
{
  is_certain = true;
  
  // Check the easy case that p is one of the subcurve endpoints.
  if (p.is_same(_ps) || p.is_same(_pt))
    return true;

  // Compare the parameter of p with the parameters of the endpoints.
  Originator_iterator  p_org = p.originator (_curve);
  CGAL_assertion (p_org != p.originators_end());
  
  Originator_iterator  s_org = _ps.originator (_curve);
  CGAL_assertion (s_org != _ps.originators_end());
  
  Originator_iterator  t_org = _pt.originator (_curve);
  CGAL_assertion (t_org != _pt.originators_end());
  
  bool      can_refine_p = ! p.is_exact();
  bool      can_refine_s = ! _ps.is_exact();
  bool      can_refine_t = ! _pt.is_exact();
  
  while (can_refine_p || can_refine_s || can_refine_t)
  {
    bool  p_lt_ps = (CGAL::compare (p_org->point_bound().t_max,
                                    s_org->point_bound().t_min) == SMALLER);
    bool  p_gt_ps = (CGAL::compare (p_org->point_bound().t_min,
                                    s_org->point_bound().t_max) == LARGER);
    bool  p_lt_pt = (CGAL::compare (p_org->point_bound().t_max,
                                    t_org->point_bound().t_min) == SMALLER);
    bool  p_gt_pt = (CGAL::compare (p_org->point_bound().t_min,
                                    t_org->point_bound().t_max) == LARGER);

    if ((p_gt_ps && p_lt_pt) || (p_lt_ps && p_gt_pt))
    {
      // The point p is definately in the x-range of the subcurve, as its
      // parameter is between the source and target parameters.
      return (true);
    }

    if ((p_lt_ps && p_lt_pt) || (p_gt_ps && p_gt_pt))
    {
      // The point p is definately not in the x-range of the subcurve,
      // as its parameter is smaller than both source and target parameter
      // (or greater than both of them).
      return (false);
    }
    
    // Try to refine the points.
    if (can_refine_p)
      can_refine_p = p.refine();
    
    if (can_refine_s)
      can_refine_s = _ps.refine();
    
    if (can_refine_t)
      can_refine_t = _pt.refine();
  }

  // If we reached here, we do not have a certain answer.
  is_certain = false;
  return (false);
}

// ---------------------------------------------------------------------------
// Compute a y-coordinate of a point on the x-monotone subcurve with a
// given x-coordinate.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
typename _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::Algebraic
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_get_y
        (const Rational& x0,
         Bezier_cache& cache) const
{
  // Obtain the t-values for with the x-coordinates of the supporting
  // curve equal x0.
  std::list<Algebraic>  t_vals;
  
  _curve.t_at_x (x0, std::back_inserter(t_vals));
  
  // Find a t-value that is in the range of the current curve.
  Nt_traits                                nt_traits;
  typename std::list<Algebraic>::iterator  t_iter;
  std::pair<Algebraic, Algebraic>          t_range = _t_range (cache);
  const Algebraic&                         t_src = t_range.first;
  const Algebraic&                         t_trg = t_range.second;
  Comparison_result                        res1, res2;
  
  for (t_iter = t_vals.begin(); t_iter != t_vals.end(); ++t_iter)
  {
    res1 = CGAL::compare (*t_iter, t_src);
    
    if (res1 == EQUAL)
    {
      // Return the y-coordinate of the source point:
      return (_ps.y());
    }
    
    res2 = CGAL::compare (*t_iter, t_trg);
    
    if (res2 == EQUAL)
    {
      // Return the y-coordinate of the source point:
      return (_pt.y());
    }

    if (res1 != res2)
    {
      // We found a t-value in the range of our x-monotone subcurve.
      // Use this value to compute the y-coordinate.
      return (nt_traits.evaluate_at (_curve.y_polynomial(), *t_iter) /
              nt_traits.convert (_curve.y_norm()));
    }
  }
  
  // If we reached here, x0 is not in the x-range of our subcurve.
  CGAL_assertion (false);
  return (0);
}

// ---------------------------------------------------------------------------
// Compare the slopes of the subcurve with another given Bezier subcurve at
// their given intersection point.
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_compare_slopes
        (const Self& cv,
         const Point_2& p) const
{
  // Get the originators of p.
  Originator_iterator     org1 = p.originator (_curve);
  CGAL_assertion (org1 != p.originators_end());

  Originator_iterator     org2 = p.originator (cv._curve);
  CGAL_assertion (org2 != p.originators_end());
  
  // If the point is only approximated, we can carry out a comparison using
  // an approximate number type.
  if (! p.is_exact())
  {
    // If the point is inexact, we assume it is a bounded intersection
    // point of two curves, and therefore the bounding angle these curves
    // span do not overlap.
    const Bez_point_bound&  bound1 = org1->point_bound();
    const Bez_point_bound&  bound2 = org2->point_bound();
    Bounding_traits         bound_tr;
    
    return (bound_tr.compare_slopes_of_bounded_intersection_point (bound1,
                                                                   bound2));
  }
  
  // Iddo: Implement slope comparison at a Rational point
  //       (using de Casteljau).
  
  // The slope of (X(t), Y(t)) at t0 is Y'(t)/X'(t).
  // Compute the slope of (*this).
  // Note that we take special care of the case X'(t) = 0, when the tangent
  // is vertical and its slope is +/- oo.
  CGAL_assertion (org1->has_parameter() && org2->has_parameter());
  
  Nt_traits         nt_traits;
  const Algebraic&  t1 = org1->parameter();
  Polynomial        derivX = nt_traits.derive (_curve.x_polynomial());
  Polynomial        derivY = nt_traits.derive (_curve.y_polynomial());
  Algebraic         numer1 = nt_traits.evaluate_at (derivY, t1) *
                             nt_traits.convert (_curve.x_norm());
  Algebraic         denom1 = nt_traits.evaluate_at (derivX, t1) *
                             nt_traits.convert (_curve.y_norm());
  CGAL::Sign        inf_slope1 = CGAL::ZERO;
  Algebraic         slope1;

  if (CGAL::sign (denom1) == CGAL::ZERO)
  {
    inf_slope1 = CGAL::sign (numer1);

    // If also Y'(t) = 0, we cannot perform the comparison:
    if (inf_slope1 == CGAL::ZERO)
      return (EQUAL);
  }
  else
  {
    slope1 = numer1 / denom1;
  }

  // Compute the slope of the other subcurve.
  const Algebraic&  t2 = org2->parameter();
  derivX = nt_traits.derive (cv._curve.x_polynomial());
  derivY = nt_traits.derive (cv._curve.y_polynomial());
  Algebraic         numer2 = nt_traits.evaluate_at (derivY, t2) *
                             nt_traits.convert (cv._curve.x_norm());
  Algebraic         denom2 = nt_traits.evaluate_at (derivX, t2) *
                             nt_traits.convert (cv._curve.y_norm());
  CGAL::Sign        inf_slope2 = CGAL::ZERO;
  Algebraic         slope2;

  if (CGAL::sign (denom2) == CGAL::ZERO)
  {
    inf_slope2 = CGAL::sign (numer2);

    // If also Y'(t) = 0, we cannot perform the comparison:
    if (inf_slope2 == CGAL::ZERO)
      return (EQUAL);
  }
  else
  {
    slope2 = numer2 / denom2;
  }

  // Handle the comparison when one slope (or both) is +/- oo.
  if (inf_slope1 == CGAL::POSITIVE)
    return (inf_slope2 == CGAL::POSITIVE ? EQUAL : LARGER);
  
  if (inf_slope1 == CGAL::NEGATIVE)
    return (inf_slope2 == CGAL::NEGATIVE ? EQUAL : SMALLER);

  if (inf_slope2 == CGAL::POSITIVE)
    return (SMALLER);

  if (inf_slope2 == CGAL::NEGATIVE)
    return (LARGER);

  // Compare the slopes.
  return (CGAL::compare (slope1, slope2));
}

// ---------------------------------------------------------------------------
// Get the range of t-value over which the subcurve is defined.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
std::pair<typename _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt,
                                        BndTrt>::Algebraic, 
          typename _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt,
                                        BndTrt>::Algebraic>
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_t_range
        (Bezier_cache& cache) const
{
  Originator_iterator  ps_org = _ps.originator(_curve);
  CGAL_assertion(ps_org != _ps.originators_end());
  
  Originator_iterator  pt_org = _pt.originator(_curve);
  CGAL_assertion(pt_org != _pt.originators_end());
  
  // Make sure that the two endpoints are exact.
  if (! ps_org->has_parameter())
    _ps.make_exact (cache);
  
  if (! pt_org->has_parameter())
    _pt.make_exact (cache);
  
  return (std::make_pair (ps_org->parameter(),
                          pt_org->parameter()));
}

// ---------------------------------------------------------------------------
// Compare the relative y-position of two x-monotone subcurve to the right
// (or to the left) of their intersection point, whose multiplicity is
// greater than 1.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_compare_to_side
        (const Self& cv,
         const Point_2& p,
         bool to_right,
         Bezier_cache& cache) const
{
  // Get the intersection points of the two curves from the cache. Note that
  // we make sure that the ID of this->_curve is smaller than of cv's curve ID.
  const bool             no_swap_curves = (_curve.id() < cv._curve.id());
  bool                   do_ovlp;
  const Intersect_list&  inter_list =
    (no_swap_curves ?
     (cache.intersections (_curve.id(),
                               _curve.x_polynomial(), _curve.x_norm(),
                               _curve.y_polynomial(), _curve.y_norm(),
                               cv._curve.id(),
                               cv._curve.x_polynomial(), cv._curve.x_norm(),
                               cv._curve.y_polynomial(), cv._curve.y_norm(),
                               do_ovlp)) :
     (cache.intersections (cv._curve.id(),
                               cv._curve.x_polynomial(), cv._curve.x_norm(),
                               cv._curve.y_polynomial(), cv._curve.y_norm(),
                               _curve.id(),
                               _curve.x_polynomial(), _curve.x_norm(),
                               _curve.y_polynomial(), _curve.y_norm(),
                               do_ovlp)));

  if (do_ovlp)
    return (EQUAL);
  
  // Get the parameter value for the point p.
  Originator_iterator          org = p.originator(_curve);
  
  CGAL_assertion (org != p.originators_end());
  CGAL_assertion (org->has_parameter());
  
  const Algebraic&             t0 = org->parameter();

  // Get the parameter range of the curve.
  const std::pair<Algebraic,
                  Algebraic>&  range = _t_range (cache);
  const Algebraic&             t_src = range.first;
  const Algebraic&             t_trg = range.second;

  // Find the next intersection point that lies to the right of p.
  Intersect_iter               iit;
  Algebraic                    next_t;
  Comparison_result            res = CGAL::EQUAL;
  bool                         found = false;
  
  for (iit = inter_list.begin(); iit != inter_list.end(); ++iit)
  {
    // Check if the current point lies to the right (left) of p. We do so by
    // considering its originating parameter value s (or t, if we swapped
    // the curves).
    const Algebraic&     t = (no_swap_curves ? (iit->s) : iit->t);
    
    res = CGAL::compare (t, t0);
    if ((to_right && ((_inc_to_right && res == LARGER) ||
                      (! _inc_to_right && res == SMALLER))) ||
        (! to_right && ((_inc_to_right && res == SMALLER) ||
                        (! _inc_to_right && res == LARGER))))
    {
      if (! found)
      {
        next_t = t;
        found = true;
      }
      else
      {
        // If we have already located an intersection point to the right
        // (left) of p, choose the leftmost (rightmost) of the two points.
        res = CGAL::compare (t, next_t);
        if ((to_right && ((_inc_to_right && res == SMALLER) ||
                          (! _inc_to_right && res == LARGER))) ||
            (! to_right && ((_inc_to_right && res == LARGER) ||
                            (! _inc_to_right && res == SMALLER))))
        {
          next_t = t;
        }
      }
    }
  }
  
  // If the next intersection point occurs before the right (left) endpoint
  // of the subcurve, keep it. Otherwise, take the parameter value at
  // the endpoint.
  if (found)
  {
    if (to_right == _dir_right)
      res = CGAL::compare (t_trg, next_t);
    else
      res = CGAL::compare (t_src, next_t);
  }
  
  if (! found ||
      (to_right && ((_inc_to_right && res == SMALLER) ||
                    (! _inc_to_right && res == LARGER))) ||
      (! to_right && ((_inc_to_right && res == LARGER) ||
                      (! _inc_to_right && res == SMALLER))))
  {
    next_t = ((to_right == _dir_right) ? t_trg : t_src);
  }
  
  // Find a rational value between t0 and t_next. Using this value, we
  // a point with rational coordinates on our subcurve. We also locate a point
  // on the other curve with the same x-coordinates.
  Nt_traits           nt_traits;
  const Rational&     mid_t = nt_traits.rational_in_interval (t0, next_t);
  const Rat_point_2&  q1 = _curve (mid_t);
  const Algebraic&    y2 = cv._get_y (q1.x(), cache);

  // We now just have to compare the y-coordinates of the two points we have
  // computed.
  return (CGAL::compare (nt_traits.convert (q1.y()), y2));
}

// ---------------------------------------------------------------------------
// Approximate the intersection points between the two given curves.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt,
                          BndTrt>::_approximate_intersection_points
        (const Curve_2& B1,
         const Curve_2& B2,
         std::list<Point_2>& inter_pts) const
{
  inter_pts.clear();

  // Make local copies of the control polygons.
  typename Bounding_traits::Control_point_vec  cp1, cp2;
  
  std::copy (B1.control_points_begin(), B1.control_points_end(),
             std::back_inserter (cp1));
  std::copy (B2.control_points_begin(), B2.control_points_end(),
             std::back_inserter(cp2));

  // Use the bounding traits to isolate the intersection points.
  Bounding_traits                              bound_tr;
  typename Bounding_traits::Bound_pair_lst     inter_pairs;

  bound_tr.CrvCrvInter (cp1, cp2, 
                        inter_pairs);
  
  // Construct the approximated points.
  typename Bounding_traits::Bound_pair_lst::const_iterator iter;
  
  for (iter = inter_pairs.begin(); iter != inter_pairs.end(); ++iter)
  {
    const Bez_point_bound&  bound1 = iter->bound1; 
    const Bez_point_bound&  bound2 = iter->bound2; 
    const Bez_point_bbox&   bbox = iter->bbox; 
    
    // In case it is impossible to further refine the point, stop here.
    if (! bound1.can_refine || ! bound2.can_refine)
      return (false);
    
    // Create the approximated intersection point.
    Point_2                 pt;
    
    if (bound1.point_type == Bounding_traits::Bez_point_bound::RATIONAL_PT &&
        bound2.point_type == Bounding_traits::Bez_point_bound::RATIONAL_PT)
    {
      CGAL_assertion (CGAL::compare (bound1.t_min, bound1.t_max) == EQUAL); 
      CGAL_assertion (CGAL::compare (bound2.t_min, bound2.t_max) == EQUAL); 
      Rational  t1 = bound1.t_min;
      Rational  t2 = bound2.t_min;

      pt = Point_2 (B1, t1);
      pt.add_originator (Originator (B2, t2));
    }
    else
    {
      pt.add_originator (Originator (B1, bound1));
      pt.add_originator (Originator (B2, bound2));
    }
    pt.set_bbox (bbox);
    
    inter_pts.push_back (pt);
  }

  // The approximation process ended OK.
  return (true);
}

// ---------------------------------------------------------------------------
// Compute the intersections with the given subcurve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_intersect
    (const Self& cv,
     Intersection_map& inter_map,
     Bezier_cache& cache,
     std::vector<Intersection_point_2>& ipts,
     Self& /* ovlp_cv */) const
{
  CGAL_precondition (_curve.id() < cv._curve.id());

  ipts.clear();

  // Construct the pair of curve IDs and look for it in the intersection map.
  Curve_pair                 curve_pair (_curve.id(), cv._curve.id());
  Intersection_map_iterator  map_iter = inter_map.find (curve_pair);
  std::list<Point_2>         inter_pts;
  bool                       app_ok = true;

  if (map_iter != inter_map.end())
  {
    // Get the intersection points between the two supporting curves as stored
    // in the map.
    inter_pts = map_iter->second;
  }
  else
  {
    // Approximate the intersection points and store them in the map.
    app_ok = _approximate_intersection_points (_curve,
                                               cv._curve,
                                               inter_pts);

    if (app_ok)
      inter_map[curve_pair] = inter_pts;
  }

  // Try to approximate the intersection points.
  bool                in_range1, in_range2;
  bool                correct_res;

  if (app_ok)
  {
    // If the approximation went OK, then we know that we just have simple
    // intersection points (with multiplicity 1). We go over the points
    // and report the ones lying in the parameter ranges of both curves.
    typename std::list<Point_2>::iterator  pit;
    
    for (pit = inter_pts.begin(); pit != inter_pts.end(); ++pit)
    {
      // Check if the point is in the range of this curve - first using
      // its parameter bounds, and if we fail we perform an exact check.
      in_range1 = _is_in_range (*pit, correct_res);
      
      if (! correct_res)
      {
        if (! pit->is_exact())
          pit->make_exact (cache);
        
        Originator_iterator  p_org = pit->originator (_curve);
        CGAL_assertion (p_org != pit->originators_end());
        
        in_range1 = _is_in_range (p_org->parameter(), cache);
      }
      
      if (! in_range1)
        continue;

      // Check if the point is in the range of the other curve - first using
      // its parameter bounds, and if we fail we perform an exact check.
      in_range2 = cv._is_in_range (*pit, correct_res);
      
      if (! correct_res)
      {
        if (! pit->is_exact())
          pit->make_exact (cache);
        
        Originator_iterator  p_org = pit->originator (cv._curve);
        CGAL_assertion (p_org != pit->originators_end());
        
        in_range2 = cv._is_in_range (p_org->parameter(), cache);
      }
      
      if (in_range1 && in_range2)
      {
        // The point lies within the parameter range of both curves, so we
        // report it as a valid intersection point with multiplicity 1.
        ipts.push_back (Intersection_point_2 (*pit, 1));
      }
    }

    // Since the apporximation went fine we cannot possible have an overlap:
    return (false);
  }

  // We did not succeed in isolate the approximate intersection points.
  // We therefore resort to the exact procedure and exactly compute them.
  bool                   do_ovlp;
  const Intersect_list&  inter_list =
    cache.intersections (_curve.id(),
                             _curve.x_polynomial(), _curve.x_norm(),
                             _curve.y_polynomial(), _curve.y_norm(),
                             cv._curve.id(),
                             cv._curve.x_polynomial(), cv._curve.x_norm(),
                             cv._curve.y_polynomial(), cv._curve.y_norm(),
                             do_ovlp);

  if (do_ovlp)
  {
    // Ron: TODO -- Handle overlapping curves ...
    CGAL_assertion (false);
    return (true);
  }

  // Go over the points and report the ones lying in the parameter ranges
  // of both curves.
  Intersect_iter          iit;

  for (iit = inter_list.begin(); iit != inter_list.end(); ++iit)
  {
    if (_is_in_range (iit->s, cache) &&
        cv._is_in_range (iit->t, cache))
    {
      // Construct an intersection point with unknown multiplicity.
      Point_2                 pt (iit->x, iit->y);
      
      pt.add_originator (Originator (_curve, iit->s));
      pt.add_originator (Originator (cv._curve, iit->t));
      
      ipts.push_back (Intersection_point_2 (pt, 0));
    }
  }

  // Mark that there is no overlap:
  return (false);
}

// ---------------------------------------------------------------------------
// Compute the exact vertical position of the point p with respect to the
// curve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_x_monotone_2<RatKer, AlgKer, NtTrt, BndTrt>::_exact_vertical_position
    (const Point_2& p,
     bool force_exact) const
{
  // If it is a rational point, obtain its rational reprsentation.
  Rat_point_2              rat_p;

  if (p.is_rational())
    rat_p = (Rat_point_2) p;

  // Get a rational approximation of the parameter values at the endpoints.
  Nt_traits                nt_traits;
  Originator_iterator      ps_org = _ps.originator(_curve);
  CGAL_assertion (ps_org != _ps.originators_end());
  
  Originator_iterator      pt_org = _pt.originator(_curve);
  CGAL_assertion (pt_org != _pt.originators_end());

  Rational                 my_t_min;
  Rational                 my_t_max;

  if (CGAL::compare (ps_org->point_bound().t_max,
                     pt_org->point_bound().t_min) == SMALLER)
  {
    // In case the parameter value of the source is smaller than the target's.
    my_t_min = ps_org->point_bound().t_max;
    my_t_max = pt_org->point_bound().t_min;
  }
  else
  {
    CGAL_assertion (CGAL::compare (pt_org->point_bound().t_max,
                                   ps_org->point_bound().t_min) == SMALLER);

    // In case the parameter value of the target is smaller than the source's.
    my_t_min = pt_org->point_bound().t_max;
    my_t_max = ps_org->point_bound().t_min;
  }

  // Start the subdivision process from the entire supporting curve.
  std::list<Subcurve>      subcurves;
  Subcurve                 init_scv;
  Rational                 x_min, y_min, x_max, y_max;
  bool                     no_x_ovlp;
  Comparison_result        res_y_min, res_y_max;

  std::copy (_curve.control_points_begin(), _curve.control_points_end(),
             std::back_inserter (init_scv.control_points));
  init_scv.t_min = 0;
  init_scv.t_max = 1;
  subcurves.push_back (init_scv);

  while (! subcurves.empty())
  {
    // Go over the list of subcurves and consider only those lying in the
    // given [t_min, t_max] bound.
    typename std::list<Subcurve>::iterator  iter = subcurves.begin();
    bool                                    is_fully_in_range;

    while (iter != subcurves.end())
    {
      if (CGAL::compare (iter->t_max, my_t_min) == SMALLER ||
          CGAL::compare (iter->t_min, my_t_max) == LARGER)
      {
        // Subcurve out of bounds of the x-monotone curve we consider - erase
        // it and continue to next subcurve.
        subcurves.erase(iter++);
        continue;
      }
        
      // Construct the bounding box of the subcurve and compare it to
      // the bounding box of the point.
      iter->bbox (x_min, y_min, x_max, y_max);

      if (p.is_rational())
      {
        no_x_ovlp = (CGAL::compare (x_min, rat_p.x()) == LARGER ||
                     CGAL::compare (x_max, rat_p.x()) == SMALLER);
      }
      else
      {
        no_x_ovlp = (CGAL::compare (nt_traits.convert (x_min),
                                    p.x()) == LARGER ||
                     CGAL::compare (nt_traits.convert (x_max),
                                    p.x()) == SMALLER);
      }

      if (no_x_ovlp)
      {
        // Subcurve out of x-bounds - erase it and continue to next subcurve.
        subcurves.erase(iter++);
        continue;
      }

      // In this case, check if there is an overlap in the y-range.
      if (p.is_rational())
      {
        res_y_min = CGAL::compare (rat_p.y(), y_min);
        res_y_max = CGAL::compare (rat_p.y(), y_max);
      }
      else
      {
        res_y_min = CGAL::compare (p.y(), nt_traits.convert (y_min));
        res_y_max = CGAL::compare (p.y(), nt_traits.convert (y_max));
      }
     
      is_fully_in_range = (CGAL::compare (iter->t_min, my_t_min) != SMALLER) &&
                          (CGAL::compare (iter->t_max, my_t_max) != LARGER);

      if (res_y_min != res_y_max || ! is_fully_in_range)
      {
        // Subdivide the current subcurve and replace iter with the two
        // resulting subcurves using de Casteljau's algorithm.
        Subcurve           scv_l, scv_r;

        scv_l.t_min = iter->t_min;
        scv_r.t_max = iter->t_max;
        scv_l.t_max = scv_r.t_min = (iter->t_min + iter->t_max) / 2;

        bisect_control_polygon_2 (iter->control_points.begin(),
                                  iter->control_points.end(),
                                  std::back_inserter(scv_l.control_points),
                                  std::front_inserter(scv_r.control_points));

        subcurves.insert (iter, scv_l);
        subcurves.insert (iter, scv_r);
        subcurves.erase(iter++);

        continue;
      }

      if (res_y_min == res_y_max)
      {
        CGAL_assertion (res_y_min != EQUAL);

        // We reached a separation, as p is either strictly above or strictly
        // below the bounding box of the current subcurve.
        return (res_y_min);
      }

      // If we got here without entering one of the clauses above,
      // then iter has not been incremented yet.
      ++iter; 
    }
  }

  // We can reach here only if we do not force an exact result.
  CGAL_assertion (! force_exact);
  return (EQUAL);
}

CGAL_END_NAMESPACE

#endif
