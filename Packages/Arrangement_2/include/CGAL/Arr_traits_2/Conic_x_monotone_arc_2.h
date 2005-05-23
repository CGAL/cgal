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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_X_MONOTONE_ARC_2_H
#define CGAL_CONIC_X_MONOTONE_ARC_2_H

/*! \file
 * Header file for the _Conic_x_monotone_arc_2<Conic_arc_2> class.
 */

#include <map>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*!
 * Representation of an x-monotone conic arc.
 * The class is templated by a representation of a general bounded conic arc.
 */

template <class Conic_arc_>
class _Conic_x_monotone_arc_2 : private Conic_arc_
{
public:

  typedef Conic_arc_                              Conic_arc_2;
  typedef _Conic_x_monotone_arc_2<Conic_arc_2>    Self;
  
  typedef typename Conic_arc_2::Alg_kernel        Alg_kernel;

  typedef typename Conic_arc_2::Alg_kernel        Alg_kernel;
  typedef typename Conic_arc_2::Algebraic         Algebraic;

  typedef typename Conic_arc_2::Point_2           Point_2;
  typedef typename Conic_arc_2::Conic_point_2     Conic_point_2;

  // Type definition for the intersection points mapping.
  typedef typename Conic_point_2::Conic_id        Conic_id;
  typedef std::pair<Conic_id, Conic_id>           Conic_pair;
  typedef std::pair<Conic_point_2, unsigned int>  Intersection_point_2;
  typedef std::list<Intersection_point_2>         Intersection_list;

  /*!
   * \struct Less functor for Conic_pair.
   */
  struct Less_conic_pair
  {
    bool operator() (const Conic_pair& cp1, const Conic_pair& cp2) const
    {
      // Compare the pairs of IDs lexicographically.
      return (cp1.first < cp2.first ||
              (cp1.first == cp2.first && cp1.second < cp2.second));
    }
  };

  typedef std::map<Conic_pair,
                   Intersection_list,
                   Less_conic_pair>               Intersection_map;
  typedef typename Intersection_map::value_type   Intersection_map_entry;
  typedef typename Intersection_map::iterator     Intersection_map_iterator;

protected:

  typedef Conic_arc_2                             Base;

  typedef typename Conic_arc_2::Integer           Integer;
  typedef typename Conic_arc_2::Nt_traits         Nt_traits;

  // Bit mask for the _info field (the least significant bit is already used
  // by the vase class).
  enum
  {
    IS_VERTICAL_SEGMENT = 2,
    IS_DIRECTED_RIGHT = 4,
    DEGREE_1 = 8,
    DEGREE_2 = 16,
    DEGREE_MASK = 8 + 16,
    PLUS_SQRT_DISC_ROOT = 32,
    FACING_UP = 64,
    FACING_DOWN = 128,
    FACING_MASK = 64 + 128
  };

  Algebraic      alg_r;      // The coefficients of the supporting conic curve:
  Algebraic      alg_s;      //
  Algebraic      alg_t;      //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 ,
  Algebraic      alg_u;      //
  Algebraic      alg_v;      // converted to algebraic numbers.
  Algebraic      alg_w;      //

  Conic_id      _id;         // The ID number of the supporting conic curve.

public:

  /// \name Constrcution methods.
  //@{

  /*!
   * Default constructor.
   */
  _Conic_x_monotone_arc_2 () :
    Base (),
    _id (0)
  {}

  /*!
   * Copy constructor.
   * \param arc The copied arc.
   */
  _Conic_x_monotone_arc_2 (const Self& arc) :
    Base (arc),
    _id (arc._id)
  {
    if (_id > 0)
      _set();
    else
      _info = 0;
  }

  /*!
   * Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \param id The ID of the base arc.
   */
  _Conic_x_monotone_arc_2 (const Base& arc,
                           const Conic_id& id) :
    Base (arc),
    _id (id)
  {
    _set ();
  }
  
  /*!
   * Construct an x-monotone sub-arc from a conic arc.
   * \param arc The given (base) arc.
   * \param source The source point.
   * \param target The target point.
   * \param id The ID of the base arc.
   */
  _Conic_x_monotone_arc_2 (const Base& arc,
                           const Point_2& source, const Point_2& target,
                           const Conic_id& id) :
    Base (arc),
    _id (id)
  {
    // Set the two endpoints.
    _source = source;
    _target = target;
    
    _set();
  }

  /*!
   * Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator= (const Self& arc)
  {
    if (this == &arc)
      return (*this);

    // Copy the base arc.
    Base::operator= (arc);

    // Set the rest of the properties.
    _id = arc._id;

    if (_id > 0)
      _set();
    else
      _info = 0;

    return (*this);
  }
  //@}

  /// \name Accessing the arc properties.
  //@{

  /*! 
   * Get the coefficients of the underlying conic.
   */
  const Integer& r () const {return (_r);}
  const Integer& s () const {return (_s);}
  const Integer& t () const {return (_t);}
  const Integer& u () const {return (_u);}
  const Integer& v () const {return (_v);}
  const Integer& w () const {return (_w);}

  /*!
   * Get the arc's source.
   * \return The source point.
   */
  const Conic_point_2& source () const
  {
    return (_source);
  }

  /*!
   * Get the arc's target.
   * \return The target point.
   */
  const Conic_point_2& target () const
  {
    return (_target);
  }

  /*!
   * Get the orientation of the arc.
   * \return The orientation.
   */
  Orientation orientation () const
  {
    return (_orient);
  }

  /*!
   * Get the left endpoint of the arc.
   */
  const Conic_point_2& left () const
  {
    if ((_info & IS_DIRECTED_RIGHT) != 0)
      return (_source);
    else
      return (_target);
  }

  /*!
   * Get the right endpoint of the arc.
   */
  const Conic_point_2& right () const
  {
    if ((_info & IS_DIRECTED_RIGHT) != 0)
      return (_target);
    else
      return (_source);
  }

  /*!
   * Get a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox () const
  {
    return (Base::bbox());
  }
  //@}

  /// \name Predicates.
  //@{

  /*!
   * Check if the conic arc is a vertical segment.
   */
  bool is_vertical () const
  {
    return ((_info & IS_VERTICAL_SEGMENT) != 0);
  }

  /*!
   * Check whether the given point lies on the arc.
   * \param p The qury point.
   * \param (true) if p lies on the arc; (false) otherwise.
   */
  bool contains_point (const Conic_point_2& p) const
  {
    // First check if p lies on the supporting conic. We first check whether
    // it is one of p's generating conic curves.
    bool       p_on_conic = false;

    if (p.is_generating_conic (_id))
    {
      p_on_conic = true;
    }
    else
    {
      // Check whether p satisfies the supporting conic equation.
      p_on_conic = _is_on_supporting_conic (p.x(), p.y());

      if (p_on_conic)
      {
        // As p lies on the supporting conic of our arc, add its ID to
        // the list of generating conics for p.
        Conic_point_2&  p_non_const = const_cast<Conic_point_2&> (p);
        p_non_const.set_generating_conic (_id);
      }
    }

    if (! p_on_conic)
      return (false);

    // Check if p is between the endpoints of the arc.
    return (_is_between_endpoints (p));
  }
  //@}

  /// \name Constructing points on the arc.
  //@{

  /*!
   * Compute a point on the arc with the same x-coordiante as the given point.
   * \param p The given point.
   * \pre The arc is not vertical and p is in the x-range of the arc.
   * \return A point on the arc with the same x-coordiante as p.
   */
  Point_2 get_point_at_x (const Point_2& p) const
  {
    // Make sure that p is in the x-range of the arc.
    CGAL_precondition ((_info & IS_VERTICAL_SEGMENT) == 0);

    CGAL_precondition_code (
      Alg_kernel   ker;
    );
    CGAL_precondition (ker.compare_x_2_object() (p, left()) != SMALLER &&
                       ker.compare_x_2_object() (p, right()) != LARGER);

    // Compute the y-coordinate according to the degree of the supporting
    // conic curve.
    Nt_traits        nt_traits;
    Algebraic        y;

    if ((_info & DEGREE_MASK) == DEGREE_1)
    {
      // In case of a linear curve, the y-coordinate is a simple linear
      // expression of x(p) (note that v is not 0 as the arc is not vertical):
      //   y = -(u*x(p) + w) / v
      y = -(alg_u*p.x() + alg_w) / alg_v;
    }
    else
    {
      CGAL_assertion ((_info & DEGREE_MASK) == DEGREE_2);

      // In this case the y-coordinate is one of solutions to the quadratic
      // equation:
      //  s*y^2 + (t*x(p) + v)*y + (r*x(p)^2 + u*x(p) + w) = 0
      Algebraic  A = alg_s;
      Algebraic  B = alg_t*p.x() + alg_v;
      Algebraic  C = (alg_r*p.x() + alg_u)*p.x() + alg_w;

      if (CGAL::sign(_s) == ZERO)
      {
	// In this case A is 0 and we have a linear equation.
	CGAL_assertion (CGAL::sign (B) != ZERO);

	y = -C / B;
      }
      else
      {
	// Solve the quadratic equation.
	Algebraic  disc = B*B - 4*A*C;
    
	CGAL_assertion (CGAL::sign (disc) != NEGATIVE);

	// We take either the root involving -sqrt(disc) or +sqrt(disc)
	// based on the information flags.
	if ((_info & PLUS_SQRT_DISC_ROOT) != 0)
	{
	  y = (nt_traits.sqrt (disc) - B) / (2*A);
	}
	else
	{
	  y = -(B + nt_traits.sqrt (disc)) / (2*A);
	}
      }
    }

    // Return the computed point.
    return (Point_2 (p.x(), y));
  }

  /*!
   * Compare to arcs immediately to the right of their intersection point.
   * \param arc The compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the right of p.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_right (const Self& arc,
                                      const Conic_point_2& p) const
  {
    CGAL_precondition ((_info & IS_VERTICAL_SEGMENT) == 0 &&
                       (arc._info & IS_VERTICAL_SEGMENT) == 0);

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (_has_same_supporting_conic (arc))
    {
      if ((_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
        return (LARGER);
      else if ((_info & FACING_DOWN) != 0 && (arc._info & FACING_UP) != 0)
        return (SMALLER);

      // In this case the two arcs overlap.
      CGAL_assertion ((_info & FACING_MASK) == (arc._info & FACING_MASK));

      return (EQUAL);
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic      slope1_numer, slope1_denom;
    Algebraic      slope2_numer, slope2_denom;

    _derive_by_x_at (p, 1, slope1_numer, slope1_denom);
    arc._derive_by_x_at (p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool     is_vertical_slope1 = (CGAL::sign (slope1_denom) == ZERO);
    const bool     is_vertical_slope2 = (CGAL::sign (slope2_denom) == ZERO);

    if (!is_vertical_slope1 && !is_vertical_slope2)
    {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res = CGAL::compare (slope1_numer*slope2_denom,
                                                   slope2_numer*slope1_denom);
      
      if (slope_res != EQUAL)
        return (slope_res);

      // Use the second-order derivative.
      _derive_by_x_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_x_at (p, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare (slope1_numer*slope2_denom, 
                                 slope2_numer*slope1_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion (slope_res != EQUAL);

      return (slope_res);
    }
    else if (!is_vertical_slope2)
    {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((_info & FACING_MASK) != 0);

      if ((_info & FACING_UP) != 0)
        return (LARGER);
      else
        return (SMALLER);
    }
    else if (!is_vertical_slope1)
    {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((arc._info & FACING_MASK) != 0);

      if ((arc._info & FACING_UP) != 0)
        return (SMALLER);
      else
        return (LARGER);
    }
    else
    {
      // The two arcs have vertical slopes at p_int: 
      // First check whether one is facing up and one down. In this case the
      // comparison result is trivial.
      if ((_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
        return (LARGER);
      else if ((_info & FACING_DOWN) != 0 && (arc._info & FACING_UP) != 0)
        return (SMALLER);

      // Compute the second-order derivative by y and act according to it.
      _derive_by_y_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_y_at (p, 2, slope2_numer, slope2_denom);

      Comparison_result slope_res = CGAL::compare (slope1_numer*slope2_denom, 
                                                   slope2_numer*slope1_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      if ((_info & FACING_UP) != 0 && (arc._info & FACING_UP) != 0)
      {
        // Both are facing up.
        return ((slope_res == LARGER) ? SMALLER : LARGER);
      }
      else
      {
        // Both are facing down.
        return (slope_res);
      }
    }

    // We should never reach here:
    CGAL_assertion(false);
    return (EQUAL);
  }

  /*!
   * Compare to arcs immediately to the leftt of their intersection point.
   * \param arc The compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the left of p.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_left (const Self& arc,
				     const Conic_point_2& p) const
  {
    CGAL_precondition ((_info & IS_VERTICAL_SEGMENT) == 0 &&
                       (arc._info & IS_VERTICAL_SEGMENT) == 0);

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (_has_same_supporting_conic (arc))
    {
      if ((_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
        return (LARGER);
      else if ((_info & FACING_DOWN) != 0 && (arc._info & FACING_UP) != 0)
        return (SMALLER);

      // In this case the two arcs overlap.
      CGAL_assertion ((_info & FACING_MASK) == (arc._info & FACING_MASK));

      return (EQUAL);
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic      slope1_numer, slope1_denom;
    Algebraic      slope2_numer, slope2_denom;
    
    _derive_by_x_at (p, 1, slope1_numer, slope1_denom);
    arc._derive_by_x_at (p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool     is_vertical_slope1 = (CGAL::sign (slope1_denom) == ZERO);
    const bool     is_vertical_slope2 = (CGAL::sign (slope2_denom) == ZERO);
    
    if (!is_vertical_slope1 && !is_vertical_slope2)
    {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result  slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                                   slope1_numer*slope2_denom);
      
      if (slope_res != EQUAL)
        return (slope_res);

      // Use the second-order derivative.
      _derive_by_x_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_x_at (p, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare (slope1_numer*slope2_denom, 
                                 slope2_numer*slope1_denom);
      
      // \todo Handle higher-order derivatives:
      CGAL_assertion (slope_res != EQUAL);

      return (slope_res);
    }
    else if (!is_vertical_slope2)
    {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((_info & FACING_MASK) != 0);

      if ((_info & FACING_UP) != 0)
        return (LARGER);
      else
        return (SMALLER);
    }
    else if (!is_vertical_slope1)
    {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((arc._info & FACING_MASK) != 0);

      if ((arc._info & FACING_UP) != 0)
        return (SMALLER);
      else
        return (LARGER);
    }
    else
    {
      // The two arcs have vertical slopes at p_int: 
      // First check whether one is facing up and one down. In this case the
      // comparison result is trivial.
      if ((_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
        return (LARGER);
      else if ((_info & FACING_DOWN) != 0 && (arc._info & FACING_UP) != 0)
        return (SMALLER);

      // Compute the second-order derivative by y and act according to it.
      _derive_by_y_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_y_at (p, 2, slope2_numer, slope2_denom);

      Comparison_result  slope_res = CGAL::compare(slope2_numer*slope1_denom, 
                                                   slope1_numer*slope2_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      if ((_info & FACING_UP) != 0 && (arc._info & FACING_UP) != 0)
      {
        // Both are facing up.
        return ((slope_res == LARGER) ? SMALLER : LARGER);
      }
      else
      {
        // Both are facing down.
        return (slope_res);
      }
    }

    // We should never reach here:
    CGAL_assertion(false);
    return (EQUAL);
  }

  /*!
   * Compute the intersections between with the given arc.
   * \param arc The given intersecting arc.
   * \param inter_map Maps conic pairs to lists of their intersection points.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template<class OutputIterator>
  OutputIterator intersect (const Self& arc,
                            Intersection_map& inter_map,
                            OutputIterator oi) const
  {
    if (_has_same_supporting_conic (arc))
    {
      // Check for overlaps between the two arcs.
      Self    overlap;

      if (_compute_overlap (arc, overlap))
      {
	// There can be just a single overlap between two x-monotone arcs:
	*oi = make_object (overlap);
	oi++;
	return (oi);
      }

      // In case there is not overlap and the supporting conics are the same,
      // there cannot be any intersection points, unless the two arcs share
      // an end point.
      // Note that in this case we do not define the multiplicity of the
      // intersection points we report.
      Alg_kernel  ker;

      if (ker.equal_2_object() (left(), arc.left()))
      {
	Intersection_point_2  ip (left(), 0);
	
	*oi = make_object (ip);
	oi++;
      }
      
      if (ker.equal_2_object() (right(), arc.right()))
      {
	Intersection_point_2  ip (right(), 0);
	
	*oi = make_object (ip);
	oi++;
      }

      return (oi);
    }
    
    // Search for the pair of supporting conics in the map (the first conic
    // ID in the pair should be smaller than the second one, to guarantee
    // uniqueness).
    Conic_pair                   conic_pair;
    Intersection_map_iterator    map_iter;
    Intersection_list            inter_list;

    if (_id < arc._id)
      conic_pair = Conic_pair (_id, arc._id);
    else
      conic_pair = Conic_pair (arc._id, _id);

    map_iter = inter_map.find (conic_pair);

    if (map_iter == inter_map.end())
    {
      // In case the intersection points between the supporting conics have
      // not been computed before, compute them now and store them in the map.
      _intersect_supporting_conics (arc, inter_list);
      inter_map[conic_pair] = inter_list;
    }
    else
    {
      // Obtain the precomputed intersection points from the map.
      inter_list = (*map_iter).second;
    }

    // Go over the list of intersection points and report those that lies on
    // both x-monotone arcs.
    typename Intersection_list::const_iterator  iter;

    for (iter = inter_list.begin(); iter != inter_list.end(); ++iter)
    {
      if (_is_between_endpoints ((*iter).first) &&
          arc._is_between_endpoints ((*iter).first))
      {
        *oi = make_object (*iter);
        ++oi;
      }
    }

    return (oi);
  }
  //@}

  /// \name Constructing x-monotone arcs.
  //@{

  /*!
   * Split the arc into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the arc (not one of its endpoints).
   */
  void split (const Conic_point_2& p,
              Self& c1, Self& c2) const
  {
    // Make sure that p lies on the interior of the arc.
    CGAL_precondition_code (
      Alg_kernel   ker;
    );
    CGAL_precondition (this->contains_point (p) &&
                       ! ker.equal_2_object() (p, _source) &&
                       ! ker.equal_2_object() (p, _target));

    // Make copies of the current arc.
    c1 = *this;
    c2 = *this;

    // Assign the endpoints of the arc.
    if ((_info & IS_DIRECTED_RIGHT) != 0)
    {
      // The arc is directed from left to right, so p becomes c1's target
      // and c2's source.
      c1._target = p;
      c2._source = p;

      if (! p.is_generating_conic (_id))
      {
        c1._target.set_generating_conic (_id);
        c2._source.set_generating_conic (_id);
      }
    }
    else
    {
      // The arc is directed from right to left, so p becomes c2's target
      // and c1's source.
      c1._source = p;
      c2._target = p;

      if (! p.is_generating_conic (_id))
      {
        c1._source.set_generating_conic (_id);
        c2._target.set_generating_conic (_id);
      }
    }

    return;
  }

  /*!
   * Check whether the two arcs are equal (have the same graph).
   * \param arc The compared arc.
   * \return (true) if the two arcs have the same graph; (false) otherwise.
   */
  bool equals (const Self& arc) const
  {
    // The two arc must have the same supporting conic curves.
    if (! _has_same_supporting_conic (arc))
      return (false);

    // Check that the arc endpoints are the same.
    Alg_kernel   ker;

    if (_orient == arc._orient)
    {
      // Same orientation - the source and target points must be the same.
      return (ker.equal_2_object() (_source, arc._source) &&
              ker.equal_2_object() (_target, arc._target));
    }
    else
    {
      // Reverse orientation - the source and target points must be swapped.
      return (ker.equal_2_object() (_source, arc._target) &&
              ker.equal_2_object() (_target, arc._source));
    }
  }

  /*!
   * Check whether it is possible to merge the arc with the given arc.
   * \param arc The query arc.
   * \return (true) if it is possible to merge the two arcs;
   *         (false) otherwise.
   */
  bool can_merge_with (const Self& arc) const
  {
    // In order to merge the two arcs, they should have the same supporting
    // conic.
    if (! _has_same_supporting_conic (arc))
      return (false);

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    Alg_kernel   ker;

    return (ker.equal_2_object() (right(), arc.left()) ||
            ker.equal_2_object() (left(), arc.right()));
  }

  /*!
   * Merge the current arc with the given arc.
   * \param arc The arc to merge with.
   * \pre The two arcs are mergeable.
   */
  void merge (const Self& arc)
  {
    CGAL_precondition (this->can_merge_with (arc));

    // Check if we should extend the arc to the left or to the right.
    Alg_kernel   ker;

    if (ker.equal_2_object() (right(), arc.left()))
    {
      // Extend the arc to the right.
      if ((_info & IS_DIRECTED_RIGHT) != 0)
	_target = arc.right();
      else
	_source = arc.right();
    }
    else
    {
      CGAL_precondition (ker.equal_2_object() (left(), arc.right()));
      
      // Extend the arc to the left.
      if ((_info & IS_DIRECTED_RIGHT) != 0)
	_source = arc.left();
      else
	_target = arc.left();
    }

    return;
  }
  //@}

private:

  /// \name Auxiliary (private) functions.
  //@{

  /*!
   * Set the properties of the x-monotone conic arc (for the usage of the 
   * constructors).
   */
  void _set ()
  {
    // Convert the coefficients of the supporting conic to algebraic numbers.
    Nt_traits        nt_traits;

    alg_r = nt_traits.convert (_r);
    alg_s = nt_traits.convert (_s);
    alg_t = nt_traits.convert (_t);
    alg_u = nt_traits.convert (_u);
    alg_v = nt_traits.convert (_v);
    alg_w = nt_traits.convert (_w);

    // Set the generating conic ID for the source and target points.
    _source.set_generating_conic (_id);
    _target.set_generating_conic (_id);

    // Clear the _info bits.
    _info = 0;

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel         ker;
    Comparison_result  dir_res = ker.compare_xy_2_object() (_source, _target);

    CGAL_assertion (dir_res != EQUAL);

    if (dir_res == SMALLER)
      _info = (_info | IS_DIRECTED_RIGHT);

    // Compute the degree of the underlying conic.
    if (CGAL::sign (_r) != ZERO ||
        CGAL::sign (_s) != ZERO ||
        CGAL::sign (_t) != ZERO)
    {
      if (_orient == COLLINEAR && 
          ker.compare_x_2_object() (_source, _target) == EQUAL)
      {
        // The arc is a vertical segment:
        _info = (_info | IS_VERTICAL_SEGMENT);
      }

      _info = (_info | DEGREE_2);
    }
    else
    {
      CGAL_assertion (CGAL::sign (_u) != ZERO ||
                      CGAL::sign (_v) != ZERO);

      if (CGAL::sign (_v) == ZERO)
      {
        // The supporting curve is of the form: _u*x + _w = 0
        _info = (_info | IS_VERTICAL_SEGMENT);
      }

      _info = (_info | DEGREE_1);

      return;
    }

    // Compute a midpoint between the source and the target and get the y-value
    // of the arc at its x-coordiante.
    Point_2          p_mid = ker.construct_midpoint_2_object() (_source,
                                                                _target);
    Algebraic        ys[2];
    int              n_ys;

    n_ys = _conic_get_y_coordinates (p_mid.x(), ys);

    CGAL_assertion (n_ys != 0);

    // Check which solution lies on the x-monotone arc.
    Point_2          p_arc_mid (p_mid.x(), ys[0]);

    if (_is_strictly_between_endpoints (p_arc_mid))
    {
      // Mark that we should use the -sqrt(disc) root for points on this
      // x-monotone arc.
      _info = (_info & ~PLUS_SQRT_DISC_ROOT);
    }
    else
    {
      CGAL_assertion (n_ys == 2);
      p_arc_mid = Point_2 (p_mid.x(), ys[1]);
      CGAL_assertion (_is_strictly_between_endpoints (p_arc_mid));

      // Mark that we should use the +sqrt(disc) root for points on this
      // x-monotone arc.
      _info = (_info | PLUS_SQRT_DISC_ROOT);
    }

    if (_orient == COLLINEAR)
      return;

    // Check whether the conic is facing up or facing down:
    // Check whether the arc (which is x-monotone of degree 2) lies above or 
    // below the segement that contects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y coordinate of a point on the arc whose x
    // coordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    Comparison_result res = ker.compare_y_2_object() (p_arc_mid, p_mid);

    if (res == LARGER)
    {
      // The arc is above the connecting segment, so it is facing upwards.
      _info = (_info | FACING_UP);
    }
    else if (res == SMALLER)
    {
      // The arc is below the connecting segment, so it is facing downwards.
      _info = _info | FACING_DOWN;
    }
    
    return;
  }

  /*!
   * Check whether the given point lies on the supporting conic of the arc.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return (true) if p lies on the supporting conic; (false) otherwise.
   */
  bool _is_on_supporting_conic (const Algebraic& px,
                                const Algebraic& py) const
  {
    // Check whether p satisfies the conic equation.
    // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
    const Algebraic  val = (alg_r*px + alg_t*py + alg_u) * px +
                           (alg_s*py + alg_v) * py +
                           alg_w;

    return (CGAL::sign (val) == ZERO);
  }

  /*!
   * Check whether the two arcs have the same supporting conic.
   * \param arc The compared arc.
   * \return (true) if the two supporting conics are the same.
   */
  bool _has_same_supporting_conic (const Self& arc) const
  {
    // Check if the two arcs originate from the same conic:
    if (_id == arc._id)
      return (true);
    
    // Check whether the coefficients of the two supporting conics are equal
    // up to a constant factor.
    Integer        factor1 = 1;
    Integer        factor2 = 1;

    if (CGAL::sign (_r) != ZERO)
      factor1 = _r;
    else if (CGAL::sign (_s) != ZERO)
      factor1 = _s;
    else if (CGAL::sign (_t) != ZERO)
      factor1 = _t;
    else if (CGAL::sign (_u) != ZERO)
      factor1 = _u;
    else if (CGAL::sign (_v) != ZERO)
      factor1 = _v;
    else if (CGAL::sign (_w) != ZERO)
      factor1 = _w;

    if (CGAL::sign (arc._r) != ZERO)
      factor2 = arc._r;
    else if (CGAL::sign (arc._s) != ZERO)
      factor2 = arc._s;
    else if (CGAL::sign (arc._t) != ZERO)
      factor2 = arc._t;
    else if (CGAL::sign (arc._u) != ZERO)
      factor2 = arc._u;
    else if (CGAL::sign (arc._v) != ZERO)
      factor2 = arc._v;
    else if (CGAL::sign (arc._w) != ZERO)
      factor2 = arc._w;

    return (CGAL::compare  (_r * factor2, arc._r * factor1) == EQUAL && 
            CGAL::compare  (_s * factor2, arc._s * factor1) == EQUAL &&
            CGAL::compare  (_t * factor2, arc._t * factor1) == EQUAL &&
            CGAL::compare  (_u * factor2, arc._u * factor1) == EQUAL &&
            CGAL::compare  (_v * factor2, arc._v * factor1) == EQUAL &&
            CGAL::compare  (_w * factor2, arc._w * factor1) == EQUAL);
  }

  /*!
   * Get the i'th order derivative by x of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1 or 2).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void _derive_by_x_at (const Point_2& p, const unsigned int& i,
			Algebraic& slope_numer, Algebraic& slope_denom) const
  {
    // The derivative by x of the conic 
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2r*x + t*y + u       alpha 
    //   y' = - ---------------- = - -------
    //           2s*y + t*x + v       beta
    //
    const Algebraic  _two = 2;
    const Algebraic  sl_numer = _two*alg_r*p.x() + alg_t*p.y() + alg_u;
    const Algebraic  sl_denom = _two*alg_s*p.y() + alg_t*p.x() + alg_v;

    if (i == 1)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else
      {
        slope_numer = sl_numer;
        slope_denom = -sl_denom;
      }

      return;
    }

    // The second derivative is given by:
    //
    //             s*alpha^2 - t*alpha*beta + r*beta^2
    //   y'' = -2 -------------------------------------
    //                           beta^3
    //
    const Algebraic  sl2_numer = alg_s * sl_numer*sl_numer -
                                 alg_t * sl_numer*sl_denom +
                                 alg_r * sl_denom*sl_denom;
    const Algebraic  sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -_two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else
      {
        slope_numer = _two *sl2_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_assertion (false);
    return;
  }

  /*!
   * Get the i'th order derivative by y of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1 or 2).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void _derive_by_y_at (const Point_2& p, const int& i,
			Algebraic& slope_numer, Algebraic& slope_denom) const
  {
    // The derivative by y of the conic 
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2s*y + t*x + v     alpha 
    //   x' = - ---------------- = -------
    //           2r*x + t*y + u      beta
    //
    const Algebraic  _two = 2;
    const Algebraic  sl_numer = _two*alg_s*p.y() + alg_t*p.x() + alg_v;
    const Algebraic  sl_denom = _two*alg_r*p.x() + alg_t*p.y() + alg_u;

    if (i == 1)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else
      {
        slope_numer = sl_numer;
        slope_denom = -sl_denom;
      }

      return;
    }

    // The second derivative is given by:
    //
    //             r*alpha^2 - t*alpha*beta + s*beta^2
    //   x'' = -2 -------------------------------------
    //                           beta^3
    //
    const Algebraic  sl2_numer = alg_r * sl_numer*sl_numer -
                                 alg_t * sl_numer*sl_denom +
                                 alg_s * sl_denom*sl_denom;
    const Algebraic  sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -_two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else
      {
        slope_numer = _two *sl2_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_assertion (false);
    return;
  }

  /*!
   * Compute the overlap with a given arc, which is supposed to have the same
   * supporting conic curve as this arc.
   * \param arc The given arc.
   * \param overlap Output: The overlapping arc (if any).
   * \return Whether we found an overlap.
   */
  bool _compute_overlap (const Self& arc, Self& overlap) const
  {
    // Check if the two arcs are identical.
    if (equals (arc))
    {
      overlap = arc;
      return (true);
    }

    if (_is_strictly_between_endpoints (arc.left()))
    {
      if (_is_strictly_between_endpoints(arc.right()))
      {
        // Case 1 - *this:     +----------->     
        //            arc:       +=====>
        overlap = arc;
        return (true);
      }
      else
      {
        // Case 2 - *this:     +----------->     
        //            arc:               +=====>
	overlap = *this;

	if ((overlap._info & IS_DIRECTED_RIGHT) != 0)
	  overlap._source = arc.left();
	else
	  overlap._target = arc.left();
	  
        return (true);
      }
    }
    else if (_is_strictly_between_endpoints (arc.right()))
    {
      // Case 3 - *this:     +----------->     
      //            arc:   +=====>
      overlap = *this;

      if ((overlap._info & IS_DIRECTED_RIGHT) != 0)
	overlap._target = arc.right();
      else
	overlap._source = arc.right();

      return (true);
    }
    else if (arc._is_between_endpoints (_source) &&
             arc._is_between_endpoints (_target) &&
             (arc._is_strictly_between_endpoints(_source) ||
              arc._is_strictly_between_endpoints(_target)))
    {
      // Case 4 - *this:     +----------->     
      //            arc:   +================>
      overlap = *this;
      return (true);
    }
    
    // If we reached here, there are no overlaps:
    return (false);
  }

  /*!
   * Intersect the supporing conic curves of this arc and the given arc.
   * \param arc The arc to intersect with.
   * \param inter_list The list of intersection points.
   */
  void _intersect_supporting_conics (const Self& arc,
                                     Intersection_list& inter_list) const
  {
    // Compute the x-coordinates of the intersection points.
    Algebraic   xs[4];
    int         n_xs;

    n_xs = _compute_resultant_roots (_r, _s, _t, _u, _v, _w,
				     _info,
                                     arc._r, arc._s, arc._t,
				     arc._u, arc._v, arc._w,
				     arc._info,
                                     xs);

    // Compute the y-coordinates of the intersection points.
    Algebraic   ys[4];
    int         n_ys;

    n_ys = _compute_resultant_roots (_s, _r, _t, _v, _u, _w,
				     _info,
                                     arc._s, arc._r, arc._t,
                                     arc._v, arc._u, arc._w,
				     arc._info,
                                     ys);

    // Pair the coordinates of the intersection points. As the vectors of
    // x and y-coordinates are sorted in ascending order, we output the
    // intersection points in lexicographically ascending order.
    unsigned int  mult;
    int           i, j;

    for (i = 0; i < n_xs; i++)
    {
      for (j = 0; j < n_ys; j++)
      {
        if (_is_on_supporting_conic (xs[i], ys[j]) &&
            arc._is_on_supporting_conic (xs[i], ys[j]))
        {
          // Create the intersection point and set its generating conics.
          Conic_point_2         ip (xs[i], ys[j]);

          ip.set_generating_conic (_id);
          ip.set_generating_conic (arc._id);

          // Compute the multiplicity of the intersection point.
	  mult = _multiplicity_of_intersection_point (arc, ip);

          // Insert the intersection point to the output list.
          inter_list.push_back (Intersection_point_2 (ip, mult));
        }
      }
    }

    return;
  }

  /*!
   * Compute the roots of the resultants of the two bivariate polynomials:
   *   C1:  r1*x^2 + s1*y^2 + t1*xy + u1*x + v1*y + w1 = 0
   *   C2:  r2*x^2 + s2*y^2 + t2*xy + u2*x + v2*y + w2 = 0
   * \param info1 The information flags of the first curve.
   * \param info2 The information flags of the second curve.
   * \param xs Output: The real-valued roots of the polynomial, sorted in an
   *                   ascending order.
   * \pre xs must be a vector of size 4.
   * \return The number of distinct roots found.
   */
  int _compute_resultant_roots (const Integer& r1, const Integer& s1,
				const Integer& t1, const Integer& u1,
				const Integer& v1, const Integer& w1,
				const int& info1,
				const Integer& r2, const Integer& s2,
				const Integer& t2, const Integer& u2,
				const Integer& v2, const Integer& w2,
				const int& info2,
				Algebraic *xs) const
  {
    if ((info1 & DEGREE_MASK) == DEGREE_2 &&
        (info2 & DEGREE_MASK) == DEGREE_1)
    {
      // If necessary, swap roles between the two curves, so that the first
      // curve always has the minimal degree.
      return (_compute_resultant_roots (r2, s2, t2, u2, v2, w2, 
					info2,
					r1, s1, t1, u1, v1, w1, 
					info1,
					xs));
    }

    // Act according to the degree of the first conic curve.
    Nt_traits      nt_traits;
    const Integer  _two = 2;
    Integer        c[5];
    unsigned int   degree = 4;
    Algebraic     *xs_end;

    if ((info1 & DEGREE_MASK) == DEGREE_1)
    {
      // The first curve has no quadratic coefficients, and represents a line.
      if (CGAL::sign (v1) == ZERO)
      {
        // The first line is u1*x + w1 = 0, therefore:
        xs[0] = nt_traits.convert(-w1) / nt_traits.convert(u1);
        return (1);
      }
      
      // We can write the first curve as: y = (u1*x + w1) / v1.
      if ((info2 & DEGREE_MASK) == DEGREE_1)
      {
        // The second curve is also a line. We therefore get the linear
        // equation c[1]*x + c[0] = 0:
        c[1] = v1*u2 - u1*v2;
        c[0] = v1*w2 - w1*v2;

        if (CGAL::sign (c[1]) == ZERO)
	  // The two lines are parallel:
          return (0);

        xs[0] =  nt_traits.convert(-c[0]) /  nt_traits.convert(c[1]);
        return (1);
      }

      // We substitute this expression into the equation of the second
      // conic, and get the quadratic equation c[2]*x^2 + c[1]*x + c[0] = 0:
      c[2] = u1*u1*s2 - u1*v1*t2 + v1*v1*r2;
      c[1] = _two*u1*w1*s2 - u1*v1*v2 - v1*w1*t2 + v1*v1*u2;
      c[0] = w1*w1*s2 - v1*w1*v2 + v1*v1*w2;

      xs_end = nt_traits.solve_quadratic_equation (c[2], c[1], c[0],
						   xs);
      return (xs_end - xs);
    }

    // At this stage, both curves have degree 2. We obtain a qaurtic polynomial
    // whose roots are the x-coordinates of the intersection points.
    if (CGAL::sign (s1) == ZERO && CGAL::sign (s2) == ZERO)
    {
      // If both s1 and s2 are zero, we can write the two curves as:
      //   C1: (t1*x + v1)*y + (r1*x^2 + u1*x + w1) = 0
      //   C2: (t2*x + v2)*y + (r2*x^2 + u2*x + w2) = 0
      // By writing the resultant of these two polynomials we get a cubic
      // polynomial, whose coefficients are given by:
      c[3] = r2*t1 - r1*t2;
      c[2] = t1*u2 - t2*u1 + r2*v1 - r1*v2;
      c[1] = t1*w2 - t2*w1 + u2*v1 - u1*v2;
      c[0] = v1*w2 - v2*w1;

      degree = 3;
    }
    else
    {
      // We can write the two curves as:
      //   C1: (s1)*y^2 + (t1*x + v1)*y + (r1*x^2 + u1*x + w1) = 0
      //   C2: (s2)*y^2 + (t2*x + v2)*y + (r2*x^2 + u2*x + w2) = 0
      // By writing the resultant of these two polynomials we get a quartic
      // polynomial, whose coefficients are given by:
      c[4] = -_two*s1*s2*r1*r2 + s1*t2*t2*r1 - s1*t2*t1*r2 +
        s1*s1*r2*r2 - s2*t1*r1*t2 + s2*t1*t1*r2 + s2*s2*r1*r1;

      c[3] = -t2*r1*v1*s2 - u2*t1*t2*s1 - v2*r1*t1*s2 -
        r2*t1*v2*s1 - _two*s1*s2*r1*u2 - t2*u1*t1*s2 + u2*t1*t1*s2 -
        r2*v1*t2*s1 + u1*t2*t2*s1 + _two*v2*r1*t2*s1 + _two*u2*r2*s1*s1 + 
        _two*r2*v1*t1*s2 + _two*u1*r1*s2*s2 - _two*s1*s2*u1*r2;

      c[2] = -r2*v1*v2*s1 + u2*u2*s1*s1 + _two*w2*r2*s1*s1 +
        _two*u2*v1*t1*s2 - u2*v1*t2*s1 + w2*t1*t1*s2 - _two*s1*s2*u1*u2 - 
        w2*t1*t2*s1 + v2*v2*r1*s1 + u1*u1*s2*s2 - v2*r1*v1*s2 +
        _two*w1*r1*s2*s2 - u2*t1*v2*s1 - t2*u1*v1*s2 - _two*s1*s2*r1*w2 -
        _two*s1*s2*w1*r2 + r2*v1*v1*s2 + w1*t2*t2*s1 - v2*u1*t1*s2 -
        t2*w1*t1*s2 + _two*v2*u1*t2*s1;

      c[1] = _two*w2*u2*s1*s1 + _two*w2*v1*t1*s2 - w2*v1*t2*s1 +
        _two*v2*w1*t2*s1 + _two*w1*u1*s2*s2 - v2*u1*v1*s2 - _two*s1*s2*u1*w2 -
        v2*w1*t1*s2 + u2*v1*v1*s2 - t2*w1*v1*s2 - w2*t1*v2*s1 + 
        v2*v2*u1*s1 - u2*v1*v2*s1 - _two*s1*s2*w1*u2;

      c[0] = s2*v1*v1*w2 - s1*v2*v1*w2 - s2*v1*w1*v2 + s2*s2*w1*w1 -
        _two*s1*s2*w1*w2 + s1*w1*v2*v2 + s1*s1*w2*w2;

      degree = 4;
    }

    // Compute the roots of the resultant polynomial.
    xs_end = nt_traits.compute_polynomial_roots (c, degree,
						 xs);
    return (xs_end - xs);
  }

  /*!
   * Compute the multiplicity of an intersection point.
   * \param arc The arc to intersect with.
   * \param p The intersection point.
   * \return The multiplicity of the intersection point.
   */
  unsigned int _multiplicity_of_intersection_point (const Self& arc,
						    const Point_2& p) const
  {
    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic      slope1_numer, slope1_denom;
    Algebraic      slope2_numer, slope2_denom;
    
    _derive_by_x_at (p, 1, slope1_numer, slope1_denom);
    arc._derive_by_x_at (p, 1, slope2_numer, slope2_denom);

    if (CGAL::compare (slope1_numer*slope2_denom,
		       slope2_numer*slope1_denom) != EQUAL)
    {
      // Different slopes at p - the mutiplicity of p is 1:
      return (1);
    }
    
    if (CGAL::sign (slope1_denom) != ZERO &&
	CGAL::sign (slope2_denom) != ZERO)
    {
      // The curves do not have a vertical slope at p.
      // Compare their second-order derivative by x:
      _derive_by_x_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_x_at (p, 2, slope2_numer, slope2_denom);
    }
    else
    {
      // Both curves have a vertical slope at p.
      // Compare their second-order derivative by y:
      _derive_by_y_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_y_at (p, 2, slope2_numer, slope2_denom);
    }

    if (CGAL::compare (slope1_numer*slope2_denom,
		       slope2_numer*slope1_denom) != EQUAL)
    {
      // Different curvatures at p - the mutiplicity of p is 2:
      return (2);
    }

    // If we reached here, the multiplicity of the intersection point is 3:
    return (3);
  }
  //@}

};

/*!
 * Exporter for x-monotone conic arcs.
 */
template <class Conic_arc_2>
std::ostream& operator<< (std::ostream& os, 
                          const _Conic_x_monotone_arc_2<Conic_arc_2>& arc)
{
  // Output the supporting conic curve.
  os << "{" << CGAL::to_double(arc.r()) << "*x^2 + "
     << CGAL::to_double(arc.s()) << "*y^2 + "
     << CGAL::to_double(arc.t()) << "*xy + " 
     << CGAL::to_double(arc.u()) << "*x + "
     << CGAL::to_double(arc.v()) << "*y + "
     << CGAL::to_double(arc.w()) << "}";

  // Output the endpoints.
  os << " : (" << CGAL::to_double(arc.source().x()) << "," 
     << CGAL::to_double(arc.source().y()) << ") ";

  if (arc.orientation() == CLOCKWISE)
    os << "--cw-->";
  else if (arc.orientation() == COUNTERCLOCKWISE)
    os << "--ccw-->";
  else
    os << "--l-->";

  os << " (" << CGAL::to_double(arc.target().x()) << "," 
     << CGAL::to_double(arc.target().y()) << ")";

  return (os);
}

CGAL_END_NAMESPACE

#endif
