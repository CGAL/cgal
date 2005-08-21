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

#ifndef CGAL_RATIONAL_ARC_2_H
#define CGAL_RATIONAL_ARC_2_H

/*! \file
 * Header file for the _Rational_arc_2 class.
 */

#include <vector>
#include <list>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of an segment of a rational function, given as:
 *
 *        D(x)
 *   y = ------              x_l <= x <= x_r
 *        N(x)
 *
 * where D and N are polynomial with integer (or rational) coefficients.
 * The class is templated with two parameters: 
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers (defined by Nt_traits::Algebraic).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types. 
 */

template <class Alg_kernel_, class Nt_traits_>
class _Rational_arc_2
{
public:

  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  typedef _Rational_arc_2<Alg_kernel, Nt_traits>  Self;
  
  typedef typename Alg_kernel::FT                 Algebraic;
  typedef typename Alg_kernel::Point_2            Point_2;

  typedef typename Nt_traits::Integer             Integer;
  typedef typename Nt_traits::Rational            Rational;
  typedef std::vector<Rational>                   Rat_vector;
  typedef typename Nt_traits::Polynomial          Polynomial;

private:

  typedef std::pair<Point_2, unsigned int>        Intersection_point_2;

  Polynomial        pnum;       // The polynomial in the numerator.
  Polynomial        pden;       // The polynomial in the denominator.
  Point_2           left_pt;    // The left endpoint.
  Point_2           right_pt;   // The right endpoint.
  bool              valid;      // Is the arc valid.

public:

  /// \name Constrcution methods.
  //@{

  /*!
   * Default constructor.
   */
  _Rational_arc_2 () :
    valid (false)
  {}

  /*!
   * Constructor of a polynomial arc, defined by y = p(x), x_min <= x <= x_max.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param x_min The minimum of the x-range.
   * \param x_max The maximum of the x-range.
   * \pre The following should hold: x_min < x_max.
   */
  _Rational_arc_2 (const Rat_vector& pcoeffs,
		   const Algebraic& x_min, const Algebraic& x_max) :
    valid (false)
  {
    CGAL_precondition (CGAL::compare (x_min, x_max) == SMALLER);
    
    // Set the numerator polynomial.
    Nt_traits    nt_traits;
    Integer      p_factor;

    if (nt_traits.construct_polynomial (&(pcoeffs[0]),
					pcoeffs.size() - 1,
					pnum,
					p_factor))
    {
      nt_traits.scale (pnum, p_factor);
    }

    // Define the denominator to be a constant polynomial.
    Integer      denom_coeffs[1];

    denom_coeffs [0] = 1;
    pden = nt_traits.construct_polynomial (denom_coeffs, 0);
    
    // Set the endpoints.
    left_pt = Point_2 (x_min, nt_traits.evaluate_at (pnum, x_min));
    right_pt = Point_2 (x_max, nt_traits.evaluate_at (pnum, x_max));

    // Mark that the arc is valid.
    valid = true;
  }    
  
  /*!
   * Constructor of a polynomial arc, defined by y = p(x)/q(x), 
   * where: x_min <= x <= x_max.
   * \param pcoeffs The rational coefficients of the polynomial p(x).
   * \param qcoeffs The rational coefficients of the polynomial q(x).
   * \param x_min The minimum of the x-range.
   * \param x_max The maximum of the x-range.
   * \pre The following should hold: x_min < x_max,
   *      and q(x) != 0 for all x_min <= x <= x_max.
   */
  _Rational_arc_2 (const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
		   const Algebraic& x_min, const Algebraic& x_max) :
    valid (false)
  {
    CGAL_precondition (CGAL::compare (x_min, x_max) == SMALLER);
    
    // Set the numerator and denominator polynomials.
    Nt_traits    nt_traits;
    Integer      p_factor, q_factor;

    if (nt_traits.construct_polynomial (&(pcoeffs[0]),
					pcoeffs.size() - 1,
					pnum,
					p_factor))
    {
      nt_traits.scale (pnum, p_factor);
    }

    if (nt_traits.construct_polynomial (&(qcoeffs[0]),
					qcoeffs.size() - 1,
					pden,
					q_factor))
    {
      nt_traits.scale (pden, q_factor);
    }
    else
    {
      // q cannot be a zero polynomial:
      CGAL_assertion_msg (false, 
			  "zero polynomial specified as the denominator.");
      return;
    }

    // Make sure that q has no real roots between x_min and x_max.
    if (nt_traits.degree (pden) > 0)
    {
      std::list<Algebraic>                          q_roots;
      typename std::list<Algebraic>::const_iterator x_iter;
      bool                                 q_has_no_roots_in_the_interval;

      nt_traits.compute_polynomial_roots (pden,
					  std::back_inserter (q_roots));
      q_has_no_roots_in_the_interval = true;

      for (x_iter = q_roots.begin(); x_iter != q_roots.end(); ++x_iter)
      {
	if (CGAL::compare (x_min, *x_iter) != LARGER && 
	    CGAL::compare (*x_iter, x_max) != LARGER)
	{
	  q_has_no_roots_in_the_interval = false;
	  break;
	}
      }

      CGAL_precondition_msg (q_has_no_roots_in_the_interval,
			     "the rational arc cannot contain poles.");
      
      if (! q_has_no_roots_in_the_interval)
	return;
    }

    // Set the endpoints.
    left_pt = Point_2 (x_min, nt_traits.evaluate_at (pnum, x_min) /
		              nt_traits.evaluate_at (pden, x_min));
    right_pt = Point_2 (x_max, nt_traits.evaluate_at (pnum, x_max) /
                               nt_traits.evaluate_at (pden, x_max));

    // Mark that the arc is valid.
    valid = true;
  }    
  
  //@}

  /// \name Accessing the arc properties.
  //@{

  /*! Get the numerator polynomial of the underlying rational function. */
  const Polynomial& numerator () const
  {
    return (pnum);
  }

  /*! Get the denominator polynomial of the underlying rational function. */
  const Polynomial& denominator () const
  {
    return (pden);
  }

  /*! Get the left endpoint. */
  const Point_2& left () const
  {
    CGAL_precondition (valid);
    return (left_pt);
  }

  /*! Get the right endpoint. */
  const Point_2& right () const
  {
    CGAL_precondition (valid);
    return (right_pt);
  }

  /*! Check if the arc is valid. */
  bool is_valid () const
  {
    return (valid);
  }
  //@}

  /// \name Predicates.
  //@{

  /*!
   * Get the relative position of the point with respect to the rational arc.
   * \param p The query point.
   * \pre p is in the x-range of the arc.
   * \return SMALLER if the point is below the arc;
   *         LARGER if the point is above the arc;
   *         EQUAL if p lies on the arc.
   */
  Comparison_result point_position (const Point_2& p) const
  {
    // Make sure that p is in the x-range of the arc.
    CGAL_precondition (valid);
    CGAL_precondition (_is_in_x_range (p.x()));

    // Evaluate the rational function at x(p):
    Nt_traits   nt_traits;
    Algebraic   y = nt_traits.evaluate_at (pnum, p.x()) /
                    nt_traits.evaluate_at (pden, p.x());
    
    // Compare the resulting y-coordinate with y(p):
    return (CGAL::compare (p.y(), y));
  }

  /*!
   * Compare the slopes of the arc with another given arc at their given 
   * intersection point.
   * \param cv The given arc.
   * \param p The intersection point.
   * \param mult Output: The mutiplicity of the intersection point.
   * \return SMALLER if (*this) slope is less than cv's;
   *         EQUAL if the two slopes are equal;
   *         LARGER if (*this) slope is greater than cv's.
   */
  Comparison_result compare_slopes (const Self& arc,
				    const Point_2& p,
				    unsigned int& mult) const
  {
    // Make sure that p is in the x-range of both arcs.
    CGAL_precondition (valid);
    CGAL_precondition (arc.valid);
    CGAL_precondition (_is_in_x_range (p.x()) &&
		       arc._is_in_x_range (p.x()));

    // Check the case of overlapping arcs.
    if (_has_same_base (arc))
      return (EQUAL);

    // The intersection point may have a maximal multiplicity value of:
    //   max (deg(p1) + deg(q2), deg(q1) + deg(p2)).
    const Algebraic&  _x = p.x();
    Nt_traits         nt_traits;
    Polynomial        pnum1 = this->pnum;
    Polynomial        pden1 = this->pden;
    const bool        simple_poly1 = nt_traits.degree (pden1);
    Polynomial        pnum2 = arc.pnum;
    Polynomial        pden2 = arc.pden;
    const bool        simple_poly2 = nt_traits.degree (pden2);
    int               max_mult;
    Algebraic         d1, d2;
    Comparison_result res;

    max_mult = nt_traits.degree (pnum1) + nt_traits.degree (pden2);

    if (nt_traits.degree (pden1) + nt_traits.degree (pnum2) > max_mult)
      max_mult = nt_traits.degree (pden1) + nt_traits.degree (pnum2);

    for (mult = 1; mult <= static_cast<unsigned int>(max_mult); mult++)
    {
      // Compute the current derivative. Use the equation:
      //
      // (p(x) / q(x))' = (p'(x)*q(x) - p(x)*q'(x)) / q^2(x)
      if (simple_poly1)
      {
	pnum1 = nt_traits.derive(pnum1);
      }
      else
      {
	pnum1 = nt_traits.derive(pnum1)*pden1 - pnum1*nt_traits.derive(pden1);
	pden1 *= pden1;
      }

      if (simple_poly2)
      {
	pnum2 = nt_traits.derive(pnum2);
      }
      else
      {
	pnum2 = nt_traits.derive(pnum2)*pden2 - pnum2*nt_traits.derive(pden2);
	pden2 *= pden2;
      }

      // Compute the two derivative values and compare them. 
      d1 = nt_traits.evaluate_at (pnum1, _x) / 
	   nt_traits.evaluate_at (pden1, _x);
      d2 = nt_traits.evaluate_at (pnum2, _x) / 
	   nt_traits.evaluate_at (pden2, _x);

      res = CGAL::compare (d1, d2);

      // Stop here in case the derivatives are not equal.
      if (res != EQUAL)
	return (res);
    }

    // If we reached here, there is an overlap.
    return (EQUAL);
  }

  /*!
   * Check whether the two arcs are equal (have the same graph).
   * \param arc The compared arc.
   * \return (true) if the two arcs have the same graph; (false) otherwise.
   */
  bool equals (const Self& arc) const
  {
    // The two arc must have the same base rational function.
    CGAL_precondition (valid);
    CGAL_precondition (arc.valid);
    if (! _has_same_base (arc))
      return (false);

    // Check that the arc endpoints are the same.
    Alg_kernel   ker;

    return (ker.equal_2_object() (left_pt, arc.left_pt) &&
	    ker.equal_2_object() (right_pt, arc.right_pt));
  }

  /*!
   * Check whether it is possible to merge the arc with the given arc.
   * \param arc The query arc.
   * \return (true) if it is possible to merge the two arcs;
   *         (false) otherwise.
   */
  bool can_merge_with (const Self& arc) const
  {
    // In order to merge the two arcs, they should have the same base rational
    // function.
    CGAL_precondition (valid);
    CGAL_precondition (arc.valid);
    if (! _has_same_base (arc))
      return (false);

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    Alg_kernel   ker;

    return (ker.equal_2_object() (right_pt, arc.left_pt) ||
            ker.equal_2_object() (left_pt, arc.right_pt));
  }
  //@}

  /// \name Constructions of points and curves.
  //@{
  
  /*!
   * Compute the intersections with the given arc.
   * \param arc The given intersecting arc.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template<class OutputIterator>
  OutputIterator intersect (const Self& arc,
                            OutputIterator oi) const
  {
    CGAL_precondition (valid);
    CGAL_precondition (arc.valid);

    if (_has_same_base (arc))
    {
      // Let p1 be the rightmost of the two left endpoints and let p2 be the 
      // leftmost of the two right endpoints.
      Alg_kernel   ker;

      const Point_2&    p1 = 
	(ker.compare_x_2_object() (left_pt, arc.left_pt) == LARGER) ?
	left_pt : arc.left_pt;
      const Point_2&    p2 = 
	(ker.compare_x_2_object() (right_pt, arc.right_pt) == SMALLER) ?
	right_pt : arc.right_pt;
      Comparison_result res = ker.compare_x_2_object() (p1, p2);

      if (res == SMALLER)
      {
	// If the range [p1, p2] is non-trivial, we have an overlapping
	// segment.
	Self      overlap_arc (*this);

	overlap_arc.left_pt = p1;
	overlap_arc.right_pt = p2;

	*oi = make_object (overlap_arc);
	++oi;
      }
      else if (res == EQUAL)
      {
	// We have a single overlapping point:
	Intersection_point_2  ip (p1, 0);
	
	*oi = make_object (ip);
	++oi;
      }

      return (oi);
    }
    
    // We wish to find the intersection points between:
    //
    //   y = p1(x)/q1(x)    and     y = p2(x)/q2(x)
    //
    // It is clear that the x-coordinates of the intersection points are
    // the roots of the polynomial: ip(x) = p1(x)*q2(x) - p2(x)*q1(x).
    Nt_traits            nt_traits;
    Polynomial           ipoly = pnum*arc.pden - arc.pnum*pden;
    std::list<Algebraic>                           xs;
    typename std::list<Algebraic>::const_iterator  x_iter;

    nt_traits.compute_polynomial_roots (ipoly,
					std::back_inserter(xs));
    
    // Go over the x-values we obtained. For each value produce an
    // intersection point if it is contained in the x-range of both curves.
    unsigned int                     mult;

    for (x_iter = xs.begin(); x_iter != xs.end(); ++x_iter)
    {
      if (_is_in_x_range(*x_iter) && arc._is_in_x_range(*x_iter))
      {
	// Compute the intersection point and obtain its multiplicity.
	Point_2    p (*x_iter, nt_traits.evaluate_at (pnum, *x_iter) /
                               nt_traits.evaluate_at (pden, *x_iter));

	this->compare_slopes (arc, p, mult);
    
	// Output the intersection point:
	Intersection_point_2  ip (p, mult);
	
	*oi = make_object (ip);
	++oi;
      }
    }

    return (oi);
  }

  /*!
   * Split the arc into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the arc (not one of its endpoints).
   */
  void split (const Point_2& p,
              Self& c1, Self& c2) const
  {
    CGAL_precondition (valid);

    // Make sure that p lies on the interior of the arc.
    CGAL_precondition_code (
      Alg_kernel   ker;
    );
    CGAL_precondition (this->point_position(p) == EQUAL &&
                       ! ker.equal_2_object() (p, left_pt) &&
                       ! ker.equal_2_object() (p, right_pt));

    // Make copies of the current arc.
    c1 = *this;
    c2 = *this;

    c1.right_pt = p;
    c2.left_pt = p;

    return;
  }

  /*!
   * Merge the current arc with the given arc.
   * \param arc The arc to merge with.
   * \pre The two arcs are mergeable.
   */
  void merge (const Self& arc)
  {
    CGAL_precondition (valid);
    CGAL_precondition (arc.valid);
    CGAL_precondition (this->can_merge_with (arc));

    // Check if we should extend the arc to the left or to the right.
    Alg_kernel   ker;

    if (ker.equal_2_object() (right_pt, arc.left_pt))
    {
      // Extend the arc to the right.
      right_pt = arc.right_pt;
    }
    else
    {
      CGAL_precondition (ker.equal_2_object() (left_pt, arc.right_pt));

      // Extend the arc to the left.
      left_pt = arc.left_pt;
    }

    return;
  }
  //@}

private:

  /// \name Auxiliary (private) functions.
  //@{

  /*!
   * Check if the given x-value is in the x-range of the arc.
   */
  bool _is_in_x_range (const Algebraic& x) const
  {
    Comparison_result  res1 = CGAL::compare (left_pt.x(), x);
    Comparison_result  res2 = CGAL::compare (x, right_pt.x());

    return (res1 != LARGER && res2 != LARGER);
  }

  /*!
   * Check if the underlying rational fucntion is the same in the given arc.
   * \param arc The given arc.
   * \return (true) if arc's underlying rational fucntion is the same
   *         as of *this; (false) otherwise.
   */
  bool _has_same_base (const Self& arc) const
  {
    // p1(x)/q1(x) == p2(x)/q2(x) if and only if p1*q2 = p2*q1:
    return (pnum * arc.pden == pden * arc.pnum);
  }
  //@}
};

/*!
 * Exporter for rational arcs.
 */
template <class Alg_kernel, class Nt_traits>
std::ostream& operator<< (std::ostream& os, 
                          const _Rational_arc_2<Alg_kernel, Nt_traits>& arc)
{
  // Output the supporting rational function and the x-range of the arc.
  os << "y = (" << arc.numerator() << ") / (" << arc.denominator() << ") : ["
     << CGAL::to_double(arc.left().x()) << " -- " 
     << CGAL::to_double(arc.right().x()) << "]";

  return (os);
}

CGAL_END_NAMESPACE

#endif
