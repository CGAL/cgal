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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_BEZIER_POINT_2_H
#define CGAL_BEZIER_POINT_2_H

/*! \file
 * Header file for the _Bezier_point_2 class.
 */

#include <CGAL/Arr_traits_2/Bezier_curve_2.h>
#include <CGAL/Handle_for.h>
#include <list>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of a point on a Bezier curve. The point has algebraic
 * coefficients, with an additional list of originator. An originator is a
 * pair of the form <B(t), t0>, meaning that this point is obtained by
 * computing B(t0) on the curve B(t).
 */

// Forward declaration:
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_point_2;

template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_point_2_rep
{
  friend class _Bezier_point_2<Rat_kernel_, Alg_kernel_, Nt_traits_>;

public:

  typedef Rat_kernel_                             Rat_kernel;
  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  
  typedef typename Rat_kernel::Point_2            Rat_point_2;
  typedef typename Alg_kernel::Point_2            Alg_point_2;
  typedef typename Nt_traits::Rational            Rational;
  typedef typename Nt_traits::Algebraic           Algebraic;

private:

  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>              Curve_2;
  typedef std::pair<Curve_2, Algebraic>           Originator;
  typedef std::list<Originator>                   Orig_list;
  typedef typename Orig_list::const_iterator      Orig_iter;

  Algebraic         _x;           // The x-coordinate.
  Algebraic         _y;           // The y-coordinate.
  Orig_list         _origs;       // The list of originators.

public:

  /*! Default constructor. */
  _Bezier_point_2_rep () :
    _x (0),
    _y (0)
  {}

  /*!
   * Constructor with coordinates.
   */
  _Bezier_point_2_rep (const Algebraic& x, const Algebraic& y) :
    _x (x),
    _y (y)
  {}

  /*!
   * Constructor from a point with rational coordinates.
   */
  _Bezier_point_2_rep (const Rat_point_2& p)
  {
    Nt_traits   nt_traits;
    _x = nt_traits.convert (p.x());
    _y = nt_traits.convert (p.y());
  }

  /*!
   * Constructor from a point with algebraic coordinates.
   */
  _Bezier_point_2_rep (const Alg_point_2& p) :
    _x (p.x()),
    _y (p.y())
  {}

  /*!
   * Constructor given an originating curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, const Rational& t0)
  {
    // Set the point coordinates.
    Nt_traits           nt_traits;
    const Rat_point_2   p = B(t0);

    _x = nt_traits.convert (p.x());
    _y = nt_traits.convert (p.y());

    // Create the originator pair <B(t), t0>.
    _origs.push_back (Originator (B, nt_traits.convert (t0)));
  }

  /*!
   * Constructor given an originating curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, const Algebraic& t0)
  {
    // Set the point coordinates.
    const Alg_point_2   p = B(t0);

    _x = p.x();
    _y = p.y();

    // Create the originator pair <B(t), t0>.
    _origs.push_back (Originator (B, t0));
  }

};

template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_point_2 :
  public Handle_for<_Bezier_point_2_rep<Rat_kernel_,
                                        Alg_kernel_,
                                        Nt_traits_> >
{
public:

  typedef Rat_kernel_                             Rat_kernel;
  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>              Self;

private:

  typedef _Bezier_point_2_rep<Rat_kernel,
                              Alg_kernel,
                              Nt_traits>          Bpt_rep;
  typedef Handle_for<Bpt_rep>                     Bpt_handle;

public:

  typedef typename Bpt_rep::Rat_point_2           Rat_point_2;
  typedef typename Bpt_rep::Alg_point_2           Alg_point_2;
  typedef typename Bpt_rep::Rational              Rational;
  typedef typename Bpt_rep::Algebraic             Algebraic;
  typedef typename Bpt_rep::Curve_2               Curve_2;
  typedef typename Bpt_rep::Orig_iter             Originator_iterator;

  /*!
   * Default constructor.
   */
  _Bezier_point_2 () :
    Bpt_handle (Bpt_rep())
  {}

  /*!
   * Copy constructor.
   */
  _Bezier_point_2 (const Self& bpt) :
    Bpt_handle (bpt)
  {}

  /*!
   * Constructor with coordinates.
   */
  _Bezier_point_2 (const Algebraic& x, const Algebraic& y) :
    Bpt_handle (Bpt_rep (x, y))
  {}

  /*!
   * Constructor from a point with rational coordinates.
   */
  _Bezier_point_2 (const Rat_point_2& p) :
    Bpt_handle (Bpt_rep (p))
  {}

  /*!
   * Constructor from a point with algebraic coordinates.
   */
  _Bezier_point_2 (const Alg_point_2& p) :
    Bpt_handle (Bpt_rep (p))
  {}

  /*!
   * Constructor given an originating curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, const Rational& t0) :
    Bpt_handle (Bpt_rep (B, t0))
  {}

  /*!
   * Constructor given an originating curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, const Algebraic& t0) :
    Bpt_handle (Bpt_rep (B, t0))
  {}

  /*!
   * Check if the two handles refer to the same object.
   */
  bool is_same (const Self& pt) const
  {
    return (this->identical (pt));
  }

  /*!
   * Check for equality.
   */
  bool equals (const Self& pt) const
  {
    if (this->identical (pt))
      return (true);
    
    return (CGAL::compare (_rep()._x, pt._rep()._x) == EQUAL &&
            CGAL::compare (_rep()._y, pt._rep()._y) == EQUAL);
  }

  /*!
   * Get the x-coordinate.
   */
  const Algebraic& x () const
  {
    return (_rep()._x);
  }

  /*!
   * Get the y-coordinate.
   */
  const Algebraic& y () const
  {
    return (_rep()._y);
  }

  /*!
   * Check if the given curve is an originator of the point.
   * \param B A Bezier curve.
   * \param t0 Output: A t-value that satisfies B(t_0) == *this
   *                   (if one exists).
   * \return Whether the point is originated by B.
   */
  bool is_originator (const Curve_2& B, Algebraic& t0) const
  {
    // Scan the list of originators and look for B.
    typename Bpt_rep::Orig_iter     it = _rep()._origs.begin();
    typename Bpt_rep::Orig_iter     end = _rep()._origs.end();

    while (it != end)
    {
      if (B.is_same (it->first))
      {
        t0 = it->second;
        return (true);
      }

      ++it;
    }

    return (false);
  }

  /*!
   * Add the given curve to the list of originators.
   * \param B A Bezier curve.
   * \param t0 The t-value.
   * \pre t0 must be between 0 and 1.
   */
  void add_originator (const Curve_2& B, const Algebraic& t0) const
  {
    Bpt_rep&  rep = const_cast<Bpt_rep&> (_rep());

    rep._origs.push_back (typename Bpt_rep::Originator (B, t0));
    return;
  }

  /*!
   * Get the range of originators (pair of <Curve_2, Algebraic>).
   */
  Originator_iterator originators_begin () const
  {
    return (_rep()._origs.begin());
  }

  Originator_iterator originators_end () const
  {
    return (_rep()._origs.end());
  }

private:

  // Get the representation.
  inline const Bpt_rep& _rep () const
  {
    return (*(this->ptr()));
  }

  inline Bpt_rep& _rep ()
  {
    return (*(this->ptr()));
  }
};

/*!
 * Exporter for Bezier points.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Bezier_point_2<Rat_kernel, Alg_kernel, Nt_traits> & pt)
{
  os << pt.x() << ' ' << pt.y();
  return (os);
}

CGAL_END_NAMESPACE

#endif
