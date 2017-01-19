// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
//
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>
//                 Iddo Hanniel <iddoh@cs.technion.ac.il>

#ifndef CGAL_BEZIER_POINT_2_H
#define CGAL_BEZIER_POINT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Header file for the _Bezier_point_2 class.
 */

#include <CGAL/Arr_geometry_traits/Bezier_curve_2.h>
#include <CGAL/Arr_geometry_traits/Bezier_cache.h>
#include <CGAL/Handle_for.h>
#include <list>
#include <ostream>

namespace CGAL {

/*! \class _Bezier_point_2
 * Representation of a point on a Bezier curve. The point has algebraic
 * coefficients, with an additional list of originator. An originator is a
 * pair of the form <B(t), t0>, meaning that this point is obtained by
 * computing B(t0) on the curve B(t).
 */

// Forward declaration:
template <class RatKernel_, class AlgKernel_, class NtTraits_,
          class BoundingTraits_>
class _Bezier_point_2;

template <class RatKernel_, class AlgKernel_, class NtTraits_,
          class BoundingTraits_>
class _Bezier_point_2_rep
{
  friend class _Bezier_point_2<RatKernel_, AlgKernel_, NtTraits_, 
                               BoundingTraits_>;

public:

  typedef RatKernel_                              Rat_kernel;
  typedef AlgKernel_                              Alg_kernel;
  typedef NtTraits_                               Nt_traits;
  typedef BoundingTraits_                         Bounding_traits;

  typedef typename Rat_kernel::Point_2            Rat_point_2;
  typedef typename Alg_kernel::Point_2            Alg_point_2;
  typedef typename Nt_traits::Rational            Rational;
  typedef typename Nt_traits::Algebraic           Algebraic;

  typedef typename Bounding_traits::Bez_point_bound   Bez_point_bound;
  typedef typename Bounding_traits::Bez_point_bbox    Bez_point_bbox;

  typedef _Bezier_point_2_rep<Rat_kernel,
                              Alg_kernel,
                              Nt_traits,
                              Bounding_traits>    Self; 

private:

  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>        Curve_2;

  typedef _Bezier_cache<Nt_traits>                Bezier_cache;

  /*! \class Originator
   * Stores information on the original curve the Bezier point came from.
   */
  class Originator
  {
  private:

    Curve_2             _curve;     /*!< The originating curve. */
    unsigned int        _xid;       /*!< Serial number of the originating
                                         x-monotone curve. */
    Bez_point_bound     _bpb;       /*!< Bounding information for the
                                         point: bouding control polygon,
                                         point type, etc. */
    Algebraic          *p_t;        /*!< The algebraic parameter for the
                                         point (if available). */

  public:

    /*! Constructor, given an exact algebraic representation. */
    Originator (const Curve_2& c, const Algebraic& t) :
      _curve (c),
      _xid (0),
      p_t (NULL)
    {
      set_parameter (t);
    }

    /*! Constructor, given an exact algebraic representation. */
    Originator (const Curve_2& c, unsigned int xid,
                const Algebraic& t) :
      _curve (c),
      _xid (xid),
      p_t (NULL)
    {
      set_parameter (t);
    }

    /*! Constructor with bounding information and no exact representation. */
    Originator (const Curve_2& c, const Bez_point_bound& bpb) :
      _curve (c),
      _xid (0),
      _bpb (bpb),
      p_t (NULL)
    {}

    /*! Constructor with bounding information and no exact representation. */
    Originator (const Curve_2& c, unsigned int xid,
                const Bez_point_bound& bpb) :
      _curve (c),
      _xid (xid),
      _bpb (bpb),
      p_t (NULL)
    {}

    /*! Copy constructor. */
    Originator (const Originator& other) :
      _curve (other._curve),
      _xid (other._xid),
      _bpb (other._bpb),
      p_t (NULL)
    {
      // Deep copy of lazy instantiation
      if (other.p_t != NULL)
        p_t = new Algebraic (*(other.p_t));
    }

    /*! Destructor. */
    ~Originator() 
    {
      if (p_t != NULL)
        delete p_t;
    }

    /*! Assignment operator. */
    Originator& operator= (const Originator& other)
    {
      // Avoid self assignments.
      if (this == &other)
        return (*this);

      // Free memory, if necessary.
      if (p_t != NULL)
        delete p_t;
      p_t = NULL;

      // Copy the data members.
      _curve = other._curve;
      _xid = other._xid;
      _bpb = other._bpb;

      // Deep copy of lazy instantiation
      if (other.p_t != NULL)
        p_t = new Algebraic (*(other.p_t));

      return (*this);
    }

    /*! Get the originating curve. */
    const Curve_2& curve () const
    {
      return (_curve);
    }

    /*! Get the serial number of the originating x-monotone curve. */
    unsigned int xid () const
    {
      return (_xid);
    }

    /*! Get the bounding information. */
    const Bez_point_bound& point_bound () const
    {
      return (_bpb);
    }

    /*! Update the bounding information. */
    void update_point_bound (const Bez_point_bound& bpb)
    {
      _bpb = bpb;
      return;
    }

    /*! Check if the algberaic parameter is available. */
    bool has_parameter () const
    {
      return (p_t != NULL);
    }

    /*! 
     * Get the algebraic parameter.
     * \pre The parameter value is available.
     */
    const Algebraic& parameter () const
    {
      CGAL_precondition (p_t != NULL);
      return (*p_t);
    }

    /*!
     * Set the parameter value. 
     * \pre The parameter value is not yet set.
     */
    void set_parameter (const Algebraic& t)
    {
      CGAL_precondition (p_t == NULL);
      
      p_t = new Algebraic (t);

      // Update the Bez_point_bound by converting t to an interval of doubles
      // and setting _bpb accordingly.
      Nt_traits                         nt_traits;
      const std::pair<double, double>&  t_bnd = nt_traits.double_interval (t);

      _bpb.t_min = t_bnd.first;
      _bpb.t_max = t_bnd.second;

      return;
    }

    /*!
     * Set the serial number of the originating x-monotone curve.
     * \param xid the new serial number of the originating x-monotone curve.
     * \pre The current xid() is 0.
     * \pre xid is possitive.
     */
    void set_xid (unsigned int xid)
    {
      CGAL_precondition (_xid == 0);
      CGAL_precondition (xid > 0);

      _xid = xid;
      return;
    }
  };

  /*! \struct Subcurve
   * Auxilary structure for the vertical_position() function.
   */
  typedef typename Bounding_traits::Control_points   Control_points;
  typedef typename Bounding_traits::NT               BoundNT;

  struct Subcurve
  {
    Control_points   ctrl;      /*!< The control points. */
    BoundNT          t_min;     /*!< Minimal parameter value. */
    BoundNT          t_max;     /*!< Maximal parameter value. */

    /*! Constructor given control points an a t-range. */
    Subcurve (const Control_points& _ctrl,
              const BoundNT& _tmin, 
              const BoundNT& _tmax) :
      ctrl (_ctrl),
      t_min (_tmin),
      t_max (_tmax)
    {}

    /*! Constructor given a t-range. */
    Subcurve (const BoundNT& _tmin, 
              const BoundNT& _tmax) :
      t_min (_tmin),
      t_max (_tmax)
    {}
  };

  typedef std::list<Originator>                   Orig_list;
  typedef typename Orig_list::const_iterator      Orig_const_iter;
  typedef typename Orig_list::iterator            Orig_iter;

  Algebraic        *p_alg_x;   /*! The exact x-coordinate (if known). */
  Rational         *p_rat_x;   /*! The x-coordinate, in case it is rational. */
  Algebraic        *p_alg_y;   /*! The exact y-coordinate (if known). */
  Rational         *p_rat_y;   /*! The y-coordinate, in case it is rational. */
  Orig_list        _origs;     /*! The list of originators. */
  Bez_point_bbox   _bbox;      /*! A bounding box. */

public:

  /*! Default constructor. */
  _Bezier_point_2_rep () :
    p_alg_x (NULL),
    p_rat_x (NULL),
    p_alg_y (NULL),
    p_rat_y (NULL)
  {}

  /*! Copy constructor. */
  _Bezier_point_2_rep (const Self& pt) :
    p_alg_x (NULL),
    p_rat_x (NULL),
    p_alg_y (NULL),
    p_rat_y (NULL),
    _origs (pt._origs),
    _bbox (pt._bbox)
  {
    if (pt.p_alg_x != NULL)
      p_alg_x = new Algebraic (*(pt.p_alg_x));
    if (pt.p_rat_x != NULL)
      p_rat_x = new Rational (*(pt.p_rat_x));
    if (pt.p_alg_y != NULL)
      p_alg_y = new Algebraic (*(pt.p_alg_y));
    if (pt.p_rat_y != NULL)
      p_rat_y = new Rational (*(pt.p_rat_y));
  }

  /*!
   * Constructor with algebraic coordinates.
   * \param x The exact x-coordinate.
   * \param y The exact y-coordinate.
   */
  _Bezier_point_2_rep (const Algebraic& x, const Algebraic& y, bool) : 
    p_rat_x (NULL),
    p_rat_y (NULL)
  {
    p_alg_x = new Algebraic (x);
    p_alg_y = new Algebraic (y);

    // Initialize the bounding box.
    Nt_traits                         nt_traits;
    const std::pair<double, double>&  x_bnd = 
                                        nt_traits.double_interval (x);
    const std::pair<double, double>&  y_bnd = 
                                        nt_traits.double_interval (y);

    _bbox.min_x = x_bnd.first;
    _bbox.max_x = x_bnd.second;
    _bbox.min_y = y_bnd.first;
    _bbox.max_y = y_bnd.second;
  }

  /*!
   * Constructor with rational coordinates.
   * \param x The exact x-coordinate.
   * \param y The exact y-coordinate.
   */
  _Bezier_point_2_rep (const Rational& x, const Rational& y)
  {
    p_rat_x = new Rational (x);
    p_rat_y = new Rational (y);

    // Convert the rational coordinates to algebraic values.
    Nt_traits                         nt_traits;

    p_alg_x = new Algebraic (nt_traits.convert (x));
    p_alg_y = new Algebraic (nt_traits.convert (y));

    // Initialize the bounding box.
    _bbox.min_x = x;
    _bbox.max_x = x;
    _bbox.min_y = y;
    _bbox.max_y = y;
  }

  /*!
   * Constructor given an originating curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, const Rational& t0);

  /*!
   * Constructor given an x-monotone curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, unsigned int xid,
                       const Rational& t0);

  /*!
   * Constructor given an originating curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, const Algebraic& t0);

  /*!
   * Constructor given an x-monotone curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, unsigned int xid,
                       const Algebraic& t0);

  /*! Destructor. */
  ~_Bezier_point_2_rep ()
  {
    if (p_rat_x != NULL)
      delete p_rat_x;
    if (p_alg_x != NULL)
      delete p_alg_x;
    if (p_rat_y != NULL)
      delete p_rat_y;
    if (p_alg_y != NULL)
      delete p_alg_y;
  }

  /*! Assignment operator. */
  Self& operator= (const Self& pt)
  {
    if (this == &pt)
      return (*this);

    if (p_rat_x != NULL)
      delete p_rat_x;
    if (p_alg_x != NULL)
      delete p_alg_x;
    if (p_rat_y != NULL)
      delete p_rat_y;
    if (p_alg_y != NULL)
      delete p_alg_y;
    p_alg_x = p_rat_x = p_alg_y = p_rat_y = NULL;


    if (pt.p_alg_x != NULL)
      p_alg_x = new Algebraic (*(pt.p_alg_x));
    if (pt.p_rat_x != NULL)
      p_rat_x = new Rational (*(pt.p_rat_x));
    if (pt.p_alg_y != NULL)
      p_alg_y = new Algebraic (*(pt.p_alg_y));
    if (pt.p_rat_y != NULL)
      p_rat_y = new Rational (*(pt.p_rat_y));
  
    _origs = pt._origs;
    _bbox = pt._bbox;

    return (*this);
  }

  /*! Check if the point is exactly computed. */
  inline bool is_exact () const
  {
    return (p_alg_x != NULL && p_alg_y != NULL);
  }

  /*! Check if the point has rational coordinates. */
  inline bool is_rational () const
  {
    return (p_rat_x != NULL && p_rat_y != NULL);
  }

  /*!
   * Compare the x-coordinate with the coordinate of the given point.
   * \param pt The other point.
   * \param cache A cache for the vertical tangency points and the
   *              intersection points.
   * \return The comparison result;
   */
  Comparison_result compare_x (Self& pt,
                               Bezier_cache& cache);

  /*!
   * Compare the two point xy-lexicographically.
   * \param pt The other point.
   * \param cache A cache for the vertical tangency points and the
   *              intersection points.
   * \return The comparison result;
   */
  Comparison_result compare_xy (Self& pt,
                                Bezier_cache& cache);

  /*!
   * Determines the vertical position of the point with respect to an
   * x-monotone subcurve, given by its control polygon.
   * \param cp The control polygon of the subcurve.
   * \param t_min Defines the smallest parameter value of the subcurve. 
   * \param t_max Defines the largest parameter value of the subcurve. 
   * \return SMALLER if the point is located below the curve;
   *         LARGER if the point is located above the curve;
   *         EQUAL if we cannot determine its precise position.
   */  
  Comparison_result vertical_position (const Control_points& cp,
                                       const BoundNT& t_min,
                                       const BoundNT& t_max);

private:

  /*!
   * Refine the bounds of the point.
   * \return Whether it was possible to further refine the point.
   */
  bool _refine ();

  /*!
   * Make sure the originator parameters fit the bound box.
   */
  void _fit_to_bbox ();

  /*!
   * Compute the exact representation of the point.
   * \param cache A cache for the vertical tangency points and the
   *              intersection points.
   */
  void _make_exact (Bezier_cache& cache);
};

template <class RatKernel_, class AlgKernel_, class NtTraits_, 
          class BoundingTraits_>
class _Bezier_point_2 :
  public Handle_for<_Bezier_point_2_rep<RatKernel_,
                                        AlgKernel_,
                                        NtTraits_,
                                        BoundingTraits_> >
{
public:

  typedef RatKernel_                              Rat_kernel;
  typedef AlgKernel_                              Alg_kernel;
  typedef NtTraits_                               Nt_traits;
  typedef BoundingTraits_                         Bounding_traits;
  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>        Self;

private:

  typedef _Bezier_point_2_rep<Rat_kernel,
                              Alg_kernel,
                              Nt_traits,
                              Bounding_traits>    Bpt_rep;
  typedef Handle_for<Bpt_rep>                     Bpt_handle;

public:

  typedef typename Bpt_rep::Rat_point_2           Rat_point_2;
  typedef typename Bpt_rep::Alg_point_2           Alg_point_2;
  typedef typename Bpt_rep::Rational              Rational;
  typedef typename Bpt_rep::Algebraic             Algebraic;
  typedef typename Bpt_rep::Curve_2               Curve_2;
  typedef typename Bpt_rep::Originator            Originator;
  typedef typename Bpt_rep::Bezier_cache          Bezier_cache;

  typedef typename Bpt_rep::Orig_const_iter       Originator_iterator;
  typedef typename Bpt_rep::Bez_point_bound       Bez_point_bound;
  typedef typename Bpt_rep::Bez_point_bbox        Bez_point_bbox;

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
   * Constructor with algebraic coordinates (only for private use).
   */
  _Bezier_point_2 (const Algebraic& x, const Algebraic& y, bool dummy) :
    Bpt_handle (Bpt_rep (x, y, dummy))
  {}

  /*!
   * Constructor with rational coordinates.
   */
  _Bezier_point_2 (const Rational& x, const Rational& y) :
    Bpt_handle (Bpt_rep (x, y))
  {}

  /*!
   * Constructor given an originating curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, const Rational& t0) :
    Bpt_handle (Bpt_rep (B, t0))
  {}

  /*!
   * Constructor given an x-monotone curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, unsigned int xid,
                   const Rational& t0) :
    Bpt_handle (Bpt_rep (B, xid, t0))
  {}

  /*!
   * Constructor given an originating curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, const Algebraic& t0) :
    Bpt_handle (Bpt_rep (B, t0))
  {}

  /*!
   * Constructor given an x-monotone curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, unsigned int xid,
                   const Algebraic& t0) :
    Bpt_handle (Bpt_rep (B, xid, t0))
  {}

  /*!
   * Assignment operator.
   */
  Self& operator= (const Self& pt)
  {
    if (this == &pt || this->identical (pt))
      return (*this);

    Bpt_handle::operator= (pt);
    return (*this);
  }

  /*!
   * Check if the two handles refer to the same object.
   */
  bool is_same (const Self& pt) const
  {
    return (this->identical (pt));
  }

  /*!
   * Check if the point is computed in an exact manner.
   */
  bool is_exact () const
  {
    return (_rep().is_exact());
  }

  /*!
   * Check if the point has rational coordinates.
   */
  bool is_rational () const
  {
    return (_rep().is_rational());
  }

  /*!
   * Get the x-coordinate.
   * \pre The point is exactly computed.
   */
  const Algebraic& x () const
  {
    CGAL_precondition (_rep().is_exact());
    return (*(_rep().p_alg_x));
  }

  /*!
   * Get the y-coordinate.
   * \pre The point is exactly computed.
   */
  const Algebraic& y () const
  {
    CGAL_precondition (_rep().is_exact());
    return (*(_rep().p_alg_y));
  }

  /*!
   * Get the approximate coordinates.
   */
  std::pair<double, double> approximate () const
  {
    double      x, y;

    if (is_exact())
    {
      x = CGAL::to_double (*(_rep().p_alg_x));
      y = CGAL::to_double (*(_rep().p_alg_y));
    }
    else
    {
      x = CGAL::to_double ((_rep()._bbox.min_x + _rep()._bbox.max_x) / 2);
      y = CGAL::to_double ((_rep()._bbox.min_y + _rep()._bbox.max_y) / 2);
    }

    return (std::make_pair (x, y));
  }

  /*!
   * Convert to a rational point.
   * \pre The point has rational coordinates.
   */
  operator Rat_point_2 () const
  {
    CGAL_precondition (_rep().is_rational());
    
    return (Rat_point_2 (*(_rep().p_rat_x), *(_rep().p_rat_y)));
  }

  /*!
   * Refine the representation of the point.
   * \return Whether a refinement was possible.
   */
  bool refine () const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());
    
    return (rep._refine());
  }

  /*!
   * Make sure the originator parameters fit the bound box.
   */
  void fit_to_bbox () const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());
    
    return (rep._fit_to_bbox());
  }

  /*!
   * Compute the exact coordinates of the point.
   */
  void make_exact (Bezier_cache& cache) const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());

    rep._make_exact (cache);
    return;
  }

  /*!
   * Compare the x-coordinate with the coordinate of the given point.
   * \param pt The other point.
   * \param cache A cache for the vertical tangency points and the
   *              intersection points.
   * \return The comparison result;
   */
  Comparison_result compare_x (const Self& pt,
                               Bezier_cache& cache) const 
  {
    if (this->identical (pt))
      return (EQUAL);

    // Const cast since we are modifying the reps.
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());
    Bpt_rep&             rep_pt = const_cast<Bpt_rep&> (pt._rep());

    return (rep.compare_x (rep_pt, cache));
  }

  /*!
   * Compare the the two points xy-lexicographically.
   * \param pt The other point.
   * \param cache A cache for the vertical tangency points and the
   *              intersection points.
   * \return The comparison result;
   */
  Comparison_result compare_xy (const Self& pt,
                                Bezier_cache& cache) const 
  {
    if (this->identical (pt))
      return EQUAL;

    Self&                p1 = const_cast<Self&> (*this);
    Self&                p2 = const_cast<Self&> (pt);
    Comparison_result    res = p1._rep().compare_xy (p2._rep(), cache);

    if (res == EQUAL)
    {
      // If we find that two points are equal, we merge their lists of
      // originators and make them refer to the same representation.
      p1.merge_originators (p2);
      p2 = p1;

      CGAL_assertion (this->identical (pt));
    }
      
    return (res);
  }

  /*!
   * Determine if the two points are equal. 
   */
  bool equals (const Self& pt,
               Bezier_cache& cache) const
  {
    return (this->compare_xy (pt, cache) == EQUAL);
  }

  /*!
   * Determines the vertical position of the point with respect to an
   * x-monotone subcurve, given by its control polygon.
   * \param cp The control polygon of the subcurve.
   * \param t_min Defines the smallest parameter value of the subcurve. 
   * \param t_max Defines the largest parameter value of the subcurve. 
   * \return SMALLER if the point is located below the curve;
   *         LARGER if the point is located above the curve;
   *         EQUAL if we cannot determine its precise position.
   */
  Comparison_result vertical_position
      (const typename Bounding_traits::Control_points& cp,
       const typename Bounding_traits::NT& t_min,
       const typename Bounding_traits::NT& t_max) const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());
    return (rep.vertical_position(cp, t_min, t_max));
  }

  /*!
   * Get the originator of the point that is associated with the given curve.
   * \param B The Bezier curve.
   * \return An iterator pointing to the requested originator;
   *         originators_end() if B is not an originator of the point.
   */
  Originator_iterator get_originator (const Curve_2& B) const
  {
    // Scan the list of originators and look for B.
    typename Bpt_rep::Orig_const_iter     it = _rep()._origs.begin();
    typename Bpt_rep::Orig_const_iter     end = _rep()._origs.end();

    while (it != end)
    {
      if (B.is_same (it->curve()))
        return (it);

      ++it;
    }

    // If we reached here, we have not found an originator:
    return (it);
  }

  /*!
   * Get the originator of the point that is associated with the given
   * x-monotone curve.
   * \param B The Bezier curve.
   * \param xid The serial number of the x-monotone subcurve.
   * \return An iterator pointing to the requested originator;
   *         originators_end() if B is not an originator of the point.
   */
  Originator_iterator get_originator (const Curve_2& B,
                                      unsigned int xid) const
  {
    // Scan the list of originators and look for B.
    typename Bpt_rep::Orig_const_iter     it = _rep()._origs.begin();
    typename Bpt_rep::Orig_const_iter     end = _rep()._origs.end();

    while (it != end)
    {
      if (B.is_same (it->curve()))
      {
        // An x-monotone id that equals 0 means that the originator is not
        // associated with a specific x-monotone curve. Otherwise, we require
        // that the IDs match.
        if (it->xid() == 0 || it->xid() == xid)
          return (it);
      }

      ++it;
    }

    // If we reached here, we have not found an originator:
    return (it);
  }

  /*!
   * Get the range of originators.
   */
  Originator_iterator originators_begin () const
  {
    return (_rep()._origs.begin());
  }

  Originator_iterator originators_end () const
  {
    return (_rep()._origs.end());
  }

  /*!
   * Add an orinigator to the point.
   */
  void add_originator (const Originator& o) const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());

    rep._origs.push_back (o);
    return;
  }

  /*!
   * Update the xid field of the given orinigator.
   */
  void update_originator_xid (const Originator& o,
                              unsigned int xid) const
  {
    Originator&     orig = const_cast<Originator&> (o);

    orig.set_xid (xid);
    return;
  }

  /*!
   * Add the originators of the given point.
   */
  void merge_originators (const Self& pt) const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());
    Originator_iterator  org_it = pt.originators_begin();

    while (org_it != pt.originators_end())
    {
      rep._origs.push_back (typename Bpt_rep::Originator (*org_it));
      ++org_it;
    }

    return;
  }

  /*! Set the bounding box for the point. */
  void set_bbox (const Bez_point_bbox& bbox)
  {
    _rep()._bbox = bbox;
  }

  /*! Get the bounding box of the point. */
  void get_bbox (typename Bounding_traits::NT& min_x, 
                 typename Bounding_traits::NT& min_y, 
                 typename Bounding_traits::NT& max_x, 
                 typename Bounding_traits::NT& max_y) const
  {
    min_x = _rep()._bbox.min_x;
    min_y = _rep()._bbox.min_y;
    max_x = _rep()._bbox.max_x;
    max_y = _rep()._bbox.max_y;
  }

private:

  /*! Get the representation (const version). */
  inline const Bpt_rep& _rep () const
  {
    return (*(this->ptr()));
  }

  /*! Get the representation (non-const version). */
  inline Bpt_rep& _rep ()
  {
    return (*(this->ptr()));
  }

};

/*!
 * Exporter for Bezier points.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits, 
          class Bounding_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Bezier_point_2<Rat_kernel, Alg_kernel, Nt_traits,
                                  Bounding_traits> & pt)
{
  if (pt.is_exact())
  {
    os << CGAL::to_double (pt.x()) << ' ' << CGAL::to_double (pt.y());
  }
  else
  {
    typename Bounding_traits::NT   min_x, min_y, max_x, max_y;
    
    pt.get_bbox(min_x, min_y, max_x, max_y);
    os << '~' << CGAL::to_double ((min_x + max_x) / 2) 
       << " ~" << CGAL::to_double ((min_y + max_y) / 2);
  }

  return (os);
}

// ---------------------------------------------------------------------------
// Constructor given an originating curve and a rational t0 value.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_Bezier_point_2_rep
        (const Curve_2& B, const Rational& t0)
{
  // Insert an originator of a rational point.
  // Note that this constructor also takes care of the Bez_bound
  // for the originator.
  Nt_traits           nt_traits;
  Originator          org (B, nt_traits.convert (t0));
  Bez_point_bound     bound = org.point_bound();

  bound.type = Bez_point_bound::RATIONAL_PT;
  org.update_point_bound (bound);
  _origs.push_back (org);

  // Evaluate the point coordinates.
  const Rat_point_2&  p = B (t0);

  p_rat_x = new Rational (p.x());
  p_alg_x = new Algebraic (nt_traits.convert (*p_rat_x));
  p_rat_y = new Rational (p.y());
  p_alg_y = new Algebraic (nt_traits.convert (*p_rat_y));

  // Set the bounding box for this point.
  _bbox.min_x = *p_rat_x;
  _bbox.max_x = *p_rat_x;
  _bbox.min_y = *p_rat_y;
  _bbox.max_y = *p_rat_y;
}

// ---------------------------------------------------------------------------
// Constructor given an x-monotone curve and a rational t0 value.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_Bezier_point_2_rep
        (const Curve_2& B, unsigned int xid,
         const Rational& t0)
{
  // Insert an originator of a rational point.
  // Note that this constructor also takes care of the Bez_bound
  // for the originator.
  Nt_traits           nt_traits;
  Originator          org (B, xid, nt_traits.convert (t0));
  Bez_point_bound     bound = org.point_bound();

  bound.type = Bez_point_bound::RATIONAL_PT;
  org.update_point_bound (bound);
  _origs.push_back (org);

  // Evaluate the point coordinates.
  const Rat_point_2&  p = B (t0);

  p_rat_x = new Rational (p.x());
  p_alg_x = new Algebraic (nt_traits.convert (*p_rat_x));
  p_rat_y = new Rational (p.y());
  p_alg_y = new Algebraic (nt_traits.convert (*p_rat_y));

  // Set the bounding box for this point.
  _bbox.min_x = *p_rat_x;
  _bbox.max_x = *p_rat_x;
  _bbox.min_y = *p_rat_y;
  _bbox.max_y = *p_rat_y;
}

// ---------------------------------------------------------------------------
// Constructor given an originating curve and an algebraic t0 value.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_Bezier_point_2_rep
        (const Curve_2& B, const Algebraic& t0) :
  p_rat_x (NULL),
  p_rat_y (NULL)
{
  // Create the originator pair <B(t), t0>.
  // Note that this constructor also takes care of the Bez_bound
  // for the originator.
  _origs.push_back (Originator (B, t0));

  // Set the point coordinates.
  const Alg_point_2   p = B (t0);
  
  p_alg_x = new Algebraic (p.x());
  p_alg_y = new Algebraic (p.y());
  
  // Set the bounding box for this point, by converting x, y to two ranges
  // of doubles.
  Nt_traits                         nt_traits;
  const std::pair<double, double>&  x_bnd = 
                                        nt_traits.double_interval (*p_alg_x);
  const std::pair<double, double>&  y_bnd = 
                                        nt_traits.double_interval (*p_alg_y);

  _bbox.min_x = x_bnd.first;
  _bbox.max_x = x_bnd.second;
  _bbox.min_y = y_bnd.first;
  _bbox.max_y = y_bnd.second;
}

// ---------------------------------------------------------------------------
// Constructor given an x-monotone curve and an algebraic t0 value.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_Bezier_point_2_rep
        (const Curve_2& B, unsigned int xid,
         const Algebraic& t0) :
  p_rat_x (NULL),
  p_rat_y (NULL)
{
  // Create the originator pair <B(t), t0>.
  // Note that this constructor also takes care of the Bez_bound
  // for the originator.
  _origs.push_back (Originator (B, xid, t0));

  // Set the point coordinates.
  const Alg_point_2   p = B (t0);
  
  p_alg_x = new Algebraic (p.x());
  p_alg_y = new Algebraic (p.y());
  
  // Set the bounding box for this point, by converting x, y  to two ranges
  // of doubles.
  Nt_traits                         nt_traits;
  const std::pair<double, double>&  x_bnd = 
                                        nt_traits.double_interval (*p_alg_x);
  const std::pair<double, double>&  y_bnd = 
                                        nt_traits.double_interval (*p_alg_y);

  _bbox.min_x = x_bnd.first;
  _bbox.max_x = x_bnd.second;
  _bbox.min_y = y_bnd.first;
  _bbox.max_y = y_bnd.second;
}

// ---------------------------------------------------------------------------
// Compare the x-coordinate to this of the given point.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::compare_x
         (Self& pt,
          Bezier_cache& cache)
{
  // Handle rational points first.
  if (is_rational() && pt.is_rational())
  {
    return (CGAL::compare (*p_rat_x, *(pt.p_rat_x)));
  }

  // Try to handle the comparison using the x-range of the bounding boxes.
  bool    can_refine1 = true;
  bool    can_refine2 = true;

  do
  {
    // Compare the x-ranges of the bounding boxes.
    if (CGAL::compare (_bbox.max_x, pt._bbox.min_x) == SMALLER)
      return (SMALLER);
    
    if (CGAL::compare (_bbox.min_x, pt._bbox.max_x) == LARGER)
      return (LARGER);

    // Check if only one of the points is exactly represented.
    Nt_traits                         nt_traits;

    if (is_exact())
    {
      // Compare the exact x-coordinate to pt's bounding box.
      if (CGAL::compare (*p_alg_x,
                         nt_traits.convert (pt._bbox.min_x)) == SMALLER)
        return (SMALLER);
      if (CGAL::compare (*p_alg_x,
                         nt_traits.convert (pt._bbox.max_x)) == LARGER)
        return (LARGER);
    }

    if (pt.is_exact())
    {
      // Compare the bounding box to pt's exact x-coordinate.
      if (CGAL::compare (nt_traits.convert (_bbox.max_x),
                         *(pt.p_alg_x)) == SMALLER)
        return (SMALLER);
      if (CGAL::compare (nt_traits.convert (_bbox.min_x),
                         *(pt.p_alg_x)) == LARGER)
        return (LARGER);
    }

    // Try to refine the representation of the points.
    // If we cannot refine any more, compute the exact representation of the
    // point.
    if (! is_exact())
    {
      if (can_refine1)
        can_refine1 = this->_refine();
      
      if (! can_refine1)
        this->_make_exact (cache);
    }
    else
    {
      can_refine1 = false;
    }

    if (! pt.is_exact())
    {
      if (can_refine2)
        can_refine2 = pt._refine();
      
      if (! can_refine2)
        pt._make_exact (cache);
    }
    else
    {
      can_refine2 = false;
    }

  } while (can_refine1 || can_refine2);

  // If we reached here, we have an exact representation of both points, and
  // we simply have to compare their x-coordinates.
  return (CGAL::compare (*p_alg_x, *(pt.p_alg_x)));
}

// ---------------------------------------------------------------------------
// Compare the two points xy-lexicographically.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::compare_xy
         (Self& pt,
          Bezier_cache& cache)
{
  // First compare the x-coordinates.
  const Comparison_result   res_x = this->compare_x (pt,
                                                     cache);

  if (res_x != EQUAL)
    return (res_x);

  // Handle rational points separately.
  if (is_rational() && pt.is_rational())
  {
    return (CGAL::compare (*p_rat_y, *(pt.p_rat_y)));
  }

  // If we reached here, the two point should have already been computed in
  // an exact manner (otherwise we could not have determined the equality of
  // their x-coordinates).
  CGAL_assertion (is_exact() && pt.is_exact());
  return (CGAL::compare (*p_alg_y, *(pt.p_alg_y)));
}

// ---------------------------------------------------------------------------
// Determine the vertical position of the point with respect to a given curve.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
Comparison_result
_Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::vertical_position
        (const Control_points& cp,
         const BoundNT& t_min,
         const BoundNT& t_max)
{
  // Initialize a list of subcurves of the original curve.
  // We start from the entire curve.
  Bounding_traits       bound_tr;
  std::list<Subcurve>   subcurves;

  subcurves.push_back (Subcurve (cp, 0, 1));

  bool                  can_refine_pt = true;
  bool                  can_refine_cv = true;
  Bez_point_bbox        scv_bbox;

  while (can_refine_pt || can_refine_cv)
  {
    // Go over the list of subcurves and consider only those lying in the
    // given [t_min, t_max] bound.
    typename std::list<Subcurve>::iterator  iter = subcurves.begin();
    bool                                    was_split = false;
    bool                                    is_fully_in_t_range;

    while (iter != subcurves.end())
    {
      if (CGAL::compare (iter->t_max, t_min) == SMALLER ||
          CGAL::compare (iter->t_min, t_max) == LARGER)
      {
        // Subcurve out of bounds - erase it and continue to next subcurve.
        subcurves.erase(iter++);
        continue;
      }

      // Construct the bounding box of the subcurve and compare it to
      // the bounding box of the point.
      bound_tr.construct_bbox (iter->ctrl, scv_bbox);

      if (! _bbox.overlaps_x (scv_bbox))
      {
        // Subcurve out of x bounds - erase it and continue to next subcurve.
        subcurves.erase(iter++);
        continue;
      }

      is_fully_in_t_range = (CGAL::compare (iter->t_min, t_min) != SMALLER) &&
                            (CGAL::compare (iter->t_max, t_max) != LARGER);

      if (_bbox.overlaps (scv_bbox) || ! is_fully_in_t_range)
      {
        // \todo This is a special case of subdividing the curve
        //       not as part of an Originator. Think again if the can_refine
        //       and de Casteljau should not be different here!!
        can_refine_cv = bound_tr.can_refine (iter->ctrl,
                                             iter->t_min, iter->t_max);

        if (! can_refine_pt && ! can_refine_cv)
        {
          // It is not possible to refine the point or the subcurve anymore:
          return (EQUAL);
        }

        if (! can_refine_cv)
        {
          // We are not able to refine the curve. However, we keep it in the
          // list as we are able to refine the point - so in the future we
          // can compare the refined point to this curve.
          ++iter;
          continue;
        }

        // Subdivide the current subcurve and replace it with the two
        // resulting subcurves (note that we insert the two new subcurves
        // before iter and remove the current one, pointed by iter).
        const BoundNT   t_mid = (iter->t_min + iter->t_max) / 2;
        Subcurve        scv_left (iter->t_min, t_mid);
        Subcurve        scv_right (t_mid, iter->t_max);

        bisect_control_polygon_2 (iter->ctrl.begin(), iter->ctrl.end(),
                                  std::back_inserter(scv_left.ctrl),
                                  std::front_inserter(scv_right.ctrl));
        
        subcurves.insert (iter, scv_left);
        subcurves.insert (iter, scv_right);
        subcurves.erase(iter++);

        was_split = true;
        continue;
      }
      else
      {
        // We found a subcurve whose x-range contains our point, but whose
        // bounding box is disjoint from the bounding box of the point.
        // We can therefore compare the y-positions of the bounding boxes.
        CGAL_assertion (! _bbox.overlaps (scv_bbox) &&
                        _bbox.overlaps_x (scv_bbox) &&
                        is_fully_in_t_range);

        return (CGAL::compare (_bbox.max_y, scv_bbox.max_y));
      }

      // If we got here without entering one of the clauses above,
      // then iter has not been incremented yet.
      ++iter; 
    }

    // If we reached here without splitting a subcurve, then we have a
    // subcurve whose bbox (given by scv_bbox) is totally below or totally
    // above the bounding box of the point.
    if (! was_split)
      break;

    // Try to refine the bounding box of the point.
    if (! is_exact())
      can_refine_pt = _refine();

    if (! can_refine_pt && ! can_refine_cv)
    {
      // It is not possible to refine the point or the subcurve anymore:
      return (EQUAL);
    }
  }

  // If we reached here, we cannot refine any more:
  return (EQUAL);
}

// ---------------------------------------------------------------------------
// Refine the bounds of the point.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
bool _Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_refine ()
{
  // Get the first originator and consider the point type.
  Bounding_traits  bound_tr;
  Originator&      orig1 = *(_origs.begin());

  if (orig1.point_bound().type == Bez_point_bound::VERTICAL_TANGENCY_PT)
  {
    CGAL_assertion(_origs.size() == 1);

    // Refine the vertical tangency point.
    typename Bounding_traits::Vertical_tangency_point vpt (orig1.point_bound(),
                                                           _bbox);
    typename Bounding_traits::Vertical_tangency_point ref_vpt;

    bound_tr.refine_vertical_tangency_point (vpt,
                                             ref_vpt);

    if (! ref_vpt.bound.can_refine)
    {
      // Indicate that it was not possible to refine the point.
      return (false);
    }

    // Update the originator and the bounding box of the point.
    orig1.update_point_bound (ref_vpt.bound);
    _bbox = ref_vpt.bbox;

    if (ref_vpt.bound.type == Bez_point_bound::RATIONAL_PT)
    {
      // In this case the point is exactly computed.
      Nt_traits                         nt_traits;

      p_rat_x = new Rational (ref_vpt.bbox.min_x);
      p_alg_x = new Algebraic (nt_traits.convert (*p_rat_x));
      p_rat_y = new Rational (ref_vpt.bbox.min_y);
      p_alg_y = new Algebraic (nt_traits.convert (*p_rat_y));

      orig1.set_parameter (nt_traits.convert (ref_vpt.bound.t_min));

      // Make sure the point is marked as a vertical tangency point.
      ref_vpt.bound.type = Bez_point_bound::VERTICAL_TANGENCY_PT;
      orig1.update_point_bound (ref_vpt.bound);
    }

    return (true);
  }
  
  if (orig1.point_bound().type == Bez_point_bound::INTERSECTION_PT)
  {
    CGAL_assertion(_origs.size() == 2);

    // Obtain the other curve that originates the intersection point and use
    // it to refine its reprsentation.
    Orig_iter    org_it = _origs.begin();
    ++org_it;
    Originator&  orig2 = *org_it;

    typename Bounding_traits::Intersection_point  ipt (orig1.point_bound(),
                                                       orig2.point_bound(),
                                                       _bbox);
    typename Bounding_traits::Intersection_point  ref_ipt;

    bound_tr.refine_intersection_point (ipt, ref_ipt);

    if (! ref_ipt.bound1.can_refine || ! ref_ipt.bound2.can_refine)
    {
      // Indicate that it was not possible to refine the point.
      return (false);
    }

    // Update the originators and the bounding box of the point.
    orig1.update_point_bound (ref_ipt.bound1);
    orig2.update_point_bound (ref_ipt.bound2);
    _bbox = ref_ipt.bbox;

    if (ref_ipt.bound1.type == Bez_point_bound::RATIONAL_PT ||
        ref_ipt.bound2.type == Bez_point_bound::RATIONAL_PT)
    {
      // In this case the point is exactly computed.
      Nt_traits                         nt_traits;

      p_rat_x = new Rational (ref_ipt.bbox.min_x);
      p_alg_x = new Algebraic (nt_traits.convert (*p_rat_x));
      p_rat_y = new Rational (ref_ipt.bbox.min_y);
      p_alg_y = new Algebraic (nt_traits.convert (*p_rat_y));

      orig1.set_parameter (nt_traits.convert (ref_ipt.bound1.t_min));
      orig2.set_parameter (nt_traits.convert (ref_ipt.bound2.t_min));

      // Make sure the point is marked as an intersection point.
      ref_ipt.bound1.type = Bez_point_bound::INTERSECTION_PT;
      ref_ipt.bound2.type = Bez_point_bound::INTERSECTION_PT;
      orig1.update_point_bound (ref_ipt.bound1);
      orig2.update_point_bound (ref_ipt.bound2);
    }

    return (true);
  }

  // If we reached here, the point is computed in an exact manner, so there
  // is not need to refine its approximation.
  return (false);
}

// ---------------------------------------------------------------------------
// Make sure the originator parameters fit the bound box.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
void _Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_fit_to_bbox ()
{
  // Get the bounding box of the point.
  const typename Bounding_traits::NT    x_min = _bbox.min_x;
  const typename Bounding_traits::NT    x_max = _bbox.max_x;
  const typename Bounding_traits::NT    y_min = _bbox.min_y;
  const typename Bounding_traits::NT    y_max = _bbox.max_y;

  // Go over all originators, and make sure that the bounding box of the
  // control polygon of each originator is not larger than the bounding box
  // of the point.
  Bounding_traits   bound_tr;
  Bez_point_bbox    org_bbox;
  Orig_iter         org_it;
  bool              refined;
  bool              all_fit;

  do
  {
    all_fit = true;
    refined = false;

    for (org_it = _origs.begin(); org_it != _origs.end(); ++org_it)
    {
      // Skip rational points.
      if (org_it->point_bound().type == Bez_point_bound::RATIONAL_PT ||
          org_it->point_bound().ctrl.empty())
      {
        continue;
      }

      // Get the bounding box of the control polygon.
      // In case this bounding box is larger than the bounding box
      // of the point, we refine the point (hence refine the control polygon;
      // note however we keep the original bounding box of the point, in order
      // to avoid an infinitive loop here ...)
      bound_tr.construct_bbox (org_it->point_bound().ctrl,
                               org_bbox);

      if (CGAL::compare (org_bbox.min_x, x_min) == SMALLER ||
          CGAL::compare (org_bbox.max_x, x_max) == LARGER ||
          CGAL::compare (org_bbox.min_y, y_min) == SMALLER ||
          CGAL::compare (org_bbox.max_y, y_max) == LARGER)
      {
        all_fit = false;
        refined |= _refine();
      }
    }
  } while (! all_fit && refined);

  return;
}

// ---------------------------------------------------------------------------
// Compute the point in an exact manner.
//
template <class RatKer, class AlgKer, class NtTrt, class BndTrt>
void _Bezier_point_2_rep<RatKer, AlgKer, NtTrt, BndTrt>::_make_exact
        (Bezier_cache& cache)
{                
  if (is_exact())
    return;

  // Check if the point is a vertical tangency point of the originator.
  Nt_traits            nt_traits;

  if (_origs.size() == 1)
  {
    Orig_iter   org_it = _origs.begin();

    CGAL_assertion (org_it->point_bound().type ==
                    Bez_point_bound::VERTICAL_TANGENCY_PT);

    // Compute (using the cache) the vertical tangency parameters of
    // the current curve.
    const typename Bezier_cache::Vertical_tangency_list&          vt_list =
      cache.get_vertical_tangencies (org_it->curve().id(),
                                     org_it->curve().x_polynomial(),
                                     org_it->curve().x_norm());
    typename Bezier_cache::Vertical_tangency_iter                 vt_it;

    // Look for a parameter within the range of the bounding interval.
    const Algebraic&  t_min = nt_traits.convert (org_it->point_bound().t_min);
    const Algebraic&  t_max = nt_traits.convert (org_it->point_bound().t_max);

    for (vt_it = vt_list.begin(); vt_it != vt_list.end(); ++vt_it)
    {
      if (CGAL::compare (*vt_it, t_min) != SMALLER &&
          CGAL::compare (*vt_it, t_max) != LARGER)
      {
        // Update the originator.
        Originator&   orig = const_cast<Originator&> (*org_it);

        orig.set_parameter (*vt_it);
        
        // Evaluate the curve at the given parameter value.
        const Alg_point_2&   p = org_it->curve() (*vt_it);

        p_alg_x = new Algebraic (p.x());
        p_alg_y = new Algebraic (p.y());

        // Update the bounding box.
        Nt_traits                         nt_traits;
        const std::pair<double, double>&  x_bnd = 
                                        nt_traits.double_interval (*p_alg_x);
        const std::pair<double, double>&  y_bnd = 
                                        nt_traits.double_interval (*p_alg_y);

        _bbox.min_x = x_bnd.first;
        _bbox.max_x = x_bnd.second;
        _bbox.min_y = y_bnd.first;
        _bbox.max_y = y_bnd.second;

        return;
      }
    }

    // We should never reach here:
    CGAL_error();
  }

  // In this case the point is an intersection between two originating
  // curves. Compute the intersections between these curves in the parameter
  // space.
  CGAL_assertion (_origs.size() == 2);

  Orig_iter    org_it1 = _origs.begin();
  Orig_iter    org_it2 = org_it1;

  ++org_it2;

  if (org_it1->curve().id() > org_it2->curve().id())
  {
    // Make sure that org_it1 refers to a curve with smaller ID than org_it2.
    --org_it2;
    ++org_it1;
  }

  Originator&  orig1 = const_cast<Originator&> (*org_it1);
  Originator&  orig2 = const_cast<Originator&> (*org_it2);
  bool         do_ovlp;

  const typename Bezier_cache::Intersection_list&           intr_list =
    cache.get_intersections (orig1.curve().id(),
                             orig1.curve().x_polynomial(),
                             orig1.curve().x_norm(),
                             orig1.curve().y_polynomial(),
                             orig1.curve().y_norm(),
                             orig2.curve().id(),
                             orig2.curve().x_polynomial(),
                             orig2.curve().x_norm(),
                             orig2.curve().y_polynomial(),
                             orig2.curve().y_norm(),
                             do_ovlp);
  typename Bezier_cache::Intersection_iter                  intr_it;
                             
  CGAL_assertion (! do_ovlp);
                    
  // Look for a parameter pair within the ranges of the bounding intervals.
  const Algebraic      s_min = nt_traits.convert (orig1.point_bound().t_min);
  const Algebraic      s_max = nt_traits.convert (orig1.point_bound().t_max);
  const Algebraic      t_min = nt_traits.convert (orig2.point_bound().t_min);
  const Algebraic      t_max = nt_traits.convert (orig2.point_bound().t_max);

  for (intr_it = intr_list.begin(); intr_it != intr_list.end(); ++intr_it)
  {
    if (CGAL::compare (intr_it->s, s_min) != SMALLER &&
        CGAL::compare (intr_it->s, s_max) != LARGER &&
        CGAL::compare (intr_it->t, t_min) != SMALLER &&
        CGAL::compare (intr_it->t, t_max) != LARGER)
    {
      // Update the originators.
      orig1.set_parameter (intr_it->s);
      orig2.set_parameter (intr_it->t);

      // Set the exact point coordinates.
      p_alg_x = new Algebraic (intr_it->x);
      p_alg_y = new Algebraic (intr_it->y);

      // Update the bounding box.
      const std::pair<double, double>&  x_bnd = 
                                          nt_traits.double_interval (*p_alg_x);
      const std::pair<double, double>&  y_bnd = 
                                          nt_traits.double_interval (*p_alg_y);

      _bbox.min_x = x_bnd.first;
      _bbox.max_x = x_bnd.second;
      _bbox.min_y = y_bnd.first;
      _bbox.max_y = y_bnd.second;

      return;
    }
  }

  // We should never reach here:
  CGAL_error();
}

} //namespace CGAL

#endif
