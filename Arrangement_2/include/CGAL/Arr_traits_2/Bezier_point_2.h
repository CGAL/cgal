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

#ifndef CGAL_BEZIER_POINT_2_H
#define CGAL_BEZIER_POINT_2_H

/*! \file
 * Header file for the _Bezier_point_2 class.
 */

#include <CGAL/Arr_traits_2/Bezier_curve_2.h>
#include <CGAL/Arr_traits_2/Bezier_cache.h>
#include <CGAL/Handle_for.h>
#include <list>
#include <ostream>

CGAL_BEGIN_NAMESPACE

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

    Curve_2             _curve;     /*! The originating curve. */
    Bez_point_bound     _bpb;       /*! Bounding information for the
                                        point: bouding control polygon,
                                        point type, etc. */
    Algebraic          *p_t;        /*! The algebraic parameter for the
                                        point (if available). */

  public:

    /*! Constructor, given an exact algebraic representation. */
    Originator(const Curve_2& c, const Algebraic& t) :
      _curve(c),
      p_t (NULL)
    {
      set_parameter (t);
    }

    /*! Constructor with bounding information and no exact representation. */
    Originator (const Curve_2& c, const Bez_point_bound& bpb) :
      _curve (c),
      _bpb (bpb),
      p_t (NULL)
    {}

    /*! Copy constructor. */
    Originator (const Originator& other) :
      _curve (other._curve),
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

      // Update the Bez_point_bound, probably by converting t to
      // an interval of doubles and setting _bpb accordingly.
      Nt_traits                         nt_traits;
      const std::pair<double, double>&  t_bnd = nt_traits.double_interval (t);

      _bpb.t_min = t_bnd.first;
      _bpb.t_max = t_bnd.second;

      return;
    }
  };

  /*! \struct Subcurve
   * Auxilary structure for the vertical_position() function.
   */
  typedef typename Bounding_traits::Control_point_vec   Control_point_vec;
  typedef typename Bounding_traits::NT                  BoundNT;

  struct Subcurve
  {
    Control_point_vec   cp;
    BoundNT             l;
    BoundNT             r;

    /*! Constructor. */
    Subcurve (const Control_point_vec& _cp,
              const BoundNT& _l, 
              const BoundNT& _r) :
      cp (_cp),
      l (_l),
      r (_r)
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
  _Bezier_point_2_rep (const Algebraic& x, const Algebraic& y) : 
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
   * Constructor given an originating curve and a rational t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, const Rational& t0);

  /*!
   * Constructor given an originating curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2_rep (const Curve_2& B, const Algebraic& t0);

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
  Comparison_result vertical_position (const Control_point_vec& cp,
                                       const BoundNT& t_min,
                                       const BoundNT& t_max);

private:

  /*!
   * Refine the bounds of the point.
   * \return Whether it was possible to further refine the point.
   */
  bool _refine ();

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
   * Constructor with coordinates.
   */
  _Bezier_point_2 (const Algebraic& x, const Algebraic& y) :
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
   * Constructor given an originating curve and an algebraic t0 value.
   * \pre t0 must be between 0 and 1.
   */
  _Bezier_point_2 (const Curve_2& B, const Algebraic& t0) :
    Bpt_handle (Bpt_rep (B, t0))
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
   * Compare the two points xy-lexicographically.
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
      (const typename Bounding_traits::Control_point_vec& cp,
       const typename Bounding_traits::NT& t_min,
       const typename Bounding_traits::NT& t_max) const
  {
    Bpt_rep&             rep = const_cast<Bpt_rep&> (_rep());
    return (rep.vertical_position(cp, t_min, t_max));
  }

  /*!
   * Get the originator of the point that is associates with the given curve.
   * \param B The Bezier curve.
   * \return An iterator pointing to the requested originator;
   *         or originators_end() if B is not an originator of the point.
   */
  // Iddo: this is a bit const-problematic since it should return
  //       const_iterator, currently Originator_iterator is typedefed to a
  //       const_iterator.
  //       (TODO - Originator_const_iterator and Originator_iterator)
  Originator_iterator get_originator(const Curve_2& B) const
  {
    // Scan the list of originators and look for B.
    typename Bpt_rep::Orig_const_iter     it = _rep()._origs.begin();
    typename Bpt_rep::Orig_const_iter     end = _rep()._origs.end();

    while (it != end)
    {
      if (B.is_same (it->curve()))
        break;
      ++it;
    }

    return it;
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

  // Iddo: workaround the ctr problems
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
  _origs.push_back (Originator(B, t0));

  // Evaluate the point coordinates.
  const Rat_point_2&  p = B (t0);
  Nt_traits           nt_traits;

  p_rat_x = new Rational (p.x());
  p_alg_x = new Algebraic (nt_traits.convert (*p_rat_x));
  p_rat_y = new Rational (p.y());
  p_alg_y = new Algebraic (nt_traits.convert (*p_rat_y));

  // Also set the bounding box for this point, by converting x, y
  // to two ranges of doubles.
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
  
  // Also set the bounding box for this point, by converting x, y
  // to two ranges of doubles.
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
    return (CGAL::compare(*p_rat_x, *(pt.p_rat_x)));
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
    if (is_exact())
    {
      // Compare the exact x-coordinate to pt's bounding box.
      if (CGAL::compare (*p_alg_x, Algebraic (pt._bbox.min_x)) == SMALLER)
        return (SMALLER);
      if (CGAL::compare (*p_alg_x, Algebraic (pt._bbox.max_x)) == LARGER)
        return (LARGER);
    }

    if (pt.is_exact())
    {
      // Compare the bounding box to pt's exact x-coordinate.
      if (CGAL::compare (Algebraic (_bbox.max_x), *(pt.p_alg_x)) == SMALLER)
        return (SMALLER);
      if (CGAL::compare (Algebraic (_bbox.min_x), *(pt.p_alg_x)) == LARGER)
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
        (const Control_point_vec& cp,
         const BoundNT& t_min,
         const BoundNT& t_max)
{
  Bounding_traits       bound_tr;
  std::list<Subcurve>   subcurves;

  subcurves.push_back(Subcurve(cp,0,1));

  // Iddo: maybe make here a comparison with an incrementor (e.g., loop 4 
  //       times first..)
  bool                  can_refine_pt = true;
  bool                  can_refine_cv = true;
  Bez_point_bbox        scv_bbox;

  while (can_refine_pt || can_refine_cv)
  {
    // Iddo: Implement Algebraic comparisons with scv_bbox if is_exact
    //       (or use the current bbox of the exact rep as is done now).

    // Go over the list of subcurves and consider only those lying in the
    // given [t_min, t_max] bound.
    typename std::list<Subcurve>::iterator  iter = subcurves.begin();
    bool                                    was_split = false;
    bool                                    is_fully_in_t_range;

    while (iter != subcurves.end())
    {
      if (CGAL::compare (iter->r, t_min) == SMALLER ||
          CGAL::compare (iter->l, t_max) == LARGER)
      {
        // Subcurve out of bounds - erase it and continue to next subcurve.
        subcurves.erase(iter++);
        continue;
      }

      // Construct the bounding box of the subcurve and compare it to
      // the bounding box of the point.
      // Iddo: In the future it might be not a bez_point_bbox
      //       but a bbox<Rational>.
      bound_tr.cp_bbox (iter->cp, scv_bbox);

      if (! _bbox.Overlaps_x(scv_bbox))
      {
        // Subcurve out of x bounds - erase it and continue to next subcurve.
        subcurves.erase(iter++);
        continue;
      }

      is_fully_in_t_range = (CGAL::compare (iter->l, t_min) != SMALLER) &&
                            (CGAL::compare (iter->r, t_max) != LARGER);

      if (_bbox.Overlaps(scv_bbox) || ! is_fully_in_t_range)
      {
        // Iddo: This is a special case of subdividing the curve
        //       not as part of an Originator. Think again if the can_refine
        //       and de Casteljau should not be different here!!
        can_refine_cv = bound_tr.can_refine (iter->cp, iter->l, iter->r);

        if (! can_refine_pt && ! can_refine_cv)
          // It is not possible to refine the point or the subcurve anymore:
          return (EQUAL);

        if (! can_refine_cv)
        {
          ++iter;
          continue;
        }

        // Subdivide the current subcurve and replace iter with the two
        // resulting subcurves.
        Control_point_vec  left_cp, right_cp;

        bisect_control_polygon_2 (iter->cp.begin(), iter->cp.end(),
                                  std::back_inserter(left_cp),
                                  std::front_inserter(right_cp));

        //bound_tr.DeCasteljau (iter->cp, 0.5, left, right);

        // Insert two new subcurves before iter and remove iter.
        BoundNT   t_mid = BoundNT(0.5) * (iter->l + iter->r);
        
        subcurves.insert (iter, Subcurve (left_cp, iter->l, t_mid));
        subcurves.insert (iter, Subcurve (right_cp, t_mid, iter->r));
        subcurves.erase(iter++);

        was_split = true;
        continue;
      }
      else
      {
        // We found a subcurve whose x-range contains our point, but whose
        // bounding box is disjoint from the bounding box of the point.
        // We can therefore compare the y-positions of the bounding boxes.
        CGAL_assertion (! _bbox.Overlaps (scv_bbox) &&
                        _bbox.Overlaps_x (scv_bbox) &&
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
      // It is not possible to refine the point or the subcurve anymore:
      return (EQUAL);
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

  if (orig1.point_bound().point_type == Bez_point_bound::VERTICAL_TANGENCY_PT)
  {
    CGAL_assertion(_origs.size() == 1);

    // Refine the vertical tangency point.
    std::pair<Bez_point_bound, Bez_point_bbox>  refined_tang_pt;

    bound_tr.refine_tangency_point (orig1.point_bound().bounding_polyline,
                                    orig1.point_bound().t_min, 
                                    orig1.point_bound().t_max,
                                    refined_tang_pt);

    if (! refined_tang_pt.first.can_refine)
      // Indicate that it was not possible to refine the point.
      return (false);

    // Update the originator and the bounding box of the point.
    orig1.update_point_bound (refined_tang_pt.first);
    _bbox = refined_tang_pt.second;
    return (true);
  }
  
  if (orig1.point_bound().point_type == Bez_point_bound::INTERSECTION_PT)
  {
    CGAL_assertion(_origs.size() == 2);

    // Obtain the other curve that originates the intersection point and use
    // it to refine its reprsentation.
    Orig_iter    org_it = _origs.begin();
    ++org_it;
    Originator&  orig2 = *org_it;

    typename Bounding_traits::Bound_pair   intersect_pt (orig1.point_bound(),
                                                         orig2.point_bound(),
                                                         _bbox);
    typename Bounding_traits::Bound_pair   refined_pt;

    bound_tr.refine_intersection_point (intersect_pt, refined_pt);

    if (! refined_pt.bound1.can_refine || ! refined_pt.bound2.can_refine)
      // Indicate that it was not possible to refine the point.
      return (false);

    // Update the originators and the bounding box of the point.
    orig1.update_point_bound (refined_pt.bound1);
    orig2.update_point_bound (refined_pt.bound2);
    _bbox = refined_pt.bbox;
    return (true);
  }

  // We should never reach here:
  //CGAL_error();
  return (false);
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

  /* Ron: For informational purposes ...
  std::cout << "***** MAKE EXACT (" << _origs.size() << " ORIGINATORS) *****"
            << std::endl;
  */

  // Check if the point is a vertical tangency point of the originator.
  if (_origs.size() == 1)
  {
    Orig_iter   org_it = _origs.begin();

    CGAL_assertion (org_it->point_bound().point_type ==
                    Bez_point_bound::VERTICAL_TANGENCY_PT);

    // Compute (using the cache) the vertical tangency parameters of
    // the current curve.
    const typename Bezier_cache::Vertical_tangency_list&          vt_list =
      cache.get_vertical_tangencies (org_it->curve().id(),
                                     org_it->curve().x_polynomial(),
                                     org_it->curve().x_norm());
    typename Bezier_cache::Vertical_tangency_iter                 vt_it;

    // Look for a parameter within the range of the bounding interval.
    const Algebraic  t_min (org_it->point_bound().t_min);
    const Algebraic  t_max (org_it->point_bound().t_max);

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
  Nt_traits            nt_traits;
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

CGAL_END_NAMESPACE

#endif
