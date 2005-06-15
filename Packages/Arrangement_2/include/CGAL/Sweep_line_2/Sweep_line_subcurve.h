// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>,
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_SUBCURVE_H
#define CGAL_SWEEP_LINE_SUBCURVE_H


#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_traits.h>
#include <CGAL/Sweep_line_2/Sweep_line_rb_tree.h>
#include <CGAL/assertions.h>



CGAL_BEGIN_NAMESPACE




/*! @class Sweep_line_subcurve
 *
 * This is a wrapper class to Curve_2 in the traits class, that contains
 * data that is used when applying the sweep algorithm on a set of curves.
 *
 * The information contained in this class is:
 * - the curve itself
 * - two points which are the source and target of the curve. We keep 
 *   the points in order to avoid many calls to the source() and 
 *   target() methods of the traits class 
 * - an indication for the direction of the curve (source point 
 *   is left or right to the target point). 
 * - a reference point that is used when comparing the y values of 
 *   any two curves. Since the curves are inserted in to a balanced 
 *   tree, and at any given time they are sorted on the status line, 
 *   and since their order may change, depending on the position of 
 *   the status line, we need to be able to compare the curves 
 *   relative to a point that will produce a correct answer.
 * - a reference to the last event point on the curve that was already 
 *   handled and also the curve that is the portion of the original 
 *   curve that is to the right of the last event point mentioned. 
 *   This is stored to avoid unneccesary splits of the curve.
 *
 */
 
template<class SweepLineTraits_2>
class Sweep_line_subcurve
{
public:

  typedef SweepLineTraits_2                               Traits;
  typedef typename Traits::Point_2                        Point_2;
  typedef typename Traits::X_monotone_curve_2             X_monotone_curve_2;

  typedef Sweep_line_subcurve<Traits> Self;
  typedef Status_line_curve_less_functor<Traits, Self>    StatusLineCurveLess;
  typedef Red_black_tree<Self*,
                         StatusLineCurveLess,
                         CGAL_ALLOCATOR(int)>             StatusLine;
  typedef typename StatusLine::iterator                   StatusLineIter;

  typedef Sweep_line_event<Traits, Self>                  Event;




  Sweep_line_subcurve() : m_overlap_subcurve(NULL),
                          m_orig_subcurve1(NULL),
                          m_orig_subcurve2(NULL)
  {}

  Sweep_line_subcurve(const X_monotone_curve_2 &curve);

  void init(const X_monotone_curve_2 &curve)
  {
    m_lastCurve = curve;
  }

  ~Sweep_line_subcurve() {}

 
  /*!
    @return a const reference to the last intersecing curve so far
  */
  const X_monotone_curve_2 &get_last_curve() const { 
    return m_lastCurve; 
  }

  
  /*! 
    updates the last intersecting curve so far.
    @param cv a reference to the curve
  */
  void set_last_curve(const X_monotone_curve_2 &cv) { 
    m_lastCurve = cv; 
  }

  
  
  template<class SweepEvent>
  bool is_end_point(const SweepEvent* event)const
  {
    return m_right_event == (Event*)event;
  }
 

  Point_2 get_right_end() const {
   return traits()->construct_max_vertex_2_object()(m_lastCurve);
  }

  
  Point_2 get_left_end() const {
    return traits()->construct_min_vertex_2_object()(m_lastCurve);  
  }

   Event* get_left_event() const
   {
     return m_left_event;
   }

   Event* get_right_event() const
   {
     return m_right_event;
   }

   template<class SweepEvent>
   void set_left_event(SweepEvent* event)
   {
     m_left_event =(Event*)event;
   }

   template<class SweepEvent>
   void set_right_event(SweepEvent* event)
   {
     m_right_event = (Event*)event;
   }



  void set_hint(StatusLineIter hint) 
  {
    m_hint = hint;
  }

  StatusLineIter get_hint() const 
  {
    return m_hint;
  }

  
  void set_overlap_subcurve(Self* overlap_sc)
  {
    m_overlap_subcurve = overlap_sc;
  }

  Self* get_overlap_subcurve()
  {
    return m_overlap_subcurve;
  }

  void set_orig_subcurve1(Self* orig_subcurve1)
  {
    m_orig_subcurve1 = orig_subcurve1;
  }

  Self* get_orig_subcurve1()
  {
    return m_orig_subcurve1;
  }

  void set_orig_subcurve2(Self* orig_subcurve2)
  {
    m_orig_subcurve2 = orig_subcurve2;
  }

  Self* get_orig_subcurve2()
  {
    return m_orig_subcurve2;
  }

  Self* getSubcurve()
  {
    Self* ptr  = m_overlap_subcurve;
    Self* prev = ptr;
    while(ptr != NULL)
    {
      prev = ptr;
      ptr = ptr->m_overlap_subcurve;
    }
    return prev;
  }

  template<class SweepEvent>
  Self* clip(const SweepEvent* e) 
  {
    if(get_right_event() != (Event*)e)
    {
      X_monotone_curve_2 dummy;
      traits()->split_2_object()(m_lastCurve, e->get_point(), dummy,m_lastCurve);

      return this;
    }
    if(m_orig_subcurve1)
    {
      Self* res;
      if((res = m_orig_subcurve1->clip(e)) == NULL)
        return m_orig_subcurve2->clip(e);
      return res;
    }
    return NULL;  
  }



  bool is_parent(Self* parent)
  {
    Self* ptr = m_overlap_subcurve;
    while(ptr != NULL)
    {
      if(ptr == parent)
        return true;
      ptr = ptr->m_overlap_subcurve;
    }
    return false;
  }

  unsigned int overlap_depth()
  {
    if(! m_orig_subcurve1)
      return 1;
    else
      return 1 + max(m_orig_subcurve1->overlap_depth(),
                     m_orig_subcurve2->overlap_depth());
  }


#ifndef NDEBUG
  void Print() const;
#endif



  
 
  protected:
  /*! the portion of the curve to the right of the last event point 
      on the curve */
  X_monotone_curve_2 m_lastCurve;

  Event* m_left_event;

  Event* m_right_event;
  
  /*! iterator at the Y-structure for the Subcurve */
  StatusLineIter m_hint;

  //the three below memvers relevant only with overlaps

  Self *m_overlap_subcurve;

  Self *m_orig_subcurve1;

  Self *m_orig_subcurve2;

 
  protected:

  Traits* traits() const
  {
    return Sweep_line_traits<Traits>::get_traits();
  }

  private:
  unsigned int max(unsigned int a, unsigned int b)
  {
    return (a >= b ? a : b);
  }

};

template<class SweepLineTraits_2>
inline Sweep_line_subcurve<SweepLineTraits_2>::
Sweep_line_subcurve(const X_monotone_curve_2 &curve): 
  m_overlap_subcurve(NULL),
  m_orig_subcurve1(NULL),
  m_orig_subcurve2(NULL)
{ 
  m_lastCurve = curve;
}

#ifndef NDEBUG
template<class SweepLineTraits_2>
void 
Sweep_line_subcurve<SweepLineTraits_2>::
Print() const
{
  std::cout << "Curve " << this << "  (" << m_lastCurve << ") " << std::endl;
}

#endif



CGAL_END_NAMESPACE

#endif
