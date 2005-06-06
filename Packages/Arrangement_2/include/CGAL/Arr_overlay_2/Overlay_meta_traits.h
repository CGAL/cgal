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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef OVERLAY_META_TRAITS_H
#define OVERLAY_META_TRAITS_H

#include <utility>

CGAL_BEGIN_NAMESPACE


/*!
 */
template <class Halfedge_handle>
class Curve_info
{

public:
  enum Color
  {
    RED,
    BLUE,
    PURPLE  //overlap
  };

private:

  Halfedge_handle    m_red_halfedge_handle;
  Halfedge_handle    m_blue_halfedge_handle;

  //hide default constructor
  Curve_info(){}

public:

  Curve_info() : m_red_halfedge_handle(),
                 m_blue_halfedge_handle(),
  {}

  Curve_info(Halfedge_handle he, Color color)
  {
    if(color == RED)
      m_red_halfedge_handle = he;
    else
    {
      CGAL_assertion(color == BLUE);
      m_blue_halfedge_handle = he;
    }
  }

  // constructor for overlap halfedges
  Curve_info(Halfedge_handle he1, Halfedge_handle he2) :
    m_red_halfedge_handle(he1),
    m_blue_halfedge_handle(he2)
  {}

 
  Halfedge_handle get_halfedge_handle() const 
  {
    Halfedge_handle null_he;
    if(m_red_halfedge_handle != null_he)
      return m_red_halfedge_handle;

    CGAL_assertion(m_blue_halfedge_handle != null_he);
    return m_blue_halfedge_handle;
  }

  //should be called only when there's an overlap
  std::pair<Halfedge_handle, Halfedge_handle> get_pair_halfedges() const
  {
    return std::make_object(m_red_halfedge_handle, m_blue_halfedge_handle);
  }

  Color get_color() const 
  {
    Halfedge_handle null_he;
    if(m_red_halfedge_handle != null_he && m_blue_halfedge_handle == null_he)
      return RED;

    if(m_blue_halfedge_handle != null_he && m_red_halfedge_handle == null_he)
      return BLUE;

    //overlap, the PURPLE color will be returned
    CGAL_assertion(m_red_halfedge_handle != null_he && 
                   m_blue_halfedge_handle != null_he);
    return PURPLE;
  }

};



template <class Traits, class Halfedge_handle>
class Overlay_meta_traits : public Traits
{
public:

  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                  Point_2;
  typedef typename Traits::Intersect_2              Base_Intersect_2;
  typedef typename Traits::Split_2                  Base_Split_2;
  typedef Curve_info<Halfedge_handle>               Curve_info;
  typedef typename Curve_info::Color                Color;


  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2 : public Base_X_monotone_curve_2 
  {
  public:
    typedef typename Traits::X_monotone_curve_2 Base;
    typedef typename Traits::Point_2            Point_2;

    friend class Overlay_meta_traits<Traits, Halfedge_handle>;
    friend class Intersect_2;

    My_X_monotone_curve_2(): Base(),
                             m_curve_info()
    {}

    My_X_monotone_curve_2(const Base& cv): Base(cv),
                                           m_curve_info()
    {}

    My_X_monotone_curve_2(const Base& cv, const Curve_info& info):
      Base(cv),
      m_curve_info(info)
    {}

    My_X_monotone_curve_2(const Base&cv, Halfedge_handle he, Color col):
      Base(cv),
      m_curve_info(he, col)
    {}

    // should be called only when there's an overlap
    My_X_monotone_curve_2(const Base&cv, Halfedge_handle he1,
                                         Halfedge_handle he2): 
      Base(cv),
      m_curve_info(he1, he2)
    {}

    const Curve_info& get_curve_info() const { return m_curve_info; }

    void set_curve_info(const Curve_info& cv_info ) {m_curve_info = cv_info;}

    Halfedge_handle get_halfedge_handle() const
    {
      return m_curve_info.get_halfedge_handle();
    }

    Color get_color() const { return m_curve_info.get_color(); }

    
    protected:
     Curve_info    m_curve_info;
  }; // nested class My_X_monotone_curve_2 - END

  
  typedef My_X_monotone_curve_2                     X_monotone_curve_2;

 
  class Intersect_2
  {
  private:

    Base_Intersect_2      m_base_intersect;

  public:
   
    /*! Constructor. */
    Intersect_2 (const Base_Intersect_2& base) :
        m_base_intersect (base)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      if (cv1.get_color() == cv2.get_color());
        return (oi); // the curves are disjoint-interior because they
                     // are already at the same Arrangement (have same color)

      OutputIterator oi_end = m_base_intersect(cv1, cv2, oi);

      // convert objects that are associated with Base_X_monotone_curve_2 to
      // the extenede X_monotone_curve_2 
      for(; oi != oi_end; ++oi)
      {
        Base_X_monotone_curve_2 overlap_cv;
        if(CGAL::assign(overlap_cv, *oi))
        {
          Halfedge_handle        red_he;
          Halfedge_handle        blue_he;

          if(cv1.get_color() == Curve_info::RED)
          {
            red_he = cv1.get_halfedge_handle();

            // overlap can occur only between curves from a different color
            CGAL_assertion(cv2.get_color() == Curve_info::BLUE);
            blue_he = cv2.get_halfedge_handle();
          }
          else
          {
            CGAL_assertion(cv1.get_color() == Curve_info::BLUE &&
                           cv2.get_color() == Curve_info::RED);

            red_he = cv2.get_halfedge_handle();
            blue_he = cv1.get_halfedge_handle();
          }

          X_monotone_curve_2 new_overlap_cv(overlap_cv, red_he, blue_he);
          *oi = make_object(new_overlap_cv);
        }
      }
      //return past-end iterator
      return oi_end;
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2(Traits::intersect_2_object()); 
  }


  class Split_2
  {
  private:
    Base_Split_2    m_base_split;

  public:

    /*! Constructor. */
    Split_2 (const Base_Split_2& base) :
        m_base_split (base)
    {}

    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      m_base_split(cv,p,c1,c2);
      c1.set_curve_info(cv.get_curve_info());
      c2.set_curve_info(cv.get_curve_info());
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2(Traits::split_2_object());
  }


  bool  are_same_color(const X_monotone_curve_2& cv1,
                       const X_monotone_curve_2& cv2)
  {
    return  (cv1.get_color() == cv2.get_color());
  }


};


CGAL_END_NAMESPACE

#endif
