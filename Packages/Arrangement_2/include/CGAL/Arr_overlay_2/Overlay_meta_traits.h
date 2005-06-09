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


CGAL_BEGIN_NAMESPACE


/*!
 */
template <class Halfedge_handle_red, class Halfedge_handle_blue>
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

  Halfedge_handle_red     m_red_halfedge_handle;
  Halfedge_handle_blue    m_blue_halfedge_handle;


public:

  Curve_info() : m_red_halfedge_handle(),
                 m_blue_halfedge_handle()
  {}



  Curve_info(Halfedge_handle_red he1, Halfedge_handle_blue he2) :
    m_red_halfedge_handle(he1),
    m_blue_halfedge_handle(he2)
  {
    CGAL_assertion(he1 != Halfedge_handle_red() ||
                   he2 != Halfedge_handle_blue());
  }

 
  Halfedge_handle_red get_red_halfedge_handle()  const 
  { 
    return m_red_halfedge_handle;  
  }

  Halfedge_handle_blue get_blue_halfedge_handle() const 
  {
    return m_blue_halfedge_handle; 
  }


  Color get_color() const 
  {
    Halfedge_handle_red     null_red_he;
    Halfedge_handle_blue    null_blue_he;
    if(m_red_halfedge_handle != null_red_he &&
       m_blue_halfedge_handle == null_blue_he)
      return RED;

    if(m_blue_halfedge_handle != null_blue_he &&
       m_red_halfedge_handle == null_red_he)
      return BLUE;

    //overlap, the PURPLE color will be returned
    CGAL_assertion(m_red_halfedge_handle != null_red_he && 
                   m_blue_halfedge_handle != null_blue_he);
    return PURPLE;
  }

};



template <class Traits,
          class Halfedge_handle_red_,
          class Halfedge_handle_blue_>
class Overlay_meta_traits : public Traits
{
public:

  typedef Halfedge_handle_red_                      Halfedge_handle_red;
  typedef Halfedge_handle_blue_                     Halfedge_handle_blue;
  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                  Point_2;
  typedef typename Traits::Intersect_2              Base_Intersect_2;
  typedef typename Traits::Split_2                  Base_Split_2;
  
  typedef Curve_info<Halfedge_handle_red,
                     Halfedge_handle_blue>          Curve_info;
  typedef typename Curve_info::Color                Color;


  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2 : public Base_X_monotone_curve_2 
  {
  public:
    typedef typename Traits::X_monotone_curve_2     Base;
    typedef typename Traits::Point_2                Point_2;

    friend class Overlay_meta_traits<Traits,
                                     Halfedge_handle_red, 
                                     Halfedge_handle_blue>;
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

    My_X_monotone_curve_2(const Base&cv, Halfedge_handle_red  he1,
                                         Halfedge_handle_blue he2): 
      Base(cv),
      m_curve_info(he1, he2)
    {}

    const Curve_info& get_curve_info() const { return m_curve_info; }

    void set_curve_info(const Curve_info& cv_info ) {m_curve_info = cv_info;}

    Halfedge_handle_red get_red_halfedge_handle() const
    {
      return m_curve_info.get_red_halfedge_handle();
    }

    Halfedge_handle_blue get_blue_halfedge_handle() const
    {
      return m_curve_info.get_blue_halfedge_handle();
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
      if (cv1.get_color() == cv2.get_color())
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
          Halfedge_handle_red        red_he;
          Halfedge_handle_blue       blue_he;

          if(cv1.get_color() == Curve_info::RED)
          {
            red_he = cv1.get_red_halfedge_handle();

            // overlap can occur only between curves from a different color
            CGAL_assertion(cv2.get_color() == Curve_info::BLUE);
            blue_he = cv2.get_blue_halfedge_handle();
          }
          else
          {
            CGAL_assertion(cv1.get_color() == Curve_info::BLUE &&
                           cv2.get_color() == Curve_info::RED);

            red_he = cv2.get_red_halfedge_handle();
            blue_he = cv1.get_blue_halfedge_handle();
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
                       const X_monotone_curve_2& cv2) const
  {
    return  (cv1.get_color() == cv2.get_color());
  }


};


CGAL_END_NAMESPACE

#endif
