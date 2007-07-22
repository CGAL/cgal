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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_OVERLAY_META_TRAITS_H
#define CGAL_OVERLAY_META_TRAITS_H


#include <CGAL/Object.h>

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


/*!
 */
template <class Vertex_handle_red, class Vertex_handle_blue>
class Point_info
{

public:
  enum Color
  {
    RED,
    BLUE,
    PURPLE  //overlap
  };

private:

  Object     m_red_obj;
  Object     m_blue_obj;


public:

  Point_info() : m_red_obj(),
                 m_blue_obj()
  {}


  Object& get_red_object()  
  { 
    return m_red_obj;  
  }

  Object&  get_blue_object()  
  {
    return m_blue_obj; 
  }

  const Object& get_red_object() const
  { 
    return m_red_obj;  
  }

  const Object&  get_blue_object() const
  {
    return m_blue_obj; 
  }

  bool is_red_object_null() const
  {
    return m_red_obj.is_empty();
  }

  bool is_blue_object_null() const
  {
    return m_blue_obj.is_empty();
  }

  void set_red_object(const Object& obj)
  {
    m_red_obj = obj;
  }

  void set_blue_object(const Object& obj)
  {
    m_blue_obj = obj;
  }

};



template <class Traits,
          class Arrangement1,
          class Arrangement2>
class Overlay_meta_traits : public Traits
{
public:

  typedef typename Arrangement1::Halfedge_const_handle 
                                                    Halfedge_handle_red;
  typedef typename Arrangement2::Halfedge_const_handle
                                                    Halfedge_handle_blue;

  typedef typename Arrangement1::Vertex_const_handle 
                                                    Vertex_handle_red;
  typedef typename Arrangement2::Vertex_const_handle
                                                    Vertex_handle_blue;

  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                  Base_Point_2;
  typedef typename Traits::Intersect_2              Base_Intersect_2;
  typedef typename Traits::Split_2                  Base_Split_2;
  typedef typename Traits::Construct_min_vertex_2   Base_Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2   Base_Construct_max_vertex_2;
  typedef typename Traits::Compare_xy_2             Base_Compare_xy_2;
  typedef typename Traits::Compare_y_at_x_2         Base_Compare_y_at_x_2;
  typedef typename Traits::Compare_y_at_x_right_2   Base_Compare_y_at_x_right_2;
  typedef typename Traits::Compare_x_2              Base_Compare_x_2;
  typedef typename Traits::Equal_2                  Base_Equal_2;
  typedef typename Traits::Has_boundary_category    Base_has_boundary_category;
  
  typedef CGAL::Curve_info<Halfedge_handle_red,
                           Halfedge_handle_blue>    Curve_info;

  typedef CGAL::Point_info<Vertex_handle_red,
                           Vertex_handle_blue>      Point_info;

  typedef typename Curve_info::Color                Color;

private:

  Traits*    m_base_traits;



public:

  Overlay_meta_traits()
  {}

  Overlay_meta_traits(Traits* base_tr) : m_base_traits(base_tr)
  {}

  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2   
  {
  public:

    typedef  Base_X_monotone_curve_2     Base;
    typedef  Base_Point_2                Point_2;

  protected:
    
    Base_X_monotone_curve_2 m_base_cv;
    Curve_info              m_curve_info;

  public:

    My_X_monotone_curve_2(): m_base_cv(),
                             m_curve_info()
    {}

    My_X_monotone_curve_2(const Base& cv): m_base_cv(cv),
                                           m_curve_info()
    {}

    My_X_monotone_curve_2(const Base& cv, const Curve_info& info):
      m_base_cv(cv),
      m_curve_info(info)
    {}

    My_X_monotone_curve_2(const Base&cv, Halfedge_handle_red  he1,
                                         Halfedge_handle_blue he2): 
      m_base_cv(cv),
      m_curve_info(he1, he2)
    {}

    const Base_X_monotone_curve_2& base_curve() const
    {
      return (m_base_cv);
    }

    Base_X_monotone_curve_2& base_curve()
    {
      return (m_base_cv);
    }

    operator const Base_X_monotone_curve_2&()const
    {
      return (m_base_cv);
    }

    operator Base_X_monotone_curve_2&()
    {
      return (m_base_cv);
    }

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
    
  }; // nested class My_X_monotone_curve_2 - END

  
  typedef My_X_monotone_curve_2                     X_monotone_curve_2;
 
  friend std::ostream& operator<< (std::ostream& os, const My_X_monotone_curve_2 & cv)
  {
    os << cv.base_curve();
    return (os);
  }



  class My_Point_2 
  {
  public:

    typedef typename Traits::Point_2    Base;

  protected:

    Base_Point_2  m_base_pt;
    Point_info    m_info;

  public:

    My_Point_2() {}

    My_Point_2(const Base& pt) : m_base_pt(pt),
                                 m_info()
    {}

    My_Point_2(const Base& pt, const Object& red, const Object& blue) : m_base_pt(pt)
    {
      m_info.set_red_object(red);
      m_info.set_blue_object(blue);
    }

    const Base_Point_2& base_point() const
    {
      return (m_base_pt);
    }

    Base_Point_2& base_point()
    {
      return (m_base_pt);
    }

    operator const Base_Point_2&()const
    {
      return (m_base_pt);
    }

    operator  Base_Point_2&()
    {
      return (m_base_pt);
    }

   
    Object& get_red_object()  
    { 
      return m_info.get_red_object();  
    }

    Object&  get_blue_object()  
    {
      return m_info.get_blue_object(); 
    }

    const Object& get_red_object() const
    { 
      return m_info.get_red_object();  
    }

    const Object&  get_blue_object() const
    {
      return m_info.get_blue_object(); 
    }

    bool is_red_object_null() const
    {
      return m_info.is_red_object_null();
    }

    bool is_blue_object_null() const
    {
      return m_info.is_blue_object_null();
    }

    void set_red_object(const Object& obj)
    {
      m_info.set_red_object(obj);
    }

    void set_blue_object(const Object& obj)
    {
      m_info.set_blue_object(obj);
    }

  };

  friend std::ostream& operator<< (std::ostream& os, const My_Point_2 & p)
  {
    os << p.base_point();
    return (os);
  }

  typedef My_Point_2   Point_2;
 
  class Intersect_2
  {
  private:

    Traits*      m_base_tr;

  public:
   
    /*! Constructor. */
    Intersect_2 (Traits* base) : m_base_tr (base)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi)
    {
      if (cv1.get_color() == cv2.get_color())
        return (oi); // the curves are disjoint-interior because they
                     // are already at the same Arrangement (have same color)
      
      
      if(cv1.get_color() == Curve_info::PURPLE ||
         cv2.get_color() == Curve_info::PURPLE)
         return (oi);

      // Compute the intersection points between the curves. Note that if
      // cv1 and cv2 are subcruves of x-monotone curves that had intersected
      // before the current point on the status line, we may get a filter
      // failure if we send the subcurve whose left endpoint is to the left
      // of the other curve - this is because their previously computed
      // intersection point p may be equal to the this left endpoint. As many
      // traits classes start by computing the intersection between the
      // supporting curves and then check whether the result is in the x-range
      // of both subcurves, this will result in a filter failure. However, if
      // we send cv1 first, then p is obviusly not in its x-range and there is
      // no filter failure.
      //
      //              / cv1
      //             /
      //            /
      //       ----+--
      //          /
      //         /
      //      p +------------- cv2
      //              ^
      //              |
      //              status line
      //
      // Note that we do not bother with curves whose left ends are unbounded,
      // since such curved did not intersect before.
      const std::pair<Base_Point_2, unsigned int>   *base_pt;
      const Base_X_monotone_curve_2                 *overlap_cv;
      bool                                           send_cv1_first = true;
      OutputIterator                                 oi_end;

      Boundary_in_x_2 inf_in_x;
      Boundary_in_y_2 inf_in_y;
      if (inf_in_x (cv1.base_curve(), MIN_END) == NO_BOUNDARY &&
          inf_in_y (cv1.base_curve(), MIN_END) == NO_BOUNDARY &&
          inf_in_x (cv2.base_curve(), MIN_END) == NO_BOUNDARY &&
          inf_in_y (cv2.base_curve(), MIN_END) == NO_BOUNDARY)
      {
        send_cv1_first =
          (m_base_tr->compare_xy_2_object()
            (m_base_tr->construct_min_vertex_2_object()(cv1.base_curve()),
             m_base_tr->construct_min_vertex_2_object()(cv2.base_curve())) == LARGER);
      }

      if (send_cv1_first)
        oi_end = m_base_tr->intersect_2_object()(cv1.base_curve(),
                                                 cv2.base_curve(), oi);
      else
        oi_end = m_base_tr->intersect_2_object()(cv2.base_curve(),
                                                 cv1.base_curve(), oi);

      // convert objects that are associated with Base_X_monotone_curve_2 to
      // the extenede X_monotone_curve_2 
      for(; oi != oi_end; ++oi)
      {
        base_pt = object_cast<std::pair<Base_Point_2, unsigned int> >(&(*oi));

        if (base_pt != NULL)
        {
          Object red_obj , blue_obj;
          if(cv1.get_color() == Curve_info::RED)
          {
            CGAL_assertion(cv2.get_color() == Curve_info::BLUE);
            red_obj = CGAL::make_object(cv1.get_red_halfedge_handle());
            blue_obj = CGAL::make_object(cv2.get_blue_halfedge_handle());
          }
          else
          {
            CGAL_assertion(cv2.get_color() == Curve_info::RED &&
                           cv1.get_color() == Curve_info::BLUE);
            red_obj = CGAL::make_object(cv2.get_red_halfedge_handle());
            blue_obj = CGAL::make_object(cv1.get_blue_halfedge_handle());
          }

          Point_2 point_plus (base_pt->first,
                              red_obj, blue_obj); // the extended point
          *oi = CGAL::make_object(std::make_pair(point_plus, 
                                                 base_pt->second));
        }
        else
        {
          overlap_cv = object_cast<Base_X_monotone_curve_2> (&(*oi));

          if (overlap_cv != NULL)
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

            *oi = CGAL::make_object (X_monotone_curve_2 (*overlap_cv,
                                                         red_he, blue_he));
          }
        }
      }
      //return past-end iterator
      return oi_end;
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2(m_base_traits); 
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
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      m_base_split(cv.base_curve(),
                   p.base_point(),
                   c1.base_curve(),
                   c2.base_curve());
      c1.set_curve_info(cv.get_curve_info());
      c2.set_curve_info(cv.get_curve_info());
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2(m_base_traits->split_2_object());
  }


  class Construct_min_vertex_2
  {
  private:
    Base_Construct_min_vertex_2  m_base_min_v;
    Base_Equal_2                 m_base_equal;

  public:

    Construct_min_vertex_2 (const Base_Construct_min_vertex_2& base_min_v,
                            const Base_Equal_2& base_equal):
      m_base_min_v (base_min_v),
      m_base_equal (base_equal)
    {}

    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      const Base_Point_2&   base_p = m_base_min_v (cv.base_curve());
      Object                red, blue;

      if(cv.get_color() == Curve_info::RED)
      {
        red = CGAL::make_object(cv.get_red_halfedge_handle()->target());
      }
      else if(cv.get_color() == Curve_info::BLUE)
      {
        blue = CGAL::make_object(cv.get_blue_halfedge_handle()->target());
      }
      else
      {
        CGAL_assertion(cv.get_color() == Curve_info::PURPLE);

        if (! cv.get_red_halfedge_handle()->target()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.get_red_halfedge_handle()->target()->point()))
        {
          red = CGAL::make_object(cv.get_red_halfedge_handle()->target());
        }

        if (! cv.get_blue_halfedge_handle()->target()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.get_blue_halfedge_handle()->target()->point()))
        {
          blue = CGAL::make_object(cv.get_blue_halfedge_handle()->target());
        }
      }

      return (Point_2 (base_p, red, blue));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return (Construct_min_vertex_2 
            (m_base_traits->construct_min_vertex_2_object(),
             m_base_traits->equal_2_object()));
  }


  class Construct_max_vertex_2
  {
  private:
    Base_Construct_max_vertex_2  m_base_max_v;
    Base_Equal_2                 m_base_equal;

  public:

    Construct_max_vertex_2 (const Base_Construct_max_vertex_2& base_max_v,
                            const Base_Equal_2& base_equal):
      m_base_max_v (base_max_v),
      m_base_equal (base_equal)
    {}

    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      const Base_Point_2&   base_p = m_base_max_v (cv.base_curve());
      Object                red, blue;

      if(cv.get_color() == Curve_info::RED)
      {
        red = CGAL::make_object(cv.get_red_halfedge_handle()->source());
      }
      else if(cv.get_color() == Curve_info::BLUE)
      {
        blue = CGAL::make_object(cv.get_blue_halfedge_handle()->source());
      }
      else
      {
        CGAL_assertion(cv.get_color() == Curve_info::PURPLE);

        if (! cv.get_red_halfedge_handle()->source()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.get_red_halfedge_handle()->source()->point()))
        {
          red = CGAL::make_object(cv.get_red_halfedge_handle()->source());
        }
          
        if (! cv.get_blue_halfedge_handle()->source()->is_at_infinity() &&
            m_base_equal (base_p,
                          cv.get_blue_halfedge_handle()->source()->point()))
        {
          blue = CGAL::make_object(cv.get_blue_halfedge_handle()->source());
        }
      }

      return (Point_2 (base_p, red, blue));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return (Construct_max_vertex_2 
            (m_base_traits->construct_max_vertex_2_object(),
             m_base_traits->equal_2_object()));
  }


  class Compare_xy_2
  {
  private:
    Base_Compare_xy_2 m_base_cmp_xy;

  public:

    Compare_xy_2(const Base_Compare_xy_2& base):
        m_base_cmp_xy(base)
    {}



    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      Vertex_handle_red vr1;
      Vertex_handle_red vr2;

      Vertex_handle_blue vb1;
      Vertex_handle_blue vb2;

      bool did_assign_v1_red  = assign(vr1, p1.get_red_object());
      bool did_assign_v2_red  = assign(vr2, p2.get_red_object());
      bool did_assign_v1_blue = assign(vb1, p1.get_blue_object());
      bool did_assign_v2_blue = assign(vb2, p2.get_blue_object());

      if((did_assign_v1_red && did_assign_v1_blue) ||
         (did_assign_v2_red && did_assign_v2_blue))
         return (m_base_cmp_xy(p1.base_point(), p2.base_point()));

      if(did_assign_v1_red && did_assign_v2_red)
      {
        if(vr1 == vr2)
          return EQUAL;

        return (m_base_cmp_xy(p1.base_point(), p2.base_point()));
      }

      if(did_assign_v1_blue && did_assign_v2_blue)
      {
        if(vb1 == vb2)
          return EQUAL;

        return (m_base_cmp_xy(p1.base_point(), p2.base_point()));
      }
      return (m_base_cmp_xy(p1.base_point(), p2.base_point()));
    }
  };


  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2(m_base_traits->compare_xy_2_object());
  }


  /*! \class
   * The Comapre_x_2 functor.
   */
  class Compare_x_2
  {
  private:
    Base_Compare_x_2 m_base_cmp_x;

  public:
    Compare_x_2(const Base_Compare_x_2& base):
        m_base_cmp_x(base)
    {}

    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (m_base_cmp_x(p1.base_point(), p2.base_point()));
    }

    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv,
                                  Curve_end ind) const
    {
      return (_compare_point_curve_imp (p, cv, ind,
                                        Base_has_boundary_category()));
    }

    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  Curve_end ind1,
                                  const X_monotone_curve_2& cv2,
                                  Curve_end ind2) const
    {
      return (_compare_curves_imp (cv1, ind1, cv2, ind2,
                                   Base_has_boundary_category()));
    }

  private:

    Comparison_result _compare_point_curve_imp (const Point_2& p,
                                                const X_monotone_curve_2& cv,
                                                Curve_end ind,
                                                Tag_true) const
    {
      return (m_base_cmp_x (p.base_point(), cv.base_curve(), ind));
    }

    Comparison_result _compare_point_curve_imp (const Point_2& ,
                                                const X_monotone_curve_2& ,
                                                Curve_end ,
                                                Tag_false) const
    {
      return (EQUAL);
    }

    Comparison_result _compare_curves_imp (const X_monotone_curve_2& cv1, 
                                           Curve_end ind1,
                                           const X_monotone_curve_2& cv2,
                                           Curve_end ind2,
                                           Tag_true) const
    {
      return (m_base_cmp_x (cv1.base_curve(), ind1, cv2.base_curve(), ind2));
    }

    Comparison_result _compare_curves_imp (const X_monotone_curve_2& ,
                                           Curve_end,
                                           const X_monotone_curve_2& , 
                                           Curve_end ind2,
                                           Tag_false) const
    {
      return (EQUAL);
    }

  };

  Compare_x_2 compare_x_2_object () const
  {
    return (Compare_x_2 (m_base_traits->compare_x_2_object()));
  }


  class Compare_y_at_x_2
  {
  private:
    Base_Compare_y_at_x_2 m_base_cmp_y_at_x;

  public:

    Compare_y_at_x_2(const Base_Compare_y_at_x_2& base):
        m_base_cmp_y_at_x(base)
    {}



    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Comparison_result operator() (const Point_2 & p,
                                  const X_monotone_curve_2 & cv) const
    {
      return (m_base_cmp_y_at_x(p.base_point(), cv.base_curve()));
    }

     Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2, 
                                  Curve_end ind) const
    {
      // The function is implemented based on the Has_infinite category.
      // If the traits class does not support unbounded curves, we just
      // return EQUAL, as this comparison will not be invoked anyway.
      return _comp_y_at_infinity_imp (cv1, cv2, ind, 
                                      Base_has_boundary_category());
    }
    private:

    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& cv1,
                                               const X_monotone_curve_2& cv2, 
                                               Curve_end ind,
                                               Tag_true) const
    {
      return (m_base_cmp_y_at_x (cv1.base_curve(), cv2.base_curve(), ind));
    }

    Comparison_result _comp_y_at_infinity_imp (const X_monotone_curve_2& ,
                                               const X_monotone_curve_2& , 
                                               Curve_end ,
                                               Tag_false) const
    {
      return (EQUAL);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(m_base_traits->compare_y_at_x_2_object());
  }


  class Compare_y_at_x_right_2
  {
  private:
    Base_Compare_y_at_x_right_2    m_base_cmp_y_at_x_right;

  public:

    Compare_y_at_x_right_2(const Base_Compare_y_at_x_right_2& base):
        m_base_cmp_y_at_x_right(base)
    {}



    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
     Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      return (m_base_cmp_y_at_x_right(cv1.base_curve(),
                                      cv2.base_curve(),
                                      p.base_point()));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2(m_base_traits->compare_y_at_x_right_2_object());
  }

  bool  are_same_color(const X_monotone_curve_2& cv1,
                       const X_monotone_curve_2& cv2) const
  {
    return  (cv1.get_color() == cv2.get_color());
  }

  class Boundary_in_x_2
  {
  public:
    
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      return _boundary_in_x_imp(cv, ind, Base_has_boundary_category());
    }

  private:
    Boundary_type _boundary_in_x_imp(const X_monotone_curve_2& cv,
                                     Curve_end ind, Tag_true) const
    {
      Traits tr;
      return (tr.boundary_in_x_2_object() (cv.base_curve(), ind));
    }

    Boundary_type _boundary_in_x_imp(const X_monotone_curve_2& ,
                                     Curve_end , Tag_false) const
    {
      return NO_BOUNDARY;
    }
  };

  /*! Get an Boundary_in_x_2 functor object. */
  Boundary_in_x_2 boundary_in_x_2_object () const
  {
    return Boundary_in_x_2();
  } 

  class Boundary_in_y_2
  {
  public:
    
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      return _boundary_in_y_imp(cv, ind, Base_has_boundary_category());
    }

  private:
    Boundary_type _boundary_in_y_imp(const X_monotone_curve_2& cv,
                                     Curve_end ind, Tag_true) const
    {
      Traits tr;
      return (tr.boundary_in_y_2_object() (cv.base_curve(), ind));
    }

    Boundary_type _boundary_in_y_imp(const X_monotone_curve_2& ,
                                     Curve_end , Tag_false) const
    {
      return NO_BOUNDARY;
    }
  };

  /*! Get an Boundary_in_x_2 functor object. */
  Boundary_in_y_2 boundary_in_y_2_object () const
  {
    return Boundary_in_y_2();
  } 
};


CGAL_END_NAMESPACE

#endif
