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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_BATCHED_POIN_LOCATION_META_TRAITS_H
#define CGAL_ARR_BATCHED_POIN_LOCATION_META_TRAITS_H

CGAL_BEGIN_NAMESPACE

template <class Traits, class Arrangement_>
class Arr_batched_point_location_meta_traits : public Traits
{
protected:
  const Traits*    m_base_traits;
 
public:
  typedef Arrangement_    Arrangement;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;

  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                  Base_Point_2;


  // Constructor
  Arr_batched_point_location_meta_traits()
  {}

  Arr_batched_point_location_meta_traits(const Traits* tr) : m_base_traits(tr)
  {}

  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2 
  {
  public:
    typedef  Base_X_monotone_curve_2    Base;

    friend class Arr_batched_point_location_meta_traits<Traits,
                                                        Halfedge_const_handle>;

    My_X_monotone_curve_2():
      m_base_cv(),
      m_he_handle(NULL)
    {}

    My_X_monotone_curve_2(const Base& cv):
      m_base_cv(cv),
      m_he_handle(NULL)
    {}

    My_X_monotone_curve_2(const Base&cv, Halfedge_const_handle he):
      m_base_cv(cv),
      m_he_handle(he)
    {}

    Halfedge_const_handle get_halfedge_handle() const
    {
      return m_he_handle;
    }

    const Base& base_curve() const
    {
      return (m_base_cv);
    }

    Base& base_curve()
    {
      return (m_base_cv);
    }

    operator const Base&() const
    {
      return (m_base_cv);
    }

    operator Base&()
    {
      return (m_base_cv);
    }

   
    protected:
      Base_X_monotone_curve_2 m_base_cv;
      Halfedge_const_handle m_he_handle;
  }; // nested class My_X_monotone_curve_2 - END

  class My_Point_2 
  {
  public:
    typedef  Base_Point_2            Base;

    friend class Arr_batched_point_location_meta_traits<Traits,
                                                        Halfedge_const_handle>;

    My_Point_2():
      m_base_pt(),
      m_v(NULL)
    {}

    My_Point_2(const Base& pt):
      m_base_pt(pt),
      m_v(NULL)
    {}

    My_Point_2(const Base& pt, Vertex_const_handle v):
      m_base_pt(pt),
      m_v(v)
    {}

    const Base& base_point() const
    {
      return (m_base_pt);
    }

    Base& base_point()
    {
      return (m_base_pt);
    }

    operator const Base&() const
    {
      return (m_base_pt);
    }

    operator Base&()
    {
      return (m_base_pt);
    }

    Vertex_const_handle get_vertex_handle() const
    {
      return m_v;
    }
   
  protected:
    Base                   m_base_pt;
    Vertex_const_handle    m_v;
  
  }; // nested class My_Point_2 - END

public:

  typedef My_X_monotone_curve_2                     X_monotone_curve_2; 
  typedef My_Point_2                                Point_2; 

  typedef typename Traits::Construct_min_vertex_2   Base_Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2   Base_Construct_max_vertex_2;
  typedef typename Traits::Compare_xy_2             Base_Compare_xy_2;
  typedef typename Traits::Compare_y_at_x_right_2   Base_Compare_y_at_x_right_2;

  class Construct_min_vertex_2
  {
  private:
    Base_Construct_min_vertex_2 m_base_min_v;

  public:

    Construct_min_vertex_2(const Base_Construct_min_vertex_2& base):
        m_base_min_v(base)
    {}



    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      Vertex_const_handle vh = cv.get_halfedge_handle()->target();
      return (Point_2(m_base_min_v(cv.base_curve()), vh));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2(m_base_traits->construct_min_vertex_2_object());
  }


  class Construct_max_vertex_2
  {
  private:
    Base_Construct_max_vertex_2 m_base_max_v;

  public:

    Construct_max_vertex_2(const Base_Construct_max_vertex_2& base):
        m_base_max_v(base)
    {}



    /*!
     * Get the right endpoint of the x-monotone curve .
     * \param cv The curve.
     * \return The right endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) 
    {
      Vertex_const_handle vh = cv.get_halfedge_handle()->source();
      return (Point_2(m_base_max_v(cv.base_curve()), vh));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2(m_base_traits->construct_max_vertex_2_object());
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
      if(p1.get_vertex_handle() == p2.get_vertex_handle() &&
         p1.get_vertex_handle() != Vertex_const_handle())
        return EQUAL;

      return (m_base_cmp_xy(p1.base_point(), p2.base_point()));
    }
  };


  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () 
  {
    return Compare_xy_2(m_base_traits->compare_xy_2_object());
  }


  class Compare_y_at_x_right_2
  {
  private:
    Base_Compare_y_at_x_right_2 m_base_comp;

  public:

    Compare_y_at_x_right_2(const Base_Compare_y_at_x_right_2& base):
      m_base_comp(base)
    {}

    /*!
     * Compare two x-monotone curve to the right of their common left endpoint.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) 
    {
      // Both cv1 and cv2 are associated with halfedges incident to the same
      // vertex. Thus, if we go in a clockwise direction from he1 and find
      // he2 after crossing to the other side (the left side) of the vertex
      // and crossing back, he1 is below he2 (and vice-versa). 
      Halfedge_const_handle  he1 = cv1.get_halfedge_handle();
      Halfedge_const_handle  he2 = cv2.get_halfedge_handle();
      Halfedge_const_handle  he;

      CGAL_assertion (he1->target() == p.get_vertex_handle());
      CGAL_assertion (he2->target() == p.get_vertex_handle());

      he = he1;
      do
      {
        he = he->next()->twin();
        if (he == he2)
          break;
        if (he->direction() == SMALLER)
          return (SMALLER);
      } while (he != he1);

      he = he2;
      do
      {
        he = he->next()->twin();
        if (he == he1)
          break;
        if (he->direction() == SMALLER)
          return (LARGER);
      } while (he != he2);

      // If we reached here, we have to compare geometrically:
      return (m_base_comp (cv1.base_curve(), cv2.base_curve(),
                           p.base_point()));
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () 
  {
    return 
      Compare_y_at_x_right_2(m_base_traits->compare_y_at_x_right_2_object());
  }

};


CGAL_END_NAMESPACE

#endif
