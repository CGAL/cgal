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

#ifndef ARR_BATCHED_POIN_LOCATION_META_TRAITS_H
#define ARR_BATCHED_POIN_LOCATION_META_TRAITS_H


CGAL_BEGIN_NAMESPACE



template <class Traits, class _Halfedge_const_handle>
class Arr_batched_point_location_meta_traits : public Traits
{

  typedef _Halfedge_const_handle  Halfedge_const_handle;
  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2 : public Traits::X_monotone_curve_2 
  {
  public:
    typedef typename Traits::X_monotone_curve_2 Base;
    typedef typename Traits::Point_2            Point_2;

    friend class Arr_batched_point_location_meta_traits<Traits,
                                                        Halfedge_const_handle>;

    My_X_monotone_curve_2():Base(),
                            m_he_handle(NULL)
    {}

    My_X_monotone_curve_2(const Base& cv):Base(cv),
                                          m_he_handle(NULL)
    {}

    My_X_monotone_curve_2(const Base&cv, Halfedge_const_handle he):Base(cv),
                                                                   m_he_handle(he)
    {}


    Halfedge_const_handle get_halfedge_handle() const
    {
      return m_he_handle;
    }

   
    protected:
    Halfedge_const_handle m_he_handle;
  }; // nested class My_X_monotone_curve_2 - END

public:
  typedef typename Traits::Point_2                  Base_Point_2;

  class My_Point_2 : public Base_Point_2
  {
  public:

    enum Type
    {
      VERTEX = 1,
      QUERY  = 2,
      ISOLATED_VERTEX = 4
    };
  private:

    char m_type;

  public:

    My_Point_2():
        Base_Point_2(),
        m_type(VERTEX)
    {}

    My_Point_2 (Type type):
        Base_Point_2(),
        m_type(type)
    {}

    /*
    My_Point_2(const Base_Point_2& base_pt):Base_Point_2(base_pt),
                                            m_type(0x0)
    {}
    */

    My_Point_2(const Base_Point_2& base_pt, Type type) :
        Base_Point_2(base_pt),
        m_type(type)
    {}

    /*
    My_Point_2& operator=(const Base_Point_2& pt)
    {
      Base_Point_2::operator= (pt);
      m_type = 0x0;
      return *this;
    }
    */

    bool is_query() const
    {
      return ((m_type & QUERY) != 0);
    }

    bool is_isolated_vertex() const
    {
      return ((m_type & ISOLATED_VERTEX) != 0);
    }

  };





public:

  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef My_Point_2                                Point_2;
  typedef My_X_monotone_curve_2                     X_monotone_curve_2;
  typedef typename Traits::Construct_min_vertex_2   Base_Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2   Base_Construct_max_vertex_2;


  class Intersect_2
  {
  public:
   
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      return (oi);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2();
  }

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
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      return Point_2 (m_base_min_v(cv), Point_2::VERTEX);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2(Traits::construct_min_vertex_2_object());
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
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      return Point_2 (m_base_max_v(cv), Point_2::VERTEX);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2(Traits::construct_max_vertex_2_object());
  }
};


CGAL_END_NAMESPACE

#endif

