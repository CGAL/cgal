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
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_PARTIAL_VD_META_TRAITS_H
#define CGAL_PARTIAL_VD_META_TRAITS_H


CGAL_BEGIN_NAMESPACE

template <class Traits, class Arrangement_>
class Partial_vd_meta_traits : public Traits
{
public:
  typedef Arrangement_    Arrangement;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;

  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                  Base_Point_2;
  typedef typename Traits::Construct_min_vertex_2   Base_construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2   Base_construct_max_vertex_2;

  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2 : public Base_X_monotone_curve_2 
  {
  public:
    typedef  Base_X_monotone_curve_2    Base;

    friend class Partial_vd_meta_traits<Traits,
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




  class My_Point_2 : public Base_Point_2 
  {
  public:
    typedef  Base_Point_2            Base;

    friend class Partial_vd_meta_traits<Traits,
                                        Halfedge_const_handle>;

    My_Point_2(): Base(),
                  m_v(NULL)
    {}

    My_Point_2(const Base& pt): Base(pt),
                                m_v(NULL)
    {}

    My_Point_2(const Base& pt, Vertex_const_handle v): Base(pt), m_v(v)
    {}


    Vertex_const_handle get_vertex_handle() const
    {
      return m_v;
    }

   
    protected:
    Vertex_const_handle    m_v;
  }; // nested class My_Point_2 - END

  typedef My_X_monotone_curve_2                     X_monotone_curve_2;
  typedef My_Point_2                                Point_2;

  class Construct_min_vertex_2
  {
  private:
    Base_construct_min_vertex_2    m_base_cons_min;

  public:

    /*! Constructor. */
    Construct_min_vertex_2 (const Base_construct_min_vertex_2& base) :
        m_base_cons_min (base)
    {}

    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      const Base_Point_2& pt = m_base_cons_min(cv);
      CGAL_assertion_code(Traits traits);
      CGAL_assertion(traits.equal_2_object()(pt, cv.get_halfedge_handle()->source()->point()));
      // the halfedge in cv should be directed from left to right
      Point_2 p(pt, cv.get_halfedge_handle()->source());
      return p;
    }

  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return Construct_min_vertex_2(Traits::construct_min_vertex_2_object());
  }

  class Construct_max_vertex_2
  {
  private:
    Base_construct_max_vertex_2    m_base_cons_max;

  public:

    /*! Constructor. */
    Construct_max_vertex_2 (const Base_construct_max_vertex_2& base) :
        m_base_cons_max (base)
    {}

    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      const Base_Point_2& pt = m_base_cons_max(cv);
      CGAL_assertion_code(Traits traits);
      CGAL_assertion(traits.equal_2_object()(pt, cv.get_halfedge_handle()->target()->point()));
      // the halfedge in cv should be directed from left to right
      Point_2 p(pt, cv.get_halfedge_handle()->target());
      return p;
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2(Traits::construct_max_vertex_2_object());
  }


public:

};


CGAL_END_NAMESPACE

#endif
