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

#ifndef ARR_ADDITION_TRAITS_H
#define ARR_ADDITION_TRAITS_H

#include <CGAL/Object.h> 
#include <list>
#include <iterator>

CGAL_BEGIN_NAMESPACE



template <class Traits, class Arrangement_>
class Arr_addition_traits : public Traits
{
public:

  typedef Arrangement_                              Arrangement;

  typedef typename Arrangement::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement::Vertex_handle       Vertex_handle;
  typedef typename Traits::X_monotone_curve_2       Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                  Base_Point_2;
  typedef typename Traits::Curve_2                  Curve_2;
  typedef typename Traits::Intersect_2              Base_Intersect_2;
  typedef typename Traits::Split_2                  Base_Split_2;
  typedef typename Traits::Make_x_monotone_2        Base_Make_x_monotone_2;



  //Constructor
  Arr_addition_traits()
  {}

  Arr_addition_traits(Traits& tr): Traits(tr)
  {}

  // nested class My_X_monotone_curve_2
  class My_X_monotone_curve_2 : public Base_X_monotone_curve_2 
  {
  public:
    typedef typename Traits::X_monotone_curve_2 Base;
    typedef typename Traits::Point_2            Point_2;

    friend class Arr_addition_traits<Traits, Halfedge_handle>;
    friend class Intersect_2;

    My_X_monotone_curve_2():Base(),
                            m_he_handle(NULL)
    {}

    My_X_monotone_curve_2(const Base& cv):Base(cv),
                                          m_he_handle(NULL)
    {}

    My_X_monotone_curve_2(const Base& cv, Halfedge_handle he):Base(cv),
                                                              m_he_handle(he)
    {}

    My_X_monotone_curve_2& operator=(const Base& cv)
    {
      Base::operator=(cv);
      m_he_handle = Halfedge_handle();
      return (*this);
    }




    Halfedge_handle get_halfedge_handle() const
    {
      return m_he_handle;
    }

    void set_halfedge_handle(Halfedge_handle he)
    {
      m_he_handle = he;
    }



    protected:
     Halfedge_handle m_he_handle;
  }; // nested class My_X_monotone_curve_2 - END


  class My_Point_2 : public Base_Point_2 
  {
  public:
    typedef  Base_Point_2            Base;

    My_Point_2(): Base(),
                  m_v(NULL)
    {}

    My_Point_2(const Base& pt): Base(pt),
                                m_v(NULL)
    {}

    My_Point_2(const Base& pt, Vertex_handle v): Base(pt),
                                                 m_v(v)
    {}


    Vertex_handle get_vertex_handle() const
    {
      return m_v;
    }

    void set_vertex_handle(Vertex_handle v)
    {
      m_v = v;
    }

   
    protected:
    Vertex_handle    m_v;
  }; // nested class My_Point_2 - END

  
  typedef My_X_monotone_curve_2                     X_monotone_curve_2;
  typedef My_Point_2                                Point_2;

 
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
      if(cv1.get_halfedge_handle() != Halfedge_handle(NULL) &&
         cv2.get_halfedge_handle() != Halfedge_handle(NULL) )
        return (oi); // the curves are disjoint-interior because they
                     // are already at the Arrangement

      OutputIterator oi_end = m_base_intersect(cv1, cv2, oi);

      // convert objects that are associated with Base_X_monotone_curve_2 to
      // X_monotone_curve_2 
      for(; oi != oi_end; ++oi)
      {
        Base_X_monotone_curve_2 overlap_cv;
        if(CGAL::assign(overlap_cv, *oi))
        {
          Halfedge_handle he(NULL);
          if(cv1.m_he_handle != Halfedge_handle(NULL))
            he = cv1.m_he_handle;
          else
            if(cv2.m_he_handle != Halfedge_handle(NULL))
              he = cv2.m_he_handle;
          X_monotone_curve_2 new_overlap_cv(overlap_cv, he);
          *oi = make_object(new_overlap_cv);
        }
        else
        {
          std::pair<Base_Point_2, unsigned int>  pt;
          CGAL_assertion(CGAL::assign(pt, *oi));
          CGAL::assign(pt, *oi);
          *oi = make_object(std::make_pair(Point_2(pt.first), pt.second));
        }
      }
      //return past-end iterator
      return oi_end;
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () 
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
      c1.set_halfedge_handle(cv.get_halfedge_handle());
      c2.set_halfedge_handle(cv.get_halfedge_handle());
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () 
  {
    return Split_2(Traits::split_2_object());
  }



  class Make_x_monotone_2
  {
  private:
    Base_Make_x_monotone_2    m_base_make_x;
    std::list<Object>         m_objects;

  public:

    /*! Constructor. */
    Make_x_monotone_2 (const Base_Make_x_monotone_2& base) :
        m_base_make_x (base)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) 
    {
      m_base_make_x(cv, std::back_inserter(m_objects));
      for(std::list<Object>::iterator iter = m_objects.begin();
          iter != m_objects.end();
          ++iter)
      {
        Base_X_monotone_curve_2 cv;
        if(assign(cv, *iter))
        {
          *oi++ = make_object(X_monotone_curve_2(cv));
        }
        else
        {
          Base_Point_2  pt;
          CGAL_assertion(assign(pt, *iter));
          assign(pt,*iter);

          *oi++ = make_object(Point_2(pt));
        }
      }
      m_objects.clear();
      return oi;
    }

  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () 
  {
    return Make_x_monotone_2(Traits::make_x_monotone_2_object());
  }

};


CGAL_END_NAMESPACE

#endif

