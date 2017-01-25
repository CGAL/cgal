// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_TRAITS_DECORATOR_H
#define CGAL_GPS_TRAITS_DECORATOR_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <boost/mpl/assert.hpp>

namespace CGAL {

template <class Traits_, class Curve_data_, class Point_data_>
class Gps_traits_decorator
{
public:
  typedef Traits_                              Base;
  typedef typename Base::X_monotone_curve_2    Base_X_monotone_curve_2;
  typedef typename Base::Point_2               Base_Point_2;
  typedef typename Base::Multiplicity          Multiplicity;

  typedef Gps_traits_decorator<Traits_, Curve_data_, Point_data_>    Self;

  typedef Curve_data_                          Curve_data;
  typedef Point_data_                          Point_data;

  typedef typename Base::Compare_x_2               Base_Compare_x_2;
  typedef typename Base::Compare_xy_2              Base_Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2    Base_Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2    Base_Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2             Base_Is_vertical_2;
  typedef typename Base::Compare_y_at_x_2          Base_Compare_y_at_x_2;
  typedef typename Base::Compare_y_at_x_right_2    Base_Compare_y_at_x_right_2;
  typedef typename Base::Equal_2                   Base_Equal_2;
  typedef typename Base::Split_2                   Base_Split_2;
  typedef typename Base::Intersect_2               Base_Intersect_2;
  typedef typename Base::Compare_endpoints_xy_2    Base_Compare_endpoints_xy_2;
  typedef typename Base::Construct_opposite_2      Base_Construct_opposite_2;
  typedef typename Base::Has_left_category         Has_left_category;
  typedef typename Base::Has_merge_category        Has_merge_category;
  typedef typename Base::Has_do_intersect_category Has_do_intersect_category;

  typedef typename Base::Left_side_category        Left_side_category;
  typedef typename Base::Bottom_side_category      Bottom_side_category;
  typedef typename Base::Top_side_category         Top_side_category;
  typedef typename Base::Right_side_category       Right_side_category;

  // a side is either oblivious or open (unbounded)
  BOOST_MPL_ASSERT(
      (boost::mpl::or_< 
       boost::is_same< Left_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Left_side_category, Arr_open_side_tag > >
      )
  );
  BOOST_MPL_ASSERT(
      (boost::mpl::or_< 
       boost::is_same< Bottom_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Bottom_side_category, Arr_open_side_tag > >
      )
  );
  BOOST_MPL_ASSERT(
      (boost::mpl::or_< 
       boost::is_same< Top_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Top_side_category, Arr_open_side_tag > >
      )
  );
  BOOST_MPL_ASSERT(
      (boost::mpl::or_< 
       boost::is_same< Right_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Right_side_category, Arr_open_side_tag > >
      )
  );

  class Ex_point_2 
  {
  protected:
    Base_Point_2 m_base;
    Point_data   m_data;

  public:
    Ex_point_2(): m_base(),
                  m_data()
    {}

    Ex_point_2(const Base_Point_2& pt) : m_base(pt),
                                         m_data()
    {}

    Ex_point_2(const Base_Point_2& pt, const Point_data& data): m_base(pt),
                                                                m_data(data)
    {}

    const Base_Point_2& base() const
    {
      return m_base;
    }

    Base_Point_2& base()
    {
      return m_base;
    }

    operator const Base_Point_2& () const
    {
      return (m_base);
    }

    Point_data& data()
    {
      return m_data;
    }

    const Point_data& data() const
    {
      return m_data;
    }

    void set_data(const Point_data& data)
    {
      m_data = data;
    }


  };

  class Ex_x_monotone_curve_2
  {
  protected:
    Base_X_monotone_curve_2    m_base;
    Curve_data                 m_data;
    
    public:
    Ex_x_monotone_curve_2(): m_base(),
                             m_data()
    {}

    Ex_x_monotone_curve_2(const Base_X_monotone_curve_2& pt) : m_base(pt),
                                                               m_data()
    {}

    Ex_x_monotone_curve_2(const Base_X_monotone_curve_2& pt, 
                          const Curve_data& data): m_base(pt),
                                                   m_data(data)
    {}

    const Base_X_monotone_curve_2& base() const
    {
      return m_base;
    }

    Base_X_monotone_curve_2& base()
    {
      return m_base;
    }

    operator const Base_X_monotone_curve_2& () const
    {
      return (m_base);
    }

    Curve_data& data()
    {
      return m_data;
    }

    const Curve_data& data() const
    {
      return m_data;
    }

    void set_data(const Curve_data& data)
    {
      m_data = data;
    }
  };

  friend std::ostream& operator<< (std::ostream& os, const Ex_x_monotone_curve_2 & cv)
  {
    os << cv.base();
    return (os);
  }

  friend std::ostream& operator<< (std::ostream& os, const Ex_point_2 & pt)
  {
    os << pt.base();
    return (os);
  }

  typedef Ex_point_2                Point_2;
  typedef Ex_x_monotone_curve_2     X_monotone_curve_2;

  protected:

  //Data members
  const Base * m_base_tr;
  bool m_traits_owner;

public:
  
  Gps_traits_decorator() :
    m_base_tr(new Base()),
    m_traits_owner(true)
  {}

  Gps_traits_decorator(const Base & base_traits) :
    m_base_tr(&base_traits),
    m_traits_owner(false)
  {}

  ~Gps_traits_decorator()
  {
    if (m_traits_owner)
      delete m_base_tr;
  }

  class Compare_x_2
  {
  protected:

    Base_Compare_x_2 m_base;

  public:
    Compare_x_2(Base_Compare_x_2& base) : m_base(base)
    {}

    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (m_base(p1.base(), p2.base()));
    }
  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return Compare_x_2(m_base_tr->compare_x_2_object());
  }


  class Compare_xy_2
  {
  protected:

    Base_Compare_xy_2 m_base;

  public:
    Compare_xy_2(const Base_Compare_xy_2& base) : m_base(base)
    {}

    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (m_base(p1.base(), p2.base()));
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2(m_base_tr->compare_xy_2_object());
  }

  class Construct_min_vertex_2
  {
  protected:

    Base_Construct_min_vertex_2 m_base;

  public:
    Construct_min_vertex_2(const Base_Construct_min_vertex_2& base) :
      m_base(base)
    {}

    Point_2 operator() (const X_monotone_curve_2& cv) const
    {
      return Point_2(m_base(cv.base()));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2(m_base_tr->construct_min_vertex_2_object());
  }

  class Construct_max_vertex_2
  {
  protected:

    Base_Construct_max_vertex_2 m_base;

  public:
    Construct_max_vertex_2(const Base_Construct_max_vertex_2& base) :
      m_base(base)
    {}

    Point_2 operator() (const X_monotone_curve_2& cv) const
    {
      return Point_2(m_base(cv.base()));
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2(m_base_tr->construct_max_vertex_2_object());
  }


  class Is_vertical_2
  {
  protected:

    Base_Is_vertical_2 m_base;

  public:
    Is_vertical_2(const Base_Is_vertical_2& base) : m_base(base)
    {}

    bool operator() (const X_monotone_curve_2& cv) const
    {
      return (m_base(cv.base()));
    }
  };

  /*! Get a Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  {
    return Is_vertical_2(m_base_tr->is_vertical_2_object());
  }


  class Compare_y_at_x_2
  {
  protected:

    Base_Compare_y_at_x_2 m_base;

  public:
    Compare_y_at_x_2(const Base_Compare_y_at_x_2& base) : m_base(base)
    {}

    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      return m_base(p.base(), cv.base());
    }
  };

  /*! Get a compare_y_at_x_2_object functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  {
    return Compare_y_at_x_2(m_base_tr->compare_y_at_x_2_object());
  }


  class Compare_y_at_x_right_2
  {
  protected:

    Base_Compare_y_at_x_right_2 m_base;

  public:
    Compare_y_at_x_right_2(const Base_Compare_y_at_x_right_2& base) :
      m_base(base)
    {}

    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      return m_base(cv1.base(), cv2.base(), p.base());
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  {
    return Compare_y_at_x_right_2(m_base_tr->compare_y_at_x_right_2_object());
  }


  class Equal_2
  {
  protected:

    Base_Equal_2 m_base;

  public:
    Equal_2(const Base_Equal_2& base) : m_base(base)
    {}

    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      return m_base(p1.base(), p2.base());
    }
  };

  /*! Get a Equal_2 functor object. */
  Equal_2 equal_2_object() const
  {
    return Equal_2(m_base_tr->equal_2_object());
  }


  class Split_2
  {
  protected:

    Base_Split_2 m_base;

  public:
    Split_2(const Base_Split_2& base) : m_base(base)
    {}

    void operator() (const X_monotone_curve_2& cv,
                     const Point_2& p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      m_base(cv.base(), p.base(), c1.base(), c2.base());
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object() const
  {
    return Split_2(m_base_tr->split_2_object());
  }


  class Intersect_2
  {
  protected:

    Base_Intersect_2             m_base;
    Base_Compare_xy_2            m_base_cmp_xy;
    Base_Construct_min_vertex_2  m_ctr_min_v;

  public:
    Intersect_2(Base_Intersect_2&            base,
                Base_Compare_xy_2&           base_cmp_xy,
                Base_Construct_min_vertex_2& base_ctr_min_v) :
      m_base(base),
      m_base_cmp_xy(base_cmp_xy),
      m_ctr_min_v(base_ctr_min_v)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      const std::pair<Base_Point_2, Multiplicity>   *base_pt;
      const Base_X_monotone_curve_2                 *overlap_cv;

       OutputIterator oi_end;
      if(m_base_cmp_xy(m_ctr_min_v(cv1.base()),
                       m_ctr_min_v(cv2.base())) == LARGER)
        oi_end = m_base(cv1.base(), cv2.base(), oi);
      else
        oi_end = m_base(cv2.base(), cv1.base(), oi);

      // convert objects that are associated with Base_X_monotone_curve_2 to
      // the extenede X_monotone_curve_2 
      for(; oi != oi_end; ++oi)
      {
        base_pt = object_cast<std::pair<Base_Point_2, Multiplicity> >(&(*oi));

        if (base_pt != NULL)
        {
          Point_2 point_plus (base_pt->first); // the extended point
          *oi = CGAL::make_object(std::make_pair(point_plus, 
                                                 base_pt->second));
        }
        else
        {
          overlap_cv = object_cast<Base_X_monotone_curve_2> (&(*oi));
          CGAL_assertion(overlap_cv != NULL);
          *oi = CGAL::make_object (X_monotone_curve_2 (*overlap_cv));
        }
      }
      //return past-end iterator
      return oi_end;
    }    
  };

  /*! Get a Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const
  {
    return Intersect_2(m_base_tr->intersect_2_object(),
                       m_base_tr->compare_xy_2_object(),
                       m_base_tr->construct_min_vertex_2_object());
  }


  class Compare_endpoints_xy_2
  {
  protected:

    Base_Compare_endpoints_xy_2 m_base;

  public:
    Compare_endpoints_xy_2(const Base_Compare_endpoints_xy_2& base) :
      m_base(base)
    {}

    Comparison_result operator()(const X_monotone_curve_2& cv) const
    {
      return (m_base(cv));

    }
  };

  /*! Get a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  {
    return Compare_endpoints_xy_2(m_base_tr->compare_endpoints_xy_2_object());
  }


  class Construct_opposite_2
  {
  protected:

    Base_Construct_opposite_2 m_base;

  public:
    Construct_opposite_2(Base_Construct_opposite_2& base) :m_base(base)
    {}

    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv) const
    {
      return (X_monotone_curve_2(m_base(cv)));
    }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2(m_base_tr->construct_opposite_2_object());
  }

};

} //namespace CGAL

#endif
