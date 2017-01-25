// Copyright (c) 1997  Tel-Aviv University (Israel).
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

#ifndef CGAL_GPS_SIMPLIFIER_TRAITS_H
#define CGAL_GPS_SIMPLIFIER_TRAITS_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_traits_decorator.h>

namespace CGAL {

class Gps_simplifier_curve_data
{
protected:
  unsigned int m_bc;
  unsigned int m_twin_bc;
  unsigned int m_index;

public:
  Gps_simplifier_curve_data()
  {}

  Gps_simplifier_curve_data(unsigned int bc,
                            unsigned int twin_bc,
                            unsigned int index):
    m_bc(bc),
    m_twin_bc(twin_bc),
    m_index(index)
  {}

  unsigned int bc() const
  {
    return m_bc;
  }

  unsigned int twin_bc() const
  {
    return m_twin_bc;
  }

  unsigned int index() const
  {
    return m_index;
  }

  unsigned int& index()
  {
    return m_index;
  }

  unsigned int& twin_bc()
  {
    return m_twin_bc;
  }

  void set_bc(unsigned int bc)
  {
    m_bc = bc;
  }

  void set_twin_bc(unsigned int twin_bc)
  {
    m_twin_bc = twin_bc;
  }

  void set_index(unsigned int index)
  {
    m_index = index;
  }
};

struct Gps_simplifier_point_data
{
protected:
  unsigned int m_index;

public:
  Gps_simplifier_point_data()
  {}

  Gps_simplifier_point_data(unsigned int index) : m_index(index)
  {}

  unsigned int index() const
  {
    return m_index;
  }

  void set_index(unsigned int index)
  {
    m_index = index;
  }
};

template <class Traits_>
class Gps_simplifier_traits :
  public Gps_traits_decorator<Traits_,
                              Gps_simplifier_curve_data,
                              Gps_simplifier_point_data>
{
public:
  typedef Traits_    Traits;
  typedef Gps_traits_decorator<Traits_,
                               Gps_simplifier_curve_data,
                               Gps_simplifier_point_data>    Base;
  typedef Gps_simplifier_traits<Traits>                      Self;
  typedef typename Traits::X_monotone_curve_2     Base_X_monotone_curve_2;
  typedef typename Traits::Point_2                Base_Point_2;
  typedef typename Traits::Construct_min_vertex_2 Base_Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2 Base_Construct_max_vertex_2;
  typedef typename Traits::Compare_endpoints_xy_2 Base_Compare_endpoints_xy_2;
  typedef typename Traits::Compare_xy_2           Base_Compare_xy_2;
  typedef typename Traits::Compare_y_at_x_right_2 Base_Compare_y_at_x_right_2;
  typedef typename Traits::Compare_y_at_x_2       Base_Compare_y_at_x_2;
  typedef typename Traits::Intersect_2            Base_Intersect_2;
  typedef typename Traits::Split_2                Base_Split_2;

protected:
  mutable unsigned int m_pgn_size;


public:

  typedef typename Base::X_monotone_curve_2       X_monotone_curve_2;
  typedef typename Base::Point_2                  Point_2;
  typedef typename Base::Multiplicity             Multiplicity;

  typedef typename Base::Curve_data               Curve_data;
  typedef typename Base::Point_data               Point_data;

  Gps_simplifier_traits()
  {}

  Gps_simplifier_traits(const Traits & tr) : Base(tr)
  {}

  unsigned int polygon_size() const
  {
    return m_pgn_size;
  }

  void set_polygon_size(unsigned int pgn_size) const
  {
    m_pgn_size = pgn_size;
  }

  bool is_valid_index(unsigned int index) const
  {
    return (index < m_pgn_size);
  }

  unsigned int invalid_index() const
  {
    return (m_pgn_size);
  }


  class Intersect_2
  {
  private:

    Base_Intersect_2             m_base;
    Base_Compare_endpoints_xy_2  m_base_cmp_endpoints;
    Base_Compare_xy_2            m_base_cmp_xy;
    Base_Construct_min_vertex_2  m_ctr_min_v;
    const Self * m_self_tr;

  public:

    /*! Constructor. */
    Intersect_2 (const Base_Intersect_2& base,
                 const Base_Compare_endpoints_xy_2& base_cmp_endpoints,
                 const Base_Compare_xy_2& base_cmp_xy,
                 const Base_Construct_min_vertex_2& ,
                 const Self*  tr) :
      m_base(base),
      m_base_cmp_endpoints(base_cmp_endpoints),
      m_base_cmp_xy(base_cmp_xy),
      m_self_tr(tr)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      //// if the two curves are incident, do not intersect them
      //if(m_self_tr->is_valid_index(cv1.data().index()) &&
      //   m_self_tr->is_valid_index(cv2.data().index()))
      //{
      //  unsigned int index_diff =
      //    (cv1.data().index() > cv2.data().index()) ?
      //    (cv1.data().index() - cv2.data().index()):
      //    (cv2.data().index() - cv1.data().index());

      //  if(index_diff == 1 ||index_diff == m_self_tr->polygon_size() -1)
      //  {
      //    return (oi);
      //  }
      //}
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
          Point_data pt_data(m_self_tr->invalid_index());
          Point_2 point_plus (base_pt->first, pt_data); // the extended point
          *oi = CGAL::make_object(std::make_pair(point_plus,
                                                 base_pt->second));
        }
        else
        {
          overlap_cv = object_cast<Base_X_monotone_curve_2> (&(*oi));

          if (overlap_cv != NULL)
          {
            unsigned int ov_bc;
            unsigned int ov_twin_bc;
            if(m_base_cmp_endpoints(cv1) == m_base_cmp_endpoints(cv2))
            {
              // cv1 and cv2 have the same directions
              ov_bc = cv1.data().bc() + cv2.data().bc();
              ov_twin_bc = cv1.data().twin_bc() + cv2.data().twin_bc();
            }
            else
            {
              // cv1 and cv2 have opposite directions
              ov_bc = cv1.data().bc() + cv2.data().twin_bc();
              ov_twin_bc = cv1.data().twin_bc() + cv2.data().bc();
            }

            if(m_base_cmp_endpoints(*overlap_cv) != m_base_cmp_endpoints(cv1))
            {
              // overlap_cv, cv1 have opposite directions
              std::swap(ov_bc, ov_twin_bc);
            }

            Curve_data cv_data(ov_bc, ov_twin_bc, m_self_tr->invalid_index());
            *oi = CGAL::make_object (X_monotone_curve_2 (*overlap_cv, cv_data));
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
    return Intersect_2(this->m_base_tr->intersect_2_object(),
                       this->m_base_tr->compare_endpoints_xy_2_object(),
                       this->m_base_tr->compare_xy_2_object(),
                       this->m_base_tr->construct_min_vertex_2_object(),
                       this);
  }

  class Split_2
  {
  private:
    Base_Split_2      m_base_split;
    const Self * m_self_tr;

  public:

    /*! Constructor. */
    Split_2 (const Base_Split_2& base, const Self* tr) :
      m_base_split(base),
      m_self_tr(tr)
    {}

    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      m_base_split(cv.base(),
                   p.base(),
                   c1.base(),
                   c2.base());
      const Curve_data& cv_data = cv.data();
      c1.set_data(Curve_data(cv_data.bc(),
                             cv_data.twin_bc(),
                             m_self_tr->invalid_index()));

      c2.set_data(Curve_data(cv_data.bc(),
                             cv_data.twin_bc(),
                             m_self_tr->invalid_index()));
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2(this->m_base_tr->split_2_object(), this);
  }

  class Construct_min_vertex_2
  {
  private:
    Base_Construct_min_vertex_2 m_base;
    Base_Compare_endpoints_xy_2 m_base_cmp_endpoints;
    const Self * m_self_tr;

  public:

    Construct_min_vertex_2(const Base_Construct_min_vertex_2& base,
                          const Base_Compare_endpoints_xy_2& base_cmp_endpoints,
                          const Self * tr):
      m_base(base),
      m_base_cmp_endpoints(base_cmp_endpoints),
      m_self_tr(tr)
    {}

    /*!
      * Get the left endpoint of the x-monotone curve (segment).
      * \param cv The curve.
      * \return The left endpoint.
      */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      if(!m_self_tr->is_valid_index(cv.data().index()))
      {
        return Point_2 (m_base(cv.base()), m_self_tr->invalid_index());
      }

      Comparison_result res = m_base_cmp_endpoints(cv);
      Point_data pt_data;
      if(res == SMALLER)
      {
        // min vertex is the source
        pt_data.set_index(cv.data().index());
      }
      else
      {
        // min vertex is the target
        pt_data.set_index((cv.data().index() + 1) % m_self_tr->polygon_size());
      }
      return Point_2 (m_base(cv.base()), pt_data);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2
      (this->m_base_tr->construct_min_vertex_2_object(),
       this->m_base_tr->compare_endpoints_xy_2_object(),
       this);
  }


  class Construct_max_vertex_2
  {
  private:
    Base_Construct_max_vertex_2      m_base;
    Base_Compare_endpoints_xy_2      m_base_cmp_endpoints;
    const Self * m_self_tr;

  public:

    Construct_max_vertex_2(const Base_Construct_max_vertex_2& base,
                          const Base_Compare_endpoints_xy_2& base_cmp_endpoints,
                          const Self * tr):
      m_base(base),
      m_base_cmp_endpoints(base_cmp_endpoints),
      m_self_tr(tr)
    {}

    /*!
      * Get the right endpoint of the x-monotone curve (segment).
      * \param cv The curve.
      * \return The left endpoint.
      */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      if(!m_self_tr->is_valid_index(cv.data().index()))
      {
        return Point_2 (m_base(cv.base()), m_self_tr->invalid_index());
      }
      Comparison_result res = m_base_cmp_endpoints(cv);
      Point_data pt_data;
      if(res == SMALLER)
      {
        // min vertex is the target
        pt_data.set_index((cv.data().index() + 1) % m_self_tr->polygon_size());
      }
      else
      {
        // min vertex is the source
        pt_data.set_index(cv.data().index());
      }
      return Point_2 (m_base(cv.base()), pt_data);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2
      (this->m_base_tr->construct_max_vertex_2_object(),
       this->m_base_tr->compare_endpoints_xy_2_object(),
       this);
  }

  class Compare_xy_2
  {
  private:
    Base_Compare_xy_2       m_base;
    const Self * m_self_tr;

  public:
    Compare_xy_2(const Base_Compare_xy_2& base,
                 const Self * tr):
      m_base(base),
      m_self_tr(tr)
    {}


    /*!
      * Get the left endpoint of the x-monotone curve (segment).
      * \param cv The curve.
      * \return The left endpoint.
      */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      //if one of the indexes is invalid, compare p1 and p2
      if(! m_self_tr->is_valid_index(p1.data().index()) ||
        ! m_self_tr->is_valid_index(p2.data().index()))
        return (m_base(p1.base(), p2.base()));

      // if the two point has the same index, return EQUAL
      if(p1.data().index() == p2.data().index())
      {
        return EQUAL;
      }

      return (m_base(p1.base(), p2.base()));
    }
  };


  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2(this->m_base_tr->compare_xy_2_object(), this);
  }
};

} //namespace CGAL

#endif
