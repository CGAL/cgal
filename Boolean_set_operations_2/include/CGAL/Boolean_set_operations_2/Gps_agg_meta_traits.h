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

#ifndef CGAL_GPS_AGG_META_TRAITS_H
#define CGAL_GPS_AGG_META_TRAITS_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <vector>
#include <boost/mpl/assert.hpp>
#include <CGAL/Boolean_set_operations_2/Gps_traits_decorator.h>
#include <CGAL/Boolean_set_operations_2/Curve_with_halfedge.h>
#include <CGAL/Boolean_set_operations_2/Point_with_vertex.h>

namespace CGAL {

template <class Arrangement_>
class Gps_agg_curve_data : public Curve_with_halfedge<Arrangement_>
{
protected:
  typedef Arrangement_                             Arrangement;
  typedef typename Arrangement::Halfedge_handle    Halfedge_handle;
  typedef Curve_with_halfedge<Arrangement_>        Base;

  const Arrangement*  m_arr; // pointer to the arrangement containing the edge.
  unsigned int        m_bc; // the boudary counter of the halfedge with the same
                          //direction as the curve

  unsigned int      m_twin_bc;  // the boudary counter of the halfedge with the same
                                 //direction as the curve

public:

  Gps_agg_curve_data() :
    Base(),
    m_arr(NULL),
    m_bc(0),
    m_twin_bc(0)
  {}

  Gps_agg_curve_data(const Arrangement* arr,
                     Halfedge_handle he,
                     unsigned int bc,
                     unsigned int twin_bc):
    Base(he),
    m_arr(arr),
    m_bc(bc),
    m_twin_bc(twin_bc)
  {}

  unsigned int bc() const
  {
    return m_bc;
  }

  unsigned int twin_bc() const
  {
    return m_twin_bc;
  }

  unsigned int& bc() 
  {
    return m_bc;
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

  const Arrangement* arr() const
  {
    return m_arr;
  }

};


template <class Arrangement_>
class Gps_agg_meta_traits :
  public Gps_traits_decorator<typename Arrangement_::Traits_adaptor_2,
                              Gps_agg_curve_data<Arrangement_>,
                              Point_with_vertex<Arrangement_> >
{
  typedef Arrangement_                            Arrangement;
  typedef typename Arrangement::Traits_adaptor_2  Traits;

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

  typedef typename Traits::Parameter_space_in_x_2 Base_Parameter_space_in_x_2;
  typedef typename Traits::Compare_y_near_boundary_2 
                                                 Base_Compare_y_near_boundary_2;

  typedef typename Traits::Parameter_space_in_y_2 Base_Parameter_space_in_y_2;
  typedef typename Traits::Compare_x_near_boundary_2 
                                                 Base_Compare_x_near_boundary_2;


  typedef Gps_agg_meta_traits<Traits>             Self;
  typedef Gps_traits_decorator<Traits,
                               Gps_agg_curve_data<Arrangement_>,
                               Point_with_vertex<Arrangement_> > 
                                                  Base;

  public:
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Base::Point_2                      Point_2;
  typedef typename Traits::Has_left_category          Has_left_category;
  typedef typename Traits::Has_merge_category         Has_merge_category;
  typedef typename Traits::Has_do_intersect_category  Has_do_intersect_category;

  typedef typename Arrangement::Left_side_category    Left_side_category;
  typedef typename Arrangement::Bottom_side_category  Bottom_side_category;
  typedef typename Arrangement::Top_side_category     Top_side_category;
  typedef typename Arrangement::Right_side_category   Right_side_category;

  typedef typename Traits::Multiplicity               Multiplicity; 
  

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

  typedef typename Base::Curve_data               Curve_data;
  typedef typename Base::Point_data               Point_data;

  typedef typename Arrangement::Halfedge_handle   Halfedge_handle;
  typedef typename Arrangement::Vertex_handle     Vertex_handle;


  Gps_agg_meta_traits()
  {}

  Gps_agg_meta_traits(const Traits & base_tr) : Base(base_tr)
  {}

  class Intersect_2
  {
  private:
    Base_Intersect_2             m_base;
    Base_Compare_endpoints_xy_2  m_base_cmp_endpoints;
    Base_Compare_xy_2            m_base_cmp_xy;
    Base_Construct_min_vertex_2  m_base_ctr_min_v;

  public:
    /*! Constructor. */
    Intersect_2 (const Base_Intersect_2& base,
                 const Base_Compare_endpoints_xy_2& base_cmp_endpoints,
                 const Base_Compare_xy_2& base_cmp_xy,
                 const Base_Construct_min_vertex_2& base_ctr_min_v) : 
      m_base(base),
      m_base_cmp_endpoints(base_cmp_endpoints),
      m_base_cmp_xy(base_cmp_xy),
      m_base_ctr_min_v(base_ctr_min_v)
    {}

    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      if (cv1.data().arr() == cv2.data().arr())
      {
        return (oi); // the curves are disjoint-interior because they
                     // are already at the same arrangement.
      }
      
      const std::pair<Base_Point_2, Multiplicity>   *base_pt;
      const Base_X_monotone_curve_2                 *overlap_cv;
      OutputIterator oi_end;
      if(m_base_cmp_xy(m_base_ctr_min_v(cv1.base()),
                       m_base_ctr_min_v(cv2.base())) == LARGER)
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

          if (overlap_cv != NULL)
          {
            unsigned int ov_bc;
            unsigned int ov_twin_bc;
            if (m_base_cmp_endpoints(cv1) == m_base_cmp_endpoints(cv2))
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

            Curve_data cv_data(cv1.data().arr(),
                               Halfedge_handle(),
                               ov_bc, ov_twin_bc);
              *oi = CGAL::make_object (X_monotone_curve_2(*overlap_cv, cv_data));
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
                       this->m_base_tr->construct_min_vertex_2_object());
  }


  class Split_2
  {
  private:
    Base_Split_2      m_base_split;

  public:

    /*! Constructor. */
    Split_2 (const Base_Split_2& base) : m_base_split (base)
    {}

    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      m_base_split(cv.base(),
                   p.base(),
                   c1.base(),
                   c2.base());
      const Curve_data& cv_data = cv.data();
      c1.set_data(Curve_data(cv_data.arr(),
                             Halfedge_handle(),
                             cv_data.bc(),
                             cv_data.twin_bc()));
      
      c2.set_data(Curve_data(cv_data.arr(),
                             Halfedge_handle(),
                             cv_data.bc(),
                             cv_data.twin_bc()));
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2(this->m_base_tr->split_2_object());
  }


  class Construct_min_vertex_2
  {
  private:
    Base_Construct_min_vertex_2 m_base;

  public:

    Construct_min_vertex_2(const Base_Construct_min_vertex_2& base) :
        m_base(base)
    {}
  
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      if(cv.data().halfedge() == Halfedge_handle())
        return (Point_2 (m_base(cv.base())));

      CGAL_assertion
        ((Arr_halfedge_direction)cv.data().halfedge()->direction() ==
         ARR_LEFT_TO_RIGHT);
      return Point_2 (m_base(cv.base()), cv.data().halfedge()->source());
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2(this->m_base_tr->
                                  construct_min_vertex_2_object());
  }


  class Construct_max_vertex_2
  {
  private:
    Base_Construct_max_vertex_2 m_base;

  public:

    Construct_max_vertex_2(const Base_Construct_max_vertex_2& base):
        m_base(base)
    {}
  
    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      if(cv.data().halfedge() == Halfedge_handle())
        return (Point_2 (m_base(cv.base())));

      CGAL_assertion((Arr_halfedge_direction)cv.data().halfedge()->direction() ==
                     ARR_LEFT_TO_RIGHT);
      return Point_2 (m_base(cv.base()), cv.data().halfedge()->target());
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2(this->m_base_tr->
                                  construct_max_vertex_2_object());
  }


  class Compare_xy_2
  {
  private:
    Base_Compare_xy_2 m_base;

  public:

    Compare_xy_2(const Base_Compare_xy_2& base):
        m_base(base)
    {}


    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      const Point_data& inf1 = p1.data();
      const Point_data& inf2 = p2.data();

      if(inf1.vertex() == Vertex_handle() || inf2.vertex() == Vertex_handle())
      {
        return(m_base(p1.base(), p2.base()));
      }

      if(inf1.vertex() == inf2.vertex())
        return (EQUAL);

      return(m_base(p1.base(), p2.base()));
    }
  };


  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2(this->m_base_tr->compare_xy_2_object());
  }



  // left-right

  
  class Parameter_space_in_x_2
  {
  private:
    Base_Parameter_space_in_x_2 m_base;
    
  public:
    
    Parameter_space_in_x_2(const Base_Parameter_space_in_x_2& base):
      m_base(base)
    {}
    
    /*! Obtains the parameter space at the end of a curve-end along the x-axis.
     */
    Arr_parameter_space operator() (const X_monotone_curve_2 & cv, 
                                    const Arr_curve_end& end) const
    {
      return m_base(cv.base(), end);
    }
    
    /*! Obtains the parameter space for a curve along the x-axis.
     */
    Arr_parameter_space operator() (const X_monotone_curve_2 & cv) const
    {
      return m_base(cv.base());
    }

    /*! Obtains the parameter space for a point along the x-axis.
     */
    Arr_parameter_space operator() (const Point_2 & pt) const
    {
      return m_base(pt.base());
    }
  };
  
  /*! Get a Construct_min_vertex_2 functor object. */
  Parameter_space_in_x_2 parameter_space_in_x_2_object () const
  {
    return Parameter_space_in_x_2(
        this->m_base_tr->parameter_space_in_x_2_object()
    );
  }

  class Compare_y_near_boundary_2
  {
  private:
    Base_Compare_y_near_boundary_2 m_base;
    
  public:
    
    Compare_y_near_boundary_2(const Base_Compare_y_near_boundary_2& base):
      m_base(base)
    {}
    
    /*!
     * Compare the relative y-positions of two curve ends.
     */
    Comparison_result operator() (const X_monotone_curve_2 & xcv1,
                                  const X_monotone_curve_2 & xcv2, 
                                  Arr_curve_end ce) const
    {
      return m_base(xcv1, xcv2, ce);
    }
  };
  
  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object () const
  {
    return Compare_y_near_boundary_2(
        this->m_base_tr->compare_y_near_boundary_2_object()
    );
  }

  // TODO Compare_y_on_boundary_2
  // TODO Is_on_x_identification_2

  // bottom-top


  class Parameter_space_in_y_2
  {
  private:
    Base_Parameter_space_in_y_2 m_base;
    
  public:
    
    Parameter_space_in_y_2(const Base_Parameter_space_in_y_2& base):
      m_base(base)
    {}
    
    /*! Obtains the parameter space at the end of a curve-end along the y-axis.
     */
    Arr_parameter_space operator() (const X_monotone_curve_2 & cv, 
                                    const Arr_curve_end& end) const
    {
      return m_base(cv.base(), end);
    }
    
    /*! Obtains the parameter space for a curve along the x-axis.
     */
    Arr_parameter_space operator() (const X_monotone_curve_2 & cv) const
    {
      return m_base(cv.base());
    }

    /*! Obtains the parameter space for a point along the x-axis.
     */
    Arr_parameter_space operator() (const Point_2 & pt) const
    {
      return m_base(pt.base());
    }
    
  };
  
  /*! Get a Construct_min_vertex_2 functor object. */
  Parameter_space_in_y_2 parameter_space_in_y_2_object () const
  {
    return Parameter_space_in_y_2
      (this->m_base_tr->parameter_space_in_y_2_object());
  }

  class Compare_x_near_boundary_2
  {
  private:
    Base_Compare_x_near_boundary_2 m_base;
    
  public:
    
    Compare_x_near_boundary_2(const Base_Compare_x_near_boundary_2& base):
      m_base(base)
    {}
    
    /*! Compare the relative x-positions of a vertical curve and another given
     * curve end.
     */
    Comparison_result operator() (const Point_2 & p,
                                  const X_monotone_curve_2 & xcv,
                                  Arr_curve_end ce) const
    {
      return m_base(p, xcv, ce);
    }

    /*! Compare the relative x-positions of two curve ends.
     */
    Comparison_result operator() (const X_monotone_curve_2 & xcv1,
                                  Arr_curve_end ce1,
                                  const X_monotone_curve_2 & xcv2,
                                  Arr_curve_end ce2) const
    {
      return m_base(xcv1, ce1, xcv2, ce2); 
    }
  };
  
  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object () const
  {
    return Compare_x_near_boundary_2
      (this->m_base_tr->compare_x_near_boundary_2_object());
  }

  // TODO Compare_x_on_boundary_2
  // TODO Is_on_y_identification_2

};

} //namespace CGAL

#endif
