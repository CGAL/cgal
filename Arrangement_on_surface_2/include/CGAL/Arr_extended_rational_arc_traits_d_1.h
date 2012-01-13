// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

//TODO: somehow use the fact the the x-value is the same in all comparisons

#ifndef CGAL_ARR_VERTICAL_SEGMENT_TRAITS
#define CGAL_ARR_VERTICAL_SEGMENT_TRAITS

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/assertions.h>
#include <CGAL/Object.h>

#include <vector>
#include "boost/variant.hpp"

#include <CGAL/Arr_ver_support/Rational_arc_with_ver_d_1.h>

namespace CGAL {

#  define Is_line(ver) \
  (((((ver)).max_parameter_space() == ARR_TOP_BOUNDARY   ) && ((ver).min_parameter_space() == ARR_BOTTOM_BOUNDARY) ))

#  define Is_ray(ver) \
    ((((ver).max_parameter_space() == ARR_TOP_BOUNDARY) && ((ver).min_parameter_space() == ARR_INTERIOR        ))||  \
     (((ver).max_parameter_space() == ARR_INTERIOR    ) && ((ver).min_parameter_space() == ARR_BOTTOM_BOUNDARY ))  )

#  define Is_segment(ver) \
    ( ((ver).max_parameter_space() == ARR_INTERIOR ) && ((ver).min_parameter_space() == ARR_INTERIOR) )

#  define Max_bounded(ver) \
    ( ((ver).max_parameter_space() == ARR_TOP_BOUNDARY ) ? false : true )

#  define Min_bounded(ver) \
    ( ((ver).min_parameter_space() == ARR_BOTTOM_BOUNDARY ) ? false : true )
  

template <class Traits_>
class Arr_traits_with_vertical_segments
{
public:
  typedef Traits_                                   Traits;
  typedef Arr_traits_with_vertical_segments<Traits> Self;
  
  typedef typename Traits::Algebraic_kernel_d_1     Algebraic_kernel_d_1;

  typedef typename Traits::Point_2                  Point_2;
  typedef typename Traits::Multiplicity             Multiplicity;
  typedef typename Traits::Algebraic_real_1         Algebraic_real_1;
  typedef typename Traits::Rat_vector               Rat_vector;
  typedef typename Traits::Coefficient              Coefficient; 

  typedef typename Traits::Integer                  Integer;
  typedef typename Traits::Rational                 Rational; 
  typedef typename Traits::Polynomial_1             Polynomial_1; 
  typedef typename Traits::Rational_function        Rational_function;

  
  typedef typename Traits::X_monotone_curve_2       Non_vertical_x_curve_2;
  typedef typename Traits::Curve_2                  Non_vertical_curve_2;
  typedef typename Traits::Vertical_segment         Vertical_segment;

  typedef Arr_vertical_rational_arc::Rational_arc_with_ver_d_1 
    <Non_vertical_x_curve_2,Algebraic_kernel_d_1>                 X_monotone_curve_2;
  typedef Arr_vertical_rational_arc::Rational_arc_with_ver_d_1 
    <Non_vertical_curve_2,Algebraic_kernel_d_1>                   Curve_2;
  
  //Category tags:
  typedef Tag_true              Has_left_category;
  typedef Tag_true              Has_merge_category;
  typedef Tag_true              Has_do_intersect_category;

  typedef Arr_open_side_tag     Left_side_category;
  typedef Arr_open_side_tag     Top_side_category;
  typedef Arr_open_side_tag     Bottom_side_category;
  typedef Arr_open_side_tag     Right_side_category;
public:
  typedef std::vector<CGAL::Object> Object_vector;
private:
  mutable Traits _traits;
public:
  const Algebraic_kernel_d_1* algebraic_kernel_d_1() const {return _traits.algebraic_kernel_d_1();}

public:
  //------------
  //Constructors
  //------------

  //---------------------
  // Default constructor.
  Arr_traits_with_vertical_segments ()
  {}

  class Construct_point_2
  {
  private:
    Traits& _traits;
  public:
    Construct_point_2 (Traits& traits) :_traits(traits) {}
    Point_2 operator() (const Rational_function& rational_function,
                        const Algebraic_real_1& x_coordinate)
    { 
      return _traits.construct_point_2_object()(rational_function,x_coordinate);
    }
    Point_2 operator() (const Rational& x,const Rational& y)
    { 
      return _traits.construct_point_2_object()(x,y);
    }
    Point_2 operator() (const Algebraic_real_1& x,const Rational& y)
    {
      return _traits.construct_point_2_object()(x,y);
    }    
  }; //Construct_point_2

  Construct_point_2 construct_point_2_object() const {return Construct_point_2(_traits);}

  class Construct_vertical_x_curve_2
  {
  private:
     Traits& _traits;
  public:
    Construct_vertical_x_curve_2(Traits& traits) 
      :_traits(traits)
    {}

    X_monotone_curve_2 operator() (const Point_2& p1,const Point_2& p2) const
    {
      return _traits.construct_vertical_segment_object() (p1,p2);
    }
    X_monotone_curve_2 operator() (const Point_2& p) const
    {
      return _traits.construct_vertical_segment_object() (p);
    }
    X_monotone_curve_2 operator() (const Point_2& p, bool is_directed_up) const
    {
      return _traits.construct_vertical_segment_object() (p,is_directed_up);
    }
  };  //Construct_vertical_x_curve_2

  Construct_vertical_x_curve_2 construct_vertical_x_curve_2_object () const
  {
    return Construct_vertical_x_curve_2(_traits);
  }

  class Construct_vertical_curve_2
  {
  private:
    Traits& _traits;
  public:
    Construct_vertical_curve_2 (Traits& traits)
      :_traits(traits)
    {}
    Curve_2 operator() (const Point_2& p1,const Point_2& p2) const
    {
      return _traits.construct_vertical_segment_object()(p1,p2);
    }
    Curve_2 operator() (const Point_2& p) const
    {
      return _traits.construct_vertical_segment_object()(p);
    }
    Curve_2 operator() (const Point_2& p, bool is_directed_up) const
    {
      return _traits.construct_vertical_segment_object()(p,is_directed_up);
    }
    
  };  //Construct_vertical_curve_2

  Construct_vertical_curve_2 construct_vertical_curve_2_object () const
  {
    return Construct_vertical_curve_2(_traits);
  }

  class Construct_curve_2
  {
  private: 
    typedef typename Traits::Algebraic_real_1 Algebraic_real_1;
    Traits& _traits;
  public:
    Construct_curve_2 (Traits& traits) 
      :_traits(traits)
    {}
    template <class InputIterator>
    Curve_2 operator() (InputIterator begin, InputIterator end) const
    {
      return _traits.construct_curve_2_object()(begin,end); 
    }
    template <class InputIterator>
    Curve_2 operator() (InputIterator begin, InputIterator end,const Algebraic_real_1& x_s, bool dir_right) const
    {
      return _traits.construct_curve_2_object()(begin,end,x_s,dir_right); 
    }
    template <class InputIterator>
    Curve_2 operator() (InputIterator begin, InputIterator end,
                        const Algebraic_real_1& x_s, const Algebraic_real_1& x_t) const
    {
      return _traits.construct_curve_2_object()(begin,end,x_s,x_t); 
    }
    template <class InputIterator>
    Curve_2 operator() (InputIterator begin_numer, InputIterator end_numer,
                        InputIterator begin_denom, InputIterator end_denom) const 
    {
      return _traits.construct_curve_2_object()(begin_numer, end_numer,begin_denom,end_denom); 
    }
    template <class InputIterator>
    Curve_2 operator() (InputIterator begin_numer, InputIterator end_numer,
                        InputIterator begin_denom, InputIterator end_denom,
                        const Algebraic_real_1& x_s, bool dir_right) const
    {
      return _traits.construct_curve_2_object()(begin_numer, end_numer,begin_denom,end_denom,x_s,dir_right); 
    }
    template <class InputIterator>
    Curve_2 operator() (InputIterator begin_numer, InputIterator end_numer,
                        InputIterator begin_denom, InputIterator end_denom,
                        const Algebraic_real_1& x_s, const Algebraic_real_1& x_t) const
    {
      return _traits.construct_curve_2_object()(begin_numer, end_numer,begin_denom,end_denom,x_s,x_t); 
    }
  };  //Construct_rational_curve_2

  Construct_curve_2 construct_curve_2_object () const
  {
    return Construct_curve_2(_traits);
  }
  class Construct_x_monotone_curve_2
  {
  private: 
    typedef typename Traits::Algebraic_real_1 Algebraic_real_1;
    Traits& _traits;
  public:
    Construct_x_monotone_curve_2 (Traits& traits)
      :_traits(traits)
    {}
    template <class InputIterator>
    X_monotone_curve_2 operator() (InputIterator begin, InputIterator end) const
    {
      return _traits.construct_x_monotone_curve_2_object()(begin,end); 
    }
    template <class InputIterator>
    X_monotone_curve_2 operator() ( InputIterator begin, InputIterator end,
                                    const Algebraic_real_1& x_s, bool dir_right) const
    {
      return _traits.construct_x_monotone_curve_2_object()(begin,end,x_s,dir_right); 
    }
    template <class InputIterator>
    X_monotone_curve_2 operator() (InputIterator begin, InputIterator end,
                                   const Algebraic_real_1& x_s, const Algebraic_real_1& x_t) const
    {
      return _traits.construct_x_monotone_curve_2_object()(begin,end,x_s,x_t); 
    }
    template <class InputIterator>
    X_monotone_curve_2 operator() ( InputIterator begin_numer, InputIterator end_numer,
                                    InputIterator begin_denom, InputIterator end_denom) const 
    {
      return _traits.construct_x_monotone_curve_2_object()(begin_numer, end_numer,begin_denom,end_denom); 
    }
    template <class InputIterator>
    X_monotone_curve_2 operator() ( InputIterator begin_numer, InputIterator end_numer,
                                    InputIterator begin_denom, InputIterator end_denom,
                                    const Algebraic_real_1& x_s, bool dir_right) const
    {
      return _traits.construct_x_monotone_curve_2_object()(begin_numer, end_numer,begin_denom,end_denom,x_s,dir_right); 
    }
    template <class InputIterator>
    X_monotone_curve_2 operator() ( InputIterator begin_numer, InputIterator end_numer,
                                    InputIterator begin_denom, InputIterator end_denom,
                                    const Algebraic_real_1& x_s, const Algebraic_real_1& x_t) const
    {
      return _traits.construct_x_monotone_curve_2_object()(begin_numer, end_numer,begin_denom,end_denom,x_s,x_t); 
    }
  };  //Construct_rational_x_curve_2

  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2(_traits);
  }
  //------------------------
  //Functor definitions.
  //------------------------

  //---------------------------------------------------------------
  //A functor that compares the x-coordinates of two points 
  class Compare_x_2
  {
  private:
    Traits& _traits;
  public:
    Compare_x_2(Traits& traits) : _traits(traits) {}
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator() (const Point_2 & p1, const Point_2 & p2) const
    {
      return _traits.compare_x_2_object()(p1,p2);
    }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return Compare_x_2(_traits);
  }

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  class Compare_xy_2
  {
  private:
    Traits& _traits;
  public:
    Compare_xy_2 (Traits& traits) : _traits(traits) {}
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      return _traits.compare_xy_2_object() (p1,p2);
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2(_traits);
  }

  /*! A functor that obtains the left endpoint of a curve. */
  class Construct_min_vertex_2
  {    
  private:
    Traits& _traits;
  public:
    Construct_min_vertex_2(Traits& traits) : _traits(traits) {}
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2 & cv) const
    {
      return (boost::apply_visitor(Construct_min_vertex_2_visitor(_traits),cv.variant()));
    }
  private:
    class Construct_min_vertex_2_visitor
      : public boost::static_visitor <const Point_2&>
    {
    private:
      Traits& _traits;
    public:
      Construct_min_vertex_2_visitor(Traits& traits) :_traits(traits) {}
      const Point_2& operator() (const Non_vertical_x_curve_2 & cv) const
      {
        return _traits.construct_min_vertex_2_object()(cv);
      }
      const Point_2& operator() (const Vertical_segment & cv) const
      {
        return (cv.min());
      }
    };  //Construct_min_vertex_2_visitor
  };  //Construct_min_vertex_2

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2(_traits);
  }

  /*! A functor that obtains the right endpoint of a curve. */
  class Construct_max_vertex_2
  {
  private:
    Traits& _traits;
  public:
    Construct_max_vertex_2 (Traits& traits) : _traits(traits) {}
    const Point_2& operator() (const X_monotone_curve_2 & cv) const
    {
      return (boost::apply_visitor(Construct_max_vertex_2_visitor(_traits),cv.variant()));
    }
  private:
    class Construct_max_vertex_2_visitor
      : public boost::static_visitor <const Point_2&>
    {
    private:
      Traits& _traits;
    public:
      Construct_max_vertex_2_visitor(Traits& traits) : _traits(traits) {}
      const Point_2& operator() (const Non_vertical_x_curve_2 & cv) const
      {
         return _traits.construct_max_vertex_2_object()(cv);
      }
      const Point_2& operator() (const Vertical_segment & cv) const
      {
        return (cv.max());
      }
    };  //Construct_max_vertex_2_visitor
  };  //Construct_max_vertex_2

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2(_traits);
  }

  /*! A functor that checks whether a given curve is vertical. */
  class Is_vertical_2
  {
  private:
    Traits& _traits;
  public:
    Is_vertical_2(Traits& traits) : _traits(traits) {}
    bool operator() (const X_monotone_curve_2& cv) const
    {
      return (boost::apply_visitor(Is_vertical_2_visitor(_traits),cv.variant()));
    }
  private:
    class Is_vertical_2_visitor
      : public boost::static_visitor <bool>
    {
    private:
      Traits& _traits;
    public:
      Is_vertical_2_visitor(Traits& traits) : _traits(traits) {}
      bool operator() (const Non_vertical_x_curve_2 & cv) const
      {
         return _traits.is_vertical_2_object()(cv);
      }
      bool operator() (const Vertical_segment & cv) const
      {
        return true;
      }
    };  //Is_vertical_2_visitor

  };  //Is_vertical_2

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const
  {
    return Is_vertical_2(_traits);
  }

  /*! A functor that compares the y-coordinates of a point and a curve at
   * the point x-coordinate.
   */
  class Compare_y_at_x_2
  {
  private:
    Traits& _traits;    
  public:
    Compare_y_at_x_2(Traits& traits) : _traits(traits) {}
    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    //Compare_y_at_x_2() {}
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      return (boost::apply_visitor(Compare_y_at_x_2_visitor(_traits,p),cv.variant()));
    }
  private:
    class Compare_y_at_x_2_visitor
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_y_at_x_2_visitor(Traits& traits,const Point_2& p) : _traits(traits),_p(p), Base() {}
      Comparison_result operator() (const Non_vertical_x_curve_2 & cv) const
      {
         return _traits.compare_y_at_x_2_object() (_p,cv);
      }
      Comparison_result operator() (const Vertical_segment & cv) const
      {
        CGAL_precondition(_traits.compare_x_2_object()(_p,cv.max()) == CGAL::EQUAL);       

        typename Traits::Compare_xy_2 compare_xy_2 = _traits.compare_xy_2_object();
        
        if (Is_line(cv))
          return CGAL::EQUAL;
        if ((Max_bounded(cv)) && (!Min_bounded(cv)) )
          return (compare_xy_2 (_p,cv.max()) == LARGER) ? LARGER : EQUAL;
        if ((!Max_bounded(cv)) && (Min_bounded(cv)) )
          return (compare_xy_2 (_p,cv.min()) == SMALLER) ? SMALLER : EQUAL;
        //line is a segment
        Comparison_result cr_t = compare_xy_2(_p,cv.max());
        Comparison_result cr_b = compare_xy_2(_p,cv.min());
        return  (cr_t == LARGER) ? LARGER :
                (cr_b == SMALLER) ? SMALLER :
                EQUAL;
      }
    private:
      Traits& _traits;
      Point_2 _p;

    };  //Compare_y_at_x_2_visitor
  };  //Compare_y_at_x_2

  //*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(_traits);
  }

  class Compare_y_at_x_left_2
  {
  private:
    Traits& _traits;
  public:
    Compare_y_at_x_left_2(Traits& traits) : _traits(traits) {}
    Comparison_result operator() (const X_monotone_curve_2& cv1,const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
      {
        return (boost::apply_visitor(Compare_y_at_x_left_2_visitor(_traits,p),cv1.variant(),cv2.variant()));
      }
  private:
    class Compare_y_at_x_left_2_visitor
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_y_at_x_left_2_visitor(Traits& traits,const Point_2& p) 
        : _traits(traits), _p(p), Base() {}
      Comparison_result operator() (const Non_vertical_x_curve_2 & cv1,
                                    const Non_vertical_x_curve_2 & cv2) const
      {
         return _traits.compare_y_at_x_left_2_object() (cv1,cv2,_p);
      }
      Comparison_result operator() (const Vertical_segment & cv1,
                                    const Non_vertical_x_curve_2 & cv2) const
      {
        assert(false);
        return CGAL::SMALLER;
      }
      Comparison_result operator() (const Non_vertical_x_curve_2& cv1,
                                    const Vertical_segment  & cv2) const
      {
        assert(false);
        return CGAL::LARGER;
      }
      Comparison_result operator() (const Vertical_segment & cv1,
                                    const Vertical_segment & cv2) const
      {
        return CGAL::EQUAL;
      }
    private:
      Traits& _traits;
      Point_2 _p;

    };  //Compare_y_at_x_left_2_visitor
  }; //Compare_y_at_x_left_2

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const {return Compare_y_at_x_left_2(_traits);}
  /*! A functor that compares compares the y-coordinates of two curves
   * immediately to the left of their intersection point.
   */
  
  /*! A functor that checks whether two points and two curves are identical. */
  class Equal_2
  {
  private:
    Traits& _traits;
  public:
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    Equal_2 (Traits& traits) :_traits(traits) {}
    bool operator() ( const X_monotone_curve_2& cv1,
                      const X_monotone_curve_2& cv2) const
    {
      return (boost::apply_visitor(Equal_2_visitor(_traits),cv1.variant(),cv2.variant()));
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      if (&p1 == &p2)
        return (true);
      return _traits.equal_2_object() (p1,p2);
   
    }
  private:
    class Equal_2_visitor
      : public boost::static_visitor <bool>
    {
    private:
      typedef boost::static_visitor <bool> Base;      
      Traits& _traits;
    public:
      Equal_2_visitor(Traits& traits) : _traits(traits), Base() {}
      bool operator() (const Non_vertical_x_curve_2& cv1,
                       const Non_vertical_x_curve_2& cv2) const      
      {
         return _traits.equal_2_object()  (cv1,cv2);
      }
      bool operator() (const Non_vertical_x_curve_2& cv1,
                       const Vertical_segment & cv2) const      
      {
        return false; 
      }
      bool operator() (const Vertical_segment & cv1,
                       const Non_vertical_x_curve_2& cv2) const      
      {
        return false; 
      }
      bool operator() (const Vertical_segment & cv1,
                       const Vertical_segment & cv2) const      
      {
        if (&cv1 == &cv2)
          return (true);

        if (CGAL::compare(cv1.x(),cv2.x()) != CGAL::EQUAL)
          return false;

        if (Is_line(cv1) && Is_line(cv2))
          return true;

        if (Is_segment(cv1) && Is_segment(cv2))
          return (_traits.equal_2_object() (cv1.min(),cv2.min()) &&
                  _traits.equal_2_object() (cv1.max(),cv2.max()) );

        if (Max_bounded(cv1) && Max_bounded(cv2))
          return (_traits.equal_2_object() (cv1.max(),cv2.max() ));

        if (Min_bounded(cv1) && Min_bounded(cv2))
          return (_traits.equal_2_object() (cv1.min(),cv2.min() ));
        
        return false;
      }
    };  //Equal_2_visitor
  };  //Equal_2

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return Equal_2(_traits);
  }

  /*! A functor that divides a curve into continues (x-monotone) curves. */
  class Make_x_monotone_2
  {
  private:    
    Traits& _traits;
  public:
    Make_x_monotone_2(Traits& traits) : _traits(traits) {}
  template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const 
    {
      Object_vector res (boost::apply_visitor(Make_x_monotone_2_visitor(_traits),cv.variant()));
      re_cast_object_vector(res,oi);
      return (oi);
    }
  private:
    class Make_x_monotone_2_visitor
      : public boost::static_visitor < Object_vector >
    {
    private:
      typedef boost::static_visitor <Object_vector> Base;      
      Traits& _traits;
    public:
      Make_x_monotone_2_visitor(Traits& traits) : _traits(traits), Base() {}
      Object_vector operator() (const Non_vertical_curve_2& cv) const      
      {
        Object_vector vec;
        _traits.make_x_monotone_2_object() (cv,std::back_inserter(vec));
        return vec;
      }
      Object_vector operator() (const Vertical_segment & cv) const      
      {
        Object_vector vec;
        vec.push_back(make_object (Vertical_segment(cv)));
        return vec;
      }
    };  //Make_x_monotone_2_visitor

  };  //Make_x_monotone_2

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return Make_x_monotone_2(_traits);
  }

  /*! A functor that splits a curve at a point. */
  class Split_2
  {
  private:
    typedef std::pair<X_monotone_curve_2,X_monotone_curve_2> Res_type;
    Traits& _traits;
  public:
    Split_2 (Traits& traits) : _traits(traits) {}
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator() ( const X_monotone_curve_2& cv, const Point_2 & p,
                      X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
       Res_type res = boost::apply_visitor(Split_2_visitor(_traits,p),cv.variant());
       c1 = res.first;
       c2 = res.second;
    }
  private:
    class Split_2_visitor
      : public boost::static_visitor <Res_type>
    {
    public:
      typedef boost::static_visitor <Res_type> Base;
      Split_2_visitor(Traits& traits,const Point_2 & p) 
                      : _traits(traits),_p(p), Base() {}
      Res_type operator() (const Non_vertical_x_curve_2& cv) const
      {
        Non_vertical_x_curve_2 _c1,_c2;
        _traits.split_2_object() (cv,_p,_c1,_c2);
        return Res_type(_c1,_c2);
      }
      Res_type operator() (const Vertical_segment & cv) const      
      {
        typename Traits::Construct_vertical_segment 
          construct_vertical_segment = _traits.construct_vertical_segment_object();
        Vertical_segment _c1,_c2;
        if (Is_line(cv))
        {
          _c1 = construct_vertical_segment(_p,false);
          _c2 = construct_vertical_segment(_p,true);
        }
        else if (Is_segment(cv))
        {
          _c1 = construct_vertical_segment(cv.min(),_p);
          _c2 = construct_vertical_segment(_p,cv.max());
        }
        else if (Max_bounded(cv)) //min unbounded
        {
          _c1 = construct_vertical_segment(_p,false);
          _c2 = construct_vertical_segment(_p,cv.max());
        }
        else //min_bounded && max unbounded
        {
          _c1 = construct_vertical_segment(cv.min(),_p);
          _c2 = construct_vertical_segment(_p,true);
        }
        return Res_type(_c1,_c2);
      }
    private:
      Traits& _traits;
      Point_2 _p;
    };  //Split_2_visitor
  };  //Split_2

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2(_traits);
  }

  /*! A functor that computes intersections between two curves. */
  class Intersect_2
  {
  private:
    Traits& _traits;    
  public:
    /*!
     * Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    Intersect_2 (Traits& traits) : _traits(traits) {}
    template<class OutputIterator>
    OutputIterator operator() ( const X_monotone_curve_2& cv1,
                                const X_monotone_curve_2& cv2,
                                OutputIterator oi) const 
    {
      Object_vector res (boost::apply_visitor(Intersect_2_visitor(_traits),cv1.variant(),cv2.variant()));
      re_cast_object_vector(res,oi);
      return (oi);
    }
  private:
    class Intersect_2_visitor
      : public boost::static_visitor <Object_vector>
    {
    private:
      typedef boost::static_visitor <Object_vector> Base;
      Traits& _traits;
    public:
      Intersect_2_visitor(Traits& traits) 
                      : _traits(traits), Base() {}

      //intersection of two Non_vertical_x_curve_2
      Object_vector operator() (const Non_vertical_x_curve_2& cv1,
                                const Non_vertical_x_curve_2& cv2) const      
      {
        Object_vector vec;
        _traits.intersect_2_object() (cv1,cv2,std::back_inserter(vec));
        return vec;
      }
      //intersection of a Non_vertical_x_curve_2 and a Vertical_segment
      Object_vector operator() (const Vertical_segment& cv1,
                                const Non_vertical_x_curve_2& cv2) const      
      {
        Object_vector vec;
        _traits.intersect_2_object() (cv2,cv1,std::back_inserter(vec));
        return vec;
      }
      //intersection of a Non_vertical_x_curve_2 and a Vertical_segment
      Object_vector operator() (const Non_vertical_x_curve_2& cv1,
                                const Vertical_segment& cv2) const      
      {
        Object_vector vec;
        _traits.intersect_2_object() (cv1,cv2,std::back_inserter(vec));
        return vec;
      }
      //intersection of two Vertical_segment
      Object_vector operator() (const Vertical_segment& cv1,
                                const Vertical_segment& cv2) const      
      {
        Object_vector vec;
        CGAL::Object object;
        if (cv1.has_same_x(cv2) == false)
          return vec;
        if (Is_line(cv1))
          object = make_object(Vertical_segment(cv2));
        else if (Is_line(cv2))
          object = make_object(Vertical_segment(cv1));
        
        if (object.empty() == false)
        {
          vec.push_back(object);
          return vec;
        }

        if (Is_ray(cv1))
        {
          if (Is_ray(cv2))
            vec = intersect_ray_ray (cv1,cv2);
          else //cv2 is a segment
            vec = intersect_ray_segment(cv1,cv2);
        }
        //cv1 is a segment
        else if (Is_ray(cv2))
        {
          vec = intersect_ray_segment(cv2,cv1);
        }
        //cv1 & cv2 are segments
        else
          vec = intersect_segment_segment(cv1,cv2);

        return vec;
      }
    private:
      //returns the intersection of a curve with a vertical line defined at p
      //CGAL::Object wraps a pair<ArrTraits::Point_2,ArrTraits::Multiplicity>
      //post return_value.second = 1
      Object_vector intersect_ray_ray ( const Vertical_segment& ray1,
                                        const Vertical_segment& ray2) const
      {
        typename Traits::Compare_xy_2 
          compare_xy_2 = _traits.compare_xy_2_object();
        typename Traits::Construct_vertical_segment 
          construct_vertical_segment = _traits.construct_vertical_segment_object();

        CGAL_precondition (Is_ray(ray1));
        CGAL_precondition (Is_ray(ray2));

        Object_vector vec;
        CGAL::Object object;


        if (Min_bounded(ray1) && Min_bounded(ray2))
        {
          //both are directed up, take the larger minimum vertex
          Comparison_result cr = compare_xy_2(ray1.min(),ray2.min());
          object  = (cr == LARGER) ?  make_object (construct_vertical_segment(ray1.min(),true)):
                                      make_object (construct_vertical_segment(ray2.min(),true))  ;
        }
        else if (Min_bounded(ray1) && Max_bounded(ray2))
        {
          //ray1 is directed up, ray2 is directed down
          //either theu have no intersection or an intersecting segment or an intersection point
          Comparison_result cr = compare_xy_2 (ray1.min(),ray2.max());
          if (cr == CGAL::LARGER)
            return vec;
          else if (cr == CGAL::EQUAL)
            object = make_object (std::pair<Point_2,Multiplicity> (ray1.min(),0));
          else  //(cr == CGAL::SMALLER)
            object = make_object (construct_vertical_segment(ray1.min(),ray2.max()));
        }
        else if (Max_bounded(ray1) && Min_bounded(ray2))
        {
          //ray1 is directed down, ray2 is directed up
          //either theu have no intersection or an intersecting segment
          Comparison_result cr = compare_xy_2 (ray1.min(),ray2.max());
          if (cr == CGAL::SMALLER)
            return vec;
          else if (cr == CGAL::EQUAL)
            object = make_object (std::pair<Point_2,Multiplicity> (ray1.min(),0));
          else
            object = make_object (construct_vertical_segment(ray1.max(),ray2.min()));
        }
        else
        {
          //both are directed down, take the smaller maximum vertex
          Comparison_result cr = compare_xy_2 (ray1.max(),ray2.max());
          object  = (cr == SMALLER) ? make_object (construct_vertical_segment(ray1.max(),false)):
                                      make_object (construct_vertical_segment(ray2.max(),false))  ;
        }
        vec.push_back(object);
        return vec;
      }
      Object_vector intersect_ray_segment(const Vertical_segment& ray,
                                          const Vertical_segment& seg) const
      {
        CGAL_precondition (Is_ray(ray));
        CGAL_precondition (Is_segment(seg));
        
        typename Traits::Compare_xy_2 
          compare_xy_2 = _traits.compare_xy_2_object();
        typename Traits::Construct_vertical_segment 
          construct_vertical_segment = _traits.construct_vertical_segment_object();

        Object_vector vec;
        CGAL::Object object;

        if (Min_bounded(ray))
        {
          //ray is directed up
          Comparison_result cr = compare_xy_2 (ray.min(),seg.min());
          if (cr != LARGER)
            object = make_object (seg);
          else 
          {
            cr = compare_xy_2 (ray.min(),seg.max());
            if (cr == LARGER)
              return vec;
            else if (cr == CGAL::EQUAL)
              object = make_object (std::pair<Point_2,Multiplicity> (ray.min(),0));
            else
              object = make_object (construct_vertical_segment(ray.min(),seg.max()));
          }
        }
        else
        {
          //ray is directed down
          Comparison_result cr = compare_xy_2 (ray.max(),seg.max());
          if (cr != SMALLER)
            object = make_object (seg);
          else 
          {
            cr = compare_xy_2 (ray.max(),seg.min());
            if (cr == SMALLER)
              return vec;
            else if (cr == CGAL::EQUAL)
              object = make_object (std::pair<Point_2,Multiplicity> (ray.max(),0));
            else
              object = make_object (construct_vertical_segment(ray.min(),seg.max()));
          }
        }
         vec.push_back(object);
        return vec;
      }

      Object_vector intersect_segment_segment(const Vertical_segment& seg1,
                                              const Vertical_segment& seg2) const
      {
        CGAL_precondition (Is_segment(seg1));
        CGAL_precondition (Is_segment(seg2));

        typename Traits::Compare_xy_2 
          compare_xy_2 = _traits.compare_xy_2_object();
        typename Traits::Construct_vertical_segment 
          construct_vertical_segment = _traits.construct_vertical_segment_object();
        
        Object_vector vec;
        CGAL::Object object;

        Comparison_result cr = compare_xy_2 (seg1.max(),seg2.min());
        if (cr == CGAL::SMALLER)
          return vec;

        cr = typename compare_xy_2 (seg1.min(),seg2.max());
        if (cr == CGAL::LARGER)
          return vec;

        Comparison_result cr1 = compare_xy_2 (seg1.min(),seg2.min());
        Comparison_result cr2 = compare_xy_2 (seg1.max(),seg2.max());

        if (  (cr1 == CGAL::LARGER ) && 
              (cr2 == CGAL::LARGER ) )
          object = make_object (construct_vertical_segment(seg1.min(),seg2.max()));
        else if (  (cr1 == CGAL::LARGER ) && 
                   (cr2 == CGAL::SMALLER ) )
          object = make_object (construct_vertical_segment(seg1.min(),seg1.max()));
        else if (  (cr1 == CGAL::SMALLER) && 
                   (cr2 == CGAL::LARGER ) )
          object = make_object (construct_vertical_segment(seg2.min(),seg2.max()));
        else  //SMALLER,SMALLER
          object = make_object (construct_vertical_segment(seg2.min(),seg1.max()));
       
        vec.push_back(object);
        return vec;
      }
    }; //Intersect_2_visitor
  };  //Intersect_2

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2(_traits);
  }

  /*! A functor that tests whether two curves can be merged. */
  class Are_mergeable_2
  {
  private:
    Traits& _traits;
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    Are_mergeable_2 (Traits& traits) : _traits(traits) {}
    bool operator() ( const X_monotone_curve_2& cv1,
                      const X_monotone_curve_2& cv2) const
    {
      return (boost::apply_visitor(Are_mergeable_2_visitor(_traits),cv1.variant(),cv2.variant()));
    }
  private:
    class Are_mergeable_2_visitor
      : public boost::static_visitor <bool>
    {
    private:
      typedef boost::static_visitor <bool> Base;      
      Traits& _traits;
    public:
      Are_mergeable_2_visitor(Traits& traits) 
                      : _traits(traits), Base() {}
      bool operator() ( const Non_vertical_x_curve_2& cv1,
                        const Non_vertical_x_curve_2& cv2) const
      {
        return  _traits.are_mergeable_2_object() (cv1,cv2);
      }
      bool operator() ( const Vertical_segment& cv1,
                        const Non_vertical_x_curve_2& cv2) const
      {
        return false;
      }
      bool operator() ( const Non_vertical_x_curve_2& cv1,
                        const Vertical_segment& cv2) const
      {
        return false;
      }
      bool operator() ( const Vertical_segment& cv1,
                        const Vertical_segment& cv2) const
      {
        bool res = false;
        //try to merge at common maximum
        if ((Max_bounded(cv1)) && (Max_bounded(cv2)))
        {
          if (_traits.equal_2_object()(cv1.max(),cv2.max()))
            return true;
        }
        //try to merge at common minimum
        if ((Min_bounded(cv1)) && (Min_bounded(cv2)))
        {
          if (_traits.equal_2_object()(cv1.min(),cv2.min()))
            return true;
        }
        //try to merge at minimum and maximum 
        if ((Max_bounded(cv1)) && (Min_bounded(cv2)))
        {
          res = (_traits.equal_2_object()(cv1.max(),cv2.min()));
        }
        if ((Min_bounded(cv1)) && (Max_bounded(cv2)))
        {
          res = (res || (_traits.equal_2_object()(cv1.min(),cv2.max())));
        }
        return res;
      }
    };  //Are_mergeable_2_visitor
  };  //Are_mergeable_2

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2(_traits);
  }

  /*! A functor that merges two curves into one. */
  class Merge_2
  {
  private:
    Traits& _traits;
  public:
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same conic curve and share a common endpoint.
     */
    Merge_2 (Traits& traits) : _traits(traits) {}
    void operator() ( const X_monotone_curve_2& cv1,
                      const X_monotone_curve_2& cv2,
                      X_monotone_curve_2& c) const
    {
      c = boost::apply_visitor(Merge_2_visitor(_traits),cv1.variant(),cv2.variant()); 
      return;
    }
  private:
    class Merge_2_visitor
      : public boost::static_visitor <X_monotone_curve_2>
    {
    private:
      typedef boost::static_visitor <X_monotone_curve_2> Base;      
     Traits& _traits;
    public:
      Merge_2_visitor(Traits& traits) : _traits(traits), Base() {}
      X_monotone_curve_2 operator() (const Non_vertical_x_curve_2& cv1,const Non_vertical_x_curve_2& cv2) const
      {
        CGAL_precondition(_traits.are_mergeable_2_object()(cv1,cv2));
        Non_vertical_x_curve_2 _c;
        _traits.merge_2_object()(cv1,cv2,_c);
        return _c;
      }
      X_monotone_curve_2 operator() (const Non_vertical_x_curve_2& cv1,const Vertical_segment& cv2) const
      {
        CGAL_precondition(false);
        return X_monotone_curve_2();
      }
      X_monotone_curve_2 operator() (const Vertical_segment& cv1,const Non_vertical_x_curve_2& cv2) const
      {
        CGAL_precondition(false);
        return X_monotone_curve_2();
      }
      X_monotone_curve_2 operator() (const Vertical_segment& cv1,const Vertical_segment& cv2) const
      {
        typename Traits::Construct_vertical_segment 
          construct_vertical_segment = _traits.construct_vertical_segment_object();

        typename Traits::Compare_xy_2 
          compare_xy_2 = _traits.compare_xy_2_object();

        Vertical_segment _c;
        if ((Max_bounded(cv1)) && 
            (Max_bounded(cv2)) &&
            (compare_xy_2(cv1.max() , cv2.max()) == CGAL::EQUAL))
        {
          //merge at cv1.max() and cv2.max()
            if ((Min_bounded(cv1) == false) || (Min_bounded(cv2) == false))
              _c = construct_vertical_segment(cv1.max(),false); //vertical ray towards -oo from common max
            //construct vertical ray with common max and minimal min
            bool get_cv1_min = (compare_xy_2(cv1.min(),cv2.min()) == CGAL::SMALLER);
            _c = construct_vertical_segment(cv1.max(), get_cv1_min ? cv1.min() : cv2.min());
            return _c;
        }

        if ((Min_bounded(cv1)) && 
            (Min_bounded(cv2)) &&
            (compare_xy_2(cv1.min() , cv2.min()) == CGAL::EQUAL))
        {
          //merge at cv1.min() and cv2.min()
            if ((Max_bounded(cv1) == false) || (Max_bounded(cv2) == false))
              _c = construct_vertical_segment(cv1.min(),true); //vertical ray towards oo from common min
            //construct vertical ray with common min and maximal max
            bool get_cv1_max = (compare_xy_2(cv1.max(),cv2.max()) == CGAL::LARGER);
            _c = construct_vertical_segment(cv1.min(), get_cv1_max ? cv1.max() : cv2.max());
            return _c;
        }

        if ((Max_bounded(cv1)) && 
            (Min_bounded(cv2)) &&
            (compare_xy_2(cv1.max() , cv2.min()) == CGAL::EQUAL))
          {
            //merge at cv1.max() and cv2.min
            if ((Min_bounded(cv1) == false) && (Max_bounded(cv2) == false))
              _c = construct_vertical_segment(cv1.max()); //vertical line at max.x
            else if ((Min_bounded(cv1) == true) && (Max_bounded(cv2) == false))
              _c = construct_vertical_segment(cv1.min(),true); //vertical ray towards +oo from cv1.min
            else if ((Min_bounded(cv1) == false) && (Max_bounded(cv2) == true))
              _c = construct_vertical_segment(cv2.max(),false); //vertical ray towards -oo from cv2.max
            else
              _c = construct_vertical_segment(cv1.min(), cv2.max());
            return _c;
          }
        else
          {
            //merge at cv1.min() and cv2.max
            if ((Max_bounded(cv1) == false) && (Min_bounded(cv2) == false))
              _c = construct_vertical_segment(cv1.max()); //vertical line at max.x
            else if ((Max_bounded(cv1) == true) && (Min_bounded(cv2) == false))
              _c = construct_vertical_segment(cv1.max(),false); //vertical ray towards -oo from cv1.max
            else if ((Max_bounded(cv1) == false) && (Min_bounded(cv2) == true))
              _c = construct_vertical_segment(cv2.min(),true); //vertical ray towards +oo from cv2.min
            else _c = construct_vertical_segment(cv2.min(), cv1.max());
            return _c;
          }
        CGAL_postcondition(false);
        return X_monotone_curve_2();
      }
    };  //Merge_2_visitor
  };  //Merge_2

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2(_traits);
  }

  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2
  {
  private:
    Traits& _traits;
  public:
    Parameter_space_in_x_2(Traits& traits) :_traits(traits) {};
    /*! Obtains the parameter space at the end of a line along the x-axis    */
    Arr_parameter_space operator()( const X_monotone_curve_2 & xcv,
                                    Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Parameter_space_in_x_2_visitor(_traits,ce),xcv.variant()));
    }
      
    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    /*Arr_parameter_space operator()(const Point_2 p) const
    {
      return Traits_2::Parameter_space_in_x_2(p);
    }*/
  private:
    class Parameter_space_in_x_2_visitor
      : public boost::static_visitor <Arr_parameter_space>
    {
    public:
      typedef boost::static_visitor <Arr_parameter_space> Base;
      Parameter_space_in_x_2_visitor(Traits& traits,Arr_curve_end ce) 
                      :_traits(traits), _ce(ce), Base() {}
      Arr_parameter_space operator()(const Non_vertical_x_curve_2& xcv ) const
      {
        return _traits.parameter_space_in_x_2_object()(xcv,_ce);
      }
      Arr_parameter_space operator()(const Vertical_segment& xcv) const
      {
        return ARR_INTERIOR;
      }

    private:
      Traits& _traits;
      Arr_curve_end _ce;
    };  //Parameter_space_in_x_2_visitor
  };  //Parameter_space_in_x_2

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { 
    return Parameter_space_in_x_2(_traits); 
  }
  
  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 
  {
  private:
    Traits& _traits;
  public:
    Parameter_space_in_y_2 (Traits& traits) :_traits(traits) {};
    /*! Obtains the parameter space at the end of a line along the y-axis .
     * Note that if the line end coincides with a pole, then unless the line
     * coincides with the identification arc, the line end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the line coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    Arr_parameter_space operator()( const X_monotone_curve_2 & xcv,
                                    Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Parameter_space_in_y_2_visitor(_traits,ce),xcv.variant()));
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    //Arr_parameter_space operator()(const Point_2 ) const
    //{
    //  return ARR_INTERIOR;
    //}
      
  private:
    class Parameter_space_in_y_2_visitor
      : public boost::static_visitor <Arr_parameter_space>
    {
    public:
      typedef boost::static_visitor <Arr_parameter_space> Base;
      Parameter_space_in_y_2_visitor(Traits& traits,Arr_curve_end ce) 
        :_traits(traits),_ce(ce), Base() {}
      Arr_parameter_space operator()(const Non_vertical_x_curve_2& xcv ) const
      {
        return _traits.parameter_space_in_y_2_object()(xcv,_ce);
      }
      Arr_parameter_space operator()(const Vertical_segment& xcv) const
      {
        if (_ce ==  ARR_MIN_END)
          return xcv.min_parameter_space();
        else  //ce ==  ARR_MAX_END
          return xcv.max_parameter_space();
      }
    private:
      Traits& _traits;
      Arr_curve_end _ce;
    };  //Parameter_space_in_x_2_visitor

  };  //Parameter_space_in_y_2

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(_traits); }

  class Compare_y_at_x_right_2
  {
  private:
    Traits& _traits;
  public:
     /* Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */    
    Compare_y_at_x_right_2 (Traits& traits) : _traits(traits) {}
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      return (boost::apply_visitor(Compare_y_at_x_right_2_visitor(_traits,p),cv1.variant(),cv2.variant()));
    }
  private:
    class Compare_y_at_x_right_2_visitor
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_y_at_x_right_2_visitor(Traits& traits,const Point_2& p) 
        : _traits(traits), _p(p), Base() {}
      Comparison_result operator() (const Non_vertical_x_curve_2& cv1,
                                    const Non_vertical_x_curve_2& cv2) const      
      {
         return _traits.compare_y_at_x_right_2_object() (cv1,cv2,_p);
      }
      Comparison_result operator() (const Non_vertical_x_curve_2& cv1,
                                    const Vertical_segment & cv2) const      
      {
        return CGAL::SMALLER;
      }
      Comparison_result operator() (const Vertical_segment & cv1,
                                    const Non_vertical_x_curve_2& cv2) const      
      {
        return CGAL::LARGER;
      }
      Comparison_result operator() (const Vertical_segment & cv1,
                                    const Vertical_segment & cv2) const      
      {
        return CGAL::EQUAL; //test bug fix        
      }
    private:
      Traits& _traits;
      Point_2 _p;

    };  //Compare_y_at_x_right_2_visitor
  };  //Compare_y_at_x_right_2

  /*! Obtain a Compare_y_at_x_right_2functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2(_traits);
  }
  class Compare_x_near_limit_2 
  {
  private:
    Traits& _traits;
  public:
    Compare_x_near_limit_2 (Traits& traits) : _traits(traits) {}
    Comparison_result operator()( const X_monotone_curve_2& xcv1, 
                                  const X_monotone_curve_2& xcv2, 
                                  Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Compare_x_near_limit_2_visitor(_traits,ce),xcv1.variant(),xcv2.variant()));
    }
  private:
    class Compare_x_near_limit_2_visitor
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_x_near_limit_2_visitor(Traits& traits,Arr_curve_end ce) 
        :_traits(traits), _ce(ce), Base() {}
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv1,const Non_vertical_x_curve_2& xcv2 ) const
      {
        return _traits.compare_x_near_limit_2_object()(xcv1,xcv2,_ce);
      }
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv1,const Vertical_segment& xcv2 ) const
      {        
        if (_ce == ARR_MIN_END)
        {
          CGAL_precondition (xcv2.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
          return CGAL::LARGER;
        }
        else //_ce2 == ARR_MAX_END
        {
          CGAL_precondition (xcv2.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
          return CGAL::SMALLER;
        }

        return _traits.compare_x_at_limit_2_object()(xcv2.max(),xcv1,_ce);
      }
      Comparison_result operator()(const Vertical_segment& xcv1,const Non_vertical_x_curve_2& xcv2 ) const
      {
        if (_ce == ARR_MIN_END)
        {
          CGAL_precondition (xcv1.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
          return CGAL::SMALLER;
        }
        else //_ce1 == ARR_MAX_END
        {
          CGAL_precondition (xcv1.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
          return CGAL::LARGER;
        }
      }
      Comparison_result operator()(const Vertical_segment& xcv1,const Vertical_segment& xcv2 ) const
      {
        if (_ce == ARR_MIN_END)
          CGAL_precondition (xcv1.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce1 == ARR_MAX_END
          CGAL_precondition (xcv1.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
         
        if (_ce == ARR_MIN_END)
          CGAL_precondition (xcv2.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce2 == ARR_MAX_END
          CGAL_precondition (xcv2.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
 
        return _traits.compare_x_2_object() (xcv1.min(),xcv2.min());
      }

   private:
     Traits& _traits;
     Arr_curve_end _ce;
    };  //Compare_x_near_limit_2_visitor
  }; //Compare_x_near_limit_2

  Compare_x_near_limit_2 compare_x_near_limit_2_object() const 
  {
    return Compare_x_near_limit_2(_traits);
  }

  class Compare_x_at_limit_2 
  {
  private:
    Traits& _traits;
  public:
    Compare_x_at_limit_2 (Traits& traits) :_traits(traits) {}
    Comparison_result operator()( const Point_2 & p,
                                  const X_monotone_curve_2 & xcv,
                                  Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Compare_x_at_limit_2_visitor_1(_traits,ce,p),xcv.variant()));
    }
    Comparison_result operator()( const X_monotone_curve_2 & xcv1,
                                  Arr_curve_end ce1,
                                  const X_monotone_curve_2 & xcv2,
                                  Arr_curve_end ce2) const
    {
      return (boost::apply_visitor(Compare_x_at_limit_2_visitor_2(_traits,ce1,ce2),xcv1.variant(),xcv2.variant()));
    }
  private:
    class Compare_x_at_limit_2_visitor_1
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_x_at_limit_2_visitor_1(Traits& traits,Arr_curve_end ce,const Point_2 & p) 
        :_traits(traits), _ce(ce), _p(p), Base() {}
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv ) const
      {
        return _traits.compare_x_at_limit_2_object()(_p,xcv,_ce);
      }
      Comparison_result operator()(const Vertical_segment& xcv) const
      {
        if (_ce == ARR_MIN_END)
          CGAL_precondition (xcv.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce == ARR_MAX_END
          CGAL_precondition (xcv.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
        return _traits.compare_x_2_object()(_p,xcv.min());
      }
    private:
      Traits& _traits;
      Point_2 _p;
      Arr_curve_end _ce;
    };  //Compare_x_at_limit_2_visitor_1
        
    class Compare_x_at_limit_2_visitor_2
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_x_at_limit_2_visitor_2(Traits& traits,Arr_curve_end ce1,Arr_curve_end ce2) 
        :_traits(traits), _ce1(ce1), _ce2(ce2), Base() {}
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv1,const Non_vertical_x_curve_2& xcv2 ) const
      {
        return _traits.compare_x_at_limit_2_object()(xcv1,_ce1,xcv2,_ce2);
      }
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv1,const Vertical_segment& xcv2 ) const
      {
        if (_ce2 == ARR_MIN_END)
          CGAL_precondition (xcv2.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce2 == ARR_MAX_END
          CGAL_precondition (xcv2.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);

        return _traits.compare_x_at_limit_2_object()(xcv2.max(),xcv1,_ce1);
      }
      Comparison_result operator()(const Vertical_segment& xcv1,const Non_vertical_x_curve_2& xcv2 ) const
      {
        if (_ce1 == ARR_MIN_END)
          CGAL_precondition (xcv1.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce1 == ARR_MAX_END
          CGAL_precondition (xcv1.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);

        return _traits.compare_x_at_limit_2_object()(xcv1.max(),xcv2,_ce2);
      }
      Comparison_result operator()(const Vertical_segment& xcv1,const Vertical_segment& xcv2 ) const
      {
        if (_ce1 == ARR_MIN_END)
          CGAL_precondition (xcv1.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce1 == ARR_MAX_END
          CGAL_precondition (xcv1.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
         
        if (_ce2 == ARR_MIN_END)
          CGAL_precondition (xcv2.min_parameter_space() == CGAL::ARR_BOTTOM_BOUNDARY);
        else //_ce2 == ARR_MAX_END
          CGAL_precondition (xcv2.max_parameter_space() == CGAL::ARR_TOP_BOUNDARY);
 
        return _traits.compare_x_2_object() (xcv1.min(),xcv2.min());
      }
    private:
      Traits& _traits;
      Arr_curve_end _ce1,_ce2;
    };  //Compare_x_at_limit_2_visitor_2
  };  //Compare_x_at_limit_2

  /*! Obtain a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(_traits); }
  /*! A function object that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 
  {
  private:
    Traits& _traits;
  public:
    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    Compare_y_near_boundary_2 (Traits& traits) :_traits(traits) {}
    Comparison_result operator()( const X_monotone_curve_2 & xcv1,
                                  const X_monotone_curve_2 & xcv2,
                                  Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Compare_y_near_boundary_2_visitor(_traits,ce),xcv1.variant(),xcv2.variant()));
    }
  private:
    class Compare_y_near_boundary_2_visitor
      : public boost::static_visitor <Comparison_result>
    {
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Compare_y_near_boundary_2_visitor(Traits& traits,Arr_curve_end ce) 
        :_traits(traits), _ce(ce), Base() {}
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv1,const Non_vertical_x_curve_2& xcv2 ) const
      {
        return _traits.compare_y_near_boundary_2_object()(xcv1,xcv2,_ce);
      }
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv1,const Vertical_segment& xcv2 ) const
      {
        //a vertical curve is not on the boundary at x = +-oo
        CGAL_precondition (false);
        return Comparison_result();
      }
      Comparison_result operator()(const Vertical_segment& xcv1,const Non_vertical_x_curve_2& xcv2 ) const
      {
        //a vertical curve is not on the boundary at x = +-oo
        CGAL_precondition (false);
        return Comparison_result();
      }
      Comparison_result operator()(const Vertical_segment& xcv1,const Vertical_segment& xcv2 ) const
      {
        //a vertical curve is not on the boundary at x = +-oo
        CGAL_precondition (false);
        return Comparison_result();
      }
    private:
      Traits&       _traits;
      Arr_curve_end _ce;
    };  //Compare_y_near_boundary_2_visitor
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(_traits); }
  //@}
  
  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
  class Compare_endpoints_xy_2
  {
  public:

    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv)
    {
      return (boost::apply_visitor(Compare_endpoints_xy_2_visitor(),cv.variant()));
    }
    private:
    class Compare_endpoints_xy_2_visitor
      : public boost::static_visitor <X_monotone_curve_2>
    {
    public:
      Comparison_result operator()(const Non_vertical_x_curve_2& xcv ) const
      {
        return _traits.construct_opposite_2_object()(xcv);
      }
      Comparison_result operator()(const Vertical_segment& xcv) const
      {
        return (SMALLER);
      }
    };  //Compare_endpoints_xy_2_visitor
  }; //Compare_endpoints_xy_2

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  {
    return Compare_endpoints_xy_2();
  }

  class Construct_opposite_2
  {
  public:
    X_monotone_curve_2 operator() (const X_monotone_curve_2& cv) const
    {
      return (boost::apply_visitor(Construct_opposite_2_visitor(),cv.variant()));
    }
  private:
    class Construct_opposite_2_visitor
      : public boost::static_visitor <X_monotone_curve_2>
    {
    public:
      X_monotone_curve_2 operator()(const Non_vertical_x_curve_2& xcv ) const
      {
        return _traits.construct_opposite_2_object()(xcv);
      }
      X_monotone_curve_2 operator()(const Vertical_segment& xcv) const
      {
        CGAL_precondition (false);
        return Vertical_segment();
      }
    };  //Construct_opposite_2_visitor

  }; // Construct_opposite_2

  Construct_opposite_2 construct_opposite_2_object() const
  {
      return Construct_opposite_2();
  }
public:
  template<class OutputIterator>
  static void re_cast_object_vector(const Object_vector& vec,OutputIterator& oi) 
  {
    for (Object_vector::const_iterator it = vec.begin() ; it != vec.end(); ++it)
    {
      X_monotone_curve_2 curr;
      CGAL::Object obj = *it;

      Point_2 p;
      std::pair<Point_2,Multiplicity> ip;
      Non_vertical_x_curve_2 xc;
      Vertical_segment xv;

      if (assign(p, obj))
      {
        curr = p;
        *oi = make_object (curr);
      }
      else if (assign (xc,obj))
      {
        curr = xc;
        *oi = make_object (curr);
      }
      else if (assign (xv,obj))
      {
        curr = xv;
        *oi = make_object (curr);
      }
      else if (assign(ip, obj))
      {
        *oi = make_object (ip);
      }
      else
        CGAL_precondition(false);
      
      ++oi;
    }
    return;
  }
};

} //namespace CGAL {


#endif  //CGAL_ARR_RATIONAL_ARC_TRAITS_D_1_H
