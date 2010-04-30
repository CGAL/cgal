// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_ARC_2_H
#define CGAL_CIRCULAR_ARC_2_H

CGAL_BEGIN_NAMESPACE
  
template <class CircularKernel> 
class Circular_arc_2 
  : public CircularKernel::Kernel_base::Circular_arc_2
{
  typedef typename CircularKernel::RT             RT;
  typedef typename CircularKernel::FT             FT;
  typedef typename CircularKernel::Point_2        Point_2;
  typedef typename CircularKernel::Line_2         Line_2;
  typedef typename CircularKernel::Circle_2       Circle_2;
  typedef typename CircularKernel::Circular_arc_point_2
                                                Circular_arc_point_2;
  
  typedef typename CircularKernel::Kernel_base::Circular_arc_2 RCircular_arc_2; 
  // RCircular_arc_2 to avoid clash with self 
public:
  typedef  RCircular_arc_2 Rep;
  typedef  CircularKernel   R; 
  

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }


  Circular_arc_2()
    : RCircular_arc_2(typename R::Construct_circular_arc_2()())
  {}

  Circular_arc_2(const Circle_2 &c)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(c))
  {}

  // Not Documented
  Circular_arc_2(const Circle_2 &support, 
                 const Line_2 &l1, const bool b_l1,
                 const Line_2 &l2, const bool b_l2)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(support,l1,b_l1,l2,b_l2))
  {}

  // Not Documented
  Circular_arc_2(const Circle_2 &c, 
		 const Circle_2 &c1, const bool b_1,
		 const Circle_2 &c2, const bool b_2)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(c,c1,b_1,c2,b_2))
  {}

  Circular_arc_2(const Point_2 &start,
                 const Point_2 &middle,
                 const Point_2 &end)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(start, middle, end)) 
  {}
  
  Circular_arc_2(const Circle_2 &support,
                 const Circular_arc_point_2 &begin,
                 const Circular_arc_point_2 &end)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(support, begin, end)) 
  {}

  Circular_arc_2(const Point_2 &start,
                 const Point_2 &end,
		 const FT &bulge)
    : RCircular_arc_2(typename 
		      R::Construct_circular_arc_2()(start, end, bulge)) 
  {}
  
 Circular_arc_2(const RCircular_arc_2 & a)
    : RCircular_arc_2(a)
  {}


  typename Qualified_result_of    
  <typename R::Construct_circular_source_vertex_2,Circular_arc_2>::type
  //const Circular_arc_point_2 &
  source() const
  {
    return typename R::Construct_circular_source_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_circular_target_vertex_2,Circular_arc_2>::type
  //const Circular_arc_point_2 &
  target() const
  {
    return typename R::Construct_circular_target_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_circular_min_vertex_2,Circular_arc_2>::type
  //const Circular_arc_point_2 & 
  left() const
  {
    return typename R::Construct_circular_min_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_circular_max_vertex_2,Circular_arc_2>::type
  //const Circular_arc_point_2 & 
  right() const
  {
    return typename R::Construct_circular_max_vertex_2()(*this);
  }

  bool is_x_monotone() const
  {
    return typename R::Is_x_monotone_2()(*this);
  }

  bool is_y_monotone() const
  {
    return typename R::Is_y_monotone_2()(*this);
  }

	typename Qualified_result_of
  <typename R::Construct_circle_2,Circular_arc_2>::type
  // const Circle_2 & 
  supporting_circle() const
  {
    return typename R::Construct_circle_2()(*this);
  }

	typename Qualified_result_of
  <typename R::Construct_center_2,Circular_arc_2>::type
  // const Point_2 & 
  center() const
  {
    return typename R::Construct_center_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Compute_squared_radius_2, Circular_arc_2>::type
  // const FT & 
  squared_radius() const
  {
    return typename R::Compute_squared_radius_2()(*this);
  }

  Bbox_2 bbox(void) const
  {
    return typename R::Construct_bbox_2()(*this);
  }

};

  template < typename CircularKernel >
  inline
  bool
  operator==(const Circular_arc_2<CircularKernel> &p,
	     const Circular_arc_2<CircularKernel> &q)
  {
    return CircularKernel().equal_2_object()(p, q);
  }
  
  template < typename CircularKernel >
  inline
  bool
  operator!=(const Circular_arc_2<CircularKernel> &p,
	     const Circular_arc_2<CircularKernel> &q)
  {
    return ! (p == q);
  }

  template < typename CK >
  std::ostream &
  operator<<(std::ostream & os, const Circular_arc_2<CK> &a)
  {
    // The output format is :
    // - supporting circle
    // - circle c1
    // - bool b1
    // - circle c2
    // - bool b2
    return os << a.supporting_circle() << " "
	      << a.source() << " "
	      << a.target() << " ";
  }
  
  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_2<CK> &a)
  {
    typename CK::Circle_2 s;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Circular_arc_2<CK>(s, p1, p2);
    return is;
  }

template < class CK >
struct Filtered_bbox_circular_kernel_2;

template < typename CK >
class Circular_arc_2 < Filtered_bbox_circular_kernel_2 < CK > > {

	  typedef Filtered_bbox_circular_kernel_2 < CK >         BK;
    typedef Circular_arc_2< BK >                           Self;
    typedef typename BK::FT                                FT;
    typedef typename BK::RT                                RT;
    typedef typename BK::Point_2                           Point_2;
    typedef typename BK::Line_2                            Line_2;
    typedef typename BK::Circle_2                          Circle_2;
    typedef typename BK::Circular_arc_point_2              Circular_arc_point_2;
    typedef typename CK::Circular_arc_2                    Rcircular_arc_2;
    typedef typename CK::Root_of_2                         Root_of_2;

public:
    typedef BK                       R; 
    typedef Circular_arc_2<BK>       Rep;

    const Rep& rep() const
    {
      return *this;
    }

    Rep& rep()
    {
      return *this;
    }

     ///////////Construction/////////////

    Circular_arc_2(){}

    // otherwise it will lead to ambiguos definitions
    explicit Circular_arc_2(const Circle_2 &c)
    : P_arc(c),bb(NULL)
    {}

    Circular_arc_2(const Circle_2 &support, 
                 	 const Line_2 &l1, const bool b_l1,
                   const Line_2 &l2, const bool b_l2)
    : P_arc(support,l1,b_l1,l2,b_l2),bb(NULL)
    {}

    
    Circular_arc_2(const Circle_2 &c, 
	   		   const Circle_2 &c1, const bool b_1,
	   		   const Circle_2 &c2, const bool b_2)
    : P_arc(c,c1,b_1,c2,b_2),bb(NULL)
    {}

    
    Circular_arc_2(const Rcircular_arc_2 &A, const bool b,
		   const Circle_2 &ccut, const bool b_cut)
    : P_arc(A, b, ccut, b_cut),bb(NULL)
    {}


    Circular_arc_2(const Point_2 &start,
    	   const Point_2 &middle,
    	   const Point_2 &end)
    : P_arc(start, middle, end),bb(NULL)
    {}

    Circular_arc_2(const Point_2 &begin,
                             const Point_2 &end,
	                         const FT &bulge) 
    : P_arc(begin, end, bulge),bb(NULL)
    {}

    Circular_arc_2(const Circle_2 &support,
    	   const Circular_arc_point_2 &begin,
    	   const Circular_arc_point_2 &end)
    : P_arc(support, begin.point(), end.point()),bb(NULL) 
	  {}

    Circular_arc_2(const Rcircular_arc_2 &a)
    : P_arc(a),bb(NULL) 
	  {}

	  Circular_arc_2(const Circular_arc_2 &c) : P_arc(c.P_arc) 
	  {
	    if(c.bb) bb = new Bbox_2(*(c.bb));
		  else bb = NULL;	
	  }

	  ~Circular_arc_2() { if(bb) delete bb; }


		//////////Predicates//////////

		bool is_x_monotone() const
		{ return P_arc.is_x_monotone();}

		bool is_y_monotone() const
		{ return P_arc.is_y_monotone();}

		bool on_upper_part() const
		{ return P_arc.on_upper_part();}
		
		
		//////////Accessors///////////

		const Rcircular_arc_2& arc () const
			{ return P_arc ;}
  
		///Interface of the inner arc/// 

		typename Qualified_result_of<typename BK::Construct_circular_source_vertex_2,Self>::type
                source() const
			{ return typename BK::Construct_circular_source_vertex_2()(*this);}

		typename Qualified_result_of<typename BK::Construct_circular_target_vertex_2,Self>::type
                target() const
			{ return typename BK::Construct_circular_target_vertex_2()(*this);}

    typename Qualified_result_of<typename BK::Construct_circular_min_vertex_2,Self>::type
                left() const
      {
			  return typename BK::Construct_circular_min_vertex_2()(*this);
			}
	      
     typename Qualified_result_of<typename BK::Construct_circular_max_vertex_2,Self>::type
                right() const
      {
			  return typename BK::Construct_circular_max_vertex_2()(*this);
			}

		Circle_2 supporting_circle() const
			{ return P_arc.supporting_circle();}

		Point_2 center() const
			{ return P_arc.center();}

		FT squared_radius() const
			{ return P_arc.squared_radius();}
		
		Bbox_2 bbox() const
			{ 
        if(bb==NULL)
          bb=new Bbox_2(P_arc.bbox());
        return *bb;
      }
                          
			
		///Specific check used for bbox construction///
		
		bool has_no_bbox() const
		{ return (bb==NULL);}
		
		bool equal_ref(const Circular_arc_2 &c) const
    {
      return CGAL::identical(P_arc, c.P_arc);      
    }

    bool is_full() const {
			return P_arc.is_full();
	  }

    bool is_complementary_x_monotone() const {
      return P_arc.is_complementary_x_monotone();
    } 

    bool is_complementary_y_monotone() const {
	    return P_arc.is_complementary_y_monotone();
    }

    bool two_end_points_on_upper_part() const {
      return P_arc.two_end_points_on_upper_part();
    }

    bool complementary_on_upper_part() const {
      return P_arc.complementary_on_upper_part();
    }

    bool two_end_points_on_left_part() const {
      return P_arc.two_end_points_on_left_part();
    }

    bool on_left_part() const {
      return P_arc.on_left_part();
    }

    bool complementary_on_left_part() const {
      return P_arc.complementary_on_left_part();
    }

	private:

		Rcircular_arc_2 P_arc;
		mutable Bbox_2 *bb;


};



CGAL_END_NAMESPACE

#endif // CGAL_CIRCULAR_ARC_2_H
