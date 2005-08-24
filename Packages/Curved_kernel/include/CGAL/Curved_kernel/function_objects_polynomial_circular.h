// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
//                     National & Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// Authors : Monique Teillaud     <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion         <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
// 
// Partially supported by INRIA's project "CALAMATA", a
// bilateral collaboration with National Kapodistrian University of
// Athens.
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file: include/CGAL/Curved_kernel/function_objects_polynomial.h

#ifndef CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_CIRCULAR_H
#define CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_CIRCULAR_H

#include <CGAL/kernel_basic.h>
#include <CGAL/Curved_kernel/internal_functions_on_circular_arc_2.h>
#include <CGAL/Curved_kernel/internal_functions_on_line_arc_2.h>
#include <CGAL/Bbox_2.h>


namespace CGAL {
namespace CircularFunctors {


  template < class CK >
  class Compare_x_2
    : public CK::Linear_kernel::Compare_x_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Point_2 Point_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_x<CK>(p0, p1); }

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Point_2 &p1) const
    { return compare_x<CK>(p0, p1); }

    result_type
    operator() (const Point_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_x<CK>(p0, p1); }

  };


  template < class CK >
  class Compare_y_2
    : public CK::Linear_kernel::Compare_y_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Point_2 Point_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_y<CK>(p0, p1); }

    result_type
    operator() (const Point_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_y<CK>(p0, p1); }

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Point_2 &p1) const
    { return compare_y<CK>(p0, p1); }

  };

  template < class CK >
  class Compare_xy_2
    : public CK::Linear_kernel::Compare_xy_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Point_2 Point_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_xy<CK>(p0, p1); }

    result_type
    operator() (const Point_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return compare_xy<CK>(p0, p1); }

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Point_2 &p1) const
    { return compare_xy<CK>(p0, p1); }

  };

  template < class CK >
  class In_range_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Line_arc_2              Line_arc_2;

  public:
    typedef bool result_type;

    result_type
    operator()(const Circular_arc_2 &a, const Circular_arc_endpoint_2 &p) const
    { return point_in_range<CK>(a, p); }

    result_type
    operator()(const Line_arc_2 &a, const Circular_arc_endpoint_2 &p) const
    { return point_in_range<CK>(a, p); }
    
  };

  template < class CK >
  class Compare_y_to_right_2
  {
    typedef typename CK::Circular_arc_2           Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2  Circular_arc_endpoint_2;
    typedef typename CK::Line_arc_2               Line_arc_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator()(const Circular_arc_2 &a1,
               const Circular_arc_2 &a2,
               const Circular_arc_endpoint_2 &p) const
    { return compare_y_to_right<CK>(a1, a2, p); }

    result_type
    operator()(const Line_arc_2 &a1,
               const Line_arc_2 &a2,
               const Circular_arc_endpoint_2 &p) const
    { return compare_y_to_right<CK>(a1, a2, p); }

    result_type
    operator()(const Line_arc_2 &a1,
               const Circular_arc_2 &a2,
               const Circular_arc_endpoint_2 &p) const
    { return compare_y_to_right<CK>(a1, a2, p); }

    result_type
    operator()(const Circular_arc_2 &a1,
               const Line_arc_2 &a2,
               const Circular_arc_endpoint_2 &p) const
    { if (compare_y_to_right<CK>(a2, a1, p) == CGAL::LARGER)
	return CGAL::SMALLER;
      return CGAL::LARGER;
    }
  
  };
 
  template < class CK >
  class Equal_2
    : public CK::Linear_kernel::Equal_2
  {
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Line_arc_2              Line_arc_2;
  public:
    typedef bool result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p0,
                const Circular_arc_endpoint_2 &p1) const
    { return equal<CK>(p0, p1); }
    
    result_type
    operator() (const Circular_arc_2 &a0, const Circular_arc_2 &a1) const
    { return equal<CK>(a0, a1); }

    result_type
    operator() (const Line_arc_2 &a0, const Line_arc_2 &a1) const
    { return equal<CK>(a0, a1); }

    result_type
      operator() ( const Line_arc_2 &a0, const Circular_arc_2 &a1) const
    {return false;}

    result_type
      operator() ( const Circular_arc_2 &a0, const Line_arc_2 &a1) const
    {return false;}
    
   };

  template < class CK >
  class Compare_y_at_x_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Line_arc_2              Line_arc_2;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator() (const Circular_arc_endpoint_2 &p,
                const Circular_arc_2 &A1) const
    { return compare_y_at_x<CK>(p, A1); }

    result_type
    operator() (const Circular_arc_endpoint_2 &p,
                const Line_arc_2 &A1) const
    { return compare_y_at_x<CK>(p, A1); }

  };

  template < class CK >
  class Do_overlap_2
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;
    typedef typename CK::Line_arc_2     Line_arc_2;

  public:
    typedef bool result_type;

    result_type
    operator() (const Circular_arc_2 &A1, const Circular_arc_2 &A2) const
    { return do_overlap<CK>(A1, A2); }
    
    result_type
    operator() (const Line_arc_2 &A1, const Line_arc_2 &A2) const
    { return do_overlap<CK>(A1, A2); }

    result_type
      operator() (const Line_arc_2 &A1, const Circular_arc_2 &A2) const
    { return false; }

     result_type
    operator() (const Circular_arc_2 &A1, const Line_arc_2 &A2) const
    { return false; }

  };

  template < class CK >
  class Make_x_monotone_2
  {
    typedef typename CK::Circular_arc_2 Circular_arc_2;
    typedef typename CK::Line_arc_2 Line_arc_2;

  public:
    // typedef OutputIterator result_type;

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 &A, OutputIterator res)
      { return make_x_monotone<CK> (A, res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 &A, OutputIterator res)
    { *res++ = make_object(A);
      return res;
    }

  };
  
  template < class CK >
  class Construct_intersections_2
  {
    public:

    typedef typename CK::Circle_2                 Circle;
    typedef typename CK::Circular_arc_2           Circular_arc;
    typedef typename CK::Line_arc_2               Line_arc;

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Circle & c2, OutputIterator res) const
      { return construct_intersections_2<CK> (c1,c2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc & c1, const Circular_arc & c2, 
	       OutputIterator res) const
      { return construct_intersections_2<CK> (c1,c2,res); }  

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc & c1, const Line_arc & c2, 
	       OutputIterator res) const
      {	return construct_intersections_2<CK> (c1,c2,res); }  

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc & c1, const Circle & c2, 
	       OutputIterator res) const
    { return construct_intersections_2<CK> (c1,c2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Line_arc & c2, 
	       OutputIterator res) const
    { return construct_intersections_2<CK> (c2,c1,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc & c1, const Circular_arc & c2, 
	       OutputIterator res) const
    { return construct_intersections_2<CK> (c1,c2,res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc & c1, const Line_arc & c2, 
	       OutputIterator res) const
    { return construct_intersections_2<CK> (c2,c1,res); }

    
  };

  template < class CK >
  class Nearest_intersection_to_right_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;

  public:
    typedef bool result_type;

    result_type
    operator()(const Circular_arc_2 &A1, const Circular_arc_2 &A2,
               const Circular_arc_endpoint_2 &pt,
               Circular_arc_endpoint_2 &p1, Circular_arc_endpoint_2 &p2) const
    { return nearest_intersection_to_right<CK>(A1, A2, pt, p1, p2); }

  };

  template < class CK >
  class Split_2
  {
    typedef typename CK::Circular_arc_2          Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2 Circular_arc_endpoint_2;
    typedef typename CK::Line_arc_2              Line_arc_2;

  public:
    typedef void result_type;

    result_type
    operator()(const Circular_arc_2 &A, 
	       const Circular_arc_endpoint_2 &p,
	       Circular_arc_2 &ca1, Circular_arc_2 &ca2) const
    { return split<CK>(A, p, ca1, ca2); }


    result_type
    operator()(const Line_arc_2 &A, 
	       const Circular_arc_endpoint_2 &p,
	       Line_arc_2 &ca1, Line_arc_2 &ca2) const
    { return split<CK>(A, p, ca1, ca2); }

  };



  template < class CK >
  class Construct_circular_arc_2
  {

    typedef typename CK::FT                           FT;
    typedef typename CK::RT                           RT;
    typedef typename CK::Linear_kernel::Point_2       Point_2;
    typedef typename CK::Line_2                       Line_2;
    typedef typename CK::Circle_2                     Circle_2;
    typedef typename CK::Circular_arc_2               Circular_arc_2;
    typedef typename CK::Kernel_base::Circular_arc_2  RCircular_arc_2;
    typedef typename Circular_arc_2::Rep              Rep;
    typedef typename CK::Circular_arc_endpoint_2      Circular_arc_endpoint_2;

  public:
    typedef  Circular_arc_2 result_type;

    result_type
    operator()(void) const
    { return Rep(); }


    result_type
    operator()(const Circle_2 &c) const
    { return Rep(c); }

    result_type
    operator()(const Circle_2 &support,
               const Line_2 &l1, bool b1,
               const Line_2 &l2, bool b2)
    { return Rep(support,l1,b1,l2,b2); }

    result_type
    operator()(const Circle_2 &c,
               const Circle_2 &c1, bool b_1,
               const Circle_2 &c2, bool b_2)
    { return Rep(c,c1,b_1,c2,b_2); }

    result_type
    operator()(const Circular_arc_2 &A,
               bool b,
               const Circle_2 &ccut, bool b_cut)
    { return Rep(A,b,ccut,b_cut); }

    result_type
    operator()(const Point_2 &begin,
               const Point_2 &middle, 
               const Point_2 &end)
    { return Rep(begin,middle,end); }

    result_type
    operator()(const Circle_2 &support,
               const Circular_arc_endpoint_2 &source, 
               const Circular_arc_endpoint_2 &target)
    { return Rep(support,source,target); }

  };



  template < class CK >
  class Construct_line_arc_2
  {

    typedef typename CK::Linear_kernel::Point_2    Point_2;
    typedef typename CK::Line_2                    Line_2;
    typedef typename CK::Circle_2                  Circle_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;
    typedef typename CK::Segment_2                 Segment_2;
    typedef typename CK::Line_arc_2                Line_arc_2;
    typedef typename CK::Kernel_base::Line_arc_2   RLine_arc_2;
    typedef typename Line_arc_2::Rep               Rep;

  public:
    typedef Line_arc_2  result_type;

    result_type
    operator()(void) const
    { return Rep(); }

    result_type
    operator()(const Line_2 &support,
	       const Circle_2 &c1,const bool b1,
	       const Circle_2 &c2,const bool b2) const
    { return Rep(support,c1,b1,c2,b2); }

    result_type
    operator()(const Line_2 &support,
	       const Line_2 &l1,
	       const Line_2 &l2) const
    { return Rep(support,l1,l2); }

    result_type
    operator()(const Line_2 &support,
	       const Circular_arc_endpoint_2 &p1,
	       const Circular_arc_endpoint_2 &p2) const
    { return Rep(support,p1,p2); }


    result_type
    operator()(const Segment_2 &s) const
    { return Rep(s); }


    result_type
    operator()(const Point_2 &p1,
	       const Point_2 &p2) const
    { return Rep(p1,p2); }


   };



  template < class CK >
  class Construct_circular_arc_endpoint_2
  {

    typedef typename CK::Circular_arc_endpoint_2               Circular_arc_endpoint_2;
    typedef typename CK::Kernel_base::Circular_arc_endpoint_2  RCircular_arc_endpoint_2;
    typedef typename Circular_arc_endpoint_2::Rep              Rep;
    typedef typename Circular_arc_endpoint_2::Root_for_circles_2_2  Root_for_circles_2_2;

  public:
    typedef  Circular_arc_endpoint_2 result_type;

    result_type
    operator()(void) const
    { return Rep(); }

    result_type
    operator()(const Root_for_circles_2_2 & np) const
    { return Rep(np); }

  };


  template <class CK>
  class Compute_x_2: Has_qrt
  {
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;
    typedef typename CK::Root_of_2                 Root_of_2;

  public:

    typedef  Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;

    qualified_result_type operator() (const Circular_arc_endpoint_2 & a) const
    {
      return (a.rep().x());
    }
  };


  template <class CK>
  class Compute_y_2: Has_qrt
  {
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;
    typedef typename CK::Root_of_2                 Root_of_2;

  public:

    typedef Root_of_2  result_type;
    typedef const result_type &        qualified_result_type;
    
    qualified_result_type operator() (const Circular_arc_endpoint_2 & a) const
    {
      return (a.rep().y());
    }
  };


  template <class CK>
  class Construct_min_vertex_2: Has_qrt
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Line_arc_2                Line_arc_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;

  public:

    typedef Circular_arc_endpoint_2 result_type;
    typedef const result_type &        qualified_result_type;

    qualified_result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().left());
    }

    qualified_result_type operator() (const Line_arc_2 & a) const
    {
      return (a.rep().left());
    }


  };

  template <class CK>
  class Construct_max_vertex_2: Has_qrt
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Line_arc_2                Line_arc_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;

  public:

    typedef Circular_arc_endpoint_2 result_type;
    typedef const result_type &        qualified_result_type;

    qualified_result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().right());
    }

    qualified_result_type operator() (const Line_arc_2 & a) const
    {
      return (a.rep().right());
    }

  };

  template <class CK>
  class Construct_source_vertex_2: Has_qrt
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Line_arc_2                Line_arc_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;

  public:

    typedef Circular_arc_endpoint_2 result_type;
    typedef const result_type &        qualified_result_type;
    
    qualified_result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().source());
    }

    qualified_result_type operator() (const Line_arc_2 & a) const
    {
      return (a.rep().source());
    }

  };


  template <class CK>
  class Construct_target_vertex_2: Has_qrt
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Line_arc_2                Line_arc_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;

  public:

    typedef Circular_arc_endpoint_2  result_type;
    typedef const result_type &        qualified_result_type;

    qualified_result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().target());
    }

    qualified_result_type operator() (const Line_arc_2 & a) const
    {
      return (a.rep().target());
    }

  };

  template <class CK>
  class Is_x_monotone_2
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;

  public:

    typedef bool result_type;

    result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().is_x_monotone());
    }
  };

  template <class CK>
  class Is_y_monotone_2
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;

  public:

    typedef bool result_type;

    result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().is_y_monotone());
    }
  };


  template <class CK>
  class Construct_supporting_circle_2: Has_qrt
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Circle_2                  Circle_2;

  public:

    typedef  Circle_2 result_type;
    typedef const result_type &        qualified_result_type;

    qualified_result_type operator() (const Circular_arc_2 & a) const
    {
      return (a.rep().supporting_circle());
    }
  };

  template <class CK>
  class Construct_supporting_line_2: Has_qrt
  {
    typedef typename CK::Line_arc_2            Line_arc_2;
    typedef typename CK::Line_2                Line_2;
    typedef typename CK::Circle_2              Circle_2;

  public:

    typedef  Line_2                    result_type;
    typedef const result_type &        qualified_result_type;

    qualified_result_type operator() (const Line_arc_2 & a) const
    {
      return (a.rep().supporting_line());
    }
  };

  template <class CK>
  class Construct_bbox_2
  {
    typedef typename CK::Circular_arc_2            Circular_arc_2;
    typedef typename CK::Circular_arc_endpoint_2   Circular_arc_endpoint_2;
    typedef typename CK::Line_arc_2                Line_arc_2;
    typedef typename CK::Circle_2                  Circle_2;

  public:

    typedef CGAL::Bbox_2   result_type;

    result_type operator() (const Circular_arc_endpoint_2 & a) const
    {
      return a.rep().bbox();
    }

    result_type operator() (const Circular_arc_2 & a) const
    {
      return a.rep().bbox();
    }

    result_type operator() (const Line_arc_2 & a) const
    {
      return a.rep().bbox();
    }


  };







} // namespace CircularFunctors
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_FUNCTION_OBJECTS_POLYNOMIAL_CIRCULAR_H
