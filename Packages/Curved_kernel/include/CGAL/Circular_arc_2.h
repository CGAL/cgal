// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Circular_arc_2.h

#ifndef CGAL_CIRCULAR_ARC_2_H
#define CGAL_CIRCULAR_ARC_2_H

namespace CGAL {

template <class CurvedKernel> 
class Circular_arc_2 
  : public CurvedKernel::Kernel_base::Circular_arc_2
{
  typedef typename CurvedKernel::RT             RT;
  typedef typename CurvedKernel::FT             FT;
  typedef typename CurvedKernel::Point_2        Point_2;
  typedef typename CurvedKernel::Line_2         Line_2;
  typedef typename CurvedKernel::Circle_2       Circle_2;
  typedef typename CurvedKernel::Circular_arc_point_2
                                                Circular_arc_point_2;

  typedef typename CurvedKernel::Kernel_base::Circular_arc_2 RCircular_arc_2; 
  // RCircular_arc_2 to avoid clash with self 
public:
  typedef  RCircular_arc_2 Rep;
  typedef  CurvedKernel   R; 
  

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

  Circular_arc_2(const Circle_2 &support, 
                 const Line_2 &l1, const bool b_l1,
                 const Line_2 &l2, const bool b_l2)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(support,l1,b_l1,l2,b_l2))
  {}

  Circular_arc_2(const Circle_2 &c, 
		 const Circle_2 &c1, const bool b_1,
		 const Circle_2 &c2, const bool b_2)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(c,c1,b_1,c2,b_2))
    {}

  Circular_arc_2(const Circular_arc_2 &A, const bool b,
		 const Circle_2 &ccut, const bool b_cut)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(A, b, ccut, b_cut))
    {}

  Circular_arc_2(const Point_2 &start,
                 const Point_2 &middle,
                 const Point_2 &end)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(start, middle, end)) {}
  
  Circular_arc_2(const Circle_2 &support,
                 const Point_2 &begin,
                 const Point_2 &end)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(support, begin, end)) {}
  
  Circular_arc_2(const Circle_2 &support,
                 const Circular_arc_point_2 &begin,
                 const Circular_arc_point_2 &end)
    : RCircular_arc_2(typename R::Construct_circular_arc_2()(support, begin, end)) {}
 
 Circular_arc_2(const RCircular_arc_2 & a)
    : RCircular_arc_2(a)
  {}


  const Circular_arc_point_2 & source() const
  {
	return typename R::Construct_source_vertex_2()(*this);
  }

  const Circular_arc_point_2 & target() const
  {
	return typename R::Construct_target_vertex_2()(*this);
  }

  const Circular_arc_point_2 & left() const
  {
	return typename R::Construct_min_vertex_2()(*this);
  }

  const Circular_arc_point_2 & right() const
  {
	return typename R::Construct_max_vertex_2()(*this);
  }

  bool is_x_monotone() const
  {
	return typename R::Is_x_monotone_2()(*this);
  }

  bool is_y_monotone() const
  {
	return typename R::Is_y_monotone_2()(*this);
  }

  const Circle_2 & supporting_circle() const
  {
	return typename R::Construct_supporting_circle_2()(*this);
  }

  const Point_2 & center() const
  {
	return typename R::Construct_supporting_circle_2()(*this).center();
  }

  const FT & squared_radius() const
  {
	return typename R::Construct_supporting_circle_2()(*this).squared_radius();
  }

 Bbox_2 bbox(void) const
 {
       return typename R::Construct_bbox_2()(*this);
 }

};

} // namespace CGAL

#endif // CGAL_CIRCULAR_ARC_2_H
