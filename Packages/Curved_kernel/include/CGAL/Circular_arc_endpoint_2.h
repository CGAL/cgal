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

// file : include/CGAL/Circular_arc_endpoint_2.h

#ifndef CGAL_CIRCULAR_ARC_ENDPOINT_2_H
#define CGAL_CIRCULAR_ARC_ENDPOINT_2_H
namespace CGAL {

template < typename CurvedKernel >
class Circular_arc_endpoint_2
  : public CurvedKernel::Kernel_base::Circular_arc_endpoint_2
{
  typedef typename CurvedKernel::Kernel_base::Circular_arc_endpoint_2 
                                           RCircular_arc_endpoint_2;
  typedef typename CurvedKernel::Circle_2                  Circle_2;

  typedef typename CurvedKernel::Root_of_2               Root_of_2;

public:
 typedef typename CGAL::Simple_cartesian<Root_of_2>::Point_2
                                                 Numeric_point_2;
  typedef CurvedKernel   R; 
  
  Circular_arc_endpoint_2() {}

  Circular_arc_endpoint_2(const Circle_2 & c, const Circle_2 & c1, 
			  const bool b)
    : RCircular_arc_endpoint_2(c, c1, b) {}

  Circular_arc_endpoint_2(const Circle_2 & c, const Circle_2 & c1,
  			  const bool b, const Numeric_point_2 & np)
    : RCircular_arc_endpoint_2(c, c1, b, np){}
};

template < typename CurvedKernel >
inline
bool
operator==(const Circular_arc_endpoint_2<CurvedKernel> &p,
           const Circular_arc_endpoint_2<CurvedKernel> &q)
{
  return CurvedKernel().equal_2_object()(p, q);
}

template < typename CurvedKernel >
inline
bool
operator!=(const Circular_arc_endpoint_2<CurvedKernel> &p,
           const Circular_arc_endpoint_2<CurvedKernel> &q)
{
  return ! (p == q);
}

} // namespace CGAL

#endif // CGAL_CIRCULAR_ARC_ENDPOINT_2_H
