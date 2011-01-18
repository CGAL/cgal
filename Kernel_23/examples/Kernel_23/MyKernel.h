#ifndef MYKERNEL_H
#define MYKERNEL_H

#include <CGAL/Cartesian.h>
#include "MyPointC2.h"
#include "MySegmentC2.h"
#include "MyConstruct_bbox_2.h"
#include "MyConstruct_coord_iterator.h"
#include "MyConstruct_point_2.h"

// K_ is the new kernel, and K_Base is the old kernel
template < typename K_, typename K_Base >
class MyCartesian_base
  : public K_Base::template Base<K_>::Type
{
  typedef typename K_Base::template Base<K_>::Type   OldK;
public:
  typedef K_                                Kernel;
  typedef MyPointC2                         Point_2;
  typedef MySegmentC2<Kernel>               Segment_2;
  typedef MyConstruct_point_2<Kernel, OldK>       Construct_point_2;
  typedef const double*                     Cartesian_const_iterator_2;
  typedef MyConstruct_coord_iterator        Construct_cartesian_const_iterator_2;
  typedef MyConstruct_bbox_2<typename OldK::Construct_bbox_2>
                                            Construct_bbox_2;

  Construct_point_2
  construct_point_2_object() const
  { return Construct_point_2(); }

  Construct_bbox_2
  construct_bbox_2_object() const
  { return Construct_bbox_2(); }

  Construct_cartesian_const_iterator_2
  construct_cartesian_const_iterator_2_object() const
  { return Construct_cartesian_const_iterator_2(); }

  template < typename Kernel2 >
  struct Base { typedef MyCartesian_base<Kernel2, K_Base>  Type; };
};


template < typename FT_ >
struct MyKernel
  : public CGAL::Type_equality_wrapper<
                MyCartesian_base<MyKernel<FT_>, CGAL::Cartesian<FT_> >,
                MyKernel<FT_> >
{};

#endif // MYKERNEL_H
