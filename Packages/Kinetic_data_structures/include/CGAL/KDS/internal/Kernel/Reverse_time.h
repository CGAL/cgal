#ifndef CGAL_KDS_KERNEL_REVERSE_TIME_H
#define CGAL_KDS_KERNEL_REVERSE_TIME_H
#include <CGAL/KDS/basic.h>


CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class K>
class Reverse_time {
public:
  Reverse_time(const typename K::Polynomial_kernel::Negate_variable &nv): nv_(nv){}

  typedef typename K::Point_3 argument_type;
  typedef typename K::Point_3 result_type;
  
  template <class O>
  O operator()(const O &i) const {
    return i.transformed_coordinates(nv_);
  }

protected:
  typename K::Polynomial_kernel::Negate_variable nv_;
};

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
