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

// file : include/CGAL/Curved_kernel/Debug_id.h

#ifndef CGAL_CURVED_KERNEL_DEBUG_ID_H
#define CGAL_CURVED_KERNEL_DEBUG_ID_H

namespace CGAL {
namespace CGALi {

// Debug_id is a base class storing an identifier, to help recognizing
// the geometric objects in a debugging dump.
// It is empty, unless CGAL_CURVED_KERNEL_DEBUG is defined.

#ifndef CGAL_CURVED_KERNEL_DEBUG
  template <typename = void>
  struct Debug_id {
    const Debug_id & id() const { return *this; }
  };

  template < typename T >
  std::ostream &
  operator<< (std::ostream & os, const Debug_id<T> &)
  {
    return os;
  }
#else
  // I make it a template in order to avoid linking problems
  // with the static data member.
  template <typename = void>
  class Debug_id {
    static int cnt;
  public:
    int _id;
    Debug_id() : _id(cnt++) {}
    const Debug_id & id() const { return *this; }
  };

  template <typename T>
  int Debug_id<T>::cnt = 0;

  template < typename T >
  std::ostream &
  operator<< (std::ostream & os, const Debug_id<T> &a)
  {
    return os << " [[[ DEBUG_ID = " << a._id << " ]]] ";
  }
#endif

} // namespace CGALi
} // namespace CGAL

#endif // CGAL_CURVED_KERNEL_DEBUG_ID_H
