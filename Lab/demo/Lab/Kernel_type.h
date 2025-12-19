#ifndef KERNEL_TYPE_H
#define KERNEL_TYPE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

namespace kernel_type_h {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K1;
}

typedef CGAL::Mesh_3::Robust_intersection_traits_3<kernel_type_h::K1> Kernel;

#endif // KERNEL_TYPE_H
