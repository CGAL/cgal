#ifndef ARR_LINEAR_H
#define ARR_LINEAR_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Number_type = Kernel::FT;

using Traits = CGAL::Arr_linear_traits_2<Kernel>;
using Point = Traits::Point_2;
using Segment = Traits::Segment_2;
using Ray = Traits::Ray_2;
using Line = Traits::Line_2;
using X_monotone_curve = Traits::X_monotone_curve_2;

using Arrangement = CGAL::Arrangement_2<Traits>;
using Vertex_handle = Arrangement::Vertex_handle;
using Halfedge_handle = Arrangement::Halfedge_handle;
using Face_handle = Arrangement::Face_handle;
using Vertex_const_handle = Arrangement::Vertex_const_handle;
using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
using Face_const_handle = Arrangement::Face_const_handle;

#endif
