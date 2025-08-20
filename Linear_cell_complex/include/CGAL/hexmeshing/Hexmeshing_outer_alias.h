#ifndef HEXMESHING_OUTER_ALIAS_H
#define HEXMESHING_OUTER_ALIAS_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Cell_attribute_with_point.h>


namespace CGAL::Hexmeshing {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Triangle = Kernel::Triangle_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
using Segment = CGAL::Segment_3<Kernel>;
using AABB_Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<AABB_Traits>;
using Primitive_id = typename Tree::Primitive_id;
using RandomPointGenerator = CGAL::Random_points_in_cube_3<Point>;

using Side_of_mesh = CGAL::Side_of_triangle_mesh<Polyhedron, Kernel>;

using LCCTraits = CGAL::Linear_cell_complex_traits<3,Kernel>;
}



#endif