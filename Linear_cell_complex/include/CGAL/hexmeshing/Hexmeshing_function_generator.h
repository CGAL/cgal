#ifndef HEXMESHING_FUNCTION_GENERATOR_H
#define HEXMESHING_FUNCTION_GENERATOR_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_AABBTree_tools.h>
#include <CGAL/hexmeshing/Hexmeshing_set_attributes.h>
#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>
#include <CGAL/hexmeshing/Hexmeshing_function_alias.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <functional>


namespace CGAL::internal::Hexmeshing {

  auto is_volume_intersecting_poly(Tree& tree) {
    return [&](LCC& lcc, Dart_handle dart){
      return is_intersect(lcc, dart, tree);
    };
  }

  // function to be used in move_points_onto_mesh
  auto detect_intersection_with_volume_fraction(double s, size_type inner_mark, size_type set_gradient_mark) {
    return [s, inner_mark, set_gradient_mark](LCC& lcc, Dart_handle dart) -> bool {
      auto dart3 = lcc.beta<3>(dart);
      if(lcc.is_free<3>(dart)) {
        return false;
        if(!lcc.is_marked(dart, inner_mark)) return false;
      }
      else{
        if(lcc.is_marked(dart, inner_mark) == lcc.is_marked(dart3, inner_mark)) return false;
      }
      double frac1 = lcc.attribute<3>(dart)->info().fraction, frac2;
      if(lcc.is_free<3>(dart)) frac2 = 0;
      else frac2 = lcc.attribute<3>(dart3)->info().fraction;

      auto &face_attr = lcc.attribute<2>(dart)->info();
      Point p1, p2;
      if(lcc.is_free<3>(dart)) {
        p1 = face_attr.dual_edge.source();
        p2 = face_attr.dual_edge.target();
      }
      else {
        p1 = lcc.attribute<3>(dart)->info().centroid;
        p2 = lcc.attribute<3>(dart3)->info().centroid;
      }

      face_attr.intersection = p1 + ((s - frac1)/(frac2 - frac1))*(p2 - p1);
      
      if(!lcc.is_marked(dart, set_gradient_mark)) {
        __set_gradient_at_dual_node(lcc, dart);
        lcc.mark_cell<3>(dart, set_gradient_mark);
      }
      if(!lcc.is_free<3>(dart) and !lcc.is_marked(dart3, set_gradient_mark)) {
        __set_gradient_at_dual_node(lcc, dart3);
        lcc.mark_cell<3>(dart3, set_gradient_mark);
      }
      
      Vector N_1 = lcc.attribute<3>(dart)->info().gradient, N_2;
      if(lcc.is_free<3>(dart)) {
        double frac = lcc.attribute<3>(dart)->info().fraction;
        Vector v = p2 - p1;
        N_2 = -(frac/(v*v)) * v;
      }
      else {
        N_2 = lcc.attribute<3>(dart3)->info().gradient;
      }
      Vector N = N_1 + ((s - frac1)/(frac2 - frac1))*(N_2 - N_1);
      face_attr.normal = (1. / std::sqrt(N*N)) * N;

      return true;
    };
  }

  auto is_marked_volume(size_type mark) {
    return [mark](LCC& lcc, Dart_handle dart) -> bool {
      return lcc.is_marked(dart, mark);
    };
  }

  auto detect_intersection(Tree& tree, Polyhedron& poly) {
    return [&](LCC& lcc, Dart_handle dart) -> bool {
      auto &face_attr = lcc.attribute<2>(dart)->info();
      if(!tree.do_intersect(face_attr.dual_edge)) return false;
      auto intersection = tree.any_intersection(face_attr.dual_edge);

      face_attr.intersection = std::get<Point>((*intersection).first);

      const Primitive_id fd = (*intersection).second;
      face_attr.normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, poly);

      return true;
    };
  }

  auto is_inner_centroid(Tree& tree) {
    return [&](LCC& lcc, Dart_handle dart) -> bool {
      // if(is_intersect(lcc, dart, tree)) return false;
      // return !is_outside_knowing_no_intersect(lcc.point(dart), tree);
      __set_centroid(lcc, dart);
      return !is_outside_knowing_no_intersect(lcc.attribute<3>(dart)->info().centroid, tree);
    };
  }

  auto is_inner_point(Tree& tree) {
    return [&](Point p) -> bool {
      return !is_outside_knowing_no_intersect(p, tree);
    };
  }
}




#endif