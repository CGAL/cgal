
#include "Aos.h"

#include <iostream>
#include <iterator>
#include <vector>

#include <qmath.h>
#include <qvector3d.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Vector_3.h>

#include "arr_print.h"

namespace {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
  using Point = Geom_traits::Point_2;
  using Curve = Geom_traits::Curve_2;
  using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits>;
  using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;

  // Extended DCEL & Arrangement
  struct Flag
  {
    bool v;
    Flag() : v{ false } {}
    Flag(bool init) : v{ init } {}
  };

  using Ext_dcel = CGAL::Arr_extended_dcel<Geom_traits, Flag, Flag, Flag>;
  using Ext_topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits, 
                                                                      Ext_dcel>;
  using Ext_aos = CGAL::Arrangement_on_surface_2<Geom_traits, Ext_topol_traits>;



  using Dir3 = Kernel::Direction_3;
  std::ostream& operator << (std::ostream& os, const Dir3& d)
  {
    os << d.dx() << ", " << d.dy() << ", " << d.dz();
    return os;
  }

  using Approximate_point_2 = Geom_traits::Approximate_point_2;
  std::ostream& operator << (std::ostream& os, const Approximate_point_2& d)
  {
    os << d.dx() << ", " << d.dy() << ", " << d.dz();
    return os;
  }

  using Approximate_number_type = Geom_traits::Approximate_number_type;
  using Approximate_kernel = Geom_traits::Approximate_kernel;
  using Approximate_Vector_3 = CGAL::Vector_3<Approximate_kernel>;
  using Approximate_Direction_3 = Approximate_kernel::Direction_3;
  using Direction_3 = Kernel::Direction_3;


  std::ostream& operator << (std::ostream& os, const Approximate_Vector_3& v)
  {
    os << v.x() << ", " << v.y() << ", " << v.z();
    //os << v.hx() << ", " << v.hy() << ", " << v.hz() << ", " << v.hw();
    return os;
  }
}


void Aos::check(const Kml::Placemarks& placemarks)
{
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Arrangement arr(&traits);

  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  int num_counted_nodes = 0;
  int num_counted_arcs = 0;
  int num_counted_polygons = 0;
  std::vector<Curve>  xcvs;
  for (const auto& pm : placemarks)
  {
    for (const auto& lring : pm.polygons)
    {
      num_counted_polygons++;

      // convert the nodes to points on unit-sphere
      std::vector<Approximate_Vector_3> sphere_points;
      for (const auto& node : lring.nodes)
      {
        num_counted_nodes++;
        const auto p = node.get_coords_3d();
        Approximate_Vector_3  v(p.x, p.y, p.z);
        sphere_points.push_back(v);
        CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
      }

      // add curves
      int num_points = sphere_points.size();
      for (int i = 0; i < sphere_points.size()-1; i++)
      {
        num_counted_arcs++;
        const auto p1 = sphere_points[i];
        const auto p2 = sphere_points[i + 1];
        auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
                          ctr_p(p2.x(), p2.y(), p2.z()));
        xcvs.push_back(xcv);
      }
    }
  }

  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " << 
                                          arr.number_of_vertices() << std::endl;

  // add arcs
  for(auto xcv : xcvs)
    CGAL::insert_curve(arr, xcv);

  std::cout << "-------------------------------\n";
  std::cout << "num nodes = " << num_counted_nodes << std::endl;
  std::cout << "num arr vertices = " << arr.number_of_vertices() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num counted arcs = " << num_counted_arcs << std::endl;
  std::cout << "num arr edges = " << arr.number_of_edges() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num polygons = " << num_counted_polygons << std::endl;
  std::cout << "num arr faces = " << arr.number_of_faces() << std::endl;
}

void Aos::ext_check(const Kml::Placemarks& placemarks)
{
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Ext_aos arr(&traits);

  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  int num_counted_nodes = 0;
  int num_counted_arcs = 0;
  int num_counted_polygons = 0;
  std::vector<Curve>  xcvs;
  for (const auto& pm : placemarks)
  {
    for (const auto& lring : pm.polygons)
    {
      num_counted_polygons++;

      // convert the nodes to points on unit-sphere
      std::vector<Approximate_Vector_3> sphere_points;
      for (const auto& node : lring.nodes)
      {
        num_counted_nodes++;
        const auto p = node.get_coords_3d();
        Approximate_Vector_3  v(p.x, p.y, p.z);
        sphere_points.push_back(v);
        CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
      }

      // add curves
      int num_points = sphere_points.size();
      for (int i = 0; i < sphere_points.size() - 1; i++)
      {
        num_counted_arcs++;
        const auto p1 = sphere_points[i];
        const auto p2 = sphere_points[i + 1];
        auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
          ctr_p(p2.x(), p2.y(), p2.z()));
        xcvs.push_back(xcv);
      }
    }
  }

  // MARK all vertices as true
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    vit->set_data(Flag(true));
  }

  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;

  // add arcs
  for (auto xcv : xcvs)
    CGAL::insert_curve(arr, xcv);

  // extract all vertices that are ADDED when inserting the arcs!
  int num_created_vertices = 0;
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    auto& d = vit->data();
    if (vit->data().v == false)
    {
      num_created_vertices++;
      auto p = vit->point();
      const auto x = CGAL::to_double(p.dx());
      const auto y = CGAL::to_double(p.dy());
      const auto z = CGAL::to_double(p.dz());
    }
  }
  std::cout << "*** num created vertices = " << num_created_vertices << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num nodes = " << num_counted_nodes << std::endl;
  std::cout << "num arr vertices = " << arr.number_of_vertices() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num counted arcs = " << num_counted_arcs << std::endl;
  std::cout << "num arr edges = " << arr.number_of_edges() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num polygons = " << num_counted_polygons << std::endl;
  std::cout << "num arr faces = " << arr.number_of_faces() << std::endl;
}
