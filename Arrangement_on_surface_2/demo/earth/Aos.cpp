
#include "Aos.h"

#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>

#include <qmath.h>
#include <qvector3d.h>

#include <nlohmann/json.hpp>
using json = nlohmann::ordered_json;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Vector_3.h>

#include "arr_print.h"
#include "Tools.h"

namespace {
  //#define USE_EPIC

#ifdef USE_EPIC
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
#else
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
#endif

  using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
  using Point = Geom_traits::Point_2;
  using Curve = Geom_traits::Curve_2;
  using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits>;
  using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;

  // the following is from "arr_inexact_construction_segments.h":
  using Segment = Geom_traits::X_monotone_curve_2;
  using Vertex_handle = Arrangement::Vertex_handle;


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


  //---------------------------------------------------------------------------
  // below are the helper functions used to construct the arcs from KML data
  // TODO: Revisit handling of INNER & OUTER boundaries
  using Curves = std::vector<Curve>;

  // get curves for the given kml placemark
  // NOTE: this is defined here to keep the definitions local to this cpp file
  Curves  get_arcs(const Kml::Placemark& placemark)
  {
    Geom_traits traits;
    auto ctr_p = traits.construct_point_2_object();
    auto ctr_cv = traits.construct_curve_2_object();

    std::vector<Curve>  xcvs;
    for (const auto& polygon : placemark.polygons)
    {
      // colect all rings into a single list (FOR NOW!!!)
      // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
      Kml::LinearRings linear_rings;
      linear_rings.push_back(polygon.outer_boundary);
      for (const auto& inner_boundary : polygon.inner_boundaries)
        linear_rings.push_back(inner_boundary);


      // convert the nodes to points on unit-sphere
      for (const auto& lring : linear_rings)
      {
        std::vector<Approximate_Vector_3>  sphere_points;
        for (const auto& node : lring.nodes)
        {
          const auto p = node.get_coords_3d();
          Approximate_Vector_3 v(p.x, p.y, p.z);
          sphere_points.push_back(v);
        }

        // add geodesic arcs for the current LinearRing
        int num_points = sphere_points.size();
        for (int i = 0; i < num_points - 1; i++)
        {
          const auto p1 = sphere_points[i];
          const auto p2 = sphere_points[i + 1];
          xcvs.push_back(ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
            ctr_p(p2.x(), p2.y(), p2.z())));
        }
      }
    }

    return xcvs;
  }


  // this one is used by the Aos::check and Aos::ext_check functions
  int num_counted_nodes = 0;
  int num_counted_arcs = 0;
  int num_counted_polygons = 0;
  std::map<Ext_aos::Vertex_handle, Kml::Node>  vertex_node_map;

  template<typename Arr_type>
  Curves  get_arcs(const Kml::Placemarks& placemarks, Arr_type& arr)
  {
    Geom_traits traits;
    auto ctr_p = traits.construct_point_2_object();
    auto ctr_cv = traits.construct_curve_2_object();

    num_counted_nodes = 0;
    num_counted_arcs = 0;
    num_counted_polygons = 0;
    std::vector<Curve>  xcvs;
    for (const auto& pm : placemarks)
    {
      for (const auto& polygon : pm.polygons)
      {
        num_counted_polygons++;

        // colect all rings into a single list (FOR NOW!!!)
        // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
        Kml::LinearRings linear_rings;
        linear_rings.push_back(polygon.outer_boundary);
        for (const auto& inner_boundary : polygon.inner_boundaries)
          linear_rings.push_back(inner_boundary);

        // loop on outer and inner boundaries 
        for (const auto& lring : linear_rings)
        {
          // convert the nodes to points on unit-sphere
          std::vector<Approximate_Vector_3>  sphere_points;
          for (const auto& node : lring.nodes)
          {
            num_counted_nodes++;
            const auto p = node.get_coords_3d();
            Approximate_Vector_3  v(p.x, p.y, p.z);
            sphere_points.push_back(v);
            auto vh = CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
            if constexpr (std::is_same<Arr_type, Ext_aos>::value)
            {
              vertex_node_map.insert(std::make_pair(vh, node));
            }
          }

          // add curves
          int num_points = sphere_points.size();
          for (int i = 0; i < num_points - 1; i++)
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
    }
    return xcvs;
  }


  Aos::Approx_arc  get_approx_curve(Curve xcv, double error)
  {
    Geom_traits traits;
    auto approx = traits.approximate_2_object();
    std::vector<QVector3D>  approx_curve;
    {
      std::vector<Approximate_point_2> v;
      auto oi2 = approx(xcv, error, std::back_insert_iterator(v));

      for (const auto& p : v)
      {
        const QVector3D arc_point(p.dx(), p.dy(), p.dz());
        approx_curve.push_back(arc_point);
      }
    }

    return approx_curve;
  }
  Aos::Approx_arcs get_approx_curves(std::vector<Curve>& xcvs, double error)
  {
    Aos::Approx_arcs  approx_curves;
    for (const auto& xcv : xcvs)
    {
      auto approx_curve = get_approx_curve(xcv, error);
      approx_curves.push_back(std::move(approx_curve));
    }

    return approx_curves;
  }
}


Aos::Approx_arc Aos::get_approx_identification_curve(double error)
{
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  // identification curve (meridian pierced by NEGATIVE Y-AXIS)
  auto xcv = ctr_cv(ctr_p(0, 0, -1), ctr_p(0, 0, 1), Dir3(0, 1, 0));

  auto approx = traits.approximate_2_object();
  Approx_arc approx_arc;
  {
    std::vector<Approximate_point_2> v;
    auto oi2 = approx(xcv, error, std::back_insert_iterator(v));
    for (const auto& p : v)
    {
      const QVector3D arc_point(p.dx(), p.dy(), p.dz());
      approx_arc.push_back(arc_point);
    }
  }

  return approx_arc;
}

Aos::Approx_arcs Aos::get_approx_arcs(double error)
{
  Geom_traits traits;
  Arrangement arr(&traits);

  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  std::vector<Curve>  xcvs;
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0)));
  xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 0, 1)));
  xcvs.push_back(ctr_cv(ctr_p(0, 1, 0), ctr_p(0, 0, 1)));
  //xcvs.push_back(ctr_cv(ctr_p(1, 0, 0), ctr_p(0, 1, 0), Dir3(0, 0, -1)));
  //xcvs.push_back(ctr_cv(Dir3(0, 0, -1)));

  auto approx = traits.approximate_2_object();

  auto approx_arcs = get_approx_curves(xcvs, error);
  //std::vector<std::vector<QVector3D>>  arcs;
  //for (const auto& xcv : xcvs)
  //{
  //  std::vector<Approximate_point_2> v;
  //  auto oi2 = approx(xcv, error, std::back_insert_iterator(v));
  //
  //  std::vector<QVector3D> arc_points;
  //  for (const auto& p : v)
  //  {
  //    const QVector3D arc_point(p.dx(), p.dy(), p.dz());
  //    arc_points.push_back(arc_point);
  //  }
  //  arcs.push_back(std::move(arc_points));
  //}
  //std::cout << "offset count = " << m_arc_offsets.size() << std::endl;

  return approx_arcs;
}
Aos::Approx_arcs Aos::get_approx_arcs(const Kml::Placemark& placemark, double error)
{
  Geom_traits traits;
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  auto xcvs = get_arcs(placemark);

  auto approx = traits.approximate_2_object();
  std::vector<std::vector<QVector3D>>  arcs;
  for (const auto& xcv : xcvs)
  {
    std::vector<Approximate_point_2> v;
    auto oi2 = approx(xcv, error, std::back_insert_iterator(v));

    std::vector<QVector3D> arc_points;
    for (const auto& p : v)
    {
      const QVector3D arc_point(p.dx(), p.dy(), p.dz());
      arc_points.push_back(arc_point);
    }
    arcs.push_back(std::move(arc_points));
  }
  //std::cout << "offset count = " << m_arc_offsets.size() << std::endl;

  return arcs;
}

void Aos::check(const Kml::Placemarks& placemarks)
{
  Geom_traits traits;
  Arrangement arr(&traits);

  auto xcvs = get_arcs(placemarks, arr);
  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;
  // add arcs
  for (auto xcv : xcvs)
    CGAL::insert(arr, xcv);

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

std::vector<QVector3D> Aos::ext_check(const Kml::Placemarks& placemarks)
{
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Ext_aos arr(&traits);

  std::cout << "-------------------------------\n";
  std::cout << "** num arr FACES (before adding arcs) = " <<
    arr.number_of_faces() << std::endl;

  auto xcvs = get_arcs(placemarks, arr);

  // MARK all vertices as true
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    vit->set_data(Flag(true));
  }

  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;

  // add arcs
  for (auto& xcv : xcvs)
    CGAL::insert(arr, xcv);

  // extract all vertices that are ADDED when inserting the arcs!
  int num_created_vertices = 0;
  std::vector<QVector3D> created_vertices;
  auto approx = traits.approximate_2_object();
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    auto& d = vit->data();
    if (vit->data().v == false)
    {
      std::cout << "-------------------------------------\n";
      std::cout << vit->point() << std::endl;

      if (2 == vit->degree())
        ;//continue;

      if (1 == vit->degree())
      {
        auto p = vit->point();
        auto p2 = p.location();
        std::cout << "   deg-1 vertex = " << p << std::endl;
        std::cout << "   deg-1 vertex: " << std::boolalpha << vit->incident_halfedges()->target()->data().v << std::endl;
      }


      num_created_vertices++;
      auto p = vit->point();
      auto ap = approx(p);
      QVector3D new_vertex(ap.dx(), ap.dy(), ap.dz());
      new_vertex.normalize();
      std::cout << new_vertex << std::endl;
      std::cout << "degree = " << vit->degree() << std::endl;

      created_vertices.push_back(new_vertex);

      // find the arcs that are adjacent to the vertex of degree 4
      if (4 == vit->degree())
      {
        std::cout << "**************************\n DEGREE 4 VERTEX: \n";
        const auto first = vit->incident_halfedges();
        auto curr = first;
        do {
          auto tvh = curr->twin()->target();
          //std::cout << std::boolalpha << svh->data().v << " - " << tvh->data().v << std::endl;
          auto it = vertex_node_map.find(tvh);
          if (it != vertex_node_map.end())
            std::cout << std::setprecision(16) << it->second << std::endl;
          else
            std::cout << "NOT FOUND!!\n";
        } while (++curr != first);
      }

      std::cout << "\n";
    }
  }
  Kml::Node n{ 180.0, -84.71338 };
  std::cout << "Node itself = " << n.get_coords_3d() << std::endl;
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

  return created_vertices;
}


std::vector<QVector3D>  Aos::ext_check_id_based(Kml::Placemarks& placemarks)
{
  // Construct the arrangement from 12 geodesic arcs.
  Geom_traits traits;
  Ext_aos arr(&traits);

  // 
  auto nodes = Kml::generate_ids(placemarks);
  //auto nodes = Kml::generate_ids_approx(placemarks, 0.001);

  //Segment s()
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  int num_counted_arcs = 0;
  int num_counted_polygons = 0;
  // 
  std::vector<Point> points;
  std::vector<Ext_aos::Vertex_handle> vertices;
  for (const auto& node : nodes)
  {
    auto n = node.get_coords_3d();
    auto p = ctr_p(n.x, n.y, n.z);
    auto v = CGAL::insert_point(arr, p);
    points.push_back(p);
    vertices.push_back(v);
    //arr.insert_at_vertices(Segment(p, p), v, v);
  }
  std::cout << "num nodes = " << nodes.size() << std::endl;
  std::cout << "num points = " << points.size() << std::endl;
  // MARK all vertices as true
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    vit->set_data(Flag(true));
  }


  for (auto& placemark : placemarks)
  {
    for (auto& polygon : placemark.polygons)
    {
      num_counted_polygons++;

      // TO DO : ADD the outer boundaries!
      auto& ids = polygon.outer_boundary.ids;
      int num_nodes = ids.size();
      for (int i = 0; i < num_nodes - 1; ++i)
      {
        num_counted_arcs++;
        const auto nid1 = ids[i];
        const auto nid2 = ids[i + 1];
        auto p1 = points[nid1];
        auto p2 = points[nid2];
        CGAL::insert(arr, ctr_cv(p1, p2));
      }
    }
  }

  std::cout << "-------------------------------\n";
  std::cout << "num arr vertices (before adding arcs) = " <<
    arr.number_of_vertices() << std::endl;

  // extract all vertices that are ADDED when inserting the arcs!
  int num_created_vertices = 0;
  std::vector<QVector3D> created_vertices;
  auto approx = traits.approximate_2_object();
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    auto& d = vit->data();
    if (vit->data().v == false)
    {
      std::cout << "-------------------------------------\n";
      std::cout << vit->point() << std::endl;

      if (2 == vit->degree())
        ;//continue;

      if (1 == vit->degree())
      {
        auto p = vit->point();
        auto p2 = p.location();
        std::cout << "deg-1 vertex = " << p << std::endl;
        std::cout << "deg-1 vertex: " << std::boolalpha << vit->incident_halfedges()->target()->data().v << std::endl;
      }


      num_created_vertices++;
      auto p = vit->point();
      auto ap = approx(p);
      QVector3D new_vertex(ap.dx(), ap.dy(), ap.dz());
      new_vertex.normalize();
      std::cout << new_vertex << std::endl;
      std::cout << "degree = " << vit->degree() << std::endl;

      created_vertices.push_back(new_vertex);

      //// find the arcs that are adjacent to this vertex
      //const auto first = vit->incident_halfedges();
      //auto curr = first;
      //do {

      //} while (++curr != first);
      std::cout << std::endl;
    }
  }
  std::cout << "*** num created vertices = " << num_created_vertices << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num nodes = " << nodes.size() << std::endl;
  std::cout << "num arr vertices = " << arr.number_of_vertices() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num counted arcs = " << num_counted_arcs << std::endl;
  std::cout << "num arr edges = " << arr.number_of_edges() << std::endl;

  std::cout << "-------------------------------\n";
  std::cout << "num polygons = " << num_counted_polygons << std::endl;
  std::cout << "num arr faces = " << arr.number_of_faces() << std::endl;

  return created_vertices;
}


Aos::Approx_arcs  Aos::find_new_faces(Kml::Placemarks& placemarks)
{
  Geom_traits traits;
  Ext_aos arr(&traits);
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  using Face_handle = Ext_aos::Face_handle;
  auto fh = arr.faces_begin();
  fh->data().v = true;
  std::cout << "num faces = " << arr.number_of_faces() << std::endl;

  auto nodes = Kml::generate_ids(placemarks);


  //-------------------------------------------------------------------------
  // define a set of vertex-handles: use this to check if the face is 
  // obtained from the polygon definition, or if it is an additional face
  using Vertex_handle = Ext_aos::Vertex_handle;
  std::map<Vertex_handle, int>  vertex_id_map;
  std::set<std::set<int>>  all_polygon_node_ids;


  num_counted_nodes = 0;
  num_counted_arcs = 0;
  num_counted_polygons = 0;
  std::vector<Curve>  xcvs;
  for (auto& pm : placemarks)
  {
    std::cout << pm.name << std::endl;
    for (auto& polygon : pm.polygons)
    {
      num_counted_polygons++;

      // colect all rings into a single list (FOR NOW!!!)
      // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
      auto linear_rings = polygon.get_all_boundaries();
      //Kml::LinearRings linear_rings;
      //linear_rings.push_back(polygon.outer_boundary);
      //for (const auto& inner_boundary : polygon.inner_boundaries)
      //  linear_rings.push_back(inner_boundary);

      // loop on outer and inner boundaries 
      //for (auto* lring : linear_rings)
      auto* lring = &polygon.outer_boundary;
      {
        int num_faces_before = arr.number_of_faces();
        std::set<int> polygon_node_ids;

        // convert the nodes to points on unit-sphere
        std::vector<Approximate_Vector_3>  sphere_points;
        //for (const auto& node : lring->nodes)
        //std::cout << "   NUM POLYGON-NODES SIZE = " << lring->ids.size() << std::endl;
        for (int i = 0; i < lring->ids.size(); ++i)
        {
          num_counted_nodes++;
          const auto id = lring->ids[i];
          const auto& node = lring->nodes[i];
          const auto p = node.get_coords_3d();
          Approximate_Vector_3  v(p.x, p.y, p.z);
          sphere_points.push_back(v);
          auto vh = CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
          polygon_node_ids.insert(id);
          // assert node-id and vertex-handle consistency
          //{
          //  auto it = vertex_id_map.find(vh);
          //  if (vertex_id_map.cend() != it)
          //  {
          //    if (id != it->second)
          //      std::cout << "*** ERROR!!!\n";
          //  }
          //}
          vertex_id_map[vh] = id;
          vh->data().v = true;
        }
        //std::cout << "   POLYGON-NODES SET SIZE = " << polygon_node_ids.size() << std::endl;
        if (lring->ids.size() != (1 + polygon_node_ids.size()))
          std::cout << "*** ASSERTION ERROR!!!!\n";

        all_polygon_node_ids.insert(std::move(polygon_node_ids));

        // add curves
        int num_points = sphere_points.size();
        for (int i = 0; i < num_points - 1; i++)
        {
          num_counted_arcs++;
          const auto p1 = sphere_points[i];
          const auto p2 = sphere_points[i + 1];
          auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
            ctr_p(p2.x(), p2.y(), p2.z()));
          //xcvs.push_back(xcv);
          CGAL::insert(arr, xcv);
        }

        int num_faces_after = arr.number_of_faces();
        int num_new_faces = num_faces_after - num_faces_before;
      }
    }
  }


  // mark all faces as TRUE (= as existing faces)
  int num_found = 0;
  int num_not_found = 0;
  std::vector<Curve>  new_face_arcs;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh)
  {
    // skip the spherical face
    std::cout << "num outer_ccbs = " << fh->number_of_outer_ccbs() << std::endl;
    if (fh->number_of_outer_ccbs() == 0)
    {
      continue;
    }


    // construct the set of all node-ids for the current face
    std::set<int>  face_node_ids_set;
    std::vector<int>  face_node_ids;
    std::vector<Curve>  face_arcs;
    auto first = fh->outer_ccb();
    auto curr = first;
    do {
      auto vh = curr->source();
      // skip if the vertex is due to intersection with the identification curve
      if ((vh->data().v == false) && (vh->degree() == 2))
        continue;

      auto id = vertex_id_map[vh];
      face_node_ids_set.insert(id);

      face_arcs.push_back(ctr_cv(curr->source()->point(), curr->target()->point()));
    } while (++curr != first);
    //std::cout << "counted vertices = " << num_vertices << std::endl;
    //std::cout << "vertices in the set = " << polygon_node_ids.size() << std::endl;

    auto it = all_polygon_node_ids.find(face_node_ids_set);
    if (it == all_polygon_node_ids.cend())
    {
      std::cout << "NOT FOUND!!!\n";
      std::cout << "num nodes = " << face_node_ids_set.size() << std::endl;
      num_not_found++;
      new_face_arcs.insert(new_face_arcs.end(), face_arcs.begin(), face_arcs.end());
    }
    else
      num_found++;
  }
  std::cout << "num not found = " << num_not_found << std::endl;

  auto approx_arcs = get_approx_curves(new_face_arcs, 0.001);
  return approx_arcs;
}


void Aos::save_arr(Kml::Placemarks& placemarks, const std::string& file_name)
{
#ifndef USE_EPIC
  Geom_traits traits;
  Ext_aos arr(&traits);
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  using Face_handle = Ext_aos::Face_handle;
  auto fh = arr.faces_begin();
  fh->data().v = true;
  std::cout << "num faces = " << arr.number_of_faces() << std::endl;

  auto nodes = Kml::generate_ids(placemarks);


  //-------------------------------------------------------------------------
  // define a set of vertex-handles: use this to check if the face is 
  // obtained from the polygon definition, or if it is an additional face
  using Vertex_handle = Ext_aos::Vertex_handle;
  using Halfedge_handle = Ext_aos::Halfedge_handle;
  using Face_handle = Ext_aos::Face_handle;
  std::map<Vertex_handle, int>  vertex_id_map;
  std::map<std::set<int>, std::string>  all_polygon_node_ids_map;

  // map to associate the created faces with the country names
  // CAUTION: the newly created faces

  num_counted_nodes = 0;
  num_counted_arcs = 0;
  num_counted_polygons = 0;
  std::vector<Curve>  xcvs;
  for (auto& pm : placemarks)
  {
    std::cout << pm.name << std::endl;
    for (auto& polygon : pm.polygons)
    {
      num_counted_polygons++;

      // colect all rings into a single list (FOR NOW!!!)
      // TO-DO: PROCESS OUTER & INNER BOUNDARIES SEPARATELY!!!
      //auto linear_rings = polygon.get_all_boundaries();
      //Kml::LinearRings linear_rings;
      //linear_rings.push_back(polygon.outer_boundary);
      //for (const auto& inner_boundary : polygon.inner_boundaries)
      //  linear_rings.push_back(inner_boundary);

      // loop on outer and inner boundaries 
      //for (auto* lring : linear_rings)
      auto* lring = &polygon.outer_boundary;
      {
        int num_faces_before = arr.number_of_faces();
        std::set<int> polygon_node_ids;

        // convert the nodes to points on unit-sphere
        std::vector<Approximate_Vector_3>  sphere_points;
        //for (const auto& node : lring->nodes)
        //std::cout << "   NUM POLYGON-NODES SIZE = " << lring->ids.size() << std::endl;
        for (int i = 0; i < lring->ids.size(); ++i)
        {
          num_counted_nodes++;
          const auto id = lring->ids[i];
          const auto& node = lring->nodes[i];
          const auto p = node.get_coords_3d();
          Approximate_Vector_3  v(p.x, p.y, p.z);
          sphere_points.push_back(v);
          auto vh = CGAL::insert_point(arr, ctr_p(p.x, p.y, p.z));
          polygon_node_ids.insert(id);
          // assert node-id and vertex-handle consistency
          //{
          //  auto it = vertex_id_map.find(vh);
          //  if (vertex_id_map.cend() != it)
          //  {
          //    if (id != it->second)
          //      std::cout << "*** ERROR!!!\n";
          //  }
          //}
          vertex_id_map[vh] = id;
          vh->data().v = true;
        }
        //std::cout << "   POLYGON-NODES SET SIZE = " << polygon_node_ids.size() << std::endl;
        if (lring->ids.size() != (1 + polygon_node_ids.size()))
          std::cout << "*** ASSERTION ERROR!!!!\n";

        all_polygon_node_ids_map.insert(std::make_pair(
          std::move(polygon_node_ids), pm.name));

        // add curves
        int num_points = sphere_points.size();
        for (int i = 0; i < num_points - 1; i++)
        {
          num_counted_arcs++;
          const auto p1 = sphere_points[i];
          const auto p2 = sphere_points[i + 1];
          auto xcv = ctr_cv(ctr_p(p1.x(), p1.y(), p1.z()),
            ctr_p(p2.x(), p2.y(), p2.z()));
          //xcvs.push_back(xcv);
          CGAL::insert(arr, xcv);
        }

        int num_faces_after = arr.number_of_faces();
        int num_new_faces = num_faces_after - num_faces_before;
      }
    }
  }

  std::cout << "*** arr.number_of_faces = " << arr.number_of_faces() << std::endl;
  std::cout << "*** arr.number_of_halfedges = " << arr.number_of_halfedges() << std::endl;
  std::cout << "*** arr.number_of_vertices = " << arr.number_of_vertices() << std::endl;

  // DEFINE JSON OBJECT
  json js;
  auto& js_points = js["points"] = json::array();

  ////////////////////////////////////////////////////////////////////////////
  // POINTS
  // define a map from each vertex to its position in the arrangement
  //auto get_num_denum  
  using FT = typename Kernel::FT;
  //using json = nlohmann::ordered_json;
  FT ft(0);
  auto ex = ft.exact();
  CGAL::Rational_traits<decltype(ex)> rt;
  typename CGAL::Algebraic_structure_traits<decltype(ex)>::Simplify simplify;


  auto set_num_denum = [&](decltype(ex)& x, json& ratx)
  {
    simplify(x);
    std::stringstream ss_x_num;
    CGAL::IO::set_ascii_mode(ss_x_num);
    ss_x_num << rt.numerator(x);
    std::string xnum;
    ss_x_num >> xnum;
    ratx["num"] = xnum;

    std::stringstream ss_x_den;
    CGAL::IO::set_ascii_mode(ss_x_den);
    ss_x_den << rt.denominator(x);
    std::string xden;
    ss_x_den >> xden;
    ratx["den"] = xden;
  };


  using Point_ = std::decay_t<decltype(arr.vertices_begin()->point())>;
  std::map<void*, int> point_pos_map;
  std::vector<Point_>  points;
  std::map<Vertex_handle, int>  vertex_pos_map;
  for (auto vh = arr.vertices_begin(); vh != arr.vertices_end(); ++vh)
  {
    // add the vertex if not found in the map
    auto it = vertex_pos_map.find(vh);
    if (it == vertex_pos_map.end())
    {
      int new_vh_pos = vertex_pos_map.size();
      vertex_pos_map[vh] = new_vh_pos;

      // write the vertex-data to JSON object
      auto& p = vh->point();
      points.push_back(p);
      point_pos_map[&p] = new_vh_pos;
      auto dx = p.dx().exact();
      auto dy = p.dy().exact();
      auto dz = p.dz().exact();

      json jv;
      jv["location"] = p.location();
      set_num_denum(dx, jv["dx"]);
      set_num_denum(dy, jv["dy"]);
      set_num_denum(dz, jv["dz"]);
      js_points.push_back(std::move(jv));
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // CURVES
  // define a map from each curve to its position in the arrangment
  auto& js_curves = js["curves"] = json::array();
  using Ext_curve = Ext_aos::X_monotone_curve_2;
  std::map<Ext_curve*, int>  curve_pos_map;
  int num_edges = 0;
  for (auto eh = arr.edges_begin(); eh != arr.edges_end(); ++eh)
  {
    num_edges++;
    auto& xcv = eh->curve();
    auto it = curve_pos_map.find(&xcv);
    if (it == curve_pos_map.end())
    {
      int new_xcv_pos = curve_pos_map.size();
      curve_pos_map[&xcv] = new_xcv_pos;

      json je;
      auto& sv = xcv.source();
      auto& tv = xcv.target();
      auto svit = std::find(points.begin(), points.end(), sv);
      auto tvit = std::find(points.begin(), points.end(), tv);
      auto svp = std::distance(points.begin(), svit);
      auto tvp = std::distance(points.begin(), tvit);
      //std::cout << "svp= " << point_pos_map[(void*)&sv] << std::endl;
      je["source"] = svp;
      je["target"] = tvp;
      auto& je_normal = je["normal"];

      // write the vertex-data to JSON object
      auto& n = xcv.normal();
      auto dx = n.dx().exact();
      auto dy = n.dy().exact();
      auto dz = n.dz().exact();
      set_num_denum(dx, je_normal["dx"]);
      set_num_denum(dy, je_normal["dy"]);
      set_num_denum(dz, je_normal["dz"]);

      je["is_vertical"] = xcv.is_vertical();
      je["is_directed_right"] = xcv.is_directed_right();
      je["is_full"] = xcv.is_full();

      js_curves.push_back(std::move(je));
    }
  }
  std::cout << "num edges = " << num_edges << std::endl;
  std::cout << "curve map size = " << curve_pos_map.size() << std::endl;
  std::cout << "js_curves size = " << js_curves.size() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // VERTICES
  // there is a one-to-one corresponce between vertices and points
  auto& js_vertices = js["vertices"] = json::array();
  for (auto vh = arr.vertices_begin(); vh != arr.vertices_end(); ++vh)
  {
    json js_vertex;
    auto vpos = vertex_pos_map[vh];
    js_vertex["point"] = vpos;
    js_vertices.push_back(std::move(js_vertex));
  }

  ////////////////////////////////////////////////////////////////////////////
  // HALF-EDGES
  //int num_half_edges = 0;
  //auto& js_halfedges = js["halfedges"] = json::array();
  //std::map<Halfedge_handle, int>  halfedge_pos_map;
  //auto write_half_edge = [&](auto& he)
  //{
  //  auto it = halfedge_pos_map.find(&he);
  //  if (it == halfedge_pos_map.end())
  //  {
  //    auto new_he_pos = halfedge_pos_map.size();
  //    halfedge_pos_map[&he] = new_he_pos;
  //  }
  //  auto svh = he.source();
  //  auto tvh = he.target();
  //  auto& xcv = he.curve();
  //  auto svp = vertex_pos_map[svh];
  //  auto tvp = vertex_pos_map[tvh];
  //  auto xcvp = curve_pos_map[&xcv];
  //  json js_he;
  //  js_he["source"] = svp;
  //  js_he["target"] = tvp;
  //  js_he["curve"] = xcvp;
  //  js_halfedges.push_back(std::move(js_he));
  //};
  //for (auto it = arr.halfedges_begin(); it != arr.halfedges_end(); ++it)
  //{
  //  auto& he = *it;
  //  write_half_edge(he);
  //  num_half_edges++;
  //}
  //std::cout << "HALF-EDGE CHECKS:\n";
  //std::cout << "  *** num total half-edges = " << num_half_edges << std::endl;
  //std::cout << "  *** halfedge-pos-map size = " << halfedge_pos_map.size() << std::endl;


  ////////////////////////////////////////////////////////////////////////////
  // EDGES
  num_edges = 0;
  auto& js_edges = js["edges"] = json::array();
  ////using Edge_ = decltype(*arr.edges_begin());
  //std::map<void*, int>  edge_pos_map;
  std::map<void*, int>  halfedge_pos_map;
  //Edge_const_iterator    eit;
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
  {
    auto& edge = *eit;
    //auto it = edge_pos_map.find(&edge);
    //if (it == edge_pos_map.end())
    //{
    //  auto new_edge_pos = edge_pos_map.size();
    //  edge_pos_map[&edge] = new_edge_pos;
    //}
    auto& xcv = edge.curve();
    auto xcvp = curve_pos_map[&xcv];
    json js_edge;
    js_edge["curve"] = xcvp;
    js_edge["direction"] = edge.direction();
    js_edge["source"] = vertex_pos_map[edge.source()];
    js_edge["target"] = vertex_pos_map[edge.target()];
    js_edges.push_back(std::move(js_edge));
    num_edges++;

    // add the halfedge indices to the map
    int new_halfedge_index = halfedge_pos_map.size();
    auto& twin = *edge.twin();
    halfedge_pos_map[&edge] = new_halfedge_index;
    halfedge_pos_map[&twin] = new_halfedge_index + 1;
  }
  std::cout << "EDGE CHECKS:\n";
  std::cout << "  *** num edges = " << num_edges << std::endl;
  std::cout << " *** js_edges size = " << js_edges.size() << std::endl;


  ////////////////////////////////////////////////////////////////////////////
  // FACES
  // CONDITION DATA: Caspian Sea needs to be defined
  //num_half_edges = 0;
  num_edges = 0;
  int num_found = 0;
  int num_not_found = 0;
  std::map<Face_handle, std::string>  face_name_map;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh)
  {
    // skip the spherical face
    std::cout << "num outer_ccbs = " << fh->number_of_outer_ccbs() << std::endl;
    if (fh->number_of_outer_ccbs() == 0)
    {
      continue;
    }

    // construct the set of all node-ids for the current face
    std::set<int>  face_node_ids_set;
    std::vector<int>  face_node_ids;
    auto first = fh->outer_ccb();
    auto curr = first;
    do {
      //num_half_edges++;
      auto vh = curr->source();
      // skip if the vertex is due to intersection with the identification curve
      if ((vh->data().v == false) && (vh->degree() == 2))
        continue;

      auto id = vertex_id_map[vh];
      face_node_ids_set.insert(id);

      //face_arcs.push_back(ctr_cv(curr->source()->point(), curr->target()->point()));
      auto& xcv = curr->curve();
    } while (++curr != first);
    //std::cout << "counted vertices = " << num_vertices << std::endl;
    //std::cout << "vertices in the set = " << polygon_node_ids.size() << std::endl;

    std::string name;
    auto it = all_polygon_node_ids_map.find(face_node_ids_set);
    if (it == all_polygon_node_ids_map.cend())
    {
      std::cout << "NOT FOUND!!!\n";
      std::cout << "num nodes = " << face_node_ids_set.size() << std::endl;
      num_not_found++;
      name = "Caspian Sea";
    }
    else
    {
      num_found++;
      name = it->second;
    }
    face_name_map[fh] = name;
  }
  std::cout << "num not found = " << num_not_found << std::endl;


  // RECORD FACES
  json& js_faces = js["faces"] = json::array();
  auto get_ccb_json = [&](Ext_aos::Ccb_halfedge_circulator first)
  {
    json js_edges;
    auto& ccb_edge_indices = js_edges["halfedges"] = json::array();
    auto curr = first;
    do {
      auto& he = *curr;
      //auto& xcv = he.curve();
      //auto it = curve_pos_map.find(&xcv);
      //if (it == curve_pos_map.end())
      //{
      //    std::cout << "ASSERTION ERROR!!!" << std::endl;
      //}
      auto it = halfedge_pos_map.find(&he);
      if (it == halfedge_pos_map.end())
      {
        std::cout << "ASSERTION ERROR!!!" << std::endl;
      }

      auto edge_pos = it->second;
      ccb_edge_indices.push_back(edge_pos);
    } while (++curr != first);

    return js_edges;
  };

  int total_num_half_edges = 0;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh)
  {
    //// skip the spherical face
    //if (fh->number_of_outer_ccbs() == 0)
    //  continue;

    // json object for the current face
    json js_face;
    auto face_name = face_name_map[fh];
    js_face["name"] = face_name;
    js_face["is_unbounded"] = false;
    js_face["is_valid"] = true;

    // at this point we are sure that we have at least 1 outer-ccb
    auto& js_outer_ccbs = js_face["outer_ccbs"] = json::array();
    for (auto ccb = fh->outer_ccbs_begin(); ccb != fh->outer_ccbs_end(); ++ccb)
    {
      auto js_ccb = get_ccb_json(*ccb);
      js_outer_ccbs.push_back(std::move(js_ccb));
    }

    // INNER CCBS
    if (fh->number_of_inner_ccbs() > 0)
    {
      auto& js_inner_ccbs = js_face["inner_ccbs"] = json::array();
      for (auto ccb = fh->inner_ccbs_begin(); ccb != fh->inner_ccbs_end(); ++ccb)
      {
        auto js_ccb = get_ccb_json(*ccb);
        js_inner_ccbs.push_back(std::move(js_ccb));
      }
    }

    js_faces.push_back(std::move(js_face));
  }
  std::cout << "total num half-edges = " << total_num_half_edges << std::endl;

  // save the arrangment
  std::ofstream ofile(file_name);
  ofile << js.dump(2);
  ofile.close();
#endif
}


namespace
{
  enum Error_id {
    FILE_NOT_FOUND,
    ILLEGAL_EXTENSION,
    UNABLE_TO_OPEN,
    FILE_IS_EMPTY,
    INVALID_INITIAL_POLYGON,
    UNSUPPORTED,
    INVALID_OUTPUT_FILE,
    ILLEGAL_SOLUTION_FILE
  };

  struct Illegal_input : public std::logic_error {
    Illegal_input(Error_id /* err */, const std::string& msg,
      const std::string& filename) :
      std::logic_error(std::string(msg).append(" (").append(filename).
        append(")!"))
    {}

    Illegal_input(Error_id /* err */, const std::string& msg) :
      std::logic_error(std::string(msg).append("!"))
    {}
  };

  struct Input_file_missing_error : public std::logic_error {
    Input_file_missing_error(std::string& str) : std::logic_error(str) {}
  };

  //
  struct country {
    //! Constructor
    country(std::string& name) :
      m_name(std::move(name))
    {}

    std::string m_name;
  };

  /*! Read a json file.
 */
  bool read_json(const std::string& filename, nlohmann::json& data) {
    using json = nlohmann::json;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
      throw Illegal_input(UNABLE_TO_OPEN, "Cannot open file", filename);
      return false;
    }
    data = json::parse(infile);
    infile.close();
    if (data.empty()) {
      throw Illegal_input(FILE_IS_EMPTY, "File is empty", filename);
      return false;
    }
    return true;
  }

  template <typename FT>
  FT to_ft(const nlohmann::json& js_ft) {
    using Exact_type = typename FT::Exact_type;
    const std::string& js_num = js_ft["num"];
    const std::string& js_den = js_ft["den"];
    std::string str = js_num + "/" + js_den;
    Exact_type eft(str);
    return FT(eft);
  }

  template <typename Arrangement_, typename Kernel_>
  bool read_arrangement(const std::string& filename, Arrangement_& arr,
    const Kernel_& kernel) {
    using Arrangement = Arrangement_;
    using Kernel = Kernel_;

    using json = nlohmann::json;
    json data;
    auto rc = read_json(filename, data);
    if (!rc) return false;

    // points
    auto it = data.find("points");
    if (it == data.end()) {
      std::cerr << "The points item is missing " << " (" << filename << ")\n";
      return false;
    }
    const auto& js_points = it.value();

    // curves
    it = data.find("curves");
    if (it == data.end()) {
      std::cerr << "The curves item is missing " << " (" << filename << ")\n";
      return false;
    }
    const auto& js_curves = it.value();

    // vertices
    it = data.find("vertices");
    if (it == data.end()) {
      std::cerr << "The vertices item is missing " << " (" << filename << ")\n";
      return false;
    }
    const auto& js_vertices = it.value();

    // edges
    it = data.find("edges");
    if (it == data.end()) {
      std::cerr << "The edges item is missing " << " (" << filename << ")\n";
      return false;
    }
    const auto& js_edges = it.value();

    // faces
    it = data.find("faces");
    if (it == data.end()) {
      std::cerr << "The faces item is missing " << " (" << filename << ")\n";
      return false;
    }
    const auto& js_faces = it.value();

    const std::size_t num_points = js_points.size();
    const std::size_t num_curves = js_curves.size();
    const std::size_t num_vertices = js_vertices.size();
    const std::size_t num_edges = js_edges.size();
    const std::size_t num_faces = js_faces.size();
    const std::size_t num_halfedges = num_edges * 2;

    if (num_points < num_vertices) {
      std::cerr << "The no. of points (" << num_points
        << ") is smaller than the no. of vertices (" << num_vertices
        << ")\n";
      return false;
    }

    if (num_curves < num_edges) {
      std::cerr << "The no. of curves (" << num_curves
        << ") is smaller than the no. of edge (" << num_edges << ")\n";
      return false;
    }

    std::cout << "# points: " << num_points << std::endl;
    std::cout << "# curves: " << num_curves << std::endl;
    std::cout << "# vertices: " << num_vertices << std::endl;
    std::cout << "# halfedges: " << num_halfedges << std::endl;
    std::cout << "# faces: " << num_faces << std::endl;

    using Point = typename Arrangement::Point_2;
    using X_monotone_curve = typename Arrangement::X_monotone_curve_2;
    using FT = typename Kernel::FT;
    using Exact_type = typename FT::Exact_type;

    std::vector<Point> points;
    points.reserve(num_points);
    for (const auto& js_pnt : js_points) {
      using Direction_3 = typename Kernel::Direction_3;
      using Location = typename Point::Location_type;
      auto location = static_cast<Location>(js_pnt["location"]);
      auto dx = to_ft<FT>(js_pnt["dx"]);
      auto dy = to_ft<FT>(js_pnt["dy"]);
      auto dz = to_ft<FT>(js_pnt["dz"]);
      Direction_3 dir(dx, dy, dz);
      Point pnt(dir, location);
      points.push_back(pnt);
    }

    std::vector<X_monotone_curve> xcurves;
    xcurves.reserve(num_curves);
    for (const auto& js_xcv : js_curves) {
      using Direction_3 = typename Kernel::Direction_3;
      std::size_t src_id = js_xcv["source"];
      std::size_t trg_id = js_xcv["target"];
      const auto& js_normal = js_xcv["normal"];
      auto dx = to_ft<FT>(js_normal["dx"]);
      auto dy = to_ft<FT>(js_normal["dy"]);
      auto dz = to_ft<FT>(js_normal["dz"]);
      Direction_3 normal(dx, dy, dz);
      bool is_vert = js_xcv["is_vertical"];
      bool is_directed_right = js_xcv["is_directed_right"];
      bool is_full = js_xcv["is_full"];
      const auto& src = points[src_id];
      const auto& trg = points[trg_id];
      X_monotone_curve xcv(src, trg, normal, is_vert, is_directed_right, is_full);
      xcurves.push_back(xcv);
    }

    using Arr_accessor = CGAL::Arr_accessor<Arrangement>;
    Arr_accessor arr_access(arr);
    arr_access.clear_all();

    // Vertices
    std::cout << "Vertices\n";
    using DVertex = typename Arr_accessor::Dcel_vertex;
    std::vector<DVertex*> vertices(num_vertices);
    size_t k = 0;
    for (const auto& js_vertex : js_vertices) {
      std::size_t point_id = js_vertex["point"];
      const auto& point = points[point_id];
      CGAL::Arr_parameter_space ps_x, ps_y;
      switch (point.location()) {
      case Point::NO_BOUNDARY_LOC: ps_x = ps_y = CGAL::INTERIOR; break;
      case Point::MIN_BOUNDARY_LOC:
        ps_x = CGAL::INTERIOR;
        ps_y = CGAL::ARR_BOTTOM_BOUNDARY;
        break;
      case Point::MID_BOUNDARY_LOC:
        ps_x = CGAL::LEFT_BOUNDARY;
        ps_y = CGAL::INTERIOR;
        break;
      case Point::MAX_BOUNDARY_LOC:
        ps_x = CGAL::INTERIOR;
        ps_y = CGAL::ARR_TOP_BOUNDARY;
        break;
      }
      vertices[k++] = arr_access.new_vertex(&point, ps_x, ps_y);
    }

    // Halfedges
    std::cout << "Halfedges\n";
    using DHalfedge = typename Arr_accessor::Dcel_halfedge;
    std::vector<DHalfedge*> halfedges(num_halfedges);
    k = 0;
    for (const auto& js_edge : js_edges) {
      std::size_t source_id = js_edge["source"];
      std::size_t target_id = js_edge["target"];
      std::size_t curve_id = js_edge["curve"];
      int direction = js_edge["direction"];
      DVertex* src_v = vertices[source_id];
      DVertex* trg_v = vertices[target_id];
      const auto& curve = xcurves[curve_id];
      DHalfedge* new_he = arr_access.new_edge(&curve);
      trg_v->set_halfedge(new_he);
      new_he->set_vertex(trg_v);
      src_v->set_halfedge(new_he->opposite());
      new_he->opposite()->set_vertex(src_v);
      new_he->set_direction(static_cast<CGAL::Arr_halfedge_direction>(direction));
      halfedges[k++] = new_he;
      halfedges[k++] = new_he->opposite();
    }

    // Faces
    std::cout << "Faces\n";
    using DFace = typename Arr_accessor::Dcel_face;
    using DOuter_ccb = typename Arr_accessor::Dcel_outer_ccb;
    using DInner_ccb = typename Arr_accessor::Dcel_inner_ccb;
    using DIso_vert = typename Arr_accessor::Dcel_isolated_vertex;
    const bool is_unbounded(false);
    const bool is_valid(true);
    for (const auto& js_face : js_faces) {
      DFace* new_f = arr_access.new_face();

      new_f->set_unbounded(is_unbounded);
      new_f->set_fictitious(!is_valid);

      // Read the outer CCBs of the face.
      auto oit = js_face.find("outer_ccbs");
      if (oit != js_face.end()) {
        const auto& js_outer_ccbs = *oit;
        for (const auto& js_ccb : js_outer_ccbs) {
          std::cout << "Outer CCB\n";
          // Allocate a new outer CCB record and set its incident face.
          auto* new_occb = arr_access.new_outer_ccb();
          new_occb->set_face(new_f);

          // Read the current outer CCB.
          auto bit = js_ccb.find("halfedges");
          if (bit == js_ccb.end()) {
            std::cerr << "The halfedges item is missing " << " (" << filename
              << ")\n";
            return false;
          }

          const auto& js_halfedges = *bit;
          auto hit = js_halfedges.begin();
          std::size_t first_idx = *hit;
          DHalfedge* first_he = halfedges[first_idx];
          first_he->set_outer_ccb(new_occb);
          DHalfedge* prev_he = first_he;
          for (++hit; hit != js_halfedges.end(); ++hit) {
            std::size_t curr_idx = *hit;
            auto curr_he = halfedges[curr_idx];
            prev_he->set_next(curr_he);		// connect
            curr_he->set_outer_ccb(new_occb);	// set the CCB
            prev_he = curr_he;
          }
          prev_he->set_next(first_he);		// close the loop
          new_f->add_outer_ccb(new_occb, first_he);
        }
      }

      // Read the inner CCBs of the face.
      auto iit = js_face.find("inner_ccbs");
      if (iit != js_face.end()) {
        const auto& js_inner_ccbs = *iit;
        for (const auto& js_ccb : js_inner_ccbs) {
          std::cout << "Inner CCB\n";
          // Allocate a new inner CCB record and set its incident face.
          auto* new_iccb = arr_access.new_inner_ccb();
          new_iccb->set_face(new_f);

          // Read the current inner CCB.
          auto bit = js_ccb.find("halfedges");
          if (bit == js_ccb.end()) {
            std::cerr << "The halfedges item is missing " << " (" << filename
              << ")\n";
            return false;
          }

          const auto& js_halfedges = *bit;
          auto hit = js_halfedges.begin();
          std::size_t first_idx = *hit;
          DHalfedge* first_he = halfedges[first_idx];
          first_he->set_inner_ccb(new_iccb);
          DHalfedge* prev_he = first_he;
          for (++hit; hit != js_halfedges.end(); ++hit) {
            std::size_t curr_idx = *hit;
            auto curr_he = halfedges[curr_idx];
            prev_he->set_next(curr_he);		// connect
            curr_he->set_inner_ccb(new_iccb);	// set the CCB
            prev_he = curr_he;
          }
          prev_he->set_next(first_he);		// close the loop
          new_f->add_inner_ccb(new_iccb, first_he);
        }
      }

      // // Read the isolated vertices inside the face.
      // Size n_isolated_vertices =
      //   formatter.read_size("number_of_isolated_vertices");
      // if (n_isolated_vertices) {
      //   formatter.read_isolated_vertices_begin();
      //   Size k;
      //   for (k = 0; k < n_isolated_vertices; k++) {
      //     // Allocate a new isolated vertex record and set its incident face.
      //     DIso_vert* new_iso_vert = arr_access.new_isolated_vertex();
      //     new_iso_vert->set_face(new_f);
      //     // Read the current isolated vertex.
      //     std::size_t v_idx = formatter.read_vertex_index();
      //     DVertex* iso_v = m_vertices[v_idx];
      //     iso_v->set_isolated_vertex(new_iso_vert);
      //     new_f->add_isolated_vertex(new_iso_vert, iso_v);
      //   }
      // }
    }

    return true;
  }

  Kernel kernel;
  Geom_traits traits;
}


Aos::Arr_handle Aos::load_arr(const std::string& file_name)
{
  Arrangement arr(&traits);;
  auto rc = read_arrangement(file_name, arr, kernel);
  if (!rc) {
    std::cerr << "Failed to load database!\n";
    return nullptr;
  }

  return &arr;
}


Aos::Arr_handle  Aos::construct(Kml::Placemarks& placemarks)
{
  Geom_traits traits;
  auto* arr = new Arrangement(&traits);

  auto xcvs = get_arcs(placemarks, *arr);
  for (auto& xcv : xcvs)
    CGAL::insert(*arr, xcv);

  return arr;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>
//#include <CGAL/Projection_traits_3.h>

#include <iostream>
#include <unordered_map>
#include <boost/property_map/property_map.hpp>

std::vector<QVector3D> Aos::get_triangles(Arr_handle arrh)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
  //  typedef CGAL::Projection_traits_3<K_epic>                         K;
  typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K>            Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              TDS;
  typedef CGAL::Exact_predicates_tag                                Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
  typedef CDT::Face_handle                                          Face_handle;
  typedef CDT::Point                                                Point;
  typedef CGAL::Polygon_2<K>                                        Polygon_2;

  auto& arr = *reinterpret_cast<Arrangement*>(arrh);

  Geom_traits traits;
  auto approx = traits.approximate_2_object();


  std::vector<std::vector<QVector3D>> all_faces;
  // loop on all faces of the arrangement
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
  {
    // skip any face with no OUTER-CCB
    if (0 == fit->number_of_outer_ccbs())
      continue;

    // COMPUTE THE CENTROID OF ALL FACE-POINTS
    std::vector<QVector3D> face_points;

    // loop on the egdes of the current outer-ccb
    auto first = fit->outer_ccb();
    auto curr = first;
    do {
      auto ap = approx(curr->source()->point());
      QVector3D p(ap.dx(), ap.dy(), ap.dz());
      p.normalize();
      face_points.push_back(p);
    } while (++curr != first);

    all_faces.push_back(std::move(face_points));
  }

  // RESULTING TRIANGLE POINTS (every 3 point => triangle)
  std::vector<QVector3D>  triangles;

  std::cout << "triangulating individual faces\n";

  // loop on all approximated faces
  for (auto& face_points : all_faces)
  {
    std::cout << "num face points = " << face_points.size() << std::endl;
    // no need to triangulate if the number of points is 3
    if (face_points.size() == 3)
    {
      triangles.insert(triangles.end(), face_points.begin(), face_points.end());
      continue;
    }


    // find the centroid of all face-points
    QVector3D centroid(0, 0, 0);
    for (const auto& fp : face_points)
      centroid += fp;
    centroid /= face_points.size();
    centroid.normalize();
    auto normal = centroid;

    K::Point_3  plane_origin(centroid.x(), centroid.y(), centroid.z());
    K::Vector_3 plane_normal(normal.x(), normal.y(), normal.z());
    K::Plane_3 plane(plane_origin, plane_normal);

    Polygon_2 polygon;

    // project all points onto the plane
    K::Point_3 origin(0, 0, 0);
    for (const auto& fp : face_points)
    {
      // define a ray through the origin and the current point
      K::Point_3 current_point(fp.x(), fp.y(), fp.z());
      K::Ray_3 ray(origin, current_point);

      auto intersection = CGAL::intersection(plane, ray);
      if (!intersection.has_value())
        std::cout << "INTERSECTION ASSERTION ERROR!!!\n";
      auto ip = boost::get<K::Point_3>(intersection.value());
      auto ip2 = plane.to_2d(ip);

      // add this to the polygon constraint
      polygon.push_back(ip2);
    }

    CDT cdt;
    cdt.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);

    std::unordered_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<Face_handle, bool> >
      in_domain(in_domain_map);

    //Mark facets that are inside the domain bounded by the polygon
    CGAL::mark_domain_in_triangulation(cdt, in_domain);

    // loop on all the triangles ("faces" in triangulation doc)
    for (Face_handle f : cdt.finite_face_handles())
    {
      // if the current triangles is not inside the polygon -> skip it
      if (false == get(in_domain, f))
        continue;

      for (int i = 0; i < 3; ++i)
      {
        auto tp = f->vertex(i)->point();
        auto tp3 = plane.to_3d(tp);
        QVector3D p3(tp3.x(), tp3.y(), tp3.z());
        p3.normalize();
        triangles.push_back(p3);
      }
    }
  }

  return triangles;
}