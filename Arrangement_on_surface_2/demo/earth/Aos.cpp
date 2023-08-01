
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
  auto xcv = ctr_cv(ctr_p(0, 0, -1), ctr_p(0, 0, 1), Dir3(0,1,0));

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
  for(auto xcv : xcvs)
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
      if(4 == vit->degree())
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
  std::cout << "num nodes = " << nodes.size()   << std::endl;
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
        CGAL::insert(arr, ctr_cv(p1,p2));
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
        for(int i=0; i< lring->ids.size(); ++i)
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

      face_arcs.push_back( ctr_cv(curr->source()->point(), curr->target()->point()));
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
  using Vertex_handle   = Ext_aos::Vertex_handle;
  using Halfedge_handle = Ext_aos::Halfedge_handle;
  using Face_handle     = Ext_aos::Face_handle;
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
      auto svp = vertex_pos_map[eh->source()];
      auto tvp = vertex_pos_map[eh->target()];
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
  int num_half_edges = 0;
  auto& js_halfedges = js["halfedges"] = json::array();
  std::map<Halfedge_handle, int>  halfedge_pos_map;
  //auto write_half_edges = [&](Ext_aos::Ccb_halfedge_circulator first)
  //{
  //  auto curr = first;
  //  do {
  //    num_half_edges++;
  //    auto& he = *curr;
  //    auto it = halfedge_pos_map.find(&he);
  //    if (it == halfedge_pos_map.end())
  //    {
  //      auto new_he_pos = halfedge_pos_map.size();
  //      halfedge_pos_map[&he] = new_he_pos;
  //    }
  //    auto svh =  he.source(); 
  //    auto tvh =  he.target();
  //    auto& xcv = he.curve();
  //    auto svp = vertex_pos_map[svh];
  //    auto tvp = vertex_pos_map[tvh];
  //    auto xcvp = curve_pos_map[&xcv];
  //    json js_he;
  //    js_he["source"] = svp;
  //    js_he["target"] = tvp;
  //    js_he["curve"]  = xcvp;
  //    js_halfedges.push_back(std::move(js_he));
  //  } while (++curr != first);
  //};
  auto write_half_edge = [&](auto& he)
  {
    auto it = halfedge_pos_map.find(&he);
    if (it == halfedge_pos_map.end())
    {
      auto new_he_pos = halfedge_pos_map.size();
      halfedge_pos_map[&he] = new_he_pos;
    }
    auto svh = he.source();
    auto tvh = he.target();
    auto& xcv = he.curve();
    auto svp = vertex_pos_map[svh];
    auto tvp = vertex_pos_map[tvh];
    auto xcvp = curve_pos_map[&xcv];
    json js_he;
    js_he["source"] = svp;
    js_he["target"] = tvp;
    js_he["curve"] = xcvp;
    js_halfedges.push_back(std::move(js_he));
  };

  for (auto it = arr.halfedges_begin(); it != arr.halfedges_end(); ++it)
  {
    auto& he = *it;
    write_half_edge(he);
    num_half_edges++;
  }
  //for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh)
  //{
  //  auto it = fh->outer_ccb();
  //  for (auto ccb = fh->outer_ccbs_begin(); ccb != fh->outer_ccbs_end(); ++ccb)
  //    write_half_edges(*ccb);  

  //  for (auto ccb = fh->inner_ccbs_begin(); ccb != fh->inner_ccbs_end(); ++ccb)
  //    write_half_edges(*ccb);
  //}
  std::cout << "HALF-EDGE CHECKS:\n";
  std::cout << "  *** num total half-edges = " << num_half_edges << std::endl;
  std::cout << "  *** halfedge-pos-map size = " << halfedge_pos_map.size() << std::endl;


  ////////////////////////////////////////////////////////////////////////////
  // FACES
  // CONDITION DATA: Caspian Sea needs to be defined
  num_half_edges = 0;
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
      num_half_edges++;
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
      json js_halfedges;
      auto& ccb_halfedge_indices = js_halfedges["halfedges"] = json::array();
      auto curr = first;
      do {
        auto& he = *curr;
        auto hep = halfedge_pos_map[&he];
        ccb_halfedge_indices.push_back(hep);
      } while (++curr != first);

      return js_halfedges;
    };
   
  int total_num_half_edges = 0;
  for (auto fh = arr.faces_begin(); fh != arr.faces_end(); ++fh)
  {
    // skip the spherical face
    if (fh->number_of_outer_ccbs() == 0)
      continue;

    // json object for the current face
    json js_face;
    auto face_name = face_name_map[fh];
    js_face["name"] = face_name;

    // at this point we are sure that we have at least 1 outer-ccb
    auto& js_outer_ccbs = js_face["outer_ccbs"] = json::array();
    for (auto ccb = fh->outer_ccbs_begin(); ccb != fh->outer_ccbs_end(); ++ccb)
    {
      auto js_ccb = get_ccb_json(*ccb);
      js_outer_ccbs.push_back(std::move(js_ccb));
    }

    // INNER CCBS
    if(fh->number_of_inner_ccbs() > 0)
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


Aos::Approx_arcs Aos::load_arr(const std::string& file_name)
{
  auto js_txt = read_file(file_name);
  auto js = json::parse(js_txt.begin(), js_txt.end());

  Geom_traits traits;
  Ext_aos arr(&traits);
  auto ctr_p = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();

  ////////////////////////////////////////////////////////////////////////////
  // POINTS
  std::vector<Point>  points;
  auto get_double_from_json = [&](json& js_val)
    {
      auto num = js_val["num"].get<std::string>();
      auto den = js_val["den"].get<std::string>();
      CGAL::Gmpq rat_x(num, den);
      return CGAL::to_double(rat_x);
    };
  auto& js_points = js["points"];
  for (auto it = js_points.begin(); it != js_points.end(); ++it)
  {
    auto& js_point = *it;
    auto loc = js_point["location"].get<int>();
    auto dx = get_double_from_json(js_point["dx"]);
    auto dy = get_double_from_json(js_point["dy"]);
    auto dz = get_double_from_json(js_point["dz"]);
    auto p = ctr_p(dx, dy, dz);
    p.set_location(static_cast<Point::Location_type>(loc));
    points.push_back(p);
  }
  std::cout << "num points = " << points.size() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // CURVES
  std::vector<Curve>  curves;
  auto& js_curves = js["curves"];
  
  for (auto it = js_curves.begin(); it != js_curves.end(); ++it)
  {
    auto& js_curve = *it;
    auto psi = js_curve["source"].get<int>();
    auto pti = js_curve["target"].get<int>();
    auto& js_normal = js_curve["normal"];
    auto nx = get_double_from_json(js_normal["dx"]);
    auto ny = get_double_from_json(js_normal["dy"]);
    auto nz = get_double_from_json(js_normal["dz"]);
    auto is_vertical = js_curve["is_vertical"].get<bool>();
    auto is_directed_right = js_curve["is_directed_right"].get<bool>();
    auto is_full = js_curve["is_full"].get<bool>();
    
    auto xcv = ctr_cv(points[psi], points[pti]);
    //auto xcv = ctr_cv(points[psi], points[pti], Direction_3(nx, ny, nz));
    //xcv.set_is_vertical(is_vertical);
    //xcv.set_is_directed_right(is_directed_right);
    //xcv.set_is_full(is_full);
    curves.push_back(xcv);
  }
  std::cout << "num curves = " << curves.size() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // VERTICES
  auto& js_vertices = js["vertices"];
  using Vertex_handle = Ext_aos::Vertex_handle;
  std::vector<Vertex_handle> vertices;
  for (auto it = js_vertices.begin(); it != js_vertices.end(); ++it)
  {
    auto& js_vertex = *it;
    auto pi = js_vertex["point"].get<int>();
    auto vh = CGAL::insert_point(arr, points[pi]);
    vertices.push_back(vh);
  }
  std::cout << "num vertices = " << vertices.size() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // HALF-EDGES
  struct Halfedge
  {
    int svi, tvi, ci;
  };
  std::vector<Halfedge>  halfedges;
  auto& js_halfedges = js["halfedges"];
  for (auto it = js_halfedges.begin(); it != js_halfedges.end(); ++it)
  {
    auto& js_halfedge = *it;
    Halfedge he;
    he.svi = js_halfedge["source"].get<int>();
    he.tvi = js_halfedge["target"].get<int>();
    he.ci = js_halfedge["curve"].get<int>();
    halfedges.push_back(he);
  }
  std::cout << "halfedges = " << halfedges.size() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // FACES
  auto add_ccbs_to_arr = [&](auto& js_ccbs)
    {
      for (auto cit = js_ccbs.begin(); cit != js_ccbs.end(); ++cit)
      {
        auto& js_outer_ccb = *cit;
        auto& js_halfedges = js_outer_ccb["halfedges"];
        std::cout << "num halfedges = " << js_halfedges.size() << std::endl;
        for (auto hit = js_halfedges.begin(); hit != js_halfedges.end(); ++hit)
        {
          auto& js_halfedge = *hit;
          auto hei = js_halfedge.get<int>();
          auto xcvi = halfedges[hei].ci;
          auto& xcv = curves[xcvi];
          CGAL::insert(arr, xcv);
        }
      }
    };
  auto& js_faces = js["faces"];
  std::cout << "num faces = " << js_faces.size() << "\n";
  for (auto it = js_faces.begin(); it != js_faces.end(); ++it)
  {
    auto& js_face = *it;
    auto& js_name = js_face["name"];
    auto& js_outer_ccbs = js_face["outer_ccbs"];
    //auto& js_inner_ccbs = js_face["inner_ccbs"];
    if(0)
    {
      std::cout << std::boolalpha << "is name string = " << js_name.is_string() << std::endl;
      std::cout << "name = " << js_name.get<std::string>() << std::endl;
      auto& js_ccbs = js_outer_ccbs;
      std::cout << std::boolalpha << "ccb is array = " << js_ccbs.is_array() << "\n";
      std::cout << "num ccbs = " << js_ccbs.size() << std::endl;

    }

    add_ccbs_to_arr(js_outer_ccbs);
    //add_ccbs_to_arr(js_inner_ccbs);
  }
  std::cout << "num arr-faces = " << arr.number_of_faces() << std::endl;

  return Approx_arcs{};
}