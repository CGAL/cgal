#ifndef CGAL_DEMO_TRIANGULATE_PRIMITIVE_H
#define CGAL_DEMO_TRIANGULATE_PRIMITIVE_H

#include <CGAL/Three/Scene_item.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Projection_traits_3.h>

#include <QColor>

#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <queue>
#include <vector>

// Ensure that all the facets are triangles
// @todo just use PMP::triangulate_face()...?
// or at least mark_faces_in_domain()

template<class Mesh, typename Kernel, typename Index_type>
class FacetTriangulator
{
public:
  using Traits = Kernel;
  using Point = typename Kernel::Point_3;
  using Vector = typename Kernel::Vector_3;

  using P_traits = CGAL::Projection_traits_3<Traits>;

  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

  struct Face_info
  {
    typename boost::graph_traits<Mesh>::halfedge_descriptor e[3];
    bool is_external;
  };

  using Vb = CGAL::Triangulation_vertex_base_with_info_2<halfedge_descriptor, P_traits>;
  using Fbb = CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<P_traits, Fbb>;
  using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Itag = CGAL::Exact_predicates_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>;

  using Vertex_handle = typename CDT::Vertex_handle;
  using Face_handle = typename CDT::Face_handle;

  struct PointAndId
  {
    Point point;
    Index_type id;
    PointAndId() = default;
    PointAndId(const Point& point, const Index_type id) : point(point), id(id) { }
  };

public:
  CDT* cdt;
  CGAL::Unique_hash_map<Vertex_handle, Index_type> v2v;

public:
  // Constructor
  FacetTriangulator(face_descriptor fd,
                    const Vector& normal,
                    Mesh* poly,
                    Vector offset = Vector(0,0,0))
  {
    std::vector<PointAndId> idPoints;
    for(halfedge_descriptor he_circ : halfedges_around_face(halfedge(fd, *poly), *poly))
      idPoints.emplace_back(get(CGAL::vertex_point, *poly, source(he_circ, *poly)) + offset,
                            source(he_circ, *poly));

    if(!triangulate(idPoints, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }

  FacetTriangulator(face_descriptor fd,
                    const std::vector<Point>& more_points,
                    const Vector& normal,
                    Mesh* poly,
                    Vector offset = Vector(0,0,0))
  {
    std::vector<PointAndId> idPoints;
    for(halfedge_descriptor he_circ : halfedges_around_face(halfedge(fd, *poly), *poly))
      idPoints.emplace_back(get(CGAL::vertex_point, *poly, source(he_circ, *poly)) + offset,
                            source(he_circ, *poly));

    if(!triangulate_with_points(idPoints, more_points, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }

  FacetTriangulator(std::vector<PointAndId>& idPoints,
                    const Vector& normal)
  {
    if(!triangulate(idPoints, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }
  FacetTriangulator(std::vector<PointAndId>& idPoints,
                    const std::vector<Point>& more_points,
                    const Vector& normal)
  {
    if(!triangulate_with_points(idPoints, more_points, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }

  ~FacetTriangulator()
  {
    if(cdt)
      delete cdt;
  }

private:
  bool triangulate(std::vector<PointAndId>& idPoints,
                   const Vector& normal)
  {
    P_traits cdt_traits(normal);
    cdt = new CDT(cdt_traits);

    Vertex_handle previous, first, last_inserted;

    // Iterate the points of the facet and decide if they must be inserted in the CDT
    typename Kernel::FT x(0), y(0), z(0);

    for(const PointAndId& idPoint : idPoints)
    {
      x += idPoint.point.x();
      y += idPoint.point.y();
      z += idPoint.point.z();

      Vertex_handle vh;
      // Always insert the first point, then only insert if the distance with the previous is reasonable.
      if(first == Vertex_handle() || idPoint.point != previous->point())
      {
        vh = cdt->insert(idPoint.point);
        v2v[vh] = idPoint.id;
        if(first == Vertex_handle())
          first = vh;

        if(previous != nullptr && previous != vh)
        {
          cdt->insert_constraint(previous, vh);
          last_inserted = previous;
        }
        previous = vh;
      }
    }

    if(last_inserted == Vertex_handle())
      return false;

    if(previous != first)
      cdt->insert_constraint(previous, first);

    // sets mark is_external
    for(Face_handle f2 : cdt->all_face_handles())
      f2->info().is_external = false;

    // check if the facet is external or internal
    std::queue<typename CDT::Face_handle> face_queue;
    face_queue.push(cdt->infinite_vertex()->face());
    while(! face_queue.empty())
    {
      typename CDT::Face_handle fh = face_queue.front();
      face_queue.pop();
      if(fh->info().is_external)
        continue;

      fh->info().is_external = true;
      for(int i = 0; i <3; ++i)
      {
        if(!cdt->is_constrained(std::make_pair(fh, i)))
          face_queue.push(fh->neighbor(i));
      }
    }

    return true;
  }

  bool triangulate_with_points(std::vector<PointAndId >& idPoints,
                               const std::vector<Point>& more_points,
                               const Vector& normal)
  {
    P_traits cdt_traits(normal);
    cdt = new CDT(cdt_traits);

    // Iterate the points of the facet and decide if they must be inserted in the CDT
    Vertex_handle previous, first, last_inserted;
    for(const PointAndId& idPoint : idPoints)
    {
      Vertex_handle vh;
      // Always insert the first point, then only insert if the distance with the previous is reasonable.
      if(first == Vertex_handle() || idPoint.point != previous->point())
      {
        vh = cdt->insert(idPoint.point);
        v2v[vh] = idPoint.id;
        if(first == Vertex_handle())
          first = vh;

        if(previous != nullptr && previous != vh)
        {
          cdt->insert_constraint(previous, vh);
          last_inserted = previous;
        }
        previous = vh;
      }
    }

    if(last_inserted == Vertex_handle())
      return false;

    cdt->insert_constraint(previous, first);
    for(const Point& point : more_points)
      cdt->insert(point);

    // sets mark is_external
    for(Face_handle f2 : cdt->all_face_handles())
      f2->info().is_external = false;

    // check if the facet is external or internal
    std::queue<typename CDT::Face_handle> face_queue;
    face_queue.push(cdt->infinite_vertex()->face());
    while(!face_queue.empty())
    {
      typename CDT::Face_handle fh = face_queue.front();
      face_queue.pop();
      if(fh->info().is_external)
        continue;
      fh->info().is_external = true;
      for(int i = 0; i <3; ++i)
      {
        if(!cdt->is_constrained(std::make_pair(fh, i)))
          face_queue.push(fh->neighbor(i));
      }
    }

    return true;
  }
};

#endif // CGAL_DEMO_TRIANGULATE_PRIMITIVE_H

