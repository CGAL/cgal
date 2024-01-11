#ifndef CGAL_DEMO_TRIANGULATE_PRIMITIVE_H
#define CGAL_DEMO_TRIANGULATE_PRIMITIVE_H

#include <CGAL/Three/Scene_item.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

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

template<class Mesh, typename Kernel, typename Index_type, typename CDT_ITag = CGAL::Exact_predicates_tag>
class FacetTriangulator
{
public:
  using Traits = Kernel;
  using Point = typename Kernel::Point_3;
  using Vector = typename Kernel::Vector_3;

  using Index = Index_type;

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
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, CDT_ITag>;

  using Vertex_handle = typename CDT::Vertex_handle;
  using Face_handle = typename CDT::Face_handle;

public:
  CDT* cdt;
  CGAL::Unique_hash_map<Vertex_handle, Index> v2v;

public:
  // Constructor
  FacetTriangulator(face_descriptor fd,
                    const Vector& normal,
                    Mesh* poly,
                    Vector offset = Vector(0,0,0))
  {
    std::vector<std::pair<typename Kernel::Point_3, Index> > idPoints;
    for(halfedge_descriptor he_circ : halfedges_around_face(halfedge(fd, *poly), *poly))
      idPoints.emplace_back(std::make_pair(get(CGAL::vertex_point, *poly, source(he_circ, *poly)) + offset,
                            source(he_circ, *poly)));

    if(!triangulate(idPoints, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }

  FacetTriangulator(face_descriptor fd,
                    const std::vector<Point>& more_points,
                    const Vector& normal,
                    Mesh* poly,
                    Vector offset = Vector(0,0,0))
  {
    std::vector<std::pair<typename Kernel::Point_3, Index> > idPoints;
    for(halfedge_descriptor he_circ : halfedges_around_face(halfedge(fd, *poly), *poly))
      idPoints.emplace_back(std::make_pair(get(CGAL::vertex_point, *poly, source(he_circ, *poly)) + offset,
                            source(he_circ, *poly)));

    if(!triangulate_with_points(idPoints, more_points, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }

  FacetTriangulator(std::vector<std::pair<typename Kernel::Point_3, Index> >& idPoints,
                    const Vector& normal)
  {
    if(!triangulate(idPoints, normal))
      std::cerr << "Facet not displayed" << std::endl;
  }
  FacetTriangulator(std::vector < std::pair<typename Kernel::Point_3, Index> >& idPoints,
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

  template<class Func>
  void per_face(const Func& f) {
    for (typename CDT::Finite_faces_iterator fit = cdt->finite_faces_begin(); fit != cdt->finite_faces_end(); ++fit)
      f(*fit, v2v);
  }

private:
  bool triangulate(std::vector<std::pair<typename Kernel::Point_3, Index> >& idPoints,
                   const Vector& normal)
  {
    P_traits cdt_traits(normal);
    cdt = new CDT(cdt_traits);

    std::map<std::pair<Point, Point>, std::size_t> edge_map;
    std::vector<bool> skip(idPoints.size(), false);

    for (std::size_t i = 0; i < idPoints.size(); i++) {
      std::size_t prev = (i - 1 + idPoints.size()) % idPoints.size();
      if (idPoints[i].first < idPoints[prev].first) {
        auto it = edge_map.emplace(std::make_pair(idPoints[i].first, idPoints[prev].first), i);
        if (!it.second) {
          skip[i] = true;
          skip[it.first->second] = true;
        }
      }
      else {
        auto it = edge_map.emplace(std::make_pair(idPoints[prev].first, idPoints[i].first), i);
        if (!it.second) {
          skip[i] = true;
          skip[it.first->second] = true;
        }
      }
    }

    Vertex_handle previous, first, last_inserted;

    // Iterate the points of the facet and decide if they must be inserted in the CDT
    typename Kernel::FT x(0), y(0), z(0);

    for(std::size_t i = 0;i<idPoints.size();i++)
    {
      const std::pair<typename Kernel::Point_3, Index>& idPoint = idPoints[i];
      x += idPoint.first.x();
      y += idPoint.first.y();
      z += idPoint.first.z();

      Vertex_handle vh;
      // Always insert the first point, then only insert if the distance with the previous is reasonable.
      if(first == Vertex_handle() || idPoint.first != previous->point())
      {
        vh = cdt->insert(idPoint.first);
        v2v[vh] = idPoint.second;
        if (first == Vertex_handle()) {
          first = vh;
        }

        if(previous != nullptr && previous != vh)
        {
          if (!skip[i])
            cdt->insert_constraint(previous, vh);
          last_inserted = previous;
        }
        previous = vh;
      }
    }

    if(last_inserted == Vertex_handle())
      return false;

    if(previous != first && !skip[0])
      cdt->insert_constraint(previous, first);

    // sets mark is_external
    for(Face_handle f2 : cdt->all_face_handles())
      f2->info().is_external = false;

    std::unordered_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<Face_handle, bool> > in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation<CDT>(*cdt, in_domain);

    for (const Face_handle& fh : cdt->all_face_handles())
      fh->info().is_external = !get(in_domain, fh);

    return true;
  }

  bool triangulate_with_points(std::vector<std::pair<typename Kernel::Point_3, Index> >& idPoints,
                               const std::vector<Point>& more_points,
                               const Vector& normal)
  {
    P_traits cdt_traits(normal);
    cdt = new CDT(cdt_traits);

    std::map<std::pair<Point, Point>, std::size_t> edge_map;
    std::vector<bool> skip(idPoints.size(), false);

    for (std::size_t i = 0; i < idPoints.size(); i++) {
      std::size_t prev = (i - 1 + idPoints.size()) % idPoints.size();
      if (idPoints[i].first < idPoints[prev].first) {
        auto it = edge_map.emplace(std::make_pair(idPoints[i].first, idPoints[prev].first), i);
        if (!it.second) {
          skip[i] = true;
          skip[it.first->second] = true;
        }
      }
      else {
        auto it = edge_map.emplace(std::make_pair(idPoints[prev].first, idPoints[i].first), i);
        if (!it.second) {
          skip[i] = true;
          skip[it.first->second] = true;
        }
      }
    }

    // Iterate the points of the facet and decide if they must be inserted in the CDT
    Vertex_handle previous, first, last_inserted;
    for (std::size_t i = 0; i < idPoints.size(); i++)
    {
      const std::pair<typename Kernel::Point_3, Index>& idPoint = idPoints[i];
      Vertex_handle vh;
      // Always insert the first point, then only insert if the distance with the previous is reasonable.
      if(first == Vertex_handle() || idPoint.first != previous->point())
      {
        vh = cdt->insert(idPoint.first);
        v2v[vh] = idPoint.second;
        if(first == Vertex_handle())
          first = vh;

        if(previous != nullptr && previous != vh)
        {
          if (!skip[i])
            cdt->insert_constraint(previous, vh);
          last_inserted = previous;
        }
        previous = vh;
      }
    }

    if(last_inserted == Vertex_handle())
      return false;

    if (previous != first && !skip[0])
    cdt->insert_constraint(previous, first);

    for(const Point& point : more_points)
      cdt->insert(point);

    // sets mark is_external
    for(Face_handle f2 : cdt->all_face_handles())
      f2->info().is_external = false;

    std::unordered_map<Face_handle, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<Face_handle, bool> > in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation<CDT>(*cdt, in_domain);

    for (const Face_handle& fh : cdt->all_face_handles())
      fh->info().is_external = !get(in_domain, fh);

    return true;
  }
};

#endif // CGAL_DEMO_TRIANGULATE_PRIMITIVE_H

