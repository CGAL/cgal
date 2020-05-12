#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include "Kernel_type.h"
#include "SMesh_type.h"
#include "Scene_points_with_normal_item.h"

typedef std::array<std::size_t,3> Facet;
typedef CGAL::Point_set_with_structure<Kernel> Structuring;
typedef Kernel::Plane_3 Plane_3;

struct Construct{
  typedef std::array<std::size_t,3> Facet;
  typedef typename boost::property_map<SMesh, boost::vertex_point_t>::type VPmap;
  SMesh& mesh;
  VPmap vpmap;
  std::vector<typename boost::graph_traits<SMesh>::vertex_descriptor> vertices;

  template <typename PointRange>
  Construct(SMesh& mesh, const PointRange& points)
    : mesh(mesh)
  {
    vpmap = get(boost::vertex_point, mesh);
    for (const auto& p : points)
    {
      typename boost::graph_traits<SMesh>::vertex_descriptor v;
      v = add_vertex(mesh);
      vertices.push_back(v);
      put(vpmap, v,  p);
    }
  }

  Construct& operator=(const Facet f)
  {
    typedef typename boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
    std::vector<vertex_descriptor> facet;
    facet.resize(3);
    facet[0]=vertices[f[0]];
    facet[1]=vertices[f[1]];
    facet[2]=vertices[f[2]];
    CGAL::Euler::add_face(facet, mesh);
    return *this;
  }
  Construct&
  operator*() { return *this; }
  Construct&
  operator++() { return *this; }
  Construct
  operator++(int) { return *this; }
};

struct Radius {

  double bound;

  Radius(double bound)
    : bound(bound)
  {}

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (const AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return adv.smallest_radius_delaunay_sphere (c, index);
    }

    // If radius > bound, return infinity so that facet is not used
    double d  = 0;
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                              c->vertex((index+2)%4)->point()));
    if(d>bound) return adv.infinity();
    d = sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();

    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere (c, index);
  }
};


template <typename Structuring>
struct Priority_with_structure_coherence {

  Structuring& structuring;
  double bound;

  Priority_with_structure_coherence(Structuring& structuring,
                                    double bound)
    : structuring (structuring), bound (bound)
  {}

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // If perimeter > bound, return infinity so that facet is not used
    if (bound != 0)
      {
        double d  = 0;
        d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                  c->vertex((index+2)%4)->point()));
        if(d>bound) return adv.infinity();
        d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                                   c->vertex((index+3)%4)->point()));
        if(d>bound) return adv.infinity();
        d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                   c->vertex((index+3)%4)->point()));
        if(d>bound) return adv.infinity();
      }

    Facet f = {{ c->vertex ((index + 1) % 4)->info (),
                 c->vertex ((index + 2) % 4)->info (),
                 c->vertex ((index + 3) % 4)->info () }};

    double weight = 100. * (5 - structuring.facet_coherence (f));

    return weight * adv.smallest_radius_delaunay_sphere (c, index);
  }

};

void get_planes_from_shape_map (const Point_set& points,
                                Point_set::Property_map<int> shape_map,
                                std::vector<Plane_3>& planes)
{
  std::vector<Point_set::Index> sorted_indices;
  sorted_indices.reserve (points.size());
  std::copy (points.begin(), points.end(), std::back_inserter (sorted_indices));

  std::sort (sorted_indices.begin(), sorted_indices.end(),
             [&](const Point_set::Index& a, const Point_set::Index& b) -> bool
             {
               return shape_map[a] < shape_map[b];
             });

  std::size_t nb_planes = shape_map[sorted_indices.back()] + 1;
  planes.reserve (nb_planes);
  std::cerr << nb_planes << " found in the shape map" << std::endl;

  std::vector<Point_set::Index>::iterator begin = sorted_indices.end();
  int plane_idx = shape_map[sorted_indices.front()];
  if (plane_idx != -1)
    begin = sorted_indices.begin();

  for (std::vector<Point_set::Index>::iterator it = sorted_indices.begin();
       it != sorted_indices.end(); ++ it)
  {
    if (shape_map[*it] != plane_idx)
    {
      if (begin != sorted_indices.end())
      {
        Plane_3 plane;
        CGAL::linear_least_squares_fitting_3
          (boost::make_transform_iterator
           (begin, CGAL::Property_map_to_unary_function<Point_set::Point_map>(points.point_map())),
           boost::make_transform_iterator
           (it, CGAL::Property_map_to_unary_function<Point_set::Point_map>(points.point_map())),
           plane, CGAL::Dimension_tag<0>());
        planes.push_back (plane);
      }
      begin = it;
      plane_idx = shape_map[*it];
    }
  }
  Plane_3 plane;
  CGAL::linear_least_squares_fitting_3
    (boost::make_transform_iterator
     (begin, CGAL::Property_map_to_unary_function<Point_set::Point_map>(points.point_map())),
     boost::make_transform_iterator
     (sorted_indices.end(), CGAL::Property_map_to_unary_function<Point_set::Point_map>(points.point_map())),
     plane, CGAL::Dimension_tag<0>());
  planes.push_back (plane);
}

SMesh* advancing_front (const Point_set& points,
                        double longest_edge,
                        double radius_ratio_bound,
                        double beta,
                        bool structuring,
                        double sampling)
{
  SMesh* mesh = new SMesh;

  if (structuring) // todo
  {
    Point_set::Property_map<int> shape_map
      = points.property_map<int>("shape").first;

    typedef CGAL::Point_set_with_structure<Kernel> Structuring;
    std::vector<Plane_3> planes;

    get_planes_from_shape_map (points, shape_map, planes);
    Structuring structuring
      (points, planes,
       sampling,
       points.parameters().
       plane_map(CGAL::Identity_property_map<Plane_3>()).
       plane_index_map(shape_map));

    std::vector<Point_3> structured;
    structured.reserve (structuring.size());
    for (std::size_t i = 0; i < structuring.size(); ++ i)
      structured.push_back (structuring.point(i));
    std::cerr << structured.size() << " structured point(s) generated" << std::endl;

    Priority_with_structure_coherence<Structuring> priority
      (structuring, 10. * sampling);
    Construct construct(*mesh, structured);
    CGAL::advancing_front_surface_reconstruction(structured.begin(),
                                                 structured.end(),
                                                 construct,
                                                 priority,
                                                 radius_ratio_bound,
                                                 beta);
  }
  else
  {
    Radius filter(longest_edge);
    Construct construct (*mesh, points.points());
    CGAL::advancing_front_surface_reconstruction(points.points().begin(),
                                                 points.points().end(),
                                                 construct,
                                                 filter,
                                                 radius_ratio_bound,
                                                 beta);
  }

  return mesh;
}
