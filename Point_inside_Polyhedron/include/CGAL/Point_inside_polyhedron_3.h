#ifndef CGAL_POINT_INSIDE_POLYHEDRON_H
#define CGAL_POINT_INSIDE_POLYHEDRON_H

#include <CGAL/internal/Point_inside_Polyhedron/Ray_3_Triangle_3_traversal_traits.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/point_generators_3.h>

#include <boost/optional.hpp>

namespace CGAL {

/** 
 * Class providing the functionality of locating query points with respect to closed triangulated polyhedron.
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam Kernel a %CGAL kernel
 */
template <class Polyhedron, class Kernel>
class Point_inside_polyhedron_3 {
// typedefs
  typedef CGAL::AABB_const_polyhedron_triangle_primitive<Kernel,Polyhedron> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef typename Traits::Bounding_box Bounding_box;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Ray_3 Ray;
//members
  Kernel m_kernel;
  Tree tree;

  const static unsigned int seed = 1340818006;
public:
  /** 
   * The constructor of query object. It also processes given polyhedron for multiple queries.
   * @pre @a polyhedron.is_pure_triangle()
   * @pre @a polyhedron.is_closed()
   * @param polyhedron polyhedron on with queries are executed
   * @param kernel kernel object
   */
  Point_inside_polyhedron_3(const Polyhedron& polyhedron, const Kernel& kernel=Kernel())
    : m_kernel(kernel)
  {
    CGAL_assertion(polyhedron.is_pure_triangle());
    CGAL_assertion(polyhedron.is_closed());

    tree.insert(polyhedron.facets_begin(),polyhedron.facets_end());
    tree.build();
  }

  /**
   * Query function to determine point location.
   * @return 
   *   - CGAL::ON_BOUNDED_SIDE if the point is inside the polyhedron
   *   - CGAL::ON_BOUNDARY if the point is on polyhedron
   *   - CGAL::ON_UNBOUNDED_SIDE if the point is outside polyhedron
   */
  Bounded_side operator()(const Point& p) const
  {
    const Bounding_box& bbox = tree.bbox();

    if(   p.x() < bbox.xmin() || p.x() > bbox.xmax()
       || p.y() < bbox.ymin() || p.y() > bbox.ymax()
       || p.z() < bbox.zmin() || p.z() > bbox.zmax() )
    {
      return ON_UNBOUNDED_SIDE;
    }

    typename Kernel::Construct_ray_3 make_ray = m_kernel.construct_ray_3_object();
    typename Kernel::Construct_vector_3 make_vector = m_kernel.construct_vector_3_object();
     
    CGAL::Random rg(seed); // seed some value for make it easy to debug
    Random_points_on_sphere_3<Point> random_point(1.,rg);
    //the direction of the vertical ray depends on the position of the point in the bbox
    //in order to limit the expected number of nodes visited.

#define CAST_TO_CLOSEST
#ifdef CAST_TO_CLOSEST
    double dist[6];
    dist[0] = std::pow(p.x() - tree.bbox().xmin(), 2);
    dist[1] = std::pow(p.x() - tree.bbox().xmax(), 2);
    dist[2] = std::pow(p.y() - tree.bbox().ymin(), 2);
    dist[3] = std::pow(p.y() - tree.bbox().ymax(), 2);
    dist[4] = std::pow(p.z() - tree.bbox().zmin(), 2);
    dist[5] = std::pow(p.z() - tree.bbox().zmax(), 2);
    
    int min_i = std::min_element(&dist[0], &dist[6]) - &dist[0];

    Ray query = make_ray(p, 
      make_vector( min_i == 0 ? -1 : min_i == 1 ? 1 : 0,
                   min_i == 2 ? -1 : min_i == 3 ? 1 : 0,
                   min_i == 4 ? -1 : min_i == 5 ? 1 : 0));
    boost::optional<Bounded_side> res;
    if(min_i == 0 || min_i == 1) {
      res = is_inside_ray_tree_traversal_with_axis<Ray, 0>(query);
    } else if(min_i == 2 || min_i == 3) {
      res = is_inside_ray_tree_traversal_with_axis<Ray, 1>(query);
    } else {
      res = is_inside_ray_tree_traversal_with_axis<Ray, 2>(query);
    }
#else
    Ray query = make_ray(p, make_vector(0,0,(2*p.z() <  tree.bbox().zmax()+tree.bbox().zmin()?-1:1)));
    boost::optional<Bounded_side> res = is_inside_ray_tree_traversal<Ray,true>(query);
#endif

    while (!res){
      //retry with a random ray
      query = make_ray(p, make_vector(CGAL::ORIGIN,*random_point++));
      res = is_inside_ray_tree_traversal<Ray,false>(query);
    }
    return *res;
  }

private:
  template <class Query,bool ray_is_vertical>
  boost::optional<Bounded_side>
  is_inside_ray_tree_traversal(const Query& query) const 
  {
    std::pair<boost::logic::tribool,std::size_t> status( boost::logic::tribool(boost::logic::indeterminate), 0);
    internal::Ray_3_Triangle_3_traversal_traits<Traits,Polyhedron,Kernel,Boolean_tag<ray_is_vertical>> traversal_traits(status);
    tree.traversal(query, traversal_traits);

    if ( !boost::logic::indeterminate(status.first) )
    {
      if (status.first) {
        return (status.second&1) == 1 ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
      }
      //otherwise the point is on the facet
      return ON_BOUNDARY;
    }
    return boost::optional<Bounded_side>(); // indeterminate
  }

  // Going to be removed 
  template <class Query, unsigned int Axis>
  boost::optional<Bounded_side>
  is_inside_ray_tree_traversal_with_axis(const Query& query) const {
    std::pair<boost::logic::tribool,std::size_t> status( boost::logic::tribool(boost::logic::indeterminate), 0);
    internal::Ray_3_Triangle_3_traversal_traits_axis<Traits,Polyhedron,Kernel, 
      Axis / 2
    > traversal_traits(status);
    tree.traversal(query, traversal_traits);

    if ( !boost::logic::indeterminate(status.first) )
    {
      if (status.first) {
        return (status.second&1) == 1 ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
      }
      //otherwise the point is on the facet
      return ON_BOUNDARY;
    }
    return boost::optional<Bounded_side>(); // indeterminate
  }
};

} // namespace CGAL

#endif //CGAL_POINT_INSIDE_POLYHEDRON_H
