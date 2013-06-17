#ifndef CGAL_POINT_INSIDE_POLYHEDRON_H
#define CGAL_POINT_INSIDE_POLYHEDRON_H

#include <CGAL/internal/Point_inside_Polyhedron/Ray_3_Triangle_3_traversal_traits.h>
#include <CGAL/internal/Point_inside_Polyhedron/AABB_triangle_accessor_3_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangle_accessor_3.h>

#include <boost/optional.hpp>

namespace CGAL {

/** 
 * This class provides an efficient point location functionality with respect to a domain bounded
 * by one or several disjoint triangulated closed polyhedral (manifold) surfaces.
 * In case several polyhedral surface are provided as input, a point is said to be inside the domain
 * if an odd number of surfaces is crossed when walking from infinity to the point.
 * The implementation is based on an AABB-tree.
 * @tparam Polyhedron a triangulated polyhedral surface
 * @tparam Kernel a \cgal kernel
 * @tparam TriangleAccessor a model of the concept `TriangleAccessor_3`, with `TriangleAccessor_3::Triangle_3` being `Kernel::Triangle_3`. 
 *         If `Polyhedron` is a \cgal Polyhedron, a default is provided.
 * \todo Doc: move the concept `TriangleAccessor_3` into the "Operation on Polyhedra" package
 * \todo Code: Use this class as an implementation detail of Mesh_3's Polyhedral_mesh_domain_3
 */
template <class Polyhedron, 
          class Kernel,
          class TriangleAccessor_3 = Triangle_accessor_3<Polyhedron, typename Polyhedron::Traits> 
>
class Point_inside_polyhedron_3{
  // typedefs
  typedef CGAL::internal::AABB_triangle_accessor_3_primitive<Kernel, TriangleAccessor_3> Primitive;
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
   * Default constructor. The domain is considered to be empty.
   */
  Point_inside_polyhedron_3(const Kernel& kernel=Kernel())
  : m_kernel(kernel) 
  { }
 
  /** 
   * Constructor with one polyhedral surface. `polyhedron` must be closed and triangulated.
   */
  Point_inside_polyhedron_3(const Polyhedron& polyhedron, const Kernel& kernel=Kernel()) 
  : m_kernel(kernel)
  {
    add_polyhedron(polyhedron);
  }

  /** 
   * Constructor with several polyhedral surfaces. All the polyhedral surfaces must be closed, triangulated and disjoint.
   * \tparam InputIterator is an input iterator with `Polyhedron` or `cpp11::reference_wrapper<Polyhedron>` as value type.
   */
  template <class InputIterator>
  Point_inside_polyhedron_3(InputIterator begin, InputIterator beyond, const Kernel& kernel=Kernel()) 
  : m_kernel(kernel)
  {
    add_polyhedra(begin, beyond);
  }

  /** 
   * Builds internal AABB tree. Optional to call, since the tree is automatically built at the time of first query.
   */ 
  void build() { tree.build(); }

  /**
   * `polyhedron` is added as input
   */
  void add_polyhedron(const Polyhedron& polyhedron) 
  {
    CGAL_assertion(polyhedron.is_pure_triangle());
    CGAL_assertion(polyhedron.is_closed());

    tree.insert(polyhedron.facets_begin(),polyhedron.facets_end());
  }
 
  /**
   * The polyhedral surfaces in the range `[begin,beyond[` are added as input
   * \tparam InputIterator is an input iterator with `Polyhedron` or `cpp11::reference_wrapper<Polyhedron>` as value type.
   */
  template<class InputIterator>
  void add_polyhedra(InputIterator begin, InputIterator beyond) 
  {
    for(; begin != beyond; ++begin) {
      add_polyhedron(*begin);
    }
  }
 
  /**
   * Query function to determine point location.
   * @return 
   *   - CGAL::ON_BOUNDED_SIDE if the point is inside the polyhedron
   *   - CGAL::ON_BOUNDARY if the point is on polyhedron
   *   - CGAL::ON_UNBOUNDED_SIDE if the point is outside polyhedron
   */
  Bounded_side operator()(const Point& point) const
  {
    const Bounding_box& bbox = tree.bbox();

    if(   point.x() < bbox.xmin() || point.x() > bbox.xmax()
       || point.y() < bbox.ymin() || point.y() > bbox.ymax()
       || point.z() < bbox.zmin() || point.z() > bbox.zmax() )
    {
      return ON_UNBOUNDED_SIDE;
    }

    typename Kernel::Construct_ray_3 make_ray = m_kernel.construct_ray_3_object();
    typename Kernel::Construct_vector_3 make_vector = m_kernel.construct_vector_3_object();
     
    CGAL::Random rg(seed); // seed some value for make it easy to debug
    Random_points_on_sphere_3<Point> random_point(1.,rg);
    //the direction of the vertical ray depends on the position of the point in the bbox
    //in order to limit the expected number of nodes visited.

    Ray query = make_ray(point, make_vector(0,0,(2*point.z() <  tree.bbox().zmax()+tree.bbox().zmin()?-1:1)));
    boost::optional<Bounded_side> res = is_inside_ray_tree_traversal<Ray,true>(query);

    while (!res){
      //retry with a random ray
      query = make_ray(point, make_vector(CGAL::ORIGIN,*random_point++));
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
 
};

} // namespace CGAL

#endif //CGAL_POINT_INSIDE_POLYHEDRON_H
