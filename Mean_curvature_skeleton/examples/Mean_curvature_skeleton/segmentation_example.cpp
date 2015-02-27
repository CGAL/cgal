#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>

#include <fstream>
#include <map>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor             face_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>        Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton                            Skeleton;

typedef Skeleton::vertex_descriptor                                  Skeleton_vertex;

// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template<class ValueType>
struct Facet_with_id_pmap
    : public boost::put_get_helper<ValueType&,
             Facet_with_id_pmap<ValueType> >
{
    typedef face_descriptor key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;

    Facet_with_id_pmap(
      std::vector<ValueType>& internal_vector
    ) : internal_vector(internal_vector) { }

    reference operator[](key_type key) const
    { return internal_vector[key->id()]; }
private:
    std::vector<ValueType>& internal_vector;
};


int main()
{
  Polyhedron tmesh;
  std::ifstream input("data/161.off");

  if ( !input || !(input >> tmesh) || tmesh.empty() ) {
    std::cerr << "Cannot open data/161.off" << std::endl;
    return 1;
  }

  // extract the skeleton
  Skeleton skeleton;
  CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);

  // init the polyhedron simplex indices
  CGAL::set_halfedgeds_items_id(tmesh)

  //for each input vertex compute its distance to the skeleton
  std::vector<double> distances(num_vertices(tmesh));
  BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton) )
  {
    const Point& skel = skeleton[gv].point;
    BOOST_FOREACH(vertex_descriptor mesh_v, vertices(tmesh))
    {
      const Point& mesh_point = skeleton[gv].vertices[i]->point();
      distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel, surf));
    }
  }

  // create a property-map for sdf values
  std::vector<double> sdf_values;
  Facet_with_id_pmap<double> sdf_property_map(sdf_values);

  // compute sdf values with skeleton
  BOOST_FOREACH(face_descriptor f, faces(tmesh))
  {
    double dist = 0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(f, tmesh))
      dist+=distances[target(hd, tmesh)->id()];
    sdf_property_map[f] = dist / 3.;
  }

  // post-process the sdf values
  CGAL::sdf_values_postprocessing(tmesh, sdf_property_map);

  // create a property-map for segment-ids (it is an adaptor for this case)
  std::vector<int> segment_ids;
  Facet_with_id_pmap<double> segment_property_map(segment_ids);

  // segment the mesh using default parameters
  std::cout << "Number of segments: "
            << CGAL::segmentation_from_sdf_values(tmesh, sdf_property_map, segment_property_map) <<"\n";

  return 0;
}

