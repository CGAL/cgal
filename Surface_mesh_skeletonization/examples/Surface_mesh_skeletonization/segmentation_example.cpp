#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::face_descriptor             face_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>        Skeletonization;
typedef Skeletonization::Skeleton                                    Skeleton;

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


int main(int argc, char* argv[])
{
  std::ifstream input((argc>1)?argv[1]:"data/161.off");
  Polyhedron tmesh;
  input >> tmesh;

  // extract the skeleton
  Skeleton skeleton;
  CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);

  // init the polyhedron simplex indices
  CGAL::set_halfedgeds_items_id(tmesh);

  //for each input vertex compute its distance to the skeleton
  std::vector<double> distances(num_vertices(tmesh));
  BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton) )
  {
    const Point& skel_pt = skeleton[v].point;
    BOOST_FOREACH(vertex_descriptor mesh_v, skeleton[v].vertices)
    {
      const Point& mesh_pt = mesh_v->point();
      distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
    }
  }

  // create a property-map for sdf values
  std::vector<double> sdf_values( num_faces(tmesh) );
  Facet_with_id_pmap<double> sdf_property_map(sdf_values);

  // compute sdf values with skeleton
  BOOST_FOREACH(face_descriptor f, faces(tmesh))
  {
    double dist = 0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(f, tmesh), tmesh))
      dist+=distances[target(hd, tmesh)->id()];
    sdf_property_map[f] = dist / 3.;
  }

  // post-process the sdf values
  CGAL::sdf_values_postprocessing(tmesh, sdf_property_map);

  // create a property-map for segment-ids (it is an adaptor for this case)
  std::vector<std::size_t> segment_ids( num_faces(tmesh) );
  Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);

  // segment the mesh using default parameters
  std::cout << "Number of segments: "
            << CGAL::segmentation_from_sdf_values(tmesh, sdf_property_map, segment_property_map) <<"\n";

  return 0;
}

