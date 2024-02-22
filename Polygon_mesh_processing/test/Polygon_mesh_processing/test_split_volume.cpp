#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/core/ref.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Aff_transformation_3 Trsfrm;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

int main()
{
  Surface_mesh base_cube;
  std::ifstream input("data-coref/cube.off");
  input >> base_cube;

  Surface_mesh input_mesh;
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> vol_id_map =
    input_mesh.add_property_map<Surface_mesh::Face_index, std::size_t>().first;

  // add large cube
  input_mesh.join(base_cube);

  // add middle cube
  Surface_mesh middle_cube = base_cube;
  PMP::transform(Trsfrm(CGAL::SCALING, 0.5), middle_cube);
  PMP::transform(Trsfrm(CGAL::TRANSLATION, Kernel::Vector_3(0.25,0.25,0.25)), middle_cube);
  input_mesh.join(middle_cube);

  // add central cube
  Surface_mesh small_cube = base_cube;
  PMP::transform(Trsfrm(CGAL::SCALING, 0.1), small_cube);
  PMP::transform(Trsfrm(CGAL::TRANSLATION, Kernel::Vector_3(0.5,0.5,0.5)), small_cube);
  input_mesh.join(small_cube);

  // test using orientation
  std::vector<PMP::Volume_error_code> expected_result;
  expected_result.push_back(PMP::VALID_VOLUME);
  expected_result.push_back(PMP::INCOMPATIBLE_ORIENTATION);
  expected_result.push_back(PMP::INCOMPATIBLE_ORIENTATION);
  std::vector<PMP::Volume_error_code>  error_codes;
  std::size_t nb_vol = PMP::volume_connected_components(input_mesh, vol_id_map,
                                                        params::error_codes(boost::ref(error_codes)));
  assert(nb_vol==3);
  std::sort(error_codes.begin(), error_codes.end());
  assert( error_codes==expected_result );

  // test ignoring orientation
  nb_vol = PMP::volume_connected_components(input_mesh, vol_id_map,
                                            params::do_orientation_tests(false)
                                            .error_codes(boost::ref(error_codes)));
  assert(nb_vol==2);
  expected_result.clear();
  expected_result.resize(2, PMP::VALID_VOLUME);
  assert( error_codes==expected_result );


  // test surface component self-intersection
  {
  Surface_mesh tmp = input_mesh;
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> tmp_vol_id_map =
    tmp.add_property_map<Surface_mesh::Face_index, std::size_t>().first;
  // create a self-intersection
  Surface_mesh::Halfedge_index h = *tmp.halfedges().begin();
  h = CGAL::Euler::split_edge(h, tmp);
  tmp.point( target(h, tmp) ) = tmp.point( target( tmp.next(h), tmp) );
  CGAL::Euler::split_face(h, next( next(h, tmp), tmp), tmp);
  h = opposite(h, tmp);
  CGAL::Euler::split_face(h, next( next(h, tmp), tmp), tmp);

  nb_vol = PMP::volume_connected_components(tmp, tmp_vol_id_map,
                                            params::do_orientation_tests(false)
                                            .error_codes(boost::ref(error_codes))
                                            .do_self_intersection_tests(true));

  assert(nb_vol==2);
  expected_result.clear();
  expected_result.push_back(PMP::VALID_VOLUME);
  expected_result.push_back(PMP::SURFACE_WITH_SELF_INTERSECTIONS);
  std::sort(error_codes.begin(), error_codes.end());
  assert( error_codes==expected_result );
  }

  // test single cc
  {
  Surface_mesh tmp = base_cube;
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> tmp_vol_id_map =
    tmp.add_property_map<Surface_mesh::Face_index, std::size_t>().first;

  nb_vol = PMP::volume_connected_components(tmp, tmp_vol_id_map,
                                            params::do_orientation_tests(false)
                                            .error_codes(boost::ref(error_codes))
                                            .do_self_intersection_tests(true));

  assert(nb_vol==1);
  expected_result.clear();
  expected_result.push_back(PMP::VALID_VOLUME);
  assert( error_codes==expected_result );
  }

  // test intersection between cc
  {
  Surface_mesh tmp = input_mesh;
  // create surface intersection
  tmp.join(base_cube);
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> tmp_vol_id_map =
    tmp.add_property_map<Surface_mesh::Face_index, std::size_t>().first;

  nb_vol = PMP::volume_connected_components(tmp, tmp_vol_id_map,
                                            params::do_orientation_tests(false)
                                            .error_codes(boost::ref(error_codes))
                                            .do_self_intersection_tests(true));

  assert(nb_vol==4);
  expected_result.clear();
  expected_result.resize(4, PMP::VOLUME_INTERSECTION);
  assert( error_codes==expected_result );
  }
  {
  Surface_mesh tmp = input_mesh;
  // create surface intersection
  tmp.join(middle_cube);
  tmp.join(middle_cube);
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> fccmap =
    tmp.add_property_map<Surface_mesh::Face_index, std::size_t>("f:CC").first;

  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> tmp_vol_id_map =
    tmp.add_property_map<Surface_mesh::Face_index, std::size_t>().first;
  std::vector< std::vector<std::size_t> > nested_cc_per_cc;

  nb_vol = PMP::volume_connected_components(tmp, tmp_vol_id_map,
                                            params::do_orientation_tests(false)
                                            .error_codes(boost::ref(error_codes))
                                            .face_connected_component_map(fccmap)
                                            .do_self_intersection_tests(true)
                                            .volume_inclusions(boost::ref(nested_cc_per_cc)));

  assert(nb_vol==5);
  expected_result.clear();
  expected_result.push_back(PMP::VALID_VOLUME);
  expected_result.resize(5, PMP::VOLUME_INTERSECTION);
  std::sort(error_codes.begin(), error_codes.end());
  assert( error_codes==expected_result );
  std::size_t sum=0;
  for (int i=0; i<5; ++i)
  {
    sum+=nested_cc_per_cc[i].size();
    assert(nested_cc_per_cc[i].size()==1 ||
           nested_cc_per_cc[i].size()==3 ||
           nested_cc_per_cc[i].size()==0);
  }
  assert(sum == 0+3*1+3);
  }

  // test level 0 orientation
  {
  Surface_mesh tmp = base_cube;
  PMP::transform(Trsfrm(CGAL::TRANSLATION, Kernel::Vector_3(2,0,0)), tmp);
  tmp.join(base_cube);
  PMP::reverse_face_orientations(tmp);
  Surface_mesh::Property_map<Surface_mesh::Face_index, std::size_t> tmp_vol_id_map =
    tmp.add_property_map<Surface_mesh::Face_index, std::size_t>().first;
  // create a self-intersection
  nb_vol = PMP::volume_connected_components(tmp, tmp_vol_id_map,
                                            params::do_orientation_tests(false));

  assert(nb_vol==2);
  nb_vol = PMP::volume_connected_components(tmp, tmp_vol_id_map,
                                            params::do_orientation_tests(true));

  assert(nb_vol==1);
  }
  // debug code
  /*
    std::cout << "  found " << nb_vol << " volumes\n";
    typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
    Filtered_graph vol_mesh(sm, 0, vol_id_map);
    for(std::size_t id = 0; id < nb_vol; ++id)
    {
      if(id > 0)
        vol_mesh.set_selected_faces(id, vol_id_map);
      Surface_mesh out;
      CGAL::copy_face_graph(vol_mesh, out);
      std::ostringstream oss;
      oss << "vol_" << id <<".off";
      std::ofstream os(oss.str().data());
      os << out;
    }
  */
}
