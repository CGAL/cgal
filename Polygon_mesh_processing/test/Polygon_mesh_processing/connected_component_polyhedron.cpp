#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <iostream>
#include <fstream>
#include <cstring>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Mesh;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Mesh_with_id;

void mesh_with_id(const std::string argv1, const bool save_output)
{
  typedef boost::graph_traits<Mesh_with_id>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<Mesh_with_id>::face_descriptor face_descriptor;

  Mesh_with_id sm;
  std::ifstream in(argv1);
  in >> sm;

  int i=0;
  for(face_descriptor f : faces(sm))
    f->id() = i++;

  i=0;
  for(vertex_descriptor v : vertices(sm))
    v->id() = i++;

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  PMP::connected_component(fd, sm, std::back_inserter(cc));

  std::cerr << cc.size() << " faces in the CC of " << &*fd << std::endl;

  boost::vector_property_map<int,
    boost::property_map<Mesh_with_id, CGAL::face_index_t>::type>
      fccmap(static_cast<unsigned>(num_faces(sm)), get(CGAL::face_index,sm));

  const std::size_t num = PMP::connected_components(sm, fccmap);

  if(argv1 == CGAL::data_file_path("meshes/blobby_3cc.off"))
  {
    assert(num == 3);
  }

  std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
  const std::size_t nb_faces = num_faces(sm);

  std::vector<face_descriptor> faces_to_remove;
  std::size_t nb_to_remove = PMP::keep_large_connected_components(
                               sm, 1000,
                               CGAL::parameters::face_size_map(CGAL::Constant_property_map<face_descriptor, std::size_t>(1))
                                                .dry_run(true)
                                                .output_iterator(std::back_inserter(faces_to_remove)));

  if (argv1 == CGAL::data_file_path("meshes/blobby_3cc.off"))
  {
    assert(nb_to_remove == 1);
    assert(faces_to_remove.size() == 680);
    assert(num_faces(sm) == nb_faces);
  }

  faces_to_remove.clear();
  nb_to_remove = PMP::keep_largest_connected_components(
                   sm, 2,
                   CGAL::parameters::face_size_map(CGAL::Constant_property_map<face_descriptor, std::size_t>(1))
                                    .dry_run(true)
                                    .output_iterator(std::back_inserter(faces_to_remove)));

  if (argv1 == CGAL::data_file_path("meshes/blobby_3cc.off"))
  {
    assert(nb_to_remove == 1);
    assert(faces_to_remove.size() == 680);
    assert(num_faces(sm) == nb_faces);
  }

  nb_to_remove = PMP::keep_largest_connected_components(
                   sm, 2,
                   CGAL::parameters::face_size_map(CGAL::Constant_property_map<face_descriptor, std::size_t>(1)));

  if (argv1 == CGAL::data_file_path("meshes/blobby_3cc.off"))
  {
    assert(nb_to_remove == 1);
    assert(faces(sm).size() == 2737);
  }

  if (!save_output)
    return;

  std::ofstream ofile("blobby_2cc_id.off");
  ofile << sm << std::endl;
  ofile.close();
}

void mesh_no_id(const std::string argv1, const bool save_output)
{
  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

  Mesh sm;
  std::ifstream in(argv1);
  in >> sm;

  std::vector<face_descriptor> cc;
  face_descriptor fd = *faces(sm).first;
  PMP::connected_component(fd, sm, std::back_inserter(cc));

  std::cerr << cc.size() << " faces in the CC of " << &*fd << std::endl;
  boost::property_map<Mesh,boost::vertex_external_index_t>::type vim
    = get(boost::vertex_external_index,sm);

  boost::property_map<Mesh,boost::face_external_index_t>::type fim
    = get(boost::face_external_index,sm);

  boost::vector_property_map<int,
    boost::property_map<Mesh, boost::face_external_index_t>::type>
      fccmap(static_cast<unsigned>(num_faces(sm)), fim);

  std::size_t num = PMP::connected_components(sm, fccmap);

  if (argv1 == CGAL::data_file_path("meshes/blobby_3cc.off"))
    assert(num == 3);

  std::cerr << "The graph has " << num << " connected components (face connectivity)" << std::endl;
  //for(face_descriptor f : faces(sm)){
  //  std::cout  << &*f << " in connected component " << fccmap[f] << std::endl;
  //}

  PMP::keep_largest_connected_components(sm, 2, CGAL::parameters::vertex_index_map(vim));
  if (save_output)
    return;

  std::ofstream ofile("blobby_2cc_no_id.off");
  ofile << sm << std::endl;
  ofile.close();
}

void test_border_cases()
{
  std::cerr <<"Testing border cases\n";
  typedef boost::graph_traits<Mesh_with_id>::face_descriptor face_descriptor;
  typedef boost::graph_traits<Mesh_with_id>::vertex_descriptor vertex_descriptor;

  std::ifstream input(CGAL::data_file_path("meshes/elephant.off"));
  Mesh_with_id sm;
  input >> sm;

  std::size_t i=0;
  for(face_descriptor f : faces(sm))
    f->id() = i++;

  i=0;
  for(vertex_descriptor v : vertices(sm))
    v->id() = i++;

  boost::vector_property_map<int,
    boost::property_map<Mesh_with_id, boost::face_index_t>::type>
      fccmap(static_cast<unsigned>(num_faces(sm)), get(boost::face_index,sm));

  PMP::connected_components(sm, fccmap);
  std::size_t nb_faces=num_faces(sm);

  std::cerr <<" removing no faces\n";
  PMP::remove_connected_components(sm, std::vector<std::size_t>(), fccmap);
  assert( nb_faces == num_faces(sm) );

  std::cerr <<" removing all faces\n";
  Mesh_with_id copy = sm;
  PMP::connected_components(sm, fccmap);
  PMP::remove_connected_components(copy, std::vector<std::size_t>(1,0), fccmap);
  assert(num_vertices(copy)==0);

  std::cerr <<" keeping all faces\n";
  PMP::connected_components(sm, fccmap);
  PMP::keep_connected_components(sm, std::vector<face_descriptor>(1, *faces(sm).first));
  assert( nb_faces == num_faces(sm) );

  std::cerr <<" keeping no faces\n";
  copy = sm;
  PMP::connected_components(sm, fccmap);
  PMP::keep_connected_components(copy, std::vector<face_descriptor>());
  assert(num_vertices(copy)==0);
}

void keep_nothing(const std::string argv1)
{
  typedef boost::graph_traits<Mesh_with_id>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<Mesh_with_id>::face_descriptor face_descriptor;

  Mesh_with_id sm;
  std::ifstream in(argv1);
  if(!(in >> sm)) {
    std::cerr << "ERROR reading file: " << argv1 << std::endl;
    return;
  }
  int i=0;
  for(face_descriptor f : faces(sm))
    f->id() = i++;

  i=0;
  for(vertex_descriptor v : vertices(sm))
    v->id() = i++;

  PMP::keep_largest_connected_components(sm, 0);
  assert(num_vertices(sm) == 0);
  assert(num_edges(sm) == 0);
  assert(num_faces(sm) == 0);
}

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby_3cc.off");
  const bool save_output = (argc > 2);

  mesh_with_id(filename, save_output);
  mesh_no_id(filename, save_output);
  test_border_cases();
  keep_nothing(filename);

  return EXIT_SUCCESS;
}
