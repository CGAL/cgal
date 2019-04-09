#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/partition.h>

#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>                           K;
typedef CGAL::Surface_mesh<K::Point_3>                           SM;

int main(int argc, char** argv)
{
  std::ifstream in((argc>1) ? argv[1] : "data/blobby.off");
  int number_of_parts = (argc>2) ? atoi(argv[2]) : 8;

  if(!in) {
    std::cerr << "Error: could not read input file" << std::endl;
    return EXIT_FAILURE;
  }

  SM sm;
  CGAL::read_off(in, sm);

  // The vertex <--> partition_id property map
  typedef SM::Property_map<SM::Vertex_index, std::size_t>          Vertex_id_map;
  Vertex_id_map vertex_pid_map = sm.add_property_map<SM::Vertex_index, std::size_t>("v:pid").first;

  // The face <--> partition_id property map
  typedef SM::Property_map<SM::Face_index, std::size_t>            Face_id_map;
  Face_id_map face_pid_map = sm.add_property_map<SM::Face_index, std::size_t>("f:pid").first;

  // Partition the mesh
  CGAL::METIS::partition_dual_graph(sm, number_of_parts,
                                    CGAL::parameters::vertex_partition_id_map(vertex_pid_map)
                                                     .face_partition_id_map(face_pid_map));

  // Extract the part n°0 of the partition into a new, independent mesh
  typedef CGAL::Face_filtered_graph<SM>                            Filtered_graph;
  Filtered_graph filtered_sm(sm, 0 /*id of th part*/, face_pid_map);
    CGAL_assertion(filtered_sm.is_selection_valid());
  SM part_sm;
  CGAL::copy_face_graph(filtered_sm, part_sm);

  // Output the mesh extracted from subpart n°0
  std::ofstream out("sm_part_0.off");
  out.precision(17);
  CGAL::write_off(out, part_sm);

  // Output all the vertices that are in the part n°0
  std::ofstream outxyz("out.xyz");
  outxyz.precision(17);
  boost::graph_traits<SM>::vertex_iterator vit, ve;
  boost::tie(vit, ve) = vertices(sm);
  for(; vit!=ve; ++vit)
  {
    if(get(vertex_pid_map, *vit) == 0)
      outxyz << sm.point(*vit) << std::endl;
  }

  return EXIT_SUCCESS;
}
