#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/partition.h>

#include <fstream>
#include <iostream>

int main(int argc, char** argv)
{
  typedef CGAL::Simple_cartesian<double>                            K;

  std::ifstream in((argc>1) ? argv[1] : "data/eight.off");
  int number_of_parts = (argc>2) ? atoi(argv[2]) : 8;

  if(!in)
  {
    std::cerr << "Error: could not read input file" << std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Polyhedron_3<K>                                     PM;
  PM pm;
  in >> pm; // read the mesh

  // The face <--> partition_id property map
  typedef CGAL::dynamic_face_property_t<int>                        Face_property_tag;
  typedef boost::property_map<PM, Face_property_tag>::type          Face_id_map;
  Face_id_map partition_id_map = get(Face_property_tag(), pm);

  // Set some custom options for METIS
  idx_t options[METIS_NOPTIONS];

  // Set all options to default ahead of manually editing some values
  METIS_SetDefaultOptions(options);

  // See METIS documentation for details on the many options
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
  options[METIS_OPTION_NCUTS] = 3;
  options[METIS_OPTION_NITER] = 10;
  options[METIS_OPTION_SEED] = 12345;
  options[METIS_OPTION_MINCONN] = 1;
  options[METIS_OPTION_CONTIG] = 1;
  options[METIS_OPTION_UFACTOR] = 25;
  options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;

  // Partition the mesh and output its parts
  CGAL::METIS::partition_graph(pm, number_of_parts,
                               CGAL::parameters::vertex_index_map(get(boost::vertex_external_index, pm))
                                                .face_partition_id_map(partition_id_map)
                                                .METIS_options(&options));

  return EXIT_SUCCESS;
}
