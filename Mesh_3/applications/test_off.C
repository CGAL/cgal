#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

typedef std::map<std::pair<int,int>, std::list<int> > Edge_map;
Edge_map edge_map;

char* argv0 = "";

void  add_in_edge_map(int nface, int i, int j)
{ 
 std::pair<int,int> edge = i < j ? 
   std::make_pair( i,j) :
   std::make_pair( j,i);
 Edge_map::iterator edge_it = edge_map.find(edge);
 if (edge_it == edge_map.end()) 
   edge_it = (edge_map.insert(std::make_pair(edge, std::list<int>()))).first;
 (edge_it->second).push_back(nface);
}

void usage(std::string error = "")
{
  if( error != "" )
    std:: cerr << "Error: " << error << std::endl;
  std::cerr << "Usage:\n" 
            << argv0 << " [OFF_FILE_NAME.off]\n"
            << "    will analyse OFF_FILE_NAME.off and display the list"
            << " of non manifold edges.\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  argv0 = argv[0];

  if (argc != 2)
    usage();

  std::ifstream ifs;
  ifs.open(argv[1],std::ifstream::in); 

  std::string heading;
  ifs >> heading;
  if(!ifs || heading != "OFF")
  {
    std::stringstream error_message;
    error_message << "bad OFF file.\n"
                  << "Expected \"OFF\"\n"
                  << "Got \"" << heading << "\"\n";
    usage(error_message.str());
  }

  // lecture des 3 nombres
  int nb_points,nb_faces, ent2;  
  ifs >> nb_points;
  ifs >> nb_faces;
  ifs >> ent2;  
  std::cerr << "number of points: "<< nb_points << std::endl;
  std::cerr << "number of facets: "<< nb_faces << std::endl;

  // read vertices
  Point p;
  for ( int ip = 0; ip < nb_points; ip++) {
    ifs >> p;
  }

  // read facets
  for (int iface = 0; iface < nb_faces ; iface++) {
    int npoints;
    ifs >> npoints;
    int vp[npoints];
    int nj;
    for (int j = 0; j < npoints;  j++) {
      ifs >> nj;
      vp[j] = nj ;
    }

    // add facet in edge map
    for (int j = 0; j < npoints-1 ; j++) 
      add_in_edge_map(iface,vp[j],vp[j+1]);

    add_in_edge_map(iface, vp[npoints-1], vp[0]);
  }

  int exit_code = EXIT_SUCCESS;

  // check surface
  std::cout << "List of non manifold edges:\n";
  for(Edge_map::iterator edge_it = edge_map.begin();
      edge_it != edge_map.end();
      edge_it++)
  {
    if (edge_it->second.size() != 2) 
    {
      std::pair<int,int> edge = edge_it->first;
      std::cout << "Edge " << edge.first << " "<< edge.second  
                << '\t'<< (edge_it->second).size() << '\t';
      std::cout << "Face list ";
      std::list<int>::iterator  face_it = (edge_it->second).begin();
      for( ; face_it != (edge_it->second).end() ; face_it++)
 	std::cout << *face_it << '\t';
      std::cout << std::endl;
      exit_code = EXIT_FAILURE;
    }
  }
  return exit_code;
}
