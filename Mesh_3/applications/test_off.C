#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <set>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef  K::Point_3                                           Point;

typedef  std::map<std::pair<int,int>, std::list<int> > Edge_map;
Edge_map edge_map;

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


int main( int argc, char *argv[])
{
  if (argc != 2){
    std::cerr << "usage test_off_input filemame " << std::endl;
    return 0;
  }

  std::ifstream ifs;
  ifs.open(argv[1],std::ifstream::in); 

  std::string heading;
  ifs >> heading;
  if(!ifs || heading != "OFF")
    std::cerr << "bad OFF file. Expected \"OFF\"\n"
              << "got \"" << heading << "\"\n";

  // lecture des 3 nombres
  int nb_points,nb_faces, ent2;  
  ifs >> nb_points;
  ifs >> nb_faces;
  ifs >> ent2;  
  std::cerr<<"nombre de points : "<<nb_points<<std::endl;
  std::cerr<<"nombre de faces : "<<nb_faces<<std::endl;
  std::cerr<<"ent2 : "<<ent2<<std::endl;



// read vertices
  Point p;
  for ( int ip = 0; ip < nb_points; ip++) {
    ifs >> p;
    std::cerr << p << std::endl;
  }

 // read facets
  for ( int iface = 0; iface < nb_faces ; iface++) {
    
    int npoints;
    ifs >> npoints;
    std::cerr << iface << "\t" << npoints <<  "\t";
    int vp[npoints];
    int nj;
    for (int j = 0; j < npoints;  j++) {
      ifs >> nj;
      vp[j] = nj ;
    }
    for (int j = 0; j < npoints;  j++) std::cerr <<  vp[j] << " ";
    std::cerr << std::endl;
    
     // add facet in edge map
     for (int j = 0; j < npoints-1 ; j++) add_in_edge_map(iface,vp[j],vp[j+1]);
     add_in_edge_map(iface, vp[npoints-1], vp[0]);
    
  }
    

   // check surface
   Edge_map::iterator edge_it=edge_map.begin();
   for( ; edge_it != edge_map.end(); edge_it++){
     if (edge_it->second.size() != 2) 
       {
	std::pair<int,int> edge = edge_it->first;
 	std::cerr << "Edge " << edge.first << " "<< edge.second  
 		  << '\t'<< (edge_it->second).size() << '\t';
	std::cerr << "Face list ";
	std::list<int>::iterator  face_it = (edge_it->second).begin();
        for( ; face_it != (edge_it->second).end() ; face_it++)
 	std::cerr << *face_it << '\t';
        std::cerr << std::endl;
      }
  }
  return 1;
}

