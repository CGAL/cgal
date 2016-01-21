#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef CGAL::Seam_mesh<Mesh>                                Seam_mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;


int main(int argc, char* argv[]) 
{
  Mesh sm;
  std::ifstream in((argc>1)?argv[1]:"data/cube.off");
  in >> sm;

  std::vector<edge_descriptor> seam;  

  // split cube in two connected components
  BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(*halfedges(sm).first,sm)){
    seam.push_back(edge(hd,sm));
  }

#if 0
  //cube 
  halfedge_descriptor hd = * halfedges(sm).first;
  std::cout << "center = " << target(hd,sm) << std::endl;
  int count=0;
  BOOST_FOREACH(hd, halfedges_around_target(hd,sm)){
    std::cout << "v = " << source(hd,sm) << std::endl;
    if(count==0 || count==2){ 
      seam.push_back(edge(hd,sm));
    }
    if(seam.size() == 2)
      break;
    count++;
  }
#endif 
 
#if 0
  //cube-ouvert
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm)){
    if(is_border(hd,sm)){
      seam.push_back(edge(next(opposite(hd,sm),sm),sm));
      break;
    }
  }
#endif 

  typedef std::map<Seam_mesh::halfedge_descriptor,int> H_vertex_index_map;
  typedef boost::associative_property_map<H_vertex_index_map> H_vertex_index_pmap;
  H_vertex_index_map hvim; H_vertex_index_pmap hvipm(hvim);

  halfedge_descriptor mhd = halfedge(seam.front(),sm);    
  Seam_mesh ssm(sm, seam);

  std::vector<face_descriptor> faces;
  boost::graph_traits<Seam_mesh>::halfedge_descriptor shd(halfedge(seam.front(),sm));
  CGAL::Polygon_mesh_processing::connected_component(face(shd,ssm),
                                                     ssm,
                                                     std::back_inserter(faces));
  std::cerr << faces.size() << std::endl;
  faces.clear();
  shd = opposite(halfedge(seam.front(),sm),sm);
  CGAL::Polygon_mesh_processing::connected_component(face(shd,ssm),
                                                     ssm,
                                                     std::back_inserter(faces));
  std::cerr << faces.size() << std::endl;

  return 0;  

  BOOST_FOREACH(Seam_mesh::vertex_descriptor vd, vertices(ssm)){
    halfedge_descriptor hd = vd;
    std::cerr << vd << " has incident halfedges:" << std::endl;
    BOOST_FOREACH(Seam_mesh::halfedge_descriptor hd2, halfedges_around_target(vd,ssm)){
      std::cerr << hd2 << std::endl;
    }
  }
  std::cerr << "done"<< std::endl;


  Seam_mesh::halfedge_descriptor smhd(mhd), smhd2(opposite(mhd,sm));
  std::cout << "target = " << target(smhd,ssm) << std::endl;

  std::cout << "vertices on seam" << std::endl;
  BOOST_FOREACH(Seam_mesh::halfedge_descriptor hd, halfedges_around_face(opposite(smhd,ssm),ssm)){
    std::cout << target(hd.tmhd,sm) << std::endl;
  }

  std::cout << "vertices around target in mesh" << std::endl;
   BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(mhd,sm)){
     std::cout << source(hd,sm) << std::endl;
  }
  
  std::cout << "vertices around target in seam mesh" << std::endl;
   BOOST_FOREACH(Seam_mesh::halfedge_descriptor hd, halfedges_around_target(smhd,ssm)){
     std::cout << source(hd.tmhd,sm) << std::endl;
  }

  std::cout << "vertices around source in seam mesh" << std::endl;
  BOOST_FOREACH(Seam_mesh::halfedge_descriptor hd, halfedges_around_source(opposite(smhd,ssm),ssm)){
     std::cout << target(hd.tmhd,sm) << std::endl;
  }
  std::cout << "vertices around source in seam mesh" << std::endl;
   BOOST_FOREACH(Seam_mesh::halfedge_descriptor hd, halfedges_around_source(smhd2,ssm)){
     std::cout << target(hd.tmhd,sm) << std::endl;
  }
  std::cout << "vertices around target in seam mesh" << std::endl;
  BOOST_FOREACH(Seam_mesh::halfedge_descriptor hd, halfedges_around_target(opposite(smhd2,ssm),ssm)){
     std::cout << source(hd.tmhd,sm) << std::endl;
  }

  
  boost::property_map<Seam_mesh, CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point,ssm);
  //std::cout << get(vpm, source(hd,ssm)) << std::endl;
  


  return 0;
}

