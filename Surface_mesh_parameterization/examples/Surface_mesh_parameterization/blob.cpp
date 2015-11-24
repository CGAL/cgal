#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/parameterize.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <cstdlib>
#include <fstream>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
typedef CGAL::Seam_mesh<SurfaceMesh>        Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef boost::graph_traits<SurfaceMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor SM_vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor SM_face_descriptor;



SM_vertex_descriptor find(const Point_3& p, const SurfaceMesh& sm)
{
  BOOST_FOREACH(SM_vertex_descriptor vd, vertices(sm)){
    if(squared_distance(p, sm.point(vd)) < 0.0001){
      return vd;
    }
  }
  std::cerr << "epsilon is too small" << std::endl;
  return SM_vertex_descriptor();
}

void
print(SM_halfedge_descriptor hd, const SurfaceMesh& sm)
{
  std::cout << "2 " << sm.point(source(hd, sm)) << "  " << sm.point(target(hd, sm)) << std::endl;
}

int main(int argc, char * argv[])
{
  SurfaceMesh sm;
  std::vector<SM_edge_descriptor> seam; 
  {
    std::ifstream in((argc>1)?argv[1]:"data/blob.off");
    in >> sm;
  }

  {
    std::ifstream in((argc>2)?argv[2]:"data/blob.polylines.txt");
    std::vector<Point_3> points;
    int n;
    in >> n;
    Point_3 p;
    in >> p;
    SM_vertex_descriptor v = find(p,sm);

    for(int i=1; i < n; i++){
      in >> p;
      SM_vertex_descriptor w = find(p,sm);
      seam.push_back(edge(v,w,sm).first);
      v = w;
    }
  }
  
#if 1
  typedef CGAL::Discrete_authalic_parameterizer_3<Mesh> Parameterizer;
#else
 typedef CGAL::LSCM_parameterizer_3<Mesh,
                           CGAL::Two_vertices_parameterizer_3<Mesh> > Parameterizer;
#endif

  // The 2D points of the uv parametrisation will be written into this map  
  typedef std::map<halfedge_descriptor,Point_2> H_uv_map;
  typedef boost::associative_property_map<H_uv_map> H_uv_pmap;
  H_uv_map huvm; H_uv_pmap huvpm(huvm);
  
  // When we have seams vertices are "duplicated"
  // We give halfedges that have as target the same duplicated vertex the same index
  typedef std::map<halfedge_descriptor,int> H_vertex_index_map;
  typedef boost::associative_property_map<H_vertex_index_map> H_vertex_index_pmap;
  H_vertex_index_map hvim; H_vertex_index_pmap hvipm(hvim);

  SM_halfedge_descriptor smhd = halfedge(seam.front(),sm);    
  Mesh mesh(sm, seam, smhd, hvipm);

  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd,mesh); // a halfedge on the virtual hole


  // check that all halfedges are marked
  BOOST_FOREACH(SM_halfedge_descriptor hd, halfedges(sm)){
    assert(get(hvipm,hd) != -1);
    assert(get(hvipm,halfedge_descriptor(hd)) != -1);
  }


  Parameterizer::Error_code err = CGAL::parameterize(mesh, Parameterizer(), bhd, huvpm, hvipm);
  switch(err) {
  case Parameterizer::OK: // Success
    break;
  case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
  case Parameterizer::ERROR_NON_TRIANGULAR_MESH:   
  case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:     
  case Parameterizer::ERROR_BORDER_TOO_SHORT:    
    std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
    break;
  default: // Error
    std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
    break;
  };
    
  BOOST_FOREACH(SM_face_descriptor fd, faces(sm)){  
    SM_halfedge_descriptor hd = halfedge(fd,sm);
    std::cout << "4 " << huvm[hd].x() << " " << huvm[hd].y() << " 0 ";
    hd = next(hd,sm);
    BOOST_FOREACH(SM_halfedge_descriptor hd2, halfedges_around_face(hd,sm)){
      std::cout << huvm[hd2].x() << " " << huvm[hd2].y() << " 0 ";
    }
    std::cout << std::endl;
  }
  
  return EXIT_SUCCESS;
  

}
