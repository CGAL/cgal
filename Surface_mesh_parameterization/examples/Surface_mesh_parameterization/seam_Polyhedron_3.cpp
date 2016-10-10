#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/parameterize.h>
#include <CGAL/Unique_hash_map.h>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Polyhedron_3<Kernel>          SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor SM_vertex_descriptor;

typedef CGAL::Unique_hash_map<SM_halfedge_descriptor,Point_2> UV_uhm;
typedef CGAL::Unique_hash_map<SM_edge_descriptor,bool> Seam_edge_uhm;
typedef CGAL::Unique_hash_map<SM_vertex_descriptor,bool> Seam_vertex_uhm;

typedef boost::associative_property_map<UV_uhm> UV_pmap;
typedef boost::associative_property_map<Seam_edge_uhm> Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm> Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

/// A helper class that writes a face as a polyline in the xy-plane
struct Face2Polyline
{
  const Mesh& mesh;
  const UV_pmap uvpm;

  Face2Polyline(const Mesh& mesh, const UV_pmap& uvpm)
    : mesh(mesh), uvpm(uvpm)
  { }

  void operator()(face_descriptor fd) const
  {
    halfedge_descriptor hd = halfedge(fd,mesh);

    std::cout << "4 " << uvpm[target(hd,mesh)].x()
                      << " " << uvpm[target(hd,mesh)].y() << " 0 ";

    hd = next(hd,mesh);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd,mesh)){
      std::cout << uvpm[vd].x() << " " << uvpm[vd].y() << " 0 ";
    }

    std::cout << std::endl;
  }
};

int main(int argc, char * argv[])
{
  SurfaceMesh sm;

  std::ifstream in_mesh((argc>1)?argv[1]:"data/lion.off");
  in_mesh >> sm;

  std::vector<SM_vertex_descriptor> id2v;
  BOOST_FOREACH(SM_vertex_descriptor vd, vertices(sm)){
    id2v.push_back(vd);
  }

  // Two property maps to store the seam edges and vertices
  Seam_edge_uhm seam_edge_uhm(false);
  Seam_edge_pmap seam_edge_pm(seam_edge_uhm);

  Seam_vertex_uhm seam_vertex_uhm(false);
  Seam_vertex_pmap seam_vertex_pm(seam_vertex_uhm);

  std::ifstream in((argc>2)?argv[2]:"data/lion.selection.txt");
  int s, t;
  SM_halfedge_descriptor smhd;
  while(in >> s >> t){
    SM_vertex_descriptor svd = id2v[s], tvd = id2v[t];
    SM_edge_descriptor ed = edge(svd, tvd,sm).first;
    if(!is_border(ed,sm)){
      put(seam_edge_pm, ed, true);
      put(seam_vertex_pm, svd, true);
      put(seam_vertex_pm, tvd, true);
      if(smhd == boost::graph_traits<SurfaceMesh>::null_halfedge()){
        smhd = halfedge(edge(svd,tvd,sm).first,sm);
      }
    }
  }

  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that the uv
  // is only stored for the canonical halfedges representing a vertex

  UV_uhm uv_uhm;
  UV_pmap uv_pm(uv_uhm);

  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd,mesh); // a halfedge on the virtual border

  CGAL::parameterize(mesh, bhd, uv_pm);

  Face2Polyline f2p(mesh,uv_pm);
  // The seam may define a patch, that is we take the connected component of the border halfedge

  CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd,mesh),mesh),
                                                     mesh,
                                                     boost::make_function_output_iterator(f2p));

  return 0;
}

