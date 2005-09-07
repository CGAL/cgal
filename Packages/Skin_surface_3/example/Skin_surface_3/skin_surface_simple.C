// examples/Skin_surface_3/skin_surface_simple.C
#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Mixed_complex_builder_3.h>

typedef CGAL::Skin_surface_traits_3<>                 Skin_traits;
typedef Skin_traits::Regular_traits                   Regular_traits;
typedef Skin_traits::Regular                          Regular;
typedef Regular_traits::Weighted_point                Reg_weighted_point;
typedef Regular_traits::Bare_point                    Reg_point;
typedef Skin_traits::Simplicial                       Simplicial;

typedef CGAL::Mixed_complex_builder_3<Skin_traits>    Mixed_complex_builder;

#include <fstream>

template <class T>
void write_edges(T &t, char * filename) {
  typedef typename T::Finite_vertices_iterator      Finite_vertices_iterator;
  typedef typename T::Finite_edges_iterator         Finite_edges_iterator;
  typedef typename T::Vertex_handle                 Vertex_handle;
    
  CGAL::Unique_hash_map<Vertex_handle, int> vertices;
  int i = 0;
  
  std::ofstream out(filename);
  out << "SKEL " << t.number_of_vertices() 
      << " " << t.number_of_finite_edges() << std::endl;
  for (Finite_vertices_iterator vit = t.finite_vertices_begin();
       vit != t.finite_vertices_end(); vit++) {
    // Avoid writing weighted points:
    out << vit->point().x() << " " << vit->point().y() << " " 
        << vit->point().z() << " " << std::endl;
    vertices[vit] = i;
    i++;
  }
  for (Finite_edges_iterator eit = t.finite_edges_begin();
       eit != t.finite_edges_end(); eit++) {
    out << "2 " << vertices[eit->first->vertex(eit->second)]
        << " " << vertices[eit->first->vertex(eit->third)] 
	<< std::endl;
  }
}

int main(int argc, char *argv[]) {
  double shrink = .5;
  std::ifstream is("./data/caffeine.cin");
  
  Regular regular;
  Reg_weighted_point wp;
  double max_weight=1;
  CGAL::Bbox_3 box;
  while (is >> wp) {
    max_weight = std::max(max_weight, wp.weight());
    box = box + wp.point().bbox();
    regular.insert(wp);
  }
  // add a bounding octahedron:
  Reg_point mid((box.xmin() + box.xmax())/2,
                (box.ymin() + box.ymax())/2,
                (box.zmin() + box.zmax())/2);
  double size = 1.5*((box.xmax() - box.xmin() +
                      box.ymax() - box.ymin() +
                      box.zmax() - box.zmin())/2 + max_weight);
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x()+size,mid.y(),mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x()-size,mid.y(),mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y()+size,mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y()-size,mid.z()),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()+size),-1));
  regular.insert(
    Reg_weighted_point(Reg_point(mid.x(),mid.y(),mid.z()-size),-1));
  write_edges(regular, "triangulation.skel");
  
  Simplicial simplicial;
  Mixed_complex_builder(regular, simplicial, shrink);

  CGAL_assertion(simplicial.is_valid());
  
  write_edges(simplicial, "mixed.skel");
  return 0;
}
