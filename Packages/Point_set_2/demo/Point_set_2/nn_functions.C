#include <CGAL/Cartesian.h>
//#include <CGAL/Simple_cartesian.h>
#include <list>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>

typedef double                                             coord_type;
typedef CGAL::Cartesian<coord_type>                        Gt;
typedef CGAL::Delaunay_triangulation_2<Gt>                 Delaunay;
typedef CGAL::Delaunay_triangulation_2<Gt>::Edge_iterator  Edge_iterator;
typedef CGAL::Delaunay_triangulation_2<Gt>::Vertex_handle  Vertex_handle;

void output(CGAL::Window_stream& W, const Delaunay& PSet)
{
  W.clear();
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    CGAL::Segment_2<Gt> s= PSet.segment(*eit);
    W << s;
  }
}

int main()
{
  Delaunay PSet;

  CGAL::Window_stream W(600,500, "Finding nearest neighbor / k nearest neighbors");  

  W.init(-500,500,-400);
  W.set_grid_dist(10.0);
  W.set_grid_mode(10);
  W.display(100,100);
  
#if defined(CGAL_USE_CGAL_WINDOW)
  W.set_point_style(CGAL::disc_point);
#else
  W.set_point_style(leda_disc_point);
#endif  
  
  W.draw_text(-260,20, "Input some points; quit input with the right mouse button");
  
  CGAL::Point_2<Gt> actual;
  int i=0;
  
  while (W >> actual){
   PSet.insert(actual);
   output(W,PSet);
   i++;
  }  
  
  std::cout << i << " points were inserted !\n";

  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it, it2;
 
  W.draw_text(-450,-350, "Input a point; we display the 5 nearest neighbors ... "); 
  
  for (i=0; i<10; i++){
    L.clear();
    W >> actual;
    
    CGAL::nearest_neighbors(PSet, actual, 5, std::back_inserter(L));
    std::cout << "actual point: " << actual << "\n";
    
    W.clear();
    W.draw_text(-450,-350, "Input a point; we display the 5 nearest neighbors ... "); 
    output(W,PSet);
    W << CGAL::RED << actual << CGAL::BLACK;
    
    for (it=L.begin();it != L.end(); it++){
      W << CGAL::GREEN << (*it)->point() << CGAL::BLACK;      
      std::cout << (*it)->point() << "\n";
    }
    std::cout << "\n";
  }  

  output(W,PSet);
  W.draw_text(-450,-350, "Input a point; we perform the lookup operation ... ");  
  
  for (i=0; i<10; i++){ 
    W >> actual;
    
    Vertex_handle vh = CGAL::lookup(PSet, actual);
    if (vh != NULL) {
      W << CGAL::GREEN << vh->point() << CGAL::BLACK; 
    }
  }   

  W.read_mouse();
  
  // k nearest neighbors of all vertices ...
  L.clear();
  
  CGAL::get_vertices(PSet, std::back_inserter(L));
  for (it=L.begin();it != L.end(); it++){
    output(W,PSet);
    W.draw_text(-450,-350, "We display the 7 nearest neighbors of every vertex (the number 7 includes the vertex)  ... ");
    Vertex_handle v = *it;
    std::list<Vertex_handle> act;
    CGAL::nearest_neighbors(PSet, v, 7, std::back_inserter(act));
    std::cout << "actual point: " << (*it)->point() << "\n";
    for (it2=act.begin();it2 != act.end(); it2++){
      W << CGAL::GREEN << (*it2)->point() << CGAL::BLACK; 
      std::cout << (*it2)->point() << "\n";     
    }    
    W << CGAL::RED << (*it)->point() << CGAL::BLACK; 
    W.read_mouse();
  }
  

  return 1;
}





