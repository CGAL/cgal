#include <CGAL/Cartesian.h>
//#include <CGAL/Simple_cartesian.h>
#include <list>
#include <CGAL/Point_set_2.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>    K;
//typedef CGAL::Simple_cartesian<double>    K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_2<K>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<K,Vb,Fb> Tds;

typedef CGAL::Point_set_2<K,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K,Tds>::Vertex_handle  Vertex_handle;


void output(CGAL::Window_stream& W, const CGAL::Point_set_2<K,Tds>& PSet)
{
  W.clear();
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    CGAL::Segment_2<K> s= PSet.segment(*eit);
    W << s;
  }
}

int main()
{
  CGAL::Point_set_2<K,Tds> PSet;

  CGAL::Window_stream W(600,500, "Finding nearest neighbor / k nearest neighbors");  

  W.init(-500,500,-400);
  W.display(100,100);
  
#if defined(CGAL_USE_CGAL_WINDOW)
  W.set_point_style(CGAL::disc_point);
#else
  W.set_point_style(leda_disc_point);
#endif  
  
  W.draw_text(-260,20, "Input some points; quit input with the right mouse button");
  
  CGAL::Point_2<K> actual;
  int i=0;
  
  while (W >> actual){
   PSet.insert(actual);
   output(W,PSet);
   i++;
  }  
  
  std::cout << i << " points were inserted !\n";

  // nearest neighbor ...  
  W.draw_text(-450,-350, "Input a point; we display the nearest neighbor ... ");
  
  for (i=0; i<5; i++){
    W >> actual;
    Vertex_handle v = PSet.nearest_neighbor(actual);
    
    if (v != NULL) {
    
     CGAL::Segment_2<K> my_seg(actual,v->point());
    
     W << CGAL::RED << v->point() << CGAL::BLACK;
     W << CGAL::BLUE << my_seg << CGAL::BLACK;
    }
  }
  
  // k nearest neighbors ...
  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it;
 
  output(W,PSet);
  W.draw_text(-450,-350, "Input a point; we display the 5 nearest neighbors ... "); 
  
  for (i=0; i<5; i++){
    L.clear();
    W >> actual;
    PSet.nearest_neighbors(actual,5, std::back_inserter(L));
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

  W.read_mouse();

  return 1;
}

