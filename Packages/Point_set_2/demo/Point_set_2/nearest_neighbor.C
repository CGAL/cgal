#include <CGAL/config.h>
#include <list>
#include <CGAL/Cartesian.h>
//#include <CGAL/Simple_cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Point_set_2.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>    REP;
//typedef CGAL::Simple_cartesian<double>    REP;

typedef CGAL::Triangulation_euclidean_traits_2<REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

typedef CGAL::point_set_traits_2<REP>     TRAITS;
typedef CGAL::Point_set_2<TRAITS,Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2<TRAITS,Gt,Tds>::Vertex  Vertex;




void output(CGAL::Window_stream& W, const CGAL::Point_set_2<TRAITS,Gt,Tds>& PSet)
{
  W.clear();
  Edge e;
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    e = *eit;
    CGAL::Segment_2<REP> s= PSet.seg(e);
    W << s;
  }
}

int main()
{
  CGAL::Point_set_2<TRAITS,Gt,Tds> PSet;

  CGAL::Window_stream W(600,500, "Finding nearest neighbor / k nearest neighbors");  
  //CGAL::cgalize(W);

  W.init(-500,500,-400);
  W.display(100,100);
  
  W.draw_text(-260,20, "Input some points; quit input with the right mouse button");
  
  CGAL::Point_2<REP> actual;
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
    
     CGAL::Segment_2<REP> my_seg(actual,PSet.pos(v));
    
     W << CGAL::RED << PSet.pos(v) << CGAL::BLACK;
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
      W << CGAL::GREEN << PSet.pos(*it) << CGAL::BLACK;      
      std::cout << PSet.pos(*it) << "\n";
    }
    std::cout << "\n";
  }  

  W.read_mouse();

  return 1;
}

