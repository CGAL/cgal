//#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <list>
#include <CGAL/Point_set_2.h>
#include <CGAL/IO/Window_stream.h>

//typedef CGAL::Cartesian<double>          K;
typedef CGAL::Simple_cartesian<double>    K;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_face_base_2<K>   Fb;
typedef CGAL::Triangulation_default_data_structure_2<K,Vb,Fb> Tds;

typedef CGAL::Point_set_2<K,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K,Tds>::Vertex_handle  Vertex_handle;

typedef CGAL::Iso_rectangle_2<K>     Rectangle;
typedef CGAL::Triangle_2<K>          Triangle;


CGAL::Point_set_2<K,Tds> PS;

void output(CGAL::Window_stream& W, const CGAL::Point_set_2<K,Tds>& PSet)
{
  W.clear();
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    CGAL::Segment_2<K> s= PSet.segment(*eit);
    W << s;
  }
}

void redraw(CGAL::Window_stream* wptr)
{
  output(*wptr,PS);
}

int main()
{
  CGAL::Window_stream W(600,500,"Range search operations on a point set");  
  //CGAL::cgalize( W);

  W.init(-500,500,-400);
  W.set_redraw(redraw); 
  W.display(100,100);
  
#if defined(CGAL_USE_CGAL_WINDOW)
  W.set_point_style(CGAL::disc_point);
#else
  W.set_point_style(leda_disc_point);
#endif   
  
  W.draw_text(-260,20, "Input some points; quit input with the right mouse button");

  CGAL::Point_2<K> pnew;

  while (W >> pnew) {
    PS.insert(pnew);
    output(W,PS);    
  }
  
  std::list<Vertex_handle>::const_iterator vit;
  std::list<Vertex_handle> LV;  
  
  std::cout << "circular range search !\n";
  W << CGAL::BLACK;
  output(W,PS);
  W.draw_text(-450,-350, "Input a circle; we perform a range search (quit: right mouse button) ... ");
  
  CGAL::Circle_2<K> rc;
   
  while (W >> rc) {
   W << CGAL::BLACK;
   output(W,PS);
   W.draw_text(-450,-350, "Input a circle; we perform a range search (quit: right mouse button) ... ");
   W << rc;

   PS.range_search(rc,std::back_inserter(LV));  
   W << CGAL::RED;
  
   for(vit=LV.begin(); vit!=LV.end(); vit++){
     W << (*vit)->point();
   }
   LV.clear();
  }
 
  std::cout << "triangular range search !\n";
  W << CGAL::BLACK;
  output(W,PS);  
  W.draw_text(-450,-350, "Input a triangle; we perform a range search (quit: right mouse button) ... ");
    
  CGAL::Point_2<K> pt1,pt2,pt3,pt4;
  Triangle Tr;
  
  while (W >> Tr){
   W << CGAL::BLACK;
   output(W,PS);
   W.draw_text(-450,-350, "Input a triangle; we perform a range search (quit: right mouse button) ... ");
   W << Tr;
  
   pt1=Tr[0]; pt2=Tr[1]; pt3=Tr[2];
  
   PS.range_search(pt1,pt2,pt3,std::back_inserter(LV));
   W << CGAL::GREEN;
   for(vit=LV.begin(); vit!=LV.end(); vit++){
     W << (*vit)->point();
   }  
   LV.clear();
  }
  
  W << CGAL::BLACK;
 
  std::cout << "rectangular range search !\n";
  W << CGAL::BLACK;
  output(W,PS);
  W.draw_text(-450,-350, "Input a rectangle; we perform a range search (quit: right mouse button) ... ");
  
  Rectangle Rect;
  while (W >> Rect) {
   W << CGAL::BLACK;
   output(W,PS);
    W.draw_text(-450,-350, "Input a rectangle; we perform a range search (quit: right mouse button) ... ");
   W << Rect;
  
   pt1=Rect[3]; pt2=Rect[0]; pt3=Rect[1]; pt4=Rect[2]; 
  
   PS.range_search(pt1,pt2,pt3,pt4,std::back_inserter(LV));
   W << CGAL::YELLOW;
   for(vit=LV.begin(); vit!=LV.end(); vit++){
     W << (*vit)->point();
   }
   LV.clear();  
  }

  return 1;
}

