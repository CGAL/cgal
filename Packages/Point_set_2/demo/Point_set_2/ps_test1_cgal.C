#include <CGAL/config.h>
#include <list>

//#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Triangle_2.h>
#include <CGAL/Iso_rectangle_2.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/IO/Window_stream.h>


//typedef CGAL::Cartesian<double>          REP;
typedef CGAL::Simple_cartesian<double>    REP;

typedef CGAL::Triangulation_euclidean_traits_2<REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

typedef CGAL::Point_set_2<Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2<Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2<Gt,Tds>::Vertex  Vertex;

typedef CGAL::Iso_rectangle_2<REP>     Rectangle;
typedef CGAL::Triangle_2<REP>          Triangle;


void output(CGAL::Window_stream& W, const CGAL::Point_set_2<Gt,Tds>& PSet)
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
  CGAL::Point_set_2<Gt,Tds> PS;

  CGAL::Window_stream W(600,500,"Range search operations on a point set");  
  //CGAL::cgalize( W);

  W.init(-500,500,-400);
  W.display(100,100);
  
  W.draw_text(-260,20, "Input some points; quit input with the right mouse button");

  CGAL::Point_2<REP> pnew;

  while (W >> pnew) {
    PS.insert(pnew);
    output(W,PS);    
  }
  
/*  
  for (i=0; i<2; i++) {
    CGAL::Point_2<REP> pnew;
    W >> pnew;
    Vertex_handle v = PS.nearest_neighbor(pnew);
    PS.del(v);
    output(W,PS);    
  }
*/  
  std::list<Vertex_handle>::const_iterator vit;
  std::list<Vertex_handle> LV;  
  
  std::cout << "circular range search !\n";
  W << CGAL::BLACK;
  output(W,PS);
  W.draw_text(-450,-350, "Input a circle; we perform a range search (quit: right mouse button) ... ");
  
  CGAL::Circle_2<REP> rc;
   
  while (W >> rc) {
   W << CGAL::BLACK;
   output(W,PS);
   W.draw_text(-450,-350, "Input a circle; we perform a range search (quit: right mouse button) ... ");
   W << rc;

   PS.range_search(rc,std::back_inserter(LV));  
   W << CGAL::RED;
  
   for(vit=LV.begin(); vit!=LV.end(); vit++){
     W << PS.pos(*vit);
   }
   LV.clear();
  }
 
  std::cout << "triangular range search !\n";
  W << CGAL::BLACK;
  output(W,PS);  
  W.draw_text(-450,-350, "Input a triangle; we perform a range search (quit: right mouse button) ... ");
    
  CGAL::Point_2<REP> pt1,pt2,pt3,pt4;
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
     W << PS.pos(*vit);
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
     W << PS.pos(*vit);
   }
   LV.clear();  
  }

  return 1;
}

