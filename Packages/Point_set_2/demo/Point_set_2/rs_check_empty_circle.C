//#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/range_search_delaunay_2.h>
#include <list>

typedef double                                             coord_type;
typedef CGAL::Simple_cartesian<coord_type>                 Gt;
typedef CGAL::Delaunay_triangulation_2<Gt>                 Delaunay;
typedef CGAL::Delaunay_triangulation_2<Gt>::Edge_iterator  Edge_iterator;
typedef CGAL::Delaunay_triangulation_2<Gt>::Vertex_handle  Vertex_handle;

typedef Gt::Point_2   Point;
typedef Gt::Circle_2  Circle;

Delaunay PS;

void output(CGAL::Window_stream& W, const Delaunay& PSet)
{
  W.clear();
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    CGAL::Segment_2<Gt> s= PSet.segment(*eit);
    W << s;
  }
}

void redraw(CGAL::Window_stream* wptr)
{ output(*wptr,PS); }


class check_empty {
public:
  bool    result;
  Circle  c;
  
  check_empty(Circle cact) : result(false), c(cact) { }
  
  bool get_result() const  { return result; }
  void set_result(bool nr) { result=nr; }
  
  bool operator()(const Point& p)
  {
    return ! c.has_on_unbounded_side(p); 
  }
};

int main()
{
  CGAL::Window_stream W(600,500,"Range search operations on a point set checking emptiness of circles");  

  W.init(-500,500,-400);
  W.set_redraw(redraw);
  W.display(100,100);
  
#if defined(CGAL_USE_CGAL_WINDOW)
  W.set_point_style(CGAL::disc_point);
#else
  W.set_point_style(leda_disc_point);
#endif   
  
  W.draw_text(-260,20, "Input some points; quit input with the right mouse button");

  CGAL::Point_2<Gt> pnew;

  while (W >> pnew) {
    PS.insert(pnew);
    output(W,PS);    
  }
  
  std::list<Vertex_handle> LV;  
  
  std::cout << "circular range search !\n";
  W << CGAL::BLACK;
  output(W,PS);
  W.draw_text(-450,-350, "Input a circle; we perform a range search (quit: right mouse button) ... ");
  
  CGAL::Circle_2<Gt> rc;
   
  while (W >> rc) {
   W << CGAL::BLACK;
   output(W,PS);
   W.draw_text(-450,-350, "Input a circle; we perform an empty circle check (quit: right mouse button) ... ");
   W << CGAL::RED; W << rc; W << CGAL::BLACK;

   check_empty checker(rc);

   CGAL::range_search(PS,rc,std::back_inserter(LV),checker,true);
   
   if (checker.get_result()) std::cout << "circle not empty !\n";
   else std::cout << "circle was empty !\n";
   
   /*  
   W << CGAL::RED;
  
   for(vit=LV.begin(); vit!=LV.end(); vit++){
     W << (*vit)->point();
   }
   */
   
   LV.clear();
  }
 
  return 0;
}

