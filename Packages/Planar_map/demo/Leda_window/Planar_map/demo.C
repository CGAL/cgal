#include <iostream>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include "configuration"

#include <CGAL/basic.h>
#include <fstream>


#ifdef USE_NAIVE_POINT_LOCATION
#include <CGAL/Pm_naive_point_location.h>
#endif
#ifdef USE_WALK_POINT_LOCATION
#include <CGAL/Pm_walk_along_line_point_location.h>
#endif

#include "draw_map.h"  

//typedef CGAL::Point_2<Rep>  Point;
//typedef CGAL::Segment_2<Rep>  Segment;

typedef CGAL::Window_stream  Window_stream;

unsigned int delete_restriction=0;

void redraw(leda_window* wp, double x0, double y0, double x1, double y1) 
{ wp->flush_buffer(x0,y0,x1,y1); }



#ifdef CGAL_PM_TIMER
CGAL::Timer t_total,t_construction,t_insert,t_remove,t_locate,t_vertical,t_split,t_merge;
int n_total=0,n_insert=0,n_remove=0,n_locate=0,n_vertical=0,n_split=0,n_merge=0;
#endif

int test(Planar_map& M,char* name,int argc,char** argv)
{
  Window_stream W(400, 400);
  
  W.button("finish",10);
  W.set_redraw(redraw);
  double x0=-1.1,x1=1.1 ,y0=-1.1; 
  W.init(x0,x1,y0);   // logical window size 
  
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  W.display(leda_window::center,leda_window::center);

#ifdef CGAL_PM_TIMER
  t_construction.reset();
  t_insert.reset();
  t_remove.reset();
  t_locate.reset();
  t_vertical.reset();
  t_split.reset();
  t_merge.reset();
  n_total=0;n_insert=0;n_remove=0;n_locate=0;n_vertical=0,n_split=0,n_merge=0;
#endif
  
  if (argc>1)
    {
#ifdef CGAL_PM_TIMER
      t_construction.start();
#endif

      if (!Init(argv[1],M)) {std::cerr << "\nBad input file";return true;}
#ifdef CGAL_PM_TIMER
      t_construction.stop();
#endif
      win_border(x0,x1,y0,M); //rescale window 
      W.init(x0,x1,y0);
      if (argc>2)
        {
          Planar_map::Locate_type lt;
          std::ifstream in(argv[2]);
          if (in.bad()) {std::cerr << "\nBad locate input file";return false;}
          unsigned int i,n_pt; in >> n_pt;
          //                        bool first=true;
          Planar_map::Point* p=new Planar_map::Point[n_pt];
          
          for (i=0;i<n_pt;i++)
            {
              number_type x,y;
              in >> x >> y ;
              p[i]=Planar_map::Point(x,y);
            }
          do
            {
#ifdef CGAL_PM_TIMER
              t_locate.start();
#endif
              for(i=0;i<n_pt;i++)
                M.locate(p[i],lt);
#ifdef CGAL_PM_TIMER
              t_locate.stop();
              unsigned int k=(unsigned int)((n_locate+n_pt)/BUNDLE)-(unsigned int)(n_locate/BUNDLE);
              for (i=0;i<k;i++) std::cout << "x";
                        std::cout << std::flush;
                        n_locate+=n_pt;
#endif
                      }
#ifdef CGAL_PM_TIMER
                    while(t_locate.time()<WIDE_PRECISION*t_locate.precision());
#else
                    while(false); // don't repeat this loop
#endif
                    
                    if (delete_restriction)
                      {
                        //				Planar_map::Locate_type lt;
                        Planar_map::Halfedge_handle e;
                        unsigned int n = std::min(delete_restriction,M.number_of_halfedges()/2);
                        while (n--)
                          {
                            e=M.halfedges_begin();
#ifdef CGAL_PM_TIMER
                            t_remove.start();
#endif
                            M.remove_edge(e);
#ifdef CGAL_PM_TIMER
                            t_remove.stop();
                            if (!(n_remove%BUNDLE)) std::cout << ">" << std::flush;
                            n_remove++;
#endif
                          }
                      }
                  }
      else
        {
          W.set_redraw(redraw);
          W.set_mode(leda_src_mode);
          W.set_node_width(3);
          W.display(leda_window::center,leda_window::center);
          draw_pm(M,W);		
        }
    }
  else
    {
      W << CGAL::GREEN;
      window_input(M,W);      
      CGAL::Window_stream W2(400, 400);
      W2.init(x0,x1,y0); 
      W2.set_mode(leda_src_mode);
      W2.set_node_width(3);
      W2.button("quit",10);
      W2.display(leda_window::max,leda_window::center);
      draw_pm(M,W2);
    }
  
#ifdef CGAL_PM_TIMER
                if (delete_restriction)
                  {
                    std::cout << "\n" << name << " " << t_construction.time() << " " << 
                      t_insert.time()/n_insert << " " << 
                      t_locate.time()/n_locate << " " << 
                      t_remove.time()/n_remove << " " <<
                      t_vertical.time()/n_vertical << " cilrv" ;
                  }
                else if (argc>2)
                  {
                    std::cout << "\n" << name << " " << t_construction.time() << " " << 
                      t_insert.time()/n_insert << " " << 
                      t_locate.time()/n_locate << " " << 
                      t_vertical.time()/n_vertical << " cilv";
                  }
                else // (argc>1)
                  {
                    std::cout << "\n" << name << t_construction.time() << " " << 
                      t_insert.time()/n_insert << " " << 
                      t_locate.time()/n_locate << " cil";
                  }
#endif
                return true;
}

int main(int argc, char* argv[])
{
#ifdef CGAL_PM_TIMER
	t_total.reset();
	t_total.start();
	t_construction.reset();
	t_insert.reset();
	t_remove.reset();
	t_locate.reset();
	t_vertical.reset();
#endif
	
#ifdef USE_NAIVE_POINT_LOCATION
#ifdef USE_DEFAULT_WITHOUT_REBUILD
#error USE one point location
#endif
#ifdef USE_WALK_POINT_LOCATION
#error  USE one point location
#endif
	CGAL::Pm_naive_point_location<Planar_map> naive_pl;
	Planar_map M(&naive_pl);
	test(M,"naive",argc,argv);
        std::cout << std::endl;
#else 
#if defined(USE_WALK_POINT_LOCATION) 
#ifdef USE_NAIVE_POINT_LOCATION
#error USE one point location
#endif
#ifdef USE_DEFAULT_WITHOUT_REBUILD
#error USE one point location
#endif
	CGAL::Pm_walk_along_line_point_location<Planar_map> pl;
	Planar_map M(&pl);
	test(M,"walk_along_line",argc,argv);
        std::cout << std::endl;
#else
#if defined(USE_DEFAULT_WITHOUT_REBUILD)
#ifdef USE_NAIVE_POINT_LOCATION
#error USE one point location
#endif
#ifdef USE_WALK_POINT_LOCATION
#error  USE one point location
#endif
	CGAL::Pm_default_point_location<Planar_map> pl(false);
	Planar_map M(&pl);
	test(M,"default_without_rebuilding",argc,argv);
	std::cout << std::endl;
#else
	Planar_map M;
        test(M,"default",argc,argv);
	std::cout << std::endl;
#endif
#endif
#endif
	return 0;
}

#endif // CGAL_USE_LEDA
