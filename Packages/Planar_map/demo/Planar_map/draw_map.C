#include <iostream>

// if LEDA is not installed, a message will be issued in runtime.
#ifdef CGAL_USE_LEDA

#include <CGAL/basic.h>
#include <fstream>
#include <list>
#include "draw_map.h"
#include "read_inp.h"

/*=========================================================================
 * Start of Code
 \*=========================================================================*/


void MarkCcb (const Ccb_halfedge_circulator & b, std::list<Pm_curve>& l)
{
  Ccb_halfedge_circulator first = b, iter = b;
  do
    {
      l.push_back(iter->curve());
      iter++;
    }
  while (iter != first);
}



void draw_arrow(Pm_point p1, Pm_point p2, bool black, CGAL::Window_stream & W)
{
  if (black)
    W << CGAL::BLACK;
  else
    W << CGAL::WHITE;
  
  //float 
  number_type ar_size=(W.xmax()-W.xmin())/20 ;
  
  W << Pm_curve (p1, p2);
#ifndef USE_LEDA_RAT_KERNEL
  W << Pm_curve(p2, Pm_point (p2.x () - ar_size , p2.y () - ar_size));
  W << Pm_curve(p2, Pm_point (p2.x () + ar_size , p2.y () - ar_size));
#else
  W << Pm_curve(p2, Pm_point (p2.xcoord () - ar_size, p2.ycoord () - ar_size));
  W << Pm_curve(p2, Pm_point (p2.xcoord () + ar_size, p2.ycoord () - ar_size));
#endif
}

void Draw (CGAL::Window_stream & W , Planar_map & pm)
{
    W.set_flush( 0 );
    
    Vertex_iterator vi;
    for (vi = pm.vertices_begin (); vi != pm.vertices_end (); vi++)
    {
                W << CGAL::GREEN;
                W << (*vi).point();
    }

    Halfedge_iterator ei;
    for (ei = pm.halfedges_begin (); ei != pm.halfedges_end (); ei++)
    {
      W << CGAL::BLUE;
      W << ((*ei).curve());
    }

    W.set_flush( 1 );
    W.flush();
}

bool VerticalRayShoot (Pm_point p, Pm_point & q, bool up , Planar_map &pm)
{
  Halfedge_handle e;
  Planar_map::Locate_type lt;
  number_type y_min, y_max, y_arrow_tip;

#ifdef CGAL_PM_TIMER
  t_vertical.start();
#endif
  e=pm.vertical_ray_shoot (p,lt, up);
#ifdef CGAL_PM_TIMER
  t_vertical.stop();
  n_vertical++;
#endif
  if (lt!=Planar_map::UNBOUNDED_FACE)
  {
    Pm_point p1 = e->source()->point();
    Pm_point p2 = e->target()->point();
    // make the arrow point upwards and touch the found segment:
#ifndef USE_LEDA_RAT_KERNEL
    if (p1.x() == p2.x())
    {  
      // the found segment is vertical
      y_min = std::min (p1.y(), p2.y());
      y_max = std::max (p1.y(), p2.y());
      y_arrow_tip = (p.y() < y_min) ? y_min : y_max;
      q = Pm_point(p.x(), y_arrow_tip); 
    }
    else
    {
      y_arrow_tip = p2.y()-
        ((p2.x()-p.x())*(p2.y()-p1.y()))/(p2.x()-p1.x());
      q = Pm_point(p.x(), y_arrow_tip);
    }
#else
    if (p1.xcoord() == p2.xcoord())
    {  
      // the found segment is vertical
      y_min = std::min (p1.ycoord(), p2.ycoord());
      y_max = std::max (p1.ycoord(), p2.ycoord());
      y_arrow_tip = (p.ycoord() < y_min) ? y_min : y_max;
      q = Pm_point(p.xcoord(), y_arrow_tip); 
    }
    else
    {
      y_arrow_tip = p2.ycoord()-
        ((p2.xcoord()-p.xcoord())*(p2.ycoord()-p1.ycoord()))/
        (p2.xcoord()-p1.xcoord());
      q = Pm_point(p.xcoord(), y_arrow_tip);
    }
#endif
    return true;
  }
  else
  {
    return false;
  }
}

void FindFace (const Pm_point& p , Planar_map &pm, std::list<Pm_curve>& l)
{
  Planar_map::Locate_type lt;
#ifdef CGAL_PM_TIMER
  t_vertical.start();
#endif
  Halfedge_handle e = pm.vertical_ray_shoot (p, lt, true);
#ifdef CGAL_PM_TIMER
  t_vertical.stop();
  n_vertical++;
#endif
  Face_handle f;
  
  if (lt!=Planar_map::UNBOUNDED_FACE)
    {
      f = e->face ();
    }
  else
    {
      f = pm.unbounded_face ();
    }
  
  Ccb_halfedge_circulator ccb_circ;
  
  if (f->does_outer_ccb_exist ())
    {
      ccb_circ = f->halfedge_on_outer_ccb()->ccb();  
      MarkCcb (ccb_circ,l);
    }
  
  Planar_map::Holes_iterator iccbit;
  for (iccbit = f->holes_begin (); iccbit != f->holes_end (); ++iccbit)
    {
      MarkCcb (*iccbit,l);
    }
  
}




int draw_pm (Planar_map & pm , CGAL::Window_stream & W)
{    

  std::list<Pm_curve> l;
  Pm_point pnt (20, 20);
  Pm_point ar1, ar2;
  int mbutton = 0;
  double x, y;

  std::cerr << "Drawing the map:" << std::endl;
  Draw (W,pm);
  
  std::cerr << "1.Left button: vertical ray shooting." << std::endl;
  std::cerr << "2.Middle button: point location." << std::endl;
  std::cerr << "3.Right button: exit" << std::endl;
  
  while (mbutton != 3)
  {
    int b=W.read_mouse(x,y);
    if (b==10) return 0;
      
    mbutton = -b;
    //      pnt = Pm_point (x, y);
#ifndef USE_LEDA_RAT_KERNEL
    pnt = Pm_point(Rep::FT(x),
                 Rep::FT(y));
#else
    pnt = Pm_point(leda_rational(x),leda_rational(y));
#endif

    draw_arrow (ar1, ar2, false,W);
    if (mbutton == 1)
    {
      ar1 = pnt;
      if (VerticalRayShoot (ar1, ar2, true,pm))
        draw_arrow (ar1, ar2, true,W);
    }
    
    if (mbutton == 2)
    {
      FindFace (pnt,pm,l);
    }
      
    if (mbutton != 0)
    {
      Draw (W,pm);
      W << CGAL::RED;
      for (std::list<Pm_curve>::iterator lit=l.begin(); lit!=l.end(); ++lit)
        W << *lit;
      l.erase(l.begin(),l.end());
    }
  }
 
  return 0;
}

//-------------------------------------------------------------------
bool ReadFile(char *filename, int &num_points, Pm_point* &pnts, 
              int &num_curves, Pm_curve* &cvs )
{
  int j, k;

  std::ifstream is(filename);
  if (is.bad()) return false;
  PM_input<Traits> inp;
  is >> inp;
  is.close();

  num_points = inp.get_num_pnts();
  num_curves = inp.get_num_cvs();
  pnts = new Pm_point[num_points];
  cvs = new Pm_curve[num_curves];

  int i;
  for(i = 0; i < num_points; i++)
  {
      inp.get(i, pnts[i]);
  }

  for(i = 0; i < inp.get_num_cvs(); i++)
  {
      inp.get(i, k, j);
      cvs[i] = Pm_curve(pnts[k], pnts[j]);
  }

  return true;
  
}
//----------------------------------------------------------------

void win_border( double &x0 , double &x1 , double &y0 ,Planar_map &pm)
{
  Vertex_iterator vit = pm.vertices_begin();
#ifndef USE_LEDA_RAT_KERNEL
  x0=x1=CGAL::to_double(( vit->point() ).x());
  y0=CGAL::to_double(( vit->point() ).y());
#else
  x0=x1=vit->point().xcoordD();
  y0=vit->point().ycoordD();
#endif

  while (vit!=pm.vertices_end())
    {
#ifndef USE_LEDA_RAT_KERNEL
      if ( ((*vit).point() ).x() < x0 )
        x0 = CGAL::to_double(( (*vit).point() ).x()) ;
      if ( ( (*vit).point() ).x() > x1 )
        x1 = CGAL::to_double(( (*vit).point() ).x()) ;
      if ( ( (*vit).point() ).y() < y0 )
        y0 = CGAL::to_double(( (*vit).point() ).y()) ;
#else
      if ( vit->point().xcoordD() < x0 )
        x0 = vit->point().xcoordD() ;
      if ( vit->point().xcoordD() > x1 )
        x1 = vit->point().xcoordD() ;
      if ( vit->point().ycoordD() < y0 )
        y0 = vit->point().ycoordD() ;
#endif
      vit++;
    }

  x0=x0-(x1-x0)/2;
  x1=x1+(x1-x0)/3;
  y0=y0-(x1-x0)/4;

  if (x1<=x0)
    std::cerr << "\nIf you are trying to read an input file "
              << "(e.g. from input_files directory),"
              << "\nmake sure to define the "
              << "USE_RATIONAL flag whenever you are reading exact input "
              << "\n(i.e. *.e files), otherwise avoid using this flag."
              << "\nexample: demo input_files\\window.f\n";
  //  CGAL_postcondition(x1>x0);
}

//DEBUG
//bool Init (char *filename , Planar_map & pm, CGAL::Window_stream& W)
bool Init (char *filename , Planar_map & pm)
{
  int num_points, num_curves, i;
  Pm_point *pnts;
  Pm_curve *cvs;

#ifdef CGAL_PM_TIMER
  // ReadFile shouldn't be included in construction time.
  t_construction.stop();
#endif
  if (!ReadFile (filename, num_points, pnts, num_curves, cvs ))
    return false;

#ifdef CGAL_PM_TIMER
  t_construction.start();
#endif

  for (i = 0; i < num_curves; i++)
  {
#ifdef CGAL_PM_DEBUG
    std::cout << "inserting curve: i\n";
    std::cout << cvs[i] << std::endl;
    //     W << cvs[i] ;
#endif
#ifdef CGAL_PM_TIMER
    t_insert.start();
#endif
    pm.insert (cvs[i]);
#ifdef CGAL_PM_TIMER
    t_insert.stop();
    n_insert++;
#endif
  }

  delete[]  cvs;
  delete[]  pnts;

  return true;
}

///////////////////////////////////////////////////////////////////////

CGAL::Window_stream& operator<<(CGAL::Window_stream& os, Planar_map &M)
{
  Halfedge_iterator it = M.halfedges_begin();
  
  while(it != M.halfedges_end()){
    os << (*it).curve();
    ++it;
  }
	
  return os;
}

//function needed for window input
Vertex_handle closest_vertex(Planar_map &M, const Pm_point& p)
{
  Vertex_handle v;
  Vertex_iterator vi = M.vertices_begin();
  if (vi==M.vertices_end()) 
    return vi; 
  else v=vi;
#ifndef USE_LEDA_RAT_KERNEL
  Rep::FT d  = CGAL::squared_distance(p, (*vi).point());
  for (; vi!=M.vertices_end(); ++vi)
    {
      Rep::FT d2 = CGAL::squared_distance(p, (*vi).point());
      if(d2 < d){
        d = d2;
        v = vi;
      }
    }

#else
  for (; vi!=M.vertices_end(); ++vi)
    if(p.cmp_dist(vi->point(),v->point())<0){v = vi;}
#endif
  
  return v; 
}

void window_input(Planar_map & M, CGAL::Window_stream &W )
{
  std::cerr << "1.Left button: start or end edge at mouse position."
            << std::endl;
  std::cerr << "2.Middle button: start or end edge at closest vertex \
from mouse position"
            << std::endl;
  std::cerr << "3.Right button: remove the edge directly above the mouse \
position"
            << std::endl;

  Pm_point p;
  Pm_point first_point;
//  Vertex_handle last_vertex;
	
  bool start_flag = true;
  
  while(1) {
    double x, y;
    int b = W.get_mouse(x,y);
    if (b==10) break;
#ifndef USE_LEDA_RAT_KERNEL
    p = Pm_point(Rep::FT(x), Rep::FT(y));
#else
    p = Pm_point(leda_rational(x),leda_rational(y));
#endif
    
    if (b == MOUSE_BUTTON(1))
    {
      if (start_flag)
      {
        first_point=p;
        start_flag=false;
      }
      else 
      {
        start_flag=true;
#ifdef CGAL_PM_TIMER
        t_insert.start();
#endif
#ifdef CGAL_PM_DEBUG
        Halfedge_handle e=
#endif
              
          M.insert(Pm_curve(first_point,p));
            
#ifdef CGAL_PM_TIMER
        t_insert.stop();
        n_insert++;
#endif
#ifdef CGAL_PM_DEBUG
        Traits_wrap traits=M.get_traits();
        CGAL_postcondition(traits.point_equal(e->source()->point(),
                                              first_point));
        CGAL_postcondition(traits.point_equal(e->target()->point(),p));
#endif
      }
        
      W << M;
    }
    
    else
      if (b==MOUSE_BUTTON(2))
      {
        if (M.number_of_vertices()==0) { //an empty map do nothing
          start_flag=true;
        }
        else {  
          Vertex_handle v=closest_vertex(M,p);
          
          if (start_flag)  { 
            first_point=v->point();
            start_flag=false;
          }
          else //insert fromfirst_point to nearest
          {
#ifdef CGAL_PM_TIMER
            t_insert.start();
#endif
#ifdef CGAL_PM_DEBUG
            Halfedge_handle e=
#endif

              M.insert(Pm_curve(first_point,v->point()));

#ifdef CGAL_PM_TIMER
            t_insert.stop();
            n_insert++;
#endif
#ifdef CGAL_PM_DEBUG
            Traits_wrap traits=M.get_traits();
            CGAL_postcondition(traits.point_equal(e->source()->point(),
                                                  first_point));
            CGAL_postcondition(traits.point_equal(e->target()->point(),
                                                  v->point()));
#endif
            start_flag=true;
          }
        }
        
        W << M;
      }

      else if(b == MOUSE_BUTTON(3))
      {
        start_flag=true;
        Planar_map::Locate_type l;
        Halfedge_handle e;
#ifdef CGAL_PM_TIMER
        t_vertical.start();
#endif
        e=M.vertical_ray_shoot(p,l,true);
#ifdef CGAL_PM_TIMER
        t_vertical.stop();
        n_vertical++;
#endif
        if (l!=Planar_map::UNBOUNDED_FACE)
        {
#ifdef CGAL_PM_TIMER
          t_remove.start();
#endif
          M.remove_edge(e);
#ifdef CGAL_PM_TIMER
          t_remove.stop();
          n_remove++;
#endif
          W.clear();
          W << M;
              
#ifdef CGAL_PM_DEBUG
          std::cout << "\nremove()" << std::flush;
          M.debug();
#endif
        }  
      }
    
    if (!M.is_valid()) {
      std::cerr << "map is not valid - aborting" << std::endl;
      exit(1);
    }
  }   
}

#endif
