// if LEDA is not installed, a message will be issued in runtime by demo.C.
#ifdef CGAL_USE_LEDA

#include <CGAL/basic.h>
#include <fstream>
#include "draw_map.h"
#include "read_inp.h"

#ifndef CGAL_IO_PLANAR_MAP_WINDOW_STREAM_H
#include <CGAL/IO/Planar_map_Window_stream.h>
#endif

bool redraw_status=true;

/*=========================================================================
 * Start of Code
 \*=========================================================================*/


void MarkCcb (const Ccb_halfedge_circulator & b, X_curve_container& l)
{
  Ccb_halfedge_circulator first = b, iter = b;
  do
    {
      l.push_back(iter->curve());
      iter++;
    }
  while (iter != first);
}



void draw_arrow (CGAL::Window_stream & w, const Point& p1, const Point& p2, 
                 bool black)
{
  /*
    std::cout << "\nvoid draw_arrow (CGAL::Window_stream & w, " <<
    "const Point& " << p1 << ", const Point& " << p2 << ", bool " << 
    black << ")" << std::flush;
  */
  
  if (black)
    w << CGAL::BLACK;
  else
    w << CGAL::WHITE;
  
  //float 
  number_type ar_size=(w.xmax()-w.xmin())/20 ;
  
  w << X_curve (Segment(p1, p2));
#ifndef USE_LEDA_RAT_KERNEL
  w << X_curve (Segment(p2, Point (p2.x () - ar_size , p2.y () - ar_size)));
  w << X_curve (Segment(p2, Point (p2.x () + ar_size , p2.y () - ar_size)));
#else
  w << X_curve (p2, Point (p2.xcoord () - ar_size , p2.ycoord () - ar_size));
  w << X_curve (p2, Point (p2.xcoord () + ar_size , p2.ycoord () - ar_size));
#endif
  
  /*
    number_type y_min, y_max, y_arrow_tip;
    
    Point p1 = e->source()->point();
    Point p2 = e->target()->point();
    // make the arrow point upwards and touch the found segment:
    
    #ifndef USE_LEDA_RAT_KERNEL
    
    if (p1.x() == p2.x())
    {  
    // the found segment is vertical
    y_min = std::min (p1.y(), p2.y());
    y_max = std::max (p1.y(), p2.y());
    y_arrow_tip = (p.y() < y_min) ? y_min : y_max;
    q = Point(p.x(), y_arrow_tip); 
    }
    else
    {
    y_arrow_tip = p2.y()-
    ((p2.x()-p.x())*(p2.y()-p1.y()))/(p2.x()-p1.x());
    q = Point(p.x(), y_arrow_tip);
    }
    #else //USE_LEDA_RAT_KERNEL
    
    if (p1.xcoord() == p2.xcoord())
    {  
    // the found segment is vertical
    y_min = std::min (p1.ycoord(), p2.ycoord());
    y_max = std::max (p1.ycoord(), p2.ycoord());
    y_arrow_tip = (p.ycoord() < y_min) ? y_min : y_max;
    q = Point(p.xcoord(), y_arrow_tip); 
    }
    else
    {
    y_arrow_tip = p2.ycoord()-
    ((p2.xcoord()-p.xcoord())*(p2.ycoord()-p1.ycoord()))/(p2.xcoord()-p1.xcoord());
    q = Point(p.xcoord(), y_arrow_tip);
    }
    
    #endif //USE_LEDA_RAT_KERNEL
  */
}




void Draw (CGAL::Window_stream & w , Planar_map & m )
{
    w.set_flush( 0 );
    
    Vertex_iterator vi;
    for (vi = m.vertices_begin (); vi != m.vertices_end (); vi++)
    {
		w << CGAL::GREEN;
		w << vi->point();
    }
	
    Halfedge_iterator ei;
    for (ei = m.halfedges_begin (); ei != m.halfedges_end (); ei++)
    {
		w << CGAL::BLUE;
		w << ei->curve();
    }
	
    w.set_flush( 1 );
    w.flush();
	
}



Locate_type vertical_ray_shoot (Planar_map &m,const Point& p,Point & q, 
                                bool up,Halfedge_handle& e)
/* Returns true if the vertical ray shoot eminating from 'p' directed upwards 
   or downwards depending on 'up' intersects the planar map. 
   If this is indeed the case, then 'q' and 'h' are the first intersection 
   point and halfedge along this ray. */
{
  Locate_type lt;

#ifdef CGAL_PM_TIMER
  t_vertical.start();
#endif
  
  e=m.vertical_ray_shoot (p,lt, up);
  
#ifdef CGAL_PM_TIMER
  t_vertical.stop();
  n_vertical++;
#endif
  
  if (lt!=Planar_map::UNBOUNDED_FACE) 
    q=m.get_traits().curve_calc_point(e->curve(),p);
  return lt;
}

bool cooriented(const Traits_wrap* traits,const Halfedge_handle& h,const X_curve& cv)
  // Optimize here ...
{
  CGAL_precondition(&h->curve()==&cv);
  const Point& s=h->source()->point(),&t=h->target()->point();
  if (traits->point_is_same_x(s,t))
    {
      CGAL_precondition(
	traits->curve_get_status(cv)==Traits::CURVE_VERTICAL_UP || 
        traits->curve_get_status(cv)==Traits::CURVE_VERTICAL_DOWN);
      return (traits->point_is_higher(t,s)==
              (traits->curve_get_status(cv)==Traits::CURVE_VERTICAL_UP));
    }
  else
    {
      CGAL_precondition(
	traits->curve_get_status(cv)==Traits::CURVE_RIGHT ||
	traits->curve_get_status(cv)==Traits::CURVE_LEFT);
      return (traits->point_is_right(t,s)==
              (traits->curve_get_status(cv)==Traits::CURVE_RIGHT));
    }
}

void split_edge(Planar_map &m, const Halfedge_handle& h,
                X_curve& cv1,X_curve& cv2,const Point& p)
{
  m.get_traits().split_curve(h->curve(),cv1,cv2,p);
  m.split_edge(h,cv1,cv2);
}

void find_face (Planar_map &m, const Point& p, X_curve_container& l)
{
  Locate_type lt;
  
#ifdef CGAL_PM_TIMER
  t_locate.start();
#endif
  
  Halfedge_handle e = m.locate(p, lt);
  
#ifdef CGAL_PM_TIMER
  t_locate.stop();
  n_locate++;
#endif
  
  Face_handle f;
  
  if (lt!=Planar_map::UNBOUNDED_FACE)
    f = e->face ();
  else
    f = m.unbounded_face ();
  
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

Window& read_line_ray_segment_plus( Window& w, const Point& s, Point& t, 
                                    int& type, int& key)
{
  double x,y;
  key = 0;
  w.state = 1;
  Point p=s;
  int save_but[8];
  w.std_buttons(save_but);
  
  bool buf = w.is_buffering();
  if (buf) w.stop_buffering();
  
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  
  for(;;)
    { 
      switch(type = curve_type)
        {
#if (__LEDA__ >= 420)
        case 0:
          key = w.read_mouse_line(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
          break;
        case 1:
          key = w.read_mouse_ray(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
          break;
#else
        case 0:
          key = w.read_mouse_seg(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
          break;
        case 1:
          key = w.read_mouse_seg(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
          break;
#endif
        case 2:
          key = w.read_mouse_seg(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
          break;
        }	
      if (key == MOUSE_BUTTON(3))  
        { 
          w.state = 0;
          break; 
        }
      
      if (key == MOUSE_BUTTON(1) || key == MOUSE_BUTTON(2)) 
	if (!w.shift_key_down() || !CGAL::read(w,p).state) break; 
      
    }
  
  
  if (w.state) 
    { 
      switch(type)
	{
	case 0:
	  w.draw_line(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
	  break;
	case 1:
	  w.draw_ray(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
	  break;
	case 2:
	  w.draw_segment(CGAL_NTS to_double(p.x()),CGAL_NTS to_double(p.y()),x,y);
	  break;
	}			
    }
  
  w.set_mode(save);
  
  if (buf) w.start_buffering();
  
  w.set_buttons(save_but);
  t=Point(x,y);
  return w;
}

void redraw(leda_window* wp, double x0, double y0, double x1, double y1) 
{ 
	wp->flush_buffer(x0,y0,x1,y1); 
	CGAL_assertion(wp && mp);
	if (redraw_status)
	  {
	    if ((display_mode==OUTPUT_OPERATOR_BBOX_MODE || 
		 display_mode==WRITE_BBOX_MODE)&&
		!mp->get_bounding_box()->is_empty())
	      // draws the boundary in black.
	      {
		
#if (__LEDA__ >= 400)
		leda_line_style style = wp->set_line_style(leda_dotted);
#else
		line_style style = wp->set_line_style(leda_dotted);
#endif
		*wp << CGAL::BLACK << mp->get_traits().get_bounding_box();
		wp->set_line_style(style);
	      }
	    if (display_mode==OUTPUT_OPERATOR_MODE || 
		display_mode==OUTPUT_OPERATOR_BBOX_MODE) 
	      *wp << CGAL::GREEN << *mp;
	    else if (display_mode==WRITE_MODE || 
		     display_mode==WRITE_BBOX_MODE) 
	      CGAL::write(*wp << CGAL::GREEN,*mp);
	    
	    *wp << CGAL::BLUE;
	    for ( Point_iterator lit=point_list.begin();
		  lit!=point_list.end(); ++lit) *wp << *lit; 
	    
	    switch(action_type)
	      {
	      case INSERT_ACTION:
	      case REMOVE_ACTION:
		/*	
			case MERGE_ACTION:
			break;
		*/
	      case VERTICAL_RAY_SHOOT_ACTION:
		// draw the ray shoot as an arrow, and color the hit halfedge 
		// in red.
		if (curve_list.begin()!=curve_list.end()) 
			  draw_arrow (*wp,lastp,lastq,true); 

		  case SPLIT_ACTION:
	      case LOCATE_ACTION:
		{
		  *wp << CGAL::RED;
		if (display_mode==OUTPUT_OPERATOR_MODE || 
			display_mode==OUTPUT_OPERATOR_BBOX_MODE)
		  for ( X_curve_iterator lit=curve_list.begin(); 
			lit!=curve_list.end(); ++lit) 
		    *wp << *lit; // draw the boundary in red.
		else
		  for ( X_curve_iterator lit=curve_list.begin(); 
			lit!=curve_list.end(); ++lit) 
				CGAL::write(*wp,*lit,mp->get_traits());// draw the boundary in red.
		  *wp << CGAL::GREEN;
		}
		break;
	      default:
		std::cerr << "\nUnknown action performed";
		break;
	      }
	  }
}

/*
int draw_pm (Planar_map & m , CGAL::Window_stream & w)
{    
X_curve_container l;
Point pnt (20, 20);
Point ar1, ar2;
int button = 0;
double x, y;

	std::cerr << "Drawing the map:" << std::endl;
	Draw (w,m);
	
	 std::cerr << "1.Left button: vertical ray shooting." << std::endl;
	 std::cerr << "2.Middle button: point location." << std::endl;
	 std::cerr << "3.Right button: exit" << std::endl;
	 
	  while (button != 3)
	  {
	  int b=w.read_mouse(x,y);
	  if (b==10) return 0;
	  
	   button = -b;
	   //      pnt = Point (x, y);
	   
		#ifndef USE_LEDA_RAT_KERNEL
		
		 pnt = Point(FT(x),
		 FT(y));
		 
		  #else //USE_LEDA_RAT_KERNEL
		  
		   pnt = Point(leda_rational(x),leda_rational(y));
		   
			#endif //USE_LEDA_RAT_KERNEL
			
			 draw_arrow (ar1, ar2, false,w);
			 if (button == 1)
			 {
			 ar1 = pnt;
			 if (vertical_ray_shoot (ar1, ar2, true,m))
			 draw_arrow (ar1, ar2, true,w);
			 }
			 
			  if (button == 2)
			  {
			  find_face (pnt,m,l);
			  }
			  
			   if (button != 0)
			   {
			   Draw (w,m);
			   w << CGAL::RED;
			   for (X_curve_iterator lit=l.begin(); lit!=l.end(); ++lit)
			   w << *lit;
			   l.erase(l.begin(),l.end());
			   }
			   }
			   
				
				 return 0;
				 }
				 */
				 
				 //-------------------------------------------------------------------
				 bool ReadFile(char *filename, int &num_points, Point* &pnts, 
					 int &num_curves, X_curve* &cvs )
				 {
					 int j, k;
					 
					 std::ifstream is(filename);
					 if (is.bad()) return false;
					 PM_input<Traits> inp;
					 is >> inp;
					 is.close();
					 
					 num_points = inp.get_num_pnts();
					 num_curves = inp.get_num_cvs();
					 pnts = new Point[num_points];
					 cvs = new X_curve[num_curves];
					 
					 int i;
					 for(i = 0; i < num_points; i++)
					 {
						 inp.get(i, pnts[i]);
					 }
					 
					 for(i = 0; i < inp.get_num_cvs(); i++)
					 {
						 inp.get(i, k, j);
						 cvs[i] = X_curve(Segment(pnts[k], pnts[j]));
						 /* Explicitly assumes that the curves are segments that are determined by 
						 their two end points */
					 }
					 
					 return true;
					 
				 }
				 //----------------------------------------------------------------
				 
				 void win_border( double &x0 , double &x1 , double &y0 ,Planar_map &m)
				 {
					 Vertex_iterator vit = m.vertices_begin();
					 
#ifndef USE_LEDA_RAT_KERNEL
					 
					 x0=x1=CGAL::to_double(( vit->point() ).x());
					 y0=CGAL::to_double(( vit->point() ).y());
					 
#else // USE_LEDA_RAT_KERNEL
					 
					 x0=x1=vit->point().xcoordD();
					 y0=vit->point().ycoordD();
					 
#endif //USE_LEDA_RAT_KERNEL
					 
					 while (vit!=m.vertices_end())
					 {
						 
#ifndef USE_LEDA_RAT_KERNEL
						 
						 if ( ((*vit).point() ).x() < x0 )
							 x0 = CGAL::to_double(( (*vit).point() ).x()) ;
						 if ( ( (*vit).point() ).x() > x1 )
							 x1 = CGAL::to_double(( (*vit).point() ).x()) ;
						 if ( ( (*vit).point() ).y() < y0 )
							 y0 = CGAL::to_double(( (*vit).point() ).y()) ;
						 
#else //USE_LEDA_RAT_KERNEL
						 
						 if ( vit->point().xcoordD() < x0 )
							 x0 = vit->point().xcoordD() ;
						 if ( vit->point().xcoordD() > x1 )
							 x1 = vit->point().xcoordD() ;
						 if ( vit->point().ycoordD() < y0 )
							 y0 = vit->point().ycoordD() ;
						 
#endif //USE_LEDA_RAT_KERNEL
						 
						 vit++;
					 }
					 
					 x0=x0-(x1-x0)/2;
					 x1=x1+(x1-x0)/3;
					 y0=y0-(x1-x0)/4;
					 
					 if (x1<=x0) std::cerr << "\nIf you are trying to read an input file "
						 << "(e.g. from input_files directory),"
						 << "\nmake sure to define the "
						 << "USE_RATIONAL flag whenever you are reading exact input "
						 << "\n(i.e. *.e files), otherwise avoid using this flag."
						 << "\nexample: demo input_files\\window.f\n";
					 //  CGAL_postcondition(x1>x0);
				 }
				 
				 //DEBUG
				 //bool Init (char *filename , Planar_map & m, CGAL::Window_stream& w)
bool Init (char *filename , Planar_map & m)
{
  int num_points, num_curves, i;
  Point *pnts;
  X_curve *cvs;
  
#ifdef CGAL_PM_TIMER
  
  t_construction.stop();
  // the time spent reading the input file shouldn't be included in 
  // construction time.
  
#endif //CGAL_PM_TIMER
  
  if (!ReadFile (filename, num_points, pnts, num_curves, cvs ))
    return false;
  
#ifdef CGAL_PM_TIMER
  
  t_construction.start();
  t_insert.start();
  
#endif //CGAL_PM_TIMER
  
  for (i = 0; i < num_curves; i++)
    {
      
#ifdef CGAL_PM_DEBUG
#ifndef CGAL_PM_TIMER
      
      std::cout << "Inserting curve: " << i << " :" << cvs[i] << std::endl;
      //      w << cvs[i] ;
      
#endif //CGAL_PM_TIMER
#endif //CGAL_PM_DEBUG
      
      m.insert (cvs[i]);
      
    }
  
#ifdef CGAL_PM_TIMER
  
  t_insert.stop();
  n_insert+=num_curves;
  
#endif //CGAL_PM_TIMER
  
  delete[]  cvs;
  delete[]  pnts;
  
  return true;
  
}

///////////////////////////////////////////////////////////////////////
				 
//function needed for window input
Vertex_handle find_closest_vertex(Planar_map &m, const Point& p)
{
  Vertex_handle v;
  Vertex_iterator vi = m.vertices_begin();
  if (vi==m.vertices_end()) 
    return vi; 
  else v=vi;
#ifndef USE_LEDA_RAT_KERNEL
  FT d  = CGAL::squared_distance(p, (*vi).point());
  for (; vi!=m.vertices_end(); ++vi)
    {
      FT d2 = CGAL::squared_distance(p, (*vi).point());
      if(d2 < d){
        d = d2;
        v = vi;
      }
    }
  
#else
  for (; vi!=m.vertices_end(); ++vi)
    if(p.cmp_dist(vi->point(),v->point())<0){v = vi;}
#endif
  
  return v;
} 

/*
void window_input(Planar_map & m, CGAL::Window_stream &w )
{
  std::cerr << "1.Left button: start or end edge at mouse position."<< std::endl;
  std::cerr << "2.Middle button: start or end edge at closest vertex from mouse position" << std::endl;
  std::cerr << "3.Right button: remove the edge directly above the mouse position" << std::endl;
  
  Point p;
  Point first_point;
  //  Vertex_handle last_vertex;
  
  bool start_flag = true;
  
  while(1) {
    double x, y;
    int b = w.get_mouse(x,y);
    if (b==10) break;
#ifndef USE_LEDA_RAT_KERNEL
    p = Point(FT(x),
              FT(y));
#else
    p = Point(leda_rational(x),leda_rational(y));
#endif
    
    if (b == BOUNDED_CURVE || b == UNBOUNDED_CURVE || b == TARGET_UNBOUNDED_CURVE || b == SOURCE_UNBOUNDED_CURVE) status=b;
    else if (b == MOUSE_BUTTON(1))
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
            Bounding_box_base* bb=(Bounding_box_base*)m.get_bounding_box();
            bb->insert(X_bounded_curve(first_point,p));
            
            switch(status)
              {
              case UNBOUNDED_CURVE:
                m.insert(X_unbounded_curve(first_point,p));
                break;
              case TARGET_UNBOUNDED_CURVE:
                m.insert(X_target_unbounded_curve(first_point,p)); 
                // the first parameter is the finite end
                break;
              case SOURCE_UNBOUNDED_CURVE:
                m.insert(X_source_unbounded_curve(first_point,p));
                // the first parameter is the finite end
                break;
              case BOUNDED_CURVE:
              default:
                m.insert(X_bounded_curve(first_point,p));
                break;
              }
            
#ifdef CGAL_PM_TIMER
            t_insert.stop();
            n_insert++;
#endif
          }
        
        w << m;
      }
    
    else
      if (b==MOUSE_BUTTON(2))
        {
          if (m.number_of_vertices()==0) { //an empty map do nothing
            start_flag=true;
          }
          else {  
            Vertex_handle v=closest_vertex(m,p);
            
            if (start_flag)  { 
              first_point=v->point();
              start_flag=false;
            }
            else //insert fromfirst_point to nearest
              {
#ifdef CGAL_PM_TIMER
                t_insert.start();
#endif
                const Point& p = v->point();
                Bounding_box_base* bb=(Bounding_box_base*)m.get_bounding_box();
                bb->insert(X_bounded_curve(first_point,p));
                
                switch(status)
                  {
                  case UNBOUNDED_CURVE:
                    m.insert(X_unbounded_curve(first_point,p));
                    break;
                  case SOURCE_UNBOUNDED_CURVE:
                    m.insert(X_source_unbounded_curve(first_point,p));
                    break;
                  case TARGET_UNBOUNDED_CURVE:
                    m.insert(X_target_unbounded_curve(first_point,p));
                    break;
                  case BOUNDED_CURVE:
                  default:
                    m.insert(X_bounded_curve(first_point,p));
                    break;
                  }
                
#ifdef CGAL_PM_TIMER
                t_insert.stop();
                n_insert++;
#endif
                start_flag=true;
              }
          }
          
          w << m;
        }
    
      else if(b == MOUSE_BUTTON(3))
        {
          start_flag=true;
          Planar_map::Locate_type l;
          Halfedge_handle e;
#ifdef CGAL_PM_TIMER
          t_vertical.start();
#endif
          e=m.vertical_ray_shoot(p,l,true);
#ifdef CGAL_PM_TIMER
          t_vertical.stop();
          n_vertical++;
#endif
          if (l!=Planar_map::UNBOUNDED_FACE)
            {
#ifdef CGAL_PM_TIMER
              t_remove.start();
#endif
              m.remove_edge(e);
#ifdef CGAL_PM_TIMER
              t_remove.stop();
              n_remove++;
#endif
              w.clear();
              w << m;
              
#ifdef CGAL_PM_DEBUG
              std::cout << "\nremove()" << std::flush;
              m.debug();
#endif
            }
          
        }
    
    if (!m.is_valid()) {
      std::cerr << "map is not valid - aborting" << std::endl;
      exit(1);
    }
    
  }
  
}
*/

#endif // CGAL_USE_LEDA
