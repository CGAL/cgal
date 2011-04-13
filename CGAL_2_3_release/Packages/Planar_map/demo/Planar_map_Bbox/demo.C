#include <iostream>

#ifndef CGAL_USE_LEDA
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include "demo.h"

unsigned int delete_restriction=0;
int curve_type,action_type,mbutton,read_queue,display_mode;
Point p,q,lastp,lastq;

Window_stream w(400, 400);
X_curve_container curve_list;
Point_container point_list;
Planar_map* mp;

#ifdef CGAL_PM_TIMER
// Used for time measurements

CGAL::Timer t_total,t_construction,t_insert,t_remove,t_locate,t_vertical,t_split,t_merge;
int n_total=0,n_insert=0,n_remove=0,n_locate=0,n_vertical=0,n_split=0,n_merge=0;

#endif // CGAL_PM_TIMER

bool test(Planar_map& m,char* name,int argc,char** argv)
{
  curve_type=2;
  action_type=0;
  display_mode=0;
  //	available=0; // number of points read and not used.
  int type;
  //	double x, y,x0=-1.1,x1=1.1 ,y0=-1.1,sx,sy,tx,ty; 
  double x0=-1.1,x1=1.1 ,y0=-1.1; 
  leda_list<leda_string> curve_choice,action_choice,display_choice;
  Bounding_box_base* bb = (Bounding_box_base*)m.get_bounding_box();
  //	Locate_type lt;
  Halfedge_handle h;
  Line l;
  Ray r;
  Segment s;
  
  action_choice.append("Insert");
  action_choice.append("Remove");
  action_choice.append("Split");
  //	action_choice.append("Merge");
  action_choice.append("Locate");
  action_choice.append("Shoot");
  curve_choice.append(UNBOUNDED_X_CURVE);
  curve_choice.append(TARGET_UNBOUNDED_X_CURVE);
  curve_choice.append(BOUNDED_X_CURVE);
  display_choice.append("curve");
  display_choice.append("curve,bbox");
  w.choice_item("Curve type",curve_type,curve_choice);
  w.choice_item("Action type",action_type,action_choice);
  w.choice_item("Display mode",display_mode,display_choice);
  w.button("Refresh",REFRESH_BUTTON);
  w.button("Remove all",REMOVE_ALL_BUTTON);
  w.button("Exit",EXIT_BUTTON);
  w.set_redraw(redraw);
  w.init(x0,x1,y0);   // logical window size 
  w.set_mode(leda_src_mode);
  //	CGAL::cgalize(w);
  w.set_node_width(3);
  w.display(leda_window::center,leda_window::center);
  w << CGAL::GREEN;
  
#ifdef CGAL_PM_TIMER
  
  t_construction.reset();
  t_insert.reset();
  t_remove.reset();
  t_locate.reset();
  t_vertical.reset();
  t_split.reset();
  t_merge.reset();
  n_total=0;n_insert=0;n_remove=0;n_locate=0;n_vertical=0,n_split=0,n_merge=0;
  
#endif // CGAL_PM_TIMER
  
  if (argc>1)
    {
      
#ifdef CGAL_PM_TIMER
      t_construction.start();
#endif
      
      if (!Init(argv[1],m)) {std::cerr << "\nBad input file";return true;}
      
#ifdef CGAL_PM_TIMER
      t_construction.stop();
#endif
      
      win_border(x0,x1,y0,m); //rescale window 
      w.init(x0,x1,y0);
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
          
#ifdef CGAL_PM_TIMER
          t_locate.start();
          do
            {
#endif
              
              for(i=0;i<n_pt;i++)
                m.locate(p[i],lt);
              
#ifdef CGAL_PM_TIMER
              t_locate.stop();
              unsigned int k=(unsigned int)((n_locate+n_pt)/BUNDLE)-(unsigned int)(n_locate/BUNDLE);
              for (i=0;i<k;i++) std::cout << "x";
              std::cout << std::flush;
              n_locate+=n_pt;
            }
          while(t_locate.time()<WIDE_PRECISION*t_locate.precision());
#endif
          
          if (delete_restriction)
            {
				//				Planar_map::Locate_type lt;
              Planar_map::Halfedge_handle e;
              unsigned int n = std::min(delete_restriction,
                                        m.number_of_halfedges()/2);
              while (n--)
                {
                  e=m.halfedges_begin();
                  
#ifdef CGAL_PM_TIMER
                  t_remove.start();
#endif
                  
                  redraw_status=false;
                  m.remove_edge(e);
                  redraw_status=true;
                  
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
          w.set_redraw(redraw);
          w.set_mode(leda_src_mode);
          w.set_node_width(3);
          w.display(leda_window::center,leda_window::center);
          //			draw_pm(m,w);		
        }
    }
  w.acknowledge("Use the buttons above to dynamically update and query \
		the planar map."); 
  
  while(true)
    {
      w.redraw();
      CGAL::read_mouse_plus(w,p,mbutton);
      if (mbutton==EXIT_BUTTON || mbutton==REMOVE_ALL_BUTTON) 
        {
          redraw_status=false;
          m.clear();
          curve_list.clear();
          if (mbutton==EXIT_BUTTON) {w.redraw();break;}
          redraw_status=true;
          continue;
        }
      else if (mbutton==REFRESH_BUTTON) 
        {
          if (action_type!=LOCATE_ACTION &&
              action_type!=VERTICAL_RAY_SHOOT_ACTION)
            curve_list.clear();
          continue;
        }
      curve_list.clear();
      w.redraw();
      /*
        #ifndef USE_LEDA_RAT_KERNEL
        p=Point(FT(x),FT(y));
        #else
        p=Point(leda_rational(x),leda_rational(y));
        #endif
      */
      switch(action_type)
        {
        case INSERT_ACTION:
#ifdef old
          switch(curve_type)
            {
            case 0: // Line
              read_line(x,y,l,m,w);
				//if ((w >> l).state) m.insert(l);
				//read_line(p,l,m,w);
              break;
            case 1: // Ray
              read_ray(x,y,r,m,w);
              if (w >> r) m.insert(r);
				//read_ray(p,r,m,w);
              break;
            case 2: // Segment
              read_segment(x,y,s,m,w);
				//read_segment(p,s,m,w);
              break;
            }
#endif
          if (mbutton==MOUSE_BUTTON(2)) 
            {
              if (m.is_empty()) continue;
              if ((type=curve_type)==LINE_TYPE) continue;
              p=find_closest_vertex(m,p)->point();
            }
          read_line_ray_segment_plus(w,p,q,type,mbutton);
          switch(mbutton)
            {
            case MOUSE_BUTTON(2):
              if (m.is_empty()) break;
              if (type==LINE_TYPE||type==RAY_TYPE) break;
              q=find_closest_vertex(m,q)->point();
            case MOUSE_BUTTON(1):
              if (m.get_traits().point_is_same(p,q)) 
                {
                  std::cerr << "\nAttempt to insert null length curve";
                  break;
                }
              switch(type)
                {
                case LINE_TYPE:
                  l = Line(p,q);
                  // ensure that the segment along l between p and q is in the bounding box
                  redraw_status=false;
                  //bb->insert(Segment(p,q));
                  bb->insert(p);
                  m.insert(l);
                  redraw_status=true;
                  break;
                case RAY_TYPE:
                  r = Ray(p,q);
                  // ensure that the segment along r between p and q is in the bounding box
                  redraw_status=false;
                  //	bb->insert(Segment(p,q));
                  bb->insert(p);
                  m.insert(r);
                  redraw_status=true;
                  break;
                case SEGMENT_TYPE:
                  s = Segment(p,q);
                  redraw_status=false;
                  m.insert(s);
                  redraw_status=true;
                  break;
                }
              break;
            default:
              break;
            }
          break;
        case REMOVE_ACTION:
          if (vertical_ray_shoot(m,p,q,true,h)==Planar_map::EDGE) 
            redraw_status=false;
          m.remove_edge(h);
          redraw_status=true;
          break;
        case SPLIT_ACTION:
          if (vertical_ray_shoot(m,p,q,true,h)==Planar_map::EDGE)
            {
              X_curve cv1,cv2;
              redraw_status=false;
              split_edge(m,h,cv1,cv2,q);
              curve_list.push_back(cv1);
              curve_list.push_back(cv2);
              redraw_status=true;
            }
          break;
				/*			case MERGE_ACTION:
                                                        CGAL_warning_msg(curve_type!=MERGE_ACTION,
                                                        "Currently merge action is not supported.");
                                                        break;
				*/
        case LOCATE_ACTION:
          find_face(m,p,curve_list);
          break;
        case VERTICAL_RAY_SHOOT_ACTION:
          if (vertical_ray_shoot(m,p,q,true,h)!=Planar_map::UNBOUNDED_FACE)
            {
              curve_list.push_back(h->curve());
              CGAL_postcondition(p!=q);
              lastp=p;lastq=q;
            }
          break;
        default:
          std::cerr << "\nUnknown action performed";
          break;
          
        }
    }
  curve_list.clear();
  
#ifdef CGAL_PM_TIMER
  // Print time measurements log
  
  if (delete_restriction)
    {
      std::cout << "\n" << name << " " << t_construction.time() << " " << 
        t_insert.time()/n_insert << " " << 
        t_remove.time()/n_remove << " " <<
        t_locate.time()/n_locate << " " << 
        t_vertical.time()/n_vertical << " ConInRemLocVer" << std::flush;
    }
  else if (argc>2)
    {
      std::cout << "\n" << name << " " << t_construction.time() << " " << 
        t_insert.time()/n_insert << " " << 
        t_locate.time()/n_locate << " " << 
        t_vertical.time()/n_vertical << " ConInsLocVer" << std::flush;
    }
  else // (argc>1)
    {
      std::cout << "\n" << name << t_construction.time() << " " << 
        t_insert.time()/n_insert << " " << 
        t_locate.time()/n_locate << " ConInsLoc" << std::flush;
    }
  
#endif //CGAL_PM_TIMER
  return true;
}

int main(int argc, char* argv[])
{
	POINT_LOCATION pl POINT_LOCATION_ARGS;
	BOUNDING_BOX bbox;
        Traits tr;
	Planar_map m(tr,&pl,&bbox); mp=&m;
	test(m,POINT_LOCATION_NAME,argc,argv);
	std::cout << std::endl;
	
	return 0;
}

#endif // CGAL_USE_LEDA
