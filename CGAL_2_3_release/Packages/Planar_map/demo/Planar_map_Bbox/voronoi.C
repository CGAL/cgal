#include <iostream>

#ifndef CGAL_USE_LEDA
int main(int argc, char* argv[])
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include "voronoi.h"

unsigned int delete_restriction=0;
int curve_type,action_type,mbutton,read_queue,display_mode;
Point p,q,lastp,lastq;


Window_stream w(400, 400);
Planar_map* mp;
Triangulation t;
X_curve_container curve_list;
Point_container point_list;

void Triangulation_2_Planar_map_2(const Triangulation& T,Planar_map& m)
{
  X_curve_container curves;
  /*	Triangulation::Vertex_iterator vit, vbegin=T.vertices_begin(), 
        vend=T.vertices_end(); */
  Triangulation::Edge_iterator eit, ebegin=T.edges_begin(), eend=T.edges_end(); 
  for (eit=ebegin; eit != eend; ++eit)
    {
      CGAL::Object o = T.dual(eit);
      /*Triangulation::*/Segment s;
      /*Triangulation::*/Ray r;
      /*Triangulation::*/Line l;
      if (CGAL::assign(s,o)) {
        curves.push_back(s); 
#ifdef CGAL_PMBB_DEBUG
        w << CGAL::BLACK;
        w << s;
#endif
        continue;
      }
      if (CGAL::assign(r,o)) {
        curves.push_back(r); 
#ifdef CGAL_PMBB_DEBUG
        w << CGAL::BLACK;
        w << r;
#endif
        continue;
      }
      if (CGAL::assign(l,o)) {
        curves.push_back(l); 
#ifdef CGAL_PMBB_DEBUG
        w << CGAL::BLACK;
        w << l;
#endif
        continue;
      }
    }
  redraw_status=false;
  m.clear();
  Bounding_box_base* bb = (Bounding_box_base*)m.get_bounding_box();
  bb->insert(point_list.begin(),point_list.end());
  // Ensure that the generating points of the Voronoi diagram are inside 
  // the bounding box.
  m.insert(curves.begin(),curves.end());
  redraw_status=true;
  curves.clear();
}

#ifdef CGAL_PM_TIMER
// Used for time measurements

CGAL::Timer t_total,t_construction,t_insert,t_remove,t_locate,t_vertical,t_split,t_merge;
int n_total=0,n_insert=0,n_remove=0,n_locate=0,n_vertical=0,n_split=0,n_merge=0;

#endif // CGAL_PM_TIMER

bool test(Triangulation& t,Planar_map& m,char* name,int argc,char** argv)
{
  curve_type=3;
  action_type=0;
  //	available=0; // number of points read and not used.
  int type;
  double /*x, y,*/x0=-1.1,x1=1.1 ,y0=-1.1/*,sx,sy,tx,ty*/; 
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
  curve_choice.append("Point");
  display_choice.append("<<");
  display_choice.append("<<;bbox");
  display_choice.append("write");
  display_choice.append("write;bbox");
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
#endif
  
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
              unsigned int k=(unsigned int)((n_locate+n_pt)/BUNDLE)-
                (unsigned int)(n_locate/BUNDLE);
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
              unsigned int n = std::min(delete_restriction,m.number_of_halfedges()/2);
              while (n--)
                {
                  e=m.halfedges_begin();
                  
#ifdef CGAL_PM_TIMER
                  t_remove.start();
#endif
                  
                  m.remove_edge(e);
                  
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
		the planar map (induced by the Voronoi diagram)."); 
  
  while(true)
    {
      w.redraw();
      CGAL::read_mouse_plus(w,p,mbutton);
      if (mbutton==MOUSE_BUTTON(2)) p=find_closest_vertex(m,p)->point();
      w.redraw();
      if (mbutton==EXIT_BUTTON || mbutton==REMOVE_ALL_BUTTON) 
        {
          m.clear();
          t.clear();
          curve_list.clear();
          point_list.clear();
          if (mbutton==EXIT_BUTTON) break;
          continue;
        }
      else if (mbutton==REFRESH_BUTTON) 
        {
          curve_list.clear();
          continue;
        }
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
          switch(curve_type)
            {
            case LINE_TYPE:
            case RAY_TYPE:
            case SEGMENT_TYPE:
              read_line_ray_segment_plus(w,p,q,type,mbutton);
              switch(mbutton)
                {
                case MOUSE_BUTTON(2):
                  q=find_closest_vertex(m,q)->point();
                case MOUSE_BUTTON(1):
                  switch(type)
                    {
                    case LINE_TYPE:
                      l = Line(p,q);
                      redraw_status=false;
                      bb->insert(Segment(p,q));
                      m.insert(l);
                      redraw_status=true;
                      break;
                    case RAY_TYPE:
                      r = Ray(p,q);
                      redraw_status=false;
                      bb->insert(Segment(p,q));
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
                }
              break;
            case POINT_TYPE:
              t.insert(p);
              point_list.push_back(p);
              Triangulation_2_Planar_map_2(t,m);
              break;
            }
          break;
        case REMOVE_ACTION:
          CGAL_warning_msg(curve_type!=REMOVE_ACTION,
                           "Currently remove action is not supported.");
          break;
        case SPLIT_ACTION:
          CGAL_warning_msg(curve_type!=SPLIT_ACTION,
                           "Currently split action is not supported.");
          break;
/*			case MERGE_ACTION:
                        CGAL_warning_msg(curve_type!=MERGE_ACTION,"Currently merge action is not supported.");
                        break;*/
        case LOCATE_ACTION:
          curve_list.clear();
          find_face(m,p,curve_list);
          break;
        case VERTICAL_RAY_SHOOT_ACTION:
          curve_list.clear();
          if (vertical_ray_shoot(m,p,q,true,h)) 
            {
              curve_list.push_back(h->curve());
              lastp=p;
              lastq=q;
            }
          break;
        }
    }
  curve_list.clear();
  point_list.clear();
  t.clear();
  
#ifdef CGAL_PM_TIMER
  // Print time measurements log
  
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
#endif //CGAL_PM_TIMER

  return true;
}


class Options {
public:
  Options()
    :  file_input(false)
    
  {}
  
  char program[100];
  char fname[100];
  bool file_input;
};



void usage(char* program)
{
  std::cerr << "\nNAME\n     "
            << program << " - Triangulation of a point set\n\n";
  std::cerr << "SYNOPSIS\n     "
            << program << " [-file fname]\n";
  
  std::cerr << "\nDESCRIPTION\n"
            << "     Triangulates a point set that comes from a file or stdin.\n";
  std::cerr << "\nOPTIONS\n"
            << "     All options can be abbreviated by their first character\n\n";
}


bool
parse(int argc, char* argv[], Options &opt)
{
  strcpy(opt.program, argv[0]);
  --argc;
  argv++;
  
  while ((argc > 0) && (argv[0][0] == '-')){
    if ((!strcmp(argv[0], "-f")) || (!strcmp(argv[0], "-file"))) {
      strcpy(opt.fname, argv[1]);
      opt.file_input = true;
      argv += 2;
      argc -= 2;
    }
    else if ((!strcmp(argv[0], "-?")) ||
             (!strcmp(argv[0], "-h")) ||
             (!strcmp(argv[0], "-help"))) {
      usage(opt.program);
      return false;
    }
    else {
      std::cerr << "Unrecognized option " << argv[0] << std::endl;
      usage(opt.program);
      return false;
    }
  }
  if(argc > 0){
    std::cerr << "Unrecognized option " << argv[0] << std::endl;
    usage(opt.program);
    return false;
  }
  return true;
}





void input_from_file(Triangulation &T,
                     const Options& opt)
{
  if(! opt.file_input){
    return;
  }
  
  std::ifstream is(opt.fname);
  CGAL::set_ascii_mode(is);
  
  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  
  std::istream_iterator<Point> begin(is);
  std::istream_iterator<Point> end;
  T.insert(begin, end);
  //    Point p;
  //     for(int i=0; i<n; i++) {
  //       is >> p; 
  //       T.insert(p);
  //     }
}


int
main(int argc, char* argv[])
{
    Options opt;
    parse(argc, argv, opt);
    input_from_file(t, opt);
    
    POINT_LOCATION pl POINT_LOCATION_ARGS;
    Traits tr;
    BOUNDING_BOX bbox;
    Planar_map m(tr,&pl,&bbox); mp=&m;
    test(t,m,POINT_LOCATION_NAME,argc,argv);
    std::cout << std::endl;
    
    return 0;
}

#endif // CGAL_USE_LEDA
