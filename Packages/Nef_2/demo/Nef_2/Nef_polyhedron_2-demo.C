#include <CGAL/basic.h>
#ifdef CGAL_USE_LEDA
#include "xpms/nef.xpm"
#include <LEDA/pixmaps/button32/eye.xpm>
#include <LEDA/pixmaps/button32/draw.xpm>
#include <CGAL/LEDA_basic.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/IO/Filtered_extended_homogeneous_Window_stream.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/IO/Nef_polyhedron_2_Window_stream.h>

#define FILTERED_KERNEL
#ifndef FILTERED_KERNEL
typedef CGAL::Extended_homogeneous<leda_integer> EKernel;
#else
typedef CGAL::Filtered_extended_homogeneous<leda_integer> EKernel;
#endif

#if defined(_MSC_VER) || defined(__BORLANDC__)
#define WIN32CONFIG
#endif

typedef  CGAL::Nef_polyhedron_2<EKernel> Nef_polyhedron;
typedef  Nef_polyhedron::Point     Point;
typedef  Nef_polyhedron::Line      Line;
typedef  Nef_polyhedron::Direction Direction;
typedef  Nef_polyhedron::Object_handle Object_handle;
typedef  Nef_polyhedron::Explorer Explorer;
typedef  Nef_polyhedron::Topological_explorer TExplorer;
typedef  Explorer::Vertex_const_handle Vertex_const_handle;
typedef  Explorer::Halfedge_const_handle Halfedge_const_handle;
typedef  Explorer::Face_const_handle Face_const_handle;

#include <LEDA/panel.h>
#include <LEDA/list.h>
#include <LEDA/d_array.h>
#include <LEDA/file_panel.h>
#include <LEDA/file.h>
#include <LEDA/stream.h>

static leda_panel                       main_panel;
static leda_panel                       op_panel;
static CGAL_LEDA_SCOPE::panel_item      nef_menu_item;
static CGAL_LEDA_SCOPE::panel_item      op_item;
static leda_string                      nef1;
static leda_string                      nef2;
static CGAL_LEDA_SCOPE::menu            create_menu;
static int num;
#ifndef FILTERED_KERNEL
static leda_string  dname = "./homogeneous_data";
#else
static leda_string  dname = "./filtered_homogeneous_data";
#endif
static leda_string  fname = "none";
static leda_string  filter = "*.nef";

CGAL::Window_stream* pW;
Nef_polyhedron* pN;
leda_list<leda_string>* pML;
leda_d_array<leda_string,Nef_polyhedron>* pMH;

static leda_string stripped(leda_string s)
{ int i = s.pos(" = "); return s.head(i); }

void win_redraw_handler(leda_window*) 
{ pW->clear(); (*pW) << (*pN); }

void win_del_handler(leda_window*) 
{ leda_panel P("acknowledge");
  P.text_item("");
  P.text_item("\\bf\\blue Do you really want to quit~?");
  P.fbutton("no",0);
  P.button("yes",1);
  if (P.open(main_panel) == 1) exit(0);
}

static int open_panel(leda_panel& p)
{ p.display(main_panel,0,0);
  int res = p.read_mouse();
  p.close();
  return res;
}

static void store_new(const Nef_polyhedron& N, leda_string t)
{ leda_string k = leda_string("N%i",++num) ;
  (*pMH)[k] = (*pN) = N;
  pML->push_front(k+" = "+t);
  win_redraw_handler(pW);
}

static void update_history()
{
  nef1=pML->head();
  nef2=pML->head();
  main_panel.add_menu(nef_menu_item,(*pML));
  op_panel.add_menu(op_item,(*pML));
}

enum { EMPTY=31, FULL, HOPEN, HCLOSED, POPEN, PCLOSED };
enum { FILE_LOAD=111, FILE_SAVE };

void create(int i) 
{
  if (pML->back()=="none") pML->pop_back();
  Line l; Point p;
  std::list<Point> Lp;
  leda_point pd;
  leda_list<leda_point> Lpd;
  CGAL_LEDA_SCOPE::string_ostream sos; CGAL::set_pretty_mode(sos);
  pW->clear();

  switch (i) {
    case HOPEN: 
      pW->message("Insert Half-Space by Line");
      (*pW) >> l; sos << '(' << l << ')' << '\0';
      if ( l.is_degenerate() ) {
        (*pW).acknowledge("Please enter non-degenerate line."); 
        win_redraw_handler(pW);
      } else 
        store_new(Nef_polyhedron(l,Nef_polyhedron::EXCLUDED),sos.str());
      break;
    case HCLOSED: 
      pW->message("Insert Half-Space by Line");
      (*pW) >> l; sos << '[' << l << ']' << '\0';
      if ( l.is_degenerate() ) {
        (*pW).acknowledge("Please enter non-degenerate line.");
        win_redraw_handler(pW);
      } else 
        store_new(Nef_polyhedron(l,Nef_polyhedron::INCLUDED),sos.str());
      break;
    case POPEN: 
      pW->message("Insert Polygon by Point Sequence");
      Lpd = pW->read_polygon();
      forall(pd,Lpd) Lp.push_back(Point(pd.xcoord(),pd.ycoord()));
      sos << '[' << Lp.size() << "-gon"<< ']' << '\0';
      store_new(Nef_polyhedron(Lp.begin(),Lp.end(),
                               Nef_polyhedron::EXCLUDED),sos.str());
      break;
    case PCLOSED: 
      pW->message("Insert Polygon by Point Sequence");
      Lpd = pW->read_polygon();
      forall(pd,Lpd) Lp.push_back(Point(pd.xcoord(),pd.ycoord()));
      sos << '[' << Lp.size() << "-gon"<< ']' << '\0';
      store_new(Nef_polyhedron(Lp.begin(),Lp.end(),
                               Nef_polyhedron::INCLUDED),sos.str());
      break;
    default:
      std::cout << "created nothing\n";
  }
  sos.freeze(0);
  update_history();
  main_panel.redraw_panel();
}

static void read_file(leda_string fn) 
{ 
  std::ifstream in(fn);
  Nef_polyhedron N; in >> N;
  fn.replace_all(".nef","");
  store_new(N,fn);update_history();
}

static bool confirm_overwrite(const char* fname)
{
#if defined(_MSC_VER) || defined(__BORLANDC__)
  return true;
#else
  leda_panel P;
  P.buttons_per_line(2);
  P.text_item("");
  P.text_item(leda_string("\\bf\\blue File\\black %s\\blue exists.",fname));
  P.button("overwrite",0);
  P.button("cancel",1);
  return (P.open() == 0);
#endif
}

static void write_file(leda_string fname)
{ 
  if (is_file(fname) && !confirm_overwrite(fname)) return;
  // win_ptr->set_status_string(" Writing " + fname);
  leda_string nef = stripped(nef1);
  if (pMH->defined(nef)) {
    std::ofstream out(fname);
    out << (*pMH)[nef];
  } else error_handler(1,"Nef polyhedron "+nef+" not defined.");
}

static void file_handler(int what)
{ 
  CGAL_LEDA_SCOPE::file_panel FP(fname,dname);
  switch (what) {
  case FILE_LOAD: FP.set_load_handler(read_file);
                  break;
  case FILE_SAVE: FP.set_save_handler(write_file);
                  break;
  }
  if (filter != "") FP.set_pattern(filter,filter);
  FP.open();
}

void view(int i)
{
  leda_string nef = stripped(nef1);
  if (pMH->defined(nef)) {
    (*pN) = (*pMH)[nef]; win_redraw_handler(pW);
  } else error_handler(1,"Nef polyhedron "+nef+" not defined.");
}

void binop(int i) 
{
  leda_string op1; 
  switch (i) {
   case 30: op1="intersection"; break;
   case 31: op1="union";        break;
   case 32: op1="difference";   break;
   case 33: op1="symmdiff"; break;
   default: break;
  }
  op_panel.set_button_label(111,op1);
  open_panel(op_panel);
  op_panel.flush();
  leda_string arg1 = stripped(nef1), arg2 = stripped(nef2);
  Nef_polyhedron N1 = (*pMH)[arg1];
  Nef_polyhedron N2 = (*pMH)[arg2];
  std::ofstream log("nef-demo.log");
  if ( !log ) CGAL_assertion_msg(0,"no output log nef-demo.log");
  log << 2 << std::endl << N1 << N2 << std::endl;
  log.close();
  switch (i) {
   case 30: *pN = N1*N2; break;
   case 31: *pN = N1+N2; break;
   case 32: *pN = N1-N2; break;
   case 33: *pN = N1^N2; break;
   default: return;
  }
  leda_string descr = op1+"("+arg1+","+arg2+")";
  store_new(*pN,descr);update_history();
  win_redraw_handler(pW);
}

void unop(int i) 
{
  leda_string op, arg = stripped(nef1);
  Nef_polyhedron N = (*pMH)[arg];
  std::ofstream log("nef-demo.log");
  if ( !log ) CGAL_assertion_msg(0,"no output log nef-demo.log");
  log << 1 << std::endl << N << std::endl;
  TRACEN("writing nef to log");
  log.close();
  switch (i) {
   case 40: op="interior";   *pN = N.interior(); break;
   case 41: op="complement"; *pN = N.complement(); break;
   case 42: op="closure";    *pN = N.closure(); break;
   case 43: op="boundary";   *pN = N.boundary(); break;
   default: return;
  }
  leda_string descr = op+"("+arg+")";
  store_new(*pN,descr);update_history();
  win_redraw_handler(pW);
}

void draw(Object_handle h)
{ CGAL::PM_visualizor<TExplorer,EKernel> 
    PMV(*pW,pN->explorer(),pN->EK,
        CGAL::PM_DefColor<TExplorer>(CGAL::RED,CGAL::RED,6,6) );
  leda_drawing_mode prev = pW->set_mode(leda_xor_mode);
  Vertex_const_handle vh; Halfedge_const_handle eh; Face_const_handle fh; 
  if ( CGAL::assign(vh,h) ) PMV.draw(vh); 
  if ( CGAL::assign(eh,h) ) PMV.draw(eh);   
  if ( CGAL::assign(fh,h) ) PMV.draw(fh); 
  pW->set_mode(prev);
}



int main(int argc, char* argv[])
{
  /* Enable debugging by introducing prime into the product:
    RP 3 EP 5/59 PM 7 NefTOP 11 
    Overlayer 13 PointLoc 17 ConstrTriang 19
    SegSweep 23
  */
  SETDTHREAD(113); 
  CGAL::set_pretty_mode ( std::cerr );
  std::cerr << "using " << CGAL::pointlocationversion << std::endl;
  std::cerr << "using " << CGAL::sweepversion << std::endl;

  CGAL::Window_stream W(600,600); pW = &W;
  Nef_polyhedron N_display; pN = &N_display;
  leda_list<leda_string> ML; pML = &ML;
  leda_d_array<leda_string,Nef_polyhedron> MH; pMH = &MH;

  W.init(-CGAL::frame_default,CGAL::frame_default,-CGAL::frame_default);
  W.set_show_coordinates(true);
  W.set_grid_mode(5);
  W.set_node_width(3);
  W.set_redraw(&win_redraw_handler);
  W.display(0,0);

  std::list<Point> Lp;
  const int r=70;
  Lp.push_back(Point(-r,-r)); 
  Lp.push_back(Point(r,-r));
  Lp.push_back(Point(r,r));   
  Lp.push_back(Point(-r,r));

  store_new(Nef_polyhedron(),"empty");
  store_new(Nef_polyhedron(Nef_polyhedron::COMPLETE),"plane");
  store_new(Nef_polyhedron(Line(Point(0,0),Point(0,1)),
            Nef_polyhedron::INCLUDED),"neg y-plane (closed)");
  store_new(Nef_polyhedron(Line(Point(-2,-1),Point(2,1)),
            Nef_polyhedron::EXCLUDED),"above diagonal (open)");
  store_new(Nef_polyhedron(Lp.begin(),Lp.end(),
            Nef_polyhedron::INCLUDED),"square (closed)");

  if ( argc == 2 ) {
    std::ifstream log( argv[1] );
    if ( !log ) 
      CGAL_assertion_msg(0,leda_string("no input log ")+argv[1]);
    int n;
    log >> n;
    for (int i=0; i<n; ++i) {
      Nef_polyhedron Ni;
      log >> Ni;
      store_new(Ni,leda_string(argv[1])+leda_string("%d",i+1));
    }
  }

  nef1=ML.head();nef2=ML.head();
  nef_menu_item = main_panel.string_item("Polyhedra",nef1,ML);
  main_panel.set_window_delete_handler(win_del_handler);
  main_panel.buttons_per_line(10);
  main_panel.set_item_width(300);
  main_panel.set_frame_label("Operations on Nef Polyhedra"); 
  main_panel.set_icon_label("Nef Polyhedra");

  op_item = op_panel.string_item("Polyhedra",nef2,ML);
  op_panel.fbutton("                  ",111);
  op_panel.set_frame_label("Choose second argument");
  op_panel.set_item_width(300);

  char* p_inter  = main_panel.create_pixrect(intersection_xpm);
  char* p_union  = main_panel.create_pixrect(union_xpm);
  char* p_diff   = main_panel.create_pixrect(difference_xpm);
  char* p_exor   = main_panel.create_pixrect(exor_xpm);
  char* p_int    = main_panel.create_pixrect(interior_xpm);
  char* p_compl  = main_panel.create_pixrect(complement_xpm);
  char* p_clos   = main_panel.create_pixrect(closure_xpm);
  char* p_bound  = main_panel.create_pixrect(boundary_xpm);
  char* p_eye    = main_panel.create_pixrect(eye_xpm);
  char* p_draw   = main_panel.create_pixrect(draw_xpm);

  main_panel.button(p_eye,p_eye,"show", 10, view);
  main_panel.button(p_draw,p_draw,"new polyhedron", 11, create_menu);

  main_panel.button(p_inter,p_inter,"intersection", 30, binop);
  main_panel.button(p_union,p_union,"union", 31, binop);
  main_panel.button(p_diff,p_diff,"difference", 32, binop);
  main_panel.button(p_exor,p_exor,"symmetric difference", 33, binop);

  main_panel.button(p_int,p_int,"interior", 40,unop);
  main_panel.button(p_compl,p_compl,"complement", 41,unop);
  main_panel.button(p_clos,p_clos,"closure", 42,unop);
  main_panel.button(p_bound,p_bound,"boundary", 43,unop);

  create_menu.button("half-plane (open)", HOPEN, create);
  create_menu.button("half-plane (closed)", HCLOSED, create);
  create_menu.button("polygon (open)", POPEN, create);
  create_menu.button("polygon (closed)", PCLOSED, create);
  create_menu.button("load from disk",FILE_LOAD, file_handler);
  create_menu.button("save to disk",FILE_SAVE, file_handler);


  leda_drawing_mode dm;
  Point p_down(0,0);
  main_panel.display();
  int x0,y0,x1,y1;
  W.frame_box(x0,y0,x1,y1);
  W.display_help_text("help/nef-demo");
  Object_handle h;
  for(;;) { 
    double x, y;
    int val;
    switch( W.read_event(val, x, y) ) {
    case CGAL_LEDA_SCOPE::button_press_event: 
      p_down = Point(x,y);
      if (val == MOUSE_BUTTON(1)) { 
        std::cerr << "locating " << p_down << std::endl;
        dm = W.set_mode(leda_xor_mode); 
        W<<CGAL::GREEN<<p_down; W.set_mode(dm);
        h = N_display.locate(p_down);
        draw(h);
      }
      if (val == MOUSE_BUTTON(2)) { 
        std::cerr << "shooting down from " << p_down << std::endl;
        dm = W.set_mode(leda_xor_mode); 
        W<<CGAL::GREEN<<p_down; W.set_mode(dm);
//      h = N_display.ray_shoot(p_down,Direction(0,-1),Nef_polyhedron::NAIVE);
        h = N_display.ray_shoot(p_down,Direction(0,-1));
        draw(h);
      }
      if (val == MOUSE_BUTTON(3)) { 
        CGAL::show_triangulation = !CGAL::show_triangulation;
        win_redraw_handler(&W);
      }
      break;
    case CGAL_LEDA_SCOPE::button_release_event: 
      if (val == MOUSE_BUTTON(1))  
#ifndef WIN32CONFIG 
      { dm = W.set_mode(leda_xor_mode);
        W<<CGAL::GREEN<<p_down; W.set_mode(dm); draw(h); }
#else
      { win_redraw_handler(&W); }
#endif
      if (val == MOUSE_BUTTON(2)) 
#ifndef WIN32CONFIG 
      { dm = W.set_mode(leda_xor_mode);
        W<<CGAL::GREEN<<p_down; W.set_mode(dm); draw(h); }
#else
      { win_redraw_handler(&W); }
#endif
      break;
    case CGAL_LEDA_SCOPE::key_press_event:
      if (val ==  CGAL_LEDA_SCOPE::KEY_UP) { // ZOOM IN 
        CGAL::frame_default*=2;
        Nef_polyhedron::Extended_kernel::RT::set_R(CGAL::frame_default);
        int r = CGAL::frame_default+10;
        W.init(-r,r,-r);
        win_redraw_handler(&W);
      }
      if (val ==  CGAL_LEDA_SCOPE::KEY_DOWN) { // ZOOM OUT       
        CGAL::frame_default/=2;
        Nef_polyhedron::Extended_kernel::RT::set_R(CGAL::frame_default);
        int r = CGAL::frame_default+10;
        W.init(-r,r,-r);
        win_redraw_handler(&W);
      }
    default:
      break;
    }
  }
#if !defined(__KCC) && !defined(__BORLANDC__) 
  return 0; // never reached 
#endif
}

#else // CGAL_USE_LEDA

int main() { return 0; }

#endif // CGAL_USE_LEDA

