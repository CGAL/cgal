// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : src/GeoWin/geo_win.c
// package       : GeoWin (1.2.2)
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ============================================================================

#include <LEDA/geowin.h>
#include <LEDA/stream.h>
#include <LEDA/file_panel.h>
#include <LEDA/file.h>
#include <LEDA/set.h>
#include <LEDA/pixmaps/texture.h>
#include <LEDA/bitmaps/button21.h>

#include "geo_localdefs.h"


GEOWIN_BEGIN_NAMESPACE

void init_leda_float_default_types();
void init_leda_rat_default_types();

GeoWin* GeoWin::call_win = NULL;
geo_scene GeoWin::call_scene = NULL;
geo_scene GeoWin::call_inp_scene = NULL;

int GeoWin::projection_mode = 0;

static unsigned char* input_bits [] = {
  point_21_bits,
  line_21_bits,
  circle_21_bits
};

void gw_save_grid_paras(window* w,double& gdist, int& gmode, grid_style& gstyle)
{
  gdist= w->get_grid_dist();
  gmode= w->get_grid_mode();
  gstyle=w->get_grid_style();
}

void gw_restore_grid_paras(window* w,double gdist, int gmode, grid_style gstyle)
{
  w->set_grid_dist(gdist);
  w->set_grid_mode(gmode);
  w->set_grid_style(gstyle);
}

void mouse_pos(long& MASK,GeoWin* gw, bool b)
{
     if (b)
      {
	MASK |= A_OBJECT;
	gw->set_cursor(60);
      }
    else
      {
	MASK = 0;
	gw->set_cursor();
      }
}

void setup_font(window& w, double size)
{     
  int sz = w.real_to_pix(size);

  if (sz<8) sz=8;

  string font_name;
  font_name = string("F%d",sz);
                
  w.set_font(font_name);
}

void geowin_draw_text_object(window&  w, double x1, double x2, double y1, double y2, void* adr, const
geowin_text& gt, bool object_assoc)
{
  typedef  void (*PFCN)(double, double, double, double, void*, bool, double&, double&);

  double sz = gt.get_size();
  double xf =0, yf= 0;
  
  PFCN pf = gt.get_position_fcn();
  if (pf == NULL) { // center in bounding box
   if (object_assoc) { xf = (x1+x2)/2;  yf = (y1+y2)/2; }
  }
  // compute position ...
  else  pf(x1, x2, y1, y2, adr, object_assoc, xf, yf);
  
  
  double xo = gt.get_x_offset();
  double yo = gt.get_y_offset();

  // font parameter ...
  geowin_font_type font_t = gt.get_font_type();
  string fn_user = gt.get_user_font();
  
  int wsz = w.real_to_pix(sz);
  
  string font_name;
  switch (font_t) {
    case roman_font:  font_name = string("T%d",wsz); break;
    case bold_font:   font_name = string("B%d",wsz); break;
    case italic_font: font_name = string("I%d",wsz); break;
    case fixed_font:  font_name = string("F%d",wsz); break;
    case user_font:   font_name = fn_user; break;
  }
  w.set_font(font_name);   
  w.draw_ctext(xo + xf, yo + yf, gt.get_text(), gt.get_color()); 
}


// for ps_file output ...
void geowin_draw_text_object(window&  w, ps_file&  ps, double x1, double x2, double y1, double y2, void* adr, const
geowin_text& gt, bool object_assoc)
{
  typedef  void (*PFCN)(double, double, double, double, void*, bool, double&, double&);

  double sz = gt.get_size();
  double xf =0, yf= 0;
  
  PFCN pf = gt.get_position_fcn();
  if (pf == NULL) { // center in bounding box
   if (object_assoc) { xf = (x1+x2)/2;  yf = (y1+y2)/2; }
  }
  // compute position ...
  else  pf(x1, x2, y1, y2, adr, object_assoc, xf, yf);
  
  double xo = gt.get_x_offset();
  double yo = gt.get_y_offset();

  // font parameter ...
  geowin_font_type font_t = gt.get_font_type();
  string fn_user = gt.get_user_font();
  
  double wsz = (double) w.real_to_pix(sz);  
  
  string font_name;
  switch (font_t) {
    case roman_font:  font_name = string("Times-Roman"); break;
    case bold_font:   font_name = string("Times-Bold"); break;
    case italic_font: font_name = string("Times-Italic"); break;
    case fixed_font:  font_name = string("Helvetica"); break;
    case user_font:   font_name = fn_user; break;
  }
  ps.set_text_font(font_name);   
  double sold = ps.set_font_size(wsz/45); 
  ps.draw_ctext(xo + xf, yo + yf, gt.get_text(), gt.get_color()); 
  ps.set_font_size(sold); 
}




bool geowin_show_options2(GeoScene* scn,panel* p,GeoWin *gw,color& filling_color,
                   bool& high_light,bool fredraw,
		   list<string> fnames, int& choice, geowin_defining_points& hdp )
{
    panel* P = p;
 
    int own_but;
    
    if( !p )
      {
	P = new panel("Scene Options");
	P->set_inf(&own_but);
      }
    
    bool hl = high_light;
    int hd = (int) hdp;
    P->text_item("\\bf\\blue Focus");
    P->bool_item("highlight ?", hl);
    
    // defining points
    list<string> hd_options;
    hd_options.append(string("Show"));
    hd_options.append(string("Hide"));
    hd_options.append(string("Highlight"));
    P->choice_item(string("show/handle def. points?"),hd,hd_options);
    
    string ch_string;
    
    if (! fnames.empty()){
      fnames.push(string("Default"));
      P->text_item("\\bf\\blue Input objects");
      //P->choice_item(string("Input:"),choice,fnames);
      // set ch_string ...
      ch_string = fnames[fnames.get_item(choice)];
      
      P->string_item(string("Input:"),ch_string,fnames);
    }  
    
    bool val=scn->GeoScene::options(P);
    
    if( ! val ) return false;

    high_light = hl;
    hdp = (geowin_defining_points)hd;
    
    // new : set choice ...
    if (! fnames.empty()){
      //cout << ch_string << "\n";
      int counter = 0;
      string siter;
      forall(siter,fnames){
        if (siter==ch_string) choice=counter;
	counter++;
      }
    }
    
    if(!p) gw->redraw();
    return true;
}



#define FCN_COUNT        100

void GeoWin::init_and_set_visible(geo_scene sc)
{
 sc->init_default_values();
 set_visible(sc,true);
 redraw();
}

void GeoWin::call_keyboard_read_object()
{
 void (GeoEditor::*f)(int);
 f = & GeoEditor::scene_read_object;
 geo_scene sc = get_active_scene();
 geo_editor ed = sc ? sc->type_editor() : 0;
 if (ed) (ed->*f)(1);
}

void GeoWin::call_mouse_read_object()
{
 void (GeoEditor::*f)(int);
 f = & GeoEditor::scene_read_object;
 geo_scene sc = get_active_scene();
 geo_editor ed = sc ? sc->type_editor() : 0;
 if (ed) (ed->*f)(0);
}

void GeoWin::set_quit_handler(bool (*f)(const GeoWin& gw))
{
 quit_handler=f;
}
  
void GeoWin::set_done_handler(bool (*f)(const GeoWin& gw)) 
{
 done_handler=f;
}

bool GeoWin::done_button_was_pressed()
{ 
 if (! (button_info[GEO_CM_OK] == 2)) return false;
  
 if (done_handler != NULL){
       bool flag=done_handler(*this);
       if (! flag) return false;
 }
 
 return true; 
} 

bool GeoWin::quit_button_was_pressed() 
{ 
 int def_choice = 0;
 if (button_info[GEO_CM_QUIT] != 2) return false;
 else { // really quit ??
    if (quit_handler != NULL){
       bool flag=quit_handler(*this);
       if (! flag) return false;
    }
 
    panel P;

    P.text_item("");
    P.text_item("\\bf\\blue Really quit ?");
    P.text_item("");
    P.choice_item("Options ~/.geowinrc",def_choice,"keep","write","reset");
    P.button("yes");
    P.button("no");

    int w= open_panel(P);

    if (w==0) {
      if (def_choice != 0) 
      { if (def_choice == 2) reset_defaults();
        save_defaults();
      }   
      return true;
    }  
    else { quit_button_reset(); return false; }
 } 
} 

GeoWin::GeoWin(const char* label)
 : mouse_actions(0) , scene_name_map(0)
{
  Wp = new window(560,630,label);
  construct();  
  editables = GeoScenePrototypes::get_all_editables();
}

GeoWin::GeoWin(int w, int h, const char* label)
  : mouse_actions(0) , scene_name_map(0)
{
  Wp = new window(w,h,label);
  construct();  
  editables = GeoScenePrototypes::get_all_editables(); 
}


GeoWin::~GeoWin()
{  
  //cout << "Destruktor GeoWin!\n"; cout.flush();
  geo_scenegroup g;
  forall(g,scene_groups){
    delete g;
  }

  geo_scene sc;
  forall(sc,scenes) { 
    sc->set_visible(false);
    sc->del_objref();
  }

  is_locked = true;
  while( cur )  destroy(cur);
  
  for( int i = 0; i <  menu_fcn_count; i++ )  delete menu_fcn[i];
  
  menu_fcn_count = 0;
  if (menu_fcn) delete[] menu_fcn;
  
  for( int j = 0; j <  last_menu; j++ )  delete menus[j];
  delete[] menus;
  
  if (visib != NULL) delete visib;

#if !defined(__win32__)
  delete Wp;
  delete msg_win; 
  if (quick_access != NULL) delete quick_access;
#endif
 
  if( pin_point ) delete pin_point;

}

// *************************** Help *************************************

void GeoWin::help_user()
{ window* wp = window::get_call_window();
  int b = window::get_call_button();
  string label = wp->get_button_label(b);
  wp->display_help_text(label);
}


void GeoWin::help_about()
{
  string naeher = string("S. N\344her");
  string baesken = string("M. B\344sken");
  panel P;

P.text_item(string("\\bf GeoWin (Version %f)",version()));
P.text_item("");
P.text_item("");
P.text_item("by U. Bartuschka, " + baesken + " and " + naeher + " (1998-2000)");
P.text_item("");
P.text_item("");
P.text_item("");
P.text_item("\
GeoWin is a LEDA data type for constructing, displaying \
and manipulating geometrical objects and visualizing \
geometrical algorithms. \
Please see the LEDA User Manual for a detailed description.");
P.text_item("");
P.text_item("");
P.text_item(string("The design and implementation of GeoWin was \
supported in part by the GALIA project."));

P.text_item("");
P.button("ok");

 open_panel(P);
}

void GeoWin::help_news()
{
 panel P;
 P.text_item(string("\\bf New in Version %f",version()));
 P.text_item("");
 P.text_item("- Text objects");
 P.text_item("");
 P.text_item("- New I/O options");
 P.text_item("");
 P.text_item("- Z - coordinate handling changed");
 P.text_item("");
 P.text_item("- lots of small changes"); 
 P.button("ok");
 open_panel(P);
}


// ***********************************************************************
// **************** Operations for handling default values ************

bool GeoWin::save_defaults(string fname)
{
 ofstream out(fname);
 if (out.fail()) return false;

 out << "color              " << DEFAULT_color       << endl;
 out << "color2             " << DEFAULT_color2      << endl;
 out << "fill_color         " << DEFAULT_fill_color  << endl;
 out << "line_width         " << DEFAULT_line_width  << endl;
 out << "active_line_width  " << DEFAULT_active_line_width << endl;
 out << "line_style         " << (int)DEFAULT_line_style  << endl;
 out << "point_style        " << (int)DEFAULT_point_style << endl;

 out << "show_coords        " << show_coords << endl;
 out << "time_meas          " << time_meas   << endl;
 out << "gw_elim            " << gw_elim    << endl;
 out << "gw_solid           " << gw_solid    << endl;
 out << "gw_anim_speed      " << gw_anim_speed  << endl;
 out << "anim_steps         " << gw_anim_steps  << endl;
 
 return true;
}

bool GeoWin::read_defaults(string fname)
{
  if (! is_file(fname)) return false;
    
  ifstream in(fname);
  if (! in.good()) return false;
  
  string key;
  int w;
  color cl;
  bool b;
  
  while (in >> key)
  { if (key == "color") {in >> cl; DEFAULT_color=cl; }
    if (key == "color2"){in >> cl; DEFAULT_color2=cl; }  
    if (key == "fill_color") {in >> cl; DEFAULT_fill_color=cl; }  
    if (key == "line_width") {in >> w; DEFAULT_line_width=w; }
    if (key == "active_line_width") {in >> w; DEFAULT_active_line_width=w; }
    if (key == "line_style") {in >> w; DEFAULT_line_style=(line_style)w; }
    if (key == "point_style"){in >> w; DEFAULT_point_style=(point_style)w; }
    if (key == "show_coords"){in >> b; DEFAULT_show_coords=b; }  
    if (key == "time_meas") {in >> b; DEFAULT_time_meas=b; }    
    if (key == "gw_elim") {in >> b; DEFAULT_gw_elim=b; }    
    if (key == "gw_solid") {in >> b; DEFAULT_gw_solid=b; }     
    if (key == "gw_anim_speed") {in >> w; DEFAULT_gw_anim_speed=w; }    
    if (key == "anim_steps") {in >> w; DEFAULT_anim_steps=w; }        
   }
  return in.good();
}

bool GeoWin::save_defaults()
{ char* home = getenv("HOME");
  string fname;
  if (home) fname = string("%s/",home);
  return save_defaults(fname + ".geowinrc");
}

bool GeoWin::read_defaults()
{ 
  reset_defaults();
  char* home = getenv("HOME");
  string fname;
  if (home) fname = string("%s/",home); else return false;
  return read_defaults(fname + ".geowinrc");
}


void GeoWin::reset_defaults()
{
 // scene defaults
 DEFAULT_color=black;
 DEFAULT_color2=red;
 DEFAULT_fill_color=ivory;
 DEFAULT_line_width=1;
 DEFAULT_active_line_width=1;
 DEFAULT_line_style=solid;
 DEFAULT_point_style=cross_point;
 
 // geowin defaults
 DEFAULT_show_coords=true;
 DEFAULT_time_meas=false;
 DEFAULT_gw_elim=false;
 DEFAULT_gw_solid=false;
 DEFAULT_gw_anim_speed=10;
 DEFAULT_anim_steps=10;
 
 show_coords=true;
 time_meas=false;
 gw_elim=false;
 gw_solid=false;
 gw_anim_speed=10;
 gw_anim_steps=10;
}

// ***********************************************************************
// ***************************  Main  Operations *************************

void GeoWin::construct()
{  
  read_defaults();
  
  init_leda_rat_default_types();
  init_leda_float_default_types();
   
  update_while_object_dragging = true;

  input_mode = 0;
  menu_init=false;
  quit_handler=NULL;
  done_handler=NULL;
  algo_create_scenes = true;
  mouse_pos_check = true;

  Wp->init(-50,450,-50);
  START_XMIN=-50;
  START_XMAX=450;
  START_YMIN=-50;
  START_YMAX=get_ymax();

  Wp->set_inf(this);
  Wp->set_redraw(GeoWin::redraw_geowin);
  show_coords = DEFAULT_show_coords;        //true;
  time_meas   = DEFAULT_time_meas;          //false;
  Wp->set_show_coord_handler(GeoWin::coord_handler);
  scribble_mode=false;
  is_open   = false;
  is_locked = false;

  gw_elim = DEFAULT_gw_elim;                //false;
  gw_solid = DEFAULT_gw_solid;              //false;
  gw_anim_speed = DEFAULT_gw_anim_speed;    //10
  
  gw_d3edges = true;
  gw_d3_cs = true;
  gw_d3_zvalue = 0.0;
  
  gw_anim_steps = DEFAULT_anim_steps;       //10;
  
  msg_win = new window(500, 200, "Messages");
  msg_win->set_inf(this);

  cur_info = string("");

  menus = new menu*[last_menu];
  for( int i = 0; i < last_menu; i++ )  
    {
      menus[i] = new menu;
      menus[i]->set_inf(this);
    }
  menu_fcn = new GeoFunctions*[FCN_COUNT];
  menu_fcn_count = 0;
  user_fcn_count = 0;
  
  Wp->set_node_width(1);
  
  vis_buttons = false;
  activate_buttons = false;
  
  quick_access = NULL;
  show_quick_access_window = false;
  
  visib = NULL;
  
  show_grid = false;
  allow_moving_from_grid = false;
  grid = 0.0;
  
  cur = edit_scene = 0;  
  reset_actions();

  d3_win = NULL;

  bg_pixmap = "none";
  pin_point = 0;
  mouse_fct = 4;
  VPICS = NULL;
  
  replace_by_blanks = "";
  
  file_existence_test = true;
} 

GEOWIN_END_NAMESPACE

#if !defined(__win32__)

static void quick_access_action(int n)
{ 
  window* wp =  window::get_call_window();
#if defined GEOWIN_USE_NAMESPACE
  GEOWIN_NAMESPACE_NAME::GeoWin* gw = (GEOWIN_NAMESPACE_NAME::GeoWin*) wp->get_inf();
  GEOWIN_NAMESPACE_NAME::geo_scene it,sel=gw->get_quick_access_scene(n);
  // now get selected edit scene and activate it ...
  list<GEOWIN_NAMESPACE_NAME::geo_scene> scenes = gw->get_scenes();  
#else  
  GeoWin* gw =  (GeoWin*) wp->get_inf();
  geo_scene it,sel=gw->get_quick_access_scene(n);  
  // now get selected edit scene and activate it ...
  list<geo_scene> scenes = gw->get_scenes();
#endif   
  if (sel==NULL) return;

  forall(it,scenes){
    if (it==sel) {
      gw->activate(sel);
    }
  }
}

#endif

static void visibility_action(int n)
{
  window* wp =  window::get_call_window();
#if defined GEOWIN_USE_NAMESPACE
  GEOWIN_NAMESPACE_NAME::GeoWin* gw = (GEOWIN_NAMESPACE_NAME::GeoWin*) wp->get_inf();
#else  
  GeoWin* gw =  (GeoWin*) wp->get_inf();
#endif
  int i=0;
  int MS = sizeof(int)*8;
  if (gw->get_scenes().size() < 32) MS=gw->get_scenes().size();
  
  for(;i<MS;i++) {
    if (n & (1 << i)) (gw->visib)[i] = true; 
    else (gw->visib)[i] = false;
  }
   
  gw->change_visibility(gw->visib);
}


GEOWIN_BEGIN_NAMESPACE

geo_scene GeoWin::get_quick_access_scene(int n)
{
  list_item it = quick_access_scenes.get_item(n);
  if (it==NULL) return NULL;
  else return quick_access_scenes[it];
}

void GeoWin::init_quick_access_window()
{
  if (quick_access != NULL) { 
    quick_access->close();
    delete quick_access;
    quick_access=NULL;
  }

  if (! show_quick_access_window) return;
  

  int size = 1;
  geo_scene sc;
  quick_access_scenes.clear();
  
  if (cur) {
    forall(sc, scenes){
     if( (sc->owner).empty() ){
	quick_access_scenes.append(sc);
	size=size+30;
     }
    }
  }
  
  // Quick access window ...
  quick_access = new window(250,size+40, "Scene activation and visibility");
  quick_access->set_inf(this);
  quick_access->buttons_per_line(1);
  quick_access->set_item_height(25);
  quick_access->set_panel_bg_color(grey1);
  quick_access->set_bg_color(grey1);
  quick_access->set_button_space(5);
  
  // get active scene ...
#if !defined(__win32__)
  int number = 0;
  geo_scene sc_act = get_active_scene();
  
  forall(sc,quick_access_scenes) {
   if (sc != sc_act)
     quick_access->button(sc->get_name(),number,quick_access_action);
   else
     quick_access->button("->  "+sc->get_name(),number,quick_access_action);
        
   number++;
  } 
#endif  
  
  // scene visibility ...
   if (scenes.size() > 0) { 
    int i=0;
    VISIB_CHOICE=0;
    list<string> LS;
   
    if (visib!=NULL) delete visib;
    visib = new bool[scenes.size()];   
   
    forall(sc , scenes)  
    {
       if (sc->get_visible()) { VISIB_CHOICE = VISIB_CHOICE + (1 << i);  } 
       string sn = get_short_name(sc);
       if (sn == "") LS.append(string("%d",i));
       else LS.append(sn);
       i++;
    }   
    quick_access->choice_mult_item("Scene visibility",VISIB_CHOICE,LS,visibility_action,NULL);  
   }
 
#if !defined(__win32__) 
   quick_access->text_item("\\bf\\blue Scene activation");
#endif  
   
   quick_access->open();
}

double GeoWin::get_xmin() const
{
  return Wp->xmin();
}
  
double GeoWin::get_ymin() const
{
  return Wp->ymin();
}
  
double GeoWin::get_xmax() const
{
  return Wp->xmax();
}
    
double GeoWin::get_ymax() const
{
  return Wp->ymax();
}

GEOWIN_END_NAMESPACE


static void select_visible(int mask)
{
  window* wp =  window::get_call_window();
#if defined GEOWIN_USE_NAMESPACE
  GEOWIN_NAMESPACE_NAME::GeoWin* gw = (GEOWIN_NAMESPACE_NAME::GeoWin*) wp->get_geowin();
  GEOWIN_NAMESPACE_NAME::geo_scene it;
#else  
  GeoWin* gw =  (GeoWin*) wp->get_geowin();
  geo_scene it;
#endif  
  
  int i=0;
    
  forall(it,gw->VAR_1) {
    gw->set_visible(it, (mask&(1<<i)) !=0 ); i++; 
  }  
  i=0;
  forall(it,gw->VAR_2) {
    gw->set_visible(it, (mask&(1<<i)) !=0 ); i++; 
  }    
  gw->init_quick_access_window();
  
  gw->redraw();
}

static void select_input(int mask)
{
  window* wp =  window::get_call_window();
#if defined GEOWIN_USE_NAMESPACE
  GEOWIN_NAMESPACE_NAME::GeoWin* gw = (GEOWIN_NAMESPACE_NAME::GeoWin*) wp->get_geowin();
  GEOWIN_NAMESPACE_NAME::geo_scene cur;
#else  
  GeoWin* gw =  (GeoWin*) wp->get_geowin();
  geo_scene cur;
#endif  
  cur = gw->get_active_scene();
  if (cur==NULL) return;
  else {
    cur->set_alternative_input_fcn(mask);
  }
}

static void select_active(int mask)
{
  //cout << "select active - " << mask << "\n";
  window* wp =  window::get_call_window();
#if defined GEOWIN_USE_NAMESPACE
  GEOWIN_NAMESPACE_NAME::GeoWin* gw = (GEOWIN_NAMESPACE_NAME::GeoWin*) wp->get_geowin();
#else  
  GeoWin* gw =  (GeoWin*) wp->get_geowin();
#endif  
  list_item lit = gw->ACT_SCN.get_item(mask);
  gw->activate(gw->ACT_SCN[lit]);
}

GEOWIN_BEGIN_NAMESPACE

static void avoid_unused_warning(unsigned char**, const char**) {}

void GeoWin::display(int xc, int yc)
{
  typedef unsigned char* CHAR_PTR;

  avoid_unused_warning(xbm_button21,name_xbm_button21);

  if( is_open ) return;
  
  if (!menu_init) init_menu();
  
  // visibility buttons ...
  if (vis_buttons){
    // show a list of buttons ...
    geo_scene it;
    
    VIS_STATE = 0;
    int i = 0;
    
    forall(it,VAR_1) {
      if (get_visible(it)) { VIS_STATE = VIS_STATE + (1 << i);  } 
      i++;
    }
    i=0;
    forall(it,VAR_2) {
      if (get_visible(it)) { VIS_STATE = VIS_STATE + (1 << i);  } 
      i++;
    }       
    
    if ((!VAR_2.empty()) && (!VPIC_LIST.empty())){
     // construct VPICS
     VPICS = new CHAR_PTR[VPIC_LIST.size()];
     unsigned char* iter;
     int i = 0;
     forall(iter,VPIC_LIST) { VPICS[i]=iter; i++; }
    }
    
    if (! VAR_1.empty()) pnit1 = Wp->choice_mult_item(string("  "), VIS_STATE , VLABEL, select_visible); 
    if (! VAR_2.empty()) pnit2 = Wp->choice_mult_item(string("  "), VIS_STATE , VAR_2.size(), VWIDTH, VHEIGHT, VPICS, select_visible);
  }
  
  // activate buttons ...
  if (activate_buttons){
    if (scenes.empty()) ; // no scenes ....
    else {
      geo_scene it;
      ACT_STATE = 0;
      int i=0;
    
      forall(it,ACT_SCN){
       if (it == scenes.head()) { ACT_STATE = ACT_STATE + (1 << i); }
       i++;
     }

     if ((!ACT_SCN.empty()) && (!APIC_LIST.empty())){
      // construct ACT_PICS
      APICS = new CHAR_PTR[APIC_LIST.size()];
      unsigned char* iter;
      int i = 0;
      forall(iter,APIC_LIST) { APICS[i]=iter; i++; }
     }    
    
     acit = Wp->choice_item(string("  "), ACT_STATE, ACT_SCN.size(), VWIDTH, VHEIGHT, APICS, select_active);
    }
  }
  
  // input object choice_item ...
  
  INP_STATE=0;
  inpit = Wp->choice_item(" ",INP_STATE,3,21,21,input_bits, select_input);
  
  // set the x-position of the item ...
  //double max_x = Wp->xmax();
  //cout << "max x:" << max_x << "\n";
  //Wp->set_item_xpos(inpit,max_x-30);
  
  
  Wp->disable_item(inpit);
  Wp->set_clear_on_resize(0);

  Wp->geowin_ptr = this;
 
  Wp->open(xc,yc);
  Wp->open_status_window();
  is_open = true;
  
  if (! show_grid) grid = 10/Wp->scale();
  
  if( !scenes.empty() ) 
    {
      activate(scenes.head());
      no_scene_off();
    }
  else no_scene_on();
  make_scene_menu();
  redraw();
}

int GeoWin::set_edit_mode(geo_scene sc, int emode)
{
  int prev = sc->edit_mode;
  sc->edit_mode = emode;
  return prev;
}
  
int GeoWin::get_edit_mode(geo_scene sc) const
{ return sc->edit_mode; }

void GeoWin::set_label(geo_scene sc, string label)
{
  VAR_1.append(sc);
  VLABEL.append(label);
  vis_buttons = true;
}

void GeoWin::set_bitmap(geo_scene sc, unsigned char* bm)
{
  VAR_2.append(sc);
  VPIC_LIST.append(bm);
  vis_buttons = true;
}

void GeoWin::set_activate_bitmap(geo_scene sc, unsigned char* bm)
{
  ACT_SCN.append(sc);
  APIC_LIST.append(bm);
  activate_buttons = true;
}

int GeoWin::set_button_width(int w)
{
  int old=VWIDTH;
  VWIDTH=w;
  return old;
}

int GeoWin::set_button_height(int h)
{
  int old=VHEIGHT;
  VHEIGHT=h;
  return old;
}

void GeoWin::add_scene_buttons(const list<geo_scene>&  LS, const list<string>& LN)
{
  if ((LS.size() > 32) || (LS.size() != LN.size())) return; 
  VAR_1 = LS;
  VLABEL = LN;
  VPICS = NULL;
  
  vis_buttons = true;
}

void GeoWin::add_scene_buttons(const list<geo_scene>& LS, int w, int h, unsigned char** bm)
{
  if (LS.size() > 32) return;
  VAR_2 = LS;
  VWIDTH  = w;
  VHEIGHT = h;
  VPICS = bm;

  vis_buttons = true;  
}

list<geo_scene> GeoWin::get_visible_scenes() const
{
  list<geo_scene> LS;
  geo_scene it;
  
  forall(it,scenes)
   if (get_visible(it)) LS.append(it);
   
  return LS;
}

int GeoWin::open_panel(panel& P)
{ window& W = get_window();
  W.disable_panel();
  int but = P.open(W);
  W.enable_panel(); 
  redraw_inp_panel();
  return but;
}

int GeoWin::open_panel(file_panel& P)
{ window& W = get_window();
  W.disable_panel();
  P.open();
  W.enable_panel(); 
  redraw_inp_panel();
  return 0;
}

void GeoWin::round_to_grid(double&x, double& y)
{
  //cout << "Round to grid !\n";
  if (grid==0) return;
  
  x=grid*(int)(x/grid+((x>0)?0.5:-0.5));
  y=grid*(int)(y/grid+((y>0)?0.5:-0.5));
  
}

void GeoWin::edit()
{
  if (!is_open) display();
  float tmloop, tact;
  
  init_quick_access_window();
  
  while( !quit_button_was_pressed() && !done_button_was_pressed() )
  { 
    long MASK = 0;

    int but;
    window* w;
    double x,y,xval;
    double xreal=0,yreal=0;
    
    int event = read_event(w, but, x, y);
    if (show_grid) round_to_grid(x,y);
    
    xval = Wp->pix_to_real(1);

    if( xval >= 1 ) 
      set_rational_precision(1);
    else
      set_rational_precision((int)(1/xval)); 
      
    //cout << "    " << x << "/" << y << "  event:" << event << " window:" << w << "\n";
    tact=0;

    while ( event != button_press_event || w != Wp ) 
    {
      event = read_event(w,but,x,y);
      if (show_grid) { xreal=x; yreal=y; round_to_grid(x,y); }    
      //cout << "  " << x << "  " << y << "\n";

      if (w != Wp || event == no_event) continue;
          
      if (cur && mouse_pos_check && (tact< 0.1)) { 
	 tmloop = used_time();
         if (show_grid) { cur->mouse_at_pos(xreal, yreal, MASK); }
         else 
	 cur->mouse_at_pos(x, y, MASK);
	 tact = used_time(tmloop);
	 //if (tact>0.1) cout << "Achtung!\n";
      }
      else 
         Wp->set_cursor(-1);
     }

    switch (but) 
    { case MOUSE_BUTTON(1) : MASK |= A_LEFT;   break;
      case MOUSE_BUTTON(2) : MASK |= A_MIDDLE; break;
      case MOUSE_BUTTON(3) : MASK |= A_RIGHT;  break;
      default : continue;
     } 


    if (scribble_mode) 
    { if (but==MOUSE_BUTTON(1))
      { do { if (scribble_first) 
              { Wp->draw_pixel(x,y, black); 
                scribble_first = false; 
               }
             else 
               Wp->draw_segment(x,y,scribble_x,scribble_y, black);
             scribble_x = x; 
             scribble_y = y;
             event = Wp->get_event(but,x,y);
        } while (event != button_release_event);

        scribble_first=true;
      }
      continue;
    } 


    if (Wp->shift_key_down()) { MASK |= A_SHIFT; }
    if (Wp->ctrl_key_down())  { MASK |= A_CTRL; }
    if (Wp->alt_key_down())   { MASK |= A_ALT; }
    
    point p(x,y);
   
    geo_action func = mouse_actions[MASK|A_IMMEDIATE];

    if (func)
       MASK |= A_IMMEDIATE;
    else
      { bool dclick = false;
        bool drag = false;
  
        if ( but == MOUSE_BUTTON(1) ||
             but == MOUSE_BUTTON(2) ||
             but == MOUSE_BUTTON(3) )
          Wp->read_mouse(x,y,300,300,dclick,drag);
  
        if (drag)   { MASK |= A_DRAG;  if (show_grid && allow_moving_from_grid) p = point(xreal,yreal); }
        if (dclick) MASK |= A_DOUBLE;
  
        func = mouse_actions[MASK];
      }
      
    if (show_grid) { Wp->set_grid_dist(-Wp->get_grid_dist()); }
    
    Wp->disable_panel(false);
    if (func) func(*this, p);
    Wp->enable_panel();
    redraw_inp_panel();
    
    if (show_grid) { Wp->set_grid_dist(-Wp->get_grid_dist()); }
  }
}


bool GeoWin::edit(geo_scene sc)
{
  if( !sc ) return false;
  if( sc->gw != this )
    error_handler(1, "edit : Scene does not belong to this GeoWin!");
  
  if( !is_open) display(); 
  
  geo_scene old_edit = edit_scene;
  edit_scene = sc;
  enable_done_button();
  disable_new_button();
  disable_close_button();
  
  activate(edit_scene);
  edit();  

  bool ret = (button_info[GEO_CM_QUIT] == 2) ? false : true;
    
  if( old_edit && ret ) 
    {
      activate(old_edit);  
      enable_done_button();
    }
  else 
    {
      disable_done_button();
      enable_new_button();
      enable_close_button();
    }
  edit_scene = old_edit;

  return ret;
}


// **********************************************************************
// *************************** Window Operations ************************

void GeoWin::init(double x1, double x2, double y1)
{
  if (x2<x1){
   double tm;
   tm=x1; x1=x2; x2=tm;
  }
  if (x1==x2) x2=x2+1; 

  double new_scale = Wp->width()/(x2-x1);
  int g = (int)grid;
  if( show_grid && g*new_scale < 4 ) { g = (int)(6/new_scale); if (g>0) g=-g; }
  if( show_grid )  Wp->init(x1, x2, y1, g, false);
  else             Wp->init(x1, x2, y1, 0, false);
}

void GeoWin::init(double x1, double x2, double y1, double y2, int rand)
{
  double x_min, x_max, y_min, y_max;
  
  if( x2 < x1 ) { x_min = x2; x_max = x1; }
  else          { x_min = x1; x_max = x2; }
  if( y2 < y1 ) { y_min = y2; y_max = y1; }
  else          { y_min = y1; y_max = y2; }
  
  double dx = x_max-x_min;
  double dy = y_max-y_min;
  
  double w = (Wp->xmax()-Wp->xmin())*Wp->scale();
  double h = (Wp->ymax()-Wp->ymin())*Wp->scale();
  
  double weff = w - 2*rand;
  double heff = h - 16 - 2*rand;
  
  double new_scale = dx ? weff/dx : Wp->scale();
  
  if( (!dx && dy) ||  (dy &&  heff/dy < weff/dx ) )
    {
      new_scale = heff/dy;
      x_min -=  (w/new_scale - dx)/2;
      x_max +=  (w/new_scale - dx)/2;
      y_min -=  (18+rand)/new_scale; 
    }
  else
    {
      x_min -= rand/new_scale;
      x_max += rand/new_scale;
      y_min -= (h/new_scale - 16/new_scale - dy)/2;
    }
  
  int g = (int)grid;
  if( show_grid && g*new_scale < 4 ) { g = (int)(6/new_scale); if (g>0) g=-g; }
  
  if( show_grid )  Wp->init(x_min, x_max, y_min, g);
  else             Wp->init(x_min, x_max, y_min, 0);
}


void GeoWin::redraw_geowin(window* wp, double x0, double y0, 
                              double x1, double y1)
{
  ((GeoWin*)(wp->get_inf()))->redraw(x0, y0, x1, y1);
}

void GeoWin::coord_handler(window* win, double x, double y)
{
  BASE_WINDOW* swin = win->get_status_window();
  //cout << "coord handler !\n"; cout.flush();

  if( !swin || !swin->is_open()) return;
  
  double y0 = swin->ymin();
  double x1 = swin->xmax();
  double y1 = swin->ymax();
  
  string coord_str("%8.2f %8.2f", x, y);
  double dy = swin->text_width(coord_str) + 2;

  double cy = y1 - swin->pix_to_real(1);
  double cx = x1 - dy;
  swin->draw_filled_rectangle(cx,y0,x1,y1,ivory);
  swin->draw_text(cx,cy,coord_str);
}

void GeoWin::redraw() 
{
  if( is_open && !is_locked && !scribble_mode)
    redraw(Wp->xmin(), Wp->ymin(), Wp->xmax(), Wp->ymax());
}

GEOWIN_END_NAMESPACE


#if defined GEOWIN_USE_NAMESPACE
static int scene_compare_z(const GEOWIN_NAMESPACE_NAME::geo_scene& sc1, const GEOWIN_NAMESPACE_NAME::geo_scene& sc2)
{
  int z1 = sc1->z, z2 = sc2->z;  
  return - compare(z1,z2);
}
#else
static int scene_compare_z(const geo_scene& sc1, const geo_scene& sc2)
{
  int z1 = sc1->z, z2 = sc2->z;  
  return -compare(z1,z2);
}
#endif


GEOWIN_BEGIN_NAMESPACE

void GeoWin::redraw(double x0, double y0, double x1, double y1) 
{ 

  if( !is_open || is_locked || scribble_mode) return;
  
  drawing_mode old_mode = Wp->set_mode(src_mode);
  Wp->start_buffering();
  
  Wp->clear(x0, y0, x1, y1);
  
  geo_scene sc;

    // redraw the scenes in z - order
    // cout << "redraw: z order\n";
  
    list<geo_scene>  scenes_0;   // z==0
    list<geo_scene>  scenes_pos; // z>0
    
    forall(sc,scenes){
      int z = get_z_order(sc);
      if (z==0) scenes_0.append(sc);
      else scenes_pos.append(sc);
    }
    // redraw the scenes with z coordinate > 0 ...
    scenes_pos.sort(scene_compare_z);
    forall( sc, scenes_pos ) 
     if( sc->visible || is_active(sc)) sc->redraw(Wp, x0, y0, x1, y1);        
    
    // redraw the scenes with z coordinate 0 ...
    forall( sc, scenes_0 ) 
     if( sc->visible || is_active(sc)) sc->redraw(Wp, x0, y0, x1, y1);    
     
  // redraw current scene on top ...
  //if(cur) { cur->redraw(Wp, x0, y0, x1, y1); }

  draw_pin();
  draw_info();
  
  Wp->flush_buffer();
  Wp->stop_buffering();
  Wp->set_mode(old_mode);
}



void GeoWin::update_status_line()
{
  //cout << "update status line !\n"; cout.flush();

  if( is_open )
    {
      string str;
      if( cur )  str = cur->information();
      Wp->set_status_string(str);
    }
}

void GeoWin::redraw_and_update_status_line()
{
  redraw();
  update_status_line();
}


// ***********************************************************************
// *************************** Scene Operations **************************

void GeoWin::insert_scene(geo_scene sc)
{
  if (!sc) error_handler(1, "Can not insert (geo_scene)0 into a GeoWin!");
  
  string tmp = sc->get_name();
  if( tmp == string("") || scene_name_map[tmp] )
    set_name(sc, sc->get_base_name());
  
  if( sc->pos )
    error_handler(1,"insert_scene : Scene is already inserted into a GeoWindow.");
  sc->gw  = this;
  
  sc->pos = scenes.append(sc);
  
  if( scenes.length() == 1 ) 
    {
      activate(sc);
      no_scene_off();
      update_status_line();
    }
  make_scene_menu();
  init_quick_access_window();
  
  redraw();
}

void GeoWin::remove(geo_scene sc)
{
  if (!sc) error_handler(1, "remove : Can not remove (geo_scene)0 !");
  
  scenes.del_item(sc->pos);
  bool rd = sc->is_visible();
  sc->pos = 0;
  sc->gw  = 0;
  
  string tmp1 = sc->get_name();
  string tmp2 = sc->get_base_name();
  
  scene_name_map[tmp2]--;
  if( tmp1 == tmp2 ) tmp1 = string("%s-%d", tmp2, 1);
  scene_name_map[tmp1]=0;
 
  if ( cur == sc )
    {
      if( !scenes.empty() ) { activate(scenes.head()); }
      else                  { cur = 0; no_scene_on(); }
      update_status_line();
    }
  
  make_scene_menu();  
  init_quick_access_window();
  if( rd ) redraw();
}

void GeoWin::destroy(geo_scene sc)
{
  if (!sc) error_handler(1, "destroy : Can not destroy (geo_scene)0 !");
  //sc->clear();
  //cout << "delete scene!\n"; cout.flush();
  if (sc == edit_scene) edit_scene = NULL;
  delete sc;
  //cout << "after deleting scene!\n"; cout.flush();  
  init_quick_access_window();
}

void GeoWin::destroy_all()
{
  geo_scene sc = get_active_scene();
  while (sc != NULL) {
   destroy(sc);
   sc = get_active_scene();
  }
  init_quick_access_window();
}


geo_scene GeoWin::get_scene_with_name(string nm) const
{
  geo_scene sc;
  forall(sc, scenes) if( sc->get_name() == nm ) return sc;
  return 0;
}

string GeoWin::set_name(geo_scene sc, string nm)
{
  if( nm != string("") )
    {
      string tmp = sc->get_name();
      if( tmp != string("") ) scene_name_map[tmp]--;
      
      sc->set_base_name(nm);
      int i, c = 0;
      
      if( scene_name_map.defined(nm) )  c = scene_name_map[nm];
      
      for( i = 0; i <= c; i++)
        {
          tmp = string("%s-%d", nm, i+1);
          
          if( !scene_name_map[tmp] )  
            {
              scene_name_map[tmp] = 1;
              if( i==0 ) tmp = nm;
              scene_name_map[nm]  = c+1;
              sc->set_name(tmp);
              break;
            }
        }  
      make_scene_menu();
    }
  return sc->get_name();
}

string GeoWin::get_short_name(geo_scene sc) { return sc->short_name; }

string GeoWin::set_short_name(geo_scene sc, string sn)
{
 string prev = sc->short_name;
 sc->short_name = sn;
 return prev;
}

void GeoWin::redraw_inp_panel()
{
  if (cur) {
   string sc_name = cur->get_type_name();
  
   if ((sc_name == string("RationalPoints")) || (sc_name == string("Points")) ) 
    Wp->enable_item(inpit);
   else Wp->disable_item(inpit);
  }
  else Wp->disable_item(inpit); 
  
  Wp->redraw_panel(inpit);
}


void GeoWin::activate(geo_scene sc)
{
  if( !is_open || !sc )  return;
  
  if( cur && cur != sc ) cur->set_active(false);
  
  string sc_name = sc->get_type_name();
  //cout << sc_name << "\n";
  
  if ((sc_name == string("RationalPoints")) || (sc_name == string("Points")) ) 
    Wp->enable_item(inpit);
  else 
    Wp->disable_item(inpit);
    
  cur = sc;
  cur->set_active(true);
  
  if( cur->edit_menu_type != -1)
    {
      Wp->set_window(menu_const + edit_menu1, menus[cur->edit_menu_type]);
      Wp->enable_button(menu_const + edit_menu1); 
      edit_scene =sc;
    }
  else   Wp->disable_button(menu_const + edit_menu1);
  
  // handle choice item ...
  if (activate_buttons){
    //cout << "activate buttons !\n"; 
    int i=0;
    geo_scene it;
    forall(it,ACT_SCN) { if (it==sc) ACT_STATE=i; i++; }
    Wp->redraw_panel(acit);
  }
  INP_STATE = sc->fcn_state;
  Wp->redraw_panel(inpit);
  
  init_quick_access_window();
  
  make_algo_menu(sc); 
  redraw();
  update_status_line();
  
  return;
}

int GeoWin::get_z_order(geo_scene sc) const { return sc->z; }

int GeoWin::set_z_order(geo_scene sc,int n) {
 int old = sc->z;
 sc->z = n;
 return old;
}

geo_scene GeoWin::get_active_scene() const { return cur; }

color GeoWin::get_color(geo_scene sc) const { return sc->get_color(); }

color GeoWin::set_color(geo_scene sc, color c) { 
  color cold=sc->get_color(); sc->set_color(c); return cold; 
}

void GeoWin::set_color(geo_scenegroup gs, color c) {
  geo_scene iter;
  forall(iter,gs->group_members) set_color(iter,c);
}
  
color GeoWin::get_color2(geo_scene sc) const { return sc->get_color2(); }

color GeoWin::set_color2(geo_scene sc, color c) { 
  color cold=sc->get_color2(); sc->set_color2(c); return cold;
}

void GeoWin::set_selection_color(geo_scenegroup gs, color c) {
  geo_scene iter;
  forall(iter,gs->group_members) set_color2(iter,c);  
}


color GeoWin::get_selection_fill_color(geo_scene sc) const
{
  return sc->sel_filling_color;
}  
	  
color GeoWin::set_selection_fill_color(geo_scene sc, color c)
{
  color prev = sc->sel_filling_color;
  sc->sel_filling_color = c;
  return prev;
} 
  
line_style GeoWin::get_selection_line_style(geo_scene sc) const
{
  return sc->sel_l_style;
}
  
line_style GeoWin::set_selection_line_style(geo_scene sc, line_style l) const
{
  line_style prev = sc->sel_l_style;
  sc->sel_l_style = l;
  return prev;
}  
  
int GeoWin::get_selection_line_width(geo_scene sc) const
{
  return sc->sel_line_width;
}
  
int GeoWin::set_selection_line_width(geo_scene sc, int w) const
{   
  int prev = sc->sel_line_width;
  sc->sel_line_width = w;
  return prev;
}


color GeoWin::set_fill_color(geo_scene sc,color c) { 
  color cold=sc->filling_color; sc->filling_color=c; 
  return cold; 
}

void GeoWin::set_fill_color(geo_scenegroup gs, color c) {
  geo_scene iter;
  forall(iter,gs->group_members) set_fill_color(iter,c);  
}

color GeoWin::get_fill_color(geo_scene sc) const { return sc->filling_color; }

color GeoWin::set_text_color(geo_scene sc,color c) { 
  color cold=sc->text_color; sc->text_color=c; 
  return cold; 
}

void GeoWin::set_text_color(geo_scenegroup gs, color c) {
  geo_scene iter;
  forall(iter,gs->group_members) set_text_color(iter,c);  
}

color GeoWin::get_text_color(geo_scene sc) const { return sc->text_color; }


void GeoWin::set_line_width(geo_scenegroup gs, int w) {
  geo_scene iter;
  forall(iter,gs->group_members) set_line_width(iter,w);  
}

void GeoWin::set_active_line_width(geo_scenegroup gs, int w) {
  geo_scene iter;
  forall(iter,gs->group_members) set_active_line_width(iter,w);  
}

void GeoWin::set_line_style(geo_scenegroup gs, line_style l) {
  geo_scene iter;
  forall(iter,gs->group_members) set_line_style(iter,l);  
}

void GeoWin::set_point_style(geo_scenegroup gs, point_style p) {
  geo_scene iter;
  forall(iter,gs->group_members) set_point_style(iter,p);  
}

bool GeoWin::get_visible(geo_scene sc) const { return sc->get_visible(); }

bool GeoWin::set_visible(geo_scene sc, bool v) { 
 bool bold=sc->get_visible(); sc->set_visible(v); return bold; }
 
void GeoWin::set_visible(geo_scenegroup gs, bool v) {
  geo_scene iter;
  forall(iter,gs->group_members) set_visible(iter,v);  
}
 
void GeoWin::set_all_visible(bool v){
 geo_scene sc;
 forall(sc,scenes) sc->set_visible(v);                                                                                                                                              
}

string GeoWin::set_description(geo_scene sc, string desc) 
{ string prev = sc->description;
  sc->description = desc; 
  return prev;
}

bool GeoWin::get_highlight(geo_scene sc) const
{
  return sc->high_light;
}
  
bool GeoWin::set_highlight(geo_scene sc, bool hl)
{
  bool prev = sc->high_light;
  sc->high_light = hl;
  return prev;
}

void GeoWin::set_handle_defining_points(geo_scene sc, geowin_defining_points gdp)
{ sc->handle_defining_points = gdp; } 
  
geowin_defining_points GeoWin::get_handle_defining_points(geo_scene sc)
{ return sc->handle_defining_points; }    
  
string GeoWin::get_description(geo_scene sc) const { return sc->description; }

geowin_redraw* GeoWin::get_redraw_pt(geo_scene sc) { return sc->redraw_pt; }

GEOWIN_END_NAMESPACE

// ***********************************************************************
// *******************  Input and Output Operations  *********************

#if defined GEOWIN_USE_NAMESPACE
static GEOWIN_NAMESPACE_NAME::GeoWin* GeoWp = 0;
static GEOWIN_NAMESPACE_NAME::geo_scene sc_help = NULL;
#else
static GeoWin* GeoWp =0;
static geo_scene sc_help = NULL;
#endif

static string geowin_choose_string = "";

static void geo_read_scene(string fn)
{ GeoWp->read_scene(fn); }

static void geo_import_scene(string fn)
{ GeoWp->import_scene(fn); }

void geo_write_active_scene(string fn)
{ GeoWp->write_scene(fn); }

void geo_write_scene(string fn)
{ 
 //cout << "geo_write_scene !\n"; cout.flush();
 GeoWp->write_scene(sc_help, fn, true); 
}

void geo_write_ps(string fn)
{ GeoWp->write_ps(fn); }

void geo_write_ps_vis(string fn)
{ GeoWp->write_ps_vis(fn); }

void geo_export_scene(string fn)
{ GeoWp->export_scene(fn); }

static bool s_full_color = true;

void geo_write_screen(string fn)
{ GeoWp->redraw();
  GeoWp->get_window().screenshot(fn,s_full_color); 
}

void geo_dummy(string fn)
{}

static bool read_without_header = false;
static bool add_obj_to_scene = false;

GEOWIN_BEGIN_NAMESPACE

void GeoWin::read_scene(string filename)
{
  ifstream is(filename);
  if( !is )  {  error_handler(0, "No such file or directory");  return; }
  
  string type_name;
  string sc_name;
  string identifier;
  double version;
  
  char* t_name = new char[200];
  char* s_name = new char[200];

  geo_scene sc;
  
  is >> identifier >> version;
  
  if( !(identifier == string("GEOWIN")) )
    {
      is.close();
      
      if( edit_scene && read_without_header ) 
        {
          ifstream is2(filename);
          edit_scene->read(is2,add_obj_to_scene);
        }
      else error_handler(0, "Incorrect Fileheader.");
      return;
    }
  
  // read rest of line before type name ...
  is.getline(t_name,200);  
  // read type name
  is.getline(t_name,200);    
  type_name = string(t_name);
  // read scene name
  is.getline(s_name,200);  
  sc_name = string(s_name);

  delete t_name; delete s_name;
  
  if( edit_scene ){  
    if( type_name == edit_scene->get_type_name() )  sc = edit_scene;
    else  
      {
        type_name = string("read : File contains data for type \n \t") + 
          type_name;
        type_name += string(",\n \t but the edited scene contains data ");
        type_name += string("of type \n \t");
        type_name += edit_scene->get_type_name();
        error_handler(0, type_name.cstring());
        return;
      }
  }
  else
    {   
      sc = GeoScene::new_scene(type_name);
      if (sc == NULL) return;
      
      sc->set_base_name(sc_name);
      insert_scene(sc);
      activate(sc);
    }
  sc->read(is,add_obj_to_scene);
  redraw();
  update_status_line();
}

void GeoWin::import_scene(string filename)
{
  list<string> name, desc;
  cur->get_import_object_names_and_descriptions(name, desc);   
  
  //which importer was choosen ...    
  string iter;
  int i = 0, import_handler_nr = -1;
  forall(iter,name){
      if (geowin_choose_string==iter) { 
	import_handler_nr = i; break; 
      }
      i++;
  }  
  if (import_handler_nr == -1) return;

  // get import object ...  
  list_item it = (cur->import_objects).get_item(import_handler_nr);
  geowin_import* io = (cur->import_objects)[it];
 
  io->operator()(cur, filename);
}  

// --------------------------------------
// export and import data ....
// --------------------------------------

void GeoWin::file_export()
{
  if (cur){
   list<string> name, desc;
   cur->get_export_object_names_and_descriptions(name, desc); 
   name.push_front("default");
    
   string str1;
   str1 = string(cur->get_name()) + string(".geo");
   string str2(".");
  
   GeoWp = this;
   file_panel P(get_window(),str1, str2);
   P.set_save_handler(::geo_export_scene);
   P.set_cancel_handler(::geo_dummy); 
   P.set_pattern("Geo Files (*.geo)","*.geo");
   
   panel& pan = P.get_panel();
   P.init_panel();
 
   pan.text_item("\\bf\\blue Choose Exporter"); 
 
   geowin_choose_string = name.head();
   pan.string_item("Exporter", geowin_choose_string, name);
       
   open_panel(P);  
   GeoWp = 0; 
  } 
}

void GeoWin::file_import()
{
  if (cur){
    list<string> name, desc;
    cur->get_import_object_names_and_descriptions(name, desc); 
    
    if (name.empty()) return; // no importers ...
    
    string str1("no name");
    string str2(".");
  
    GeoWp = this;
    file_panel P(get_window(),str1, str2);
    P.set_load_handler(::geo_import_scene);
    P.set_cancel_handler(::geo_dummy);
    P.set_pattern("Geo Files (*.geo)","*.geo");
  
    panel& pan = P.get_panel();
    P.init_panel();
 
    pan.text_item("\\bf\\blue Choose Importer"); 
 
    geowin_choose_string = name.head();
    pan.string_item("Importer", geowin_choose_string, name);

    open_panel(P);
    GeoWp = 0;
  }
}

// save active scene ...
void GeoWin::file_save_active()
{
  string str1;
  if(cur) str1 = string(cur->get_name()) + string(".geo"); else return;
  string str2(".");
  
  GeoWp = this;
  file_panel P(get_window(),str1, str2);
  P.set_save_handler(::geo_write_active_scene);
  P.set_cancel_handler(::geo_dummy);
  P.set_pattern("Geo Files (*.geo)","*.geo");
  open_panel(P);  
  GeoWp = 0;
}

// save some scene ...
void GeoWin::file_save()
{
  geo_scene sc = NULL, sc_iter;
  
  //get scene ...
  list<geo_scene> scn = get_scenes();
  list<string> names;
  forall(sc_iter, scn) names.append(get_name(sc_iter));
  
  if (names.empty()) return;
  string choose = names.head();
  
  panel Psc("Choose scene");
  Psc.buttons_per_line(2);
  
  Psc.text_item("\\bf\\blue Scenes:");
  Psc.text_item("");
 
  Psc.string_item("", choose, names );
  
  Psc.fbutton("Apply", APPLY_BUTTON);
  Psc.button("Cancel", CANCEL_BUTTON);
  
  int but = open_panel(Psc);
  
  if( but == CANCEL_BUTTON ) return;
  // find scene ...
  forall(sc_iter, scn) { if (choose== get_name(sc_iter)) sc = sc_iter; }

  string str1;
  if(sc) str1 = string(sc->get_name()) + string(".geo"); else return;
  string str2(".");
  
  sc_help = sc;
  GeoWp = this;
  file_panel P(get_window(),str1, str2);
  
  P.set_save_handler(::geo_write_scene);
  P.set_cancel_handler(::geo_dummy);
  P.set_pattern("Geo Files (*.geo)","*.geo");
  open_panel(P);  
  GeoWp = 0;
}

void GeoWin::save_ps()
{
  string str1;
  if(cur) str1 = string(cur->get_name()) + string(".ps"); else return;
  string str2(".");
  
  GeoWp = this;
  file_panel P(get_window(),str1, str2);
  P.set_save_handler(::geo_write_ps);
  P.set_cancel_handler(::geo_dummy);
  P.set_pattern("PS Files (*.ps)","*.ps");
  open_panel(P);  
  GeoWp = 0;
}

void GeoWin::save_ps_vis()
{
  string str1("ps_output.ps");
  if(! cur) return;
  string str2(".");
  
  GeoWp = this;
  file_panel P(get_window(),str1, str2);
  P.set_save_handler(::geo_write_ps_vis);
  P.set_cancel_handler(::geo_dummy);
  P.set_pattern("PS Files (*.ps)","*.ps");
  open_panel(P);  
  GeoWp = 0;
}

bool GeoWin::test_file_exists(string filename)
{
  if (file_existence_test){
   if (is_file(filename)){ // file is already there ...
    window& wind = get_window();
    string cnf("File exists. Overwrite it ?");
    bool b = wind.confirm(cnf);
    return !b;
   }
  }
  return false;
}

void GeoWin::write_scene(string filename)
{
  if (!cur) return;
  
  //test, if file already exists
  if (test_file_exists(filename)) return;
  
  ofstream os(filename);
  if( !os )  { error_handler(0, "Can not save");  return;   }
  string str("GEOWIN ");
  os << str << version() << endl;
  os << cur->get_type_name() << endl;
  os << cur->get_base_name() << endl;
  cur->write(os,"");
  os.close();
}

void GeoWin::write_scene(geo_scene sc, string filename, bool test)
{
  if (test) {
    if (test_file_exists(filename)) return;
  }  
    
  ofstream os(filename);
  if( !os )  { error_handler(0, "Can not save");  return;   }
  string str("GEOWIN ");
  os << str << version() << endl;
  //cout << "write_scene get type name :\n" << flush;
  string tn = sc->get_type_name();
  
  if ((tn.length() > 4) && (tn.tail(4) == string("Base"))){ // cut end off ...
   tn = tn.head(tn.length()-4);
  }
  
  os << tn << endl;
  os << sc->get_base_name() << endl;
  sc->write(os,"");
  os.close();
}

void GeoWin::write_ps(string filename)
{
  if (!cur) return;
  if (test_file_exists(filename)) return;  
  
  window& win=get_window();

  ps_file F(filename);
  F.init(win.xmin(),win.xmax(),win.ymin());
  F.set_output_mode(colored_mode);
  F.set_point_style(cur->p_style);  
  cur->write_postscript(F);
  F.close();
}

/*
void GeoWin::write_ps_vis(string filename)
{
  if (!cur) return;
  if (test_file_exists(filename)) return;
  
  geo_scene sc;
  window& win=get_window();

  ps_file F(filename);
  F.init(win.xmin(),win.xmax(),win.ymin());
  F.set_output_mode(colored_mode); 

  // iterate through scenes...
  forall(sc,scenes) {
    F.set_point_style(sc->p_style);
    if (sc->visible && sc != cur) sc->write_postscript(F); 
  }  
  if(cur) {
    F.set_point_style(cur->p_style); 
    cur->write_postscript(F);
  }
  
  F.close();
}
*/


void GeoWin::write_ps_vis(string filename)
{
  if (!cur) return;
  if (test_file_exists(filename)) return;
  
  geo_scene sc;
  window& win=get_window();

  ps_file F(filename);
  F.init(win.xmin(),win.xmax(),win.ymin());
  F.set_output_mode(colored_mode); 
  
  // output scenes ...
  list<geo_scene>  scenes_0;   // z==0
  list<geo_scene>  scenes_pos; // z>0
    
  forall(sc,scenes){
      int z = get_z_order(sc);
      if (z==0) scenes_0.append(sc);
      else scenes_pos.append(sc);
  }
    
  // redraw the scenes with z coordinate > 0 ...
  scenes_pos.sort(scene_compare_z);
  forall( sc, scenes_pos ) {
     F.set_point_style(sc->p_style);
     if( sc->visible || is_active(sc)) sc->write_postscript(F);   
  }     
    
  // redraw the scenes with z coordinate 0 ...
  forall( sc, scenes_0 ) {
     F.set_point_style(sc->p_style);
     if( sc->visible || is_active(sc)) sc->write_postscript(F);
  }
  
  F.close();
}




void GeoWin::replace_characters(string& strh)
{
  int max = strh.length();
  int i;
  
  set<char> help_set;
  
  // init the set ...
  for (i=0;i<replace_by_blanks.length();i++){
    help_set.insert(replace_by_blanks[i]);
  }
  
  for(i=0;i<max;i++){ 
    char& act = strh[i];
    if (help_set.member(act)) act = ' '; 
  }
}

void GeoWin::export_scene(string filename)
// write to file without header
// or export using a geowin_export object
{
  if (!cur) return;
  if (test_file_exists(filename)) return;
  
  ofstream os_file(filename);
  
  if (geowin_choose_string == "default"){
    string_ostream os;
  
    if( !os_file )  { error_handler(0, "Can not save");  return;   }
    cur->write(os,"\n");
    os << "\n";
    os << ends;
  
    string strh = os.str();
    replace_characters(strh); 
  
    os_file << strh;
  }
  else { // export object ...
    list<string> name, desc;
    cur->get_export_object_names_and_descriptions(name, desc);   
      
    string iter;
    int i = 0;
    int export_handler_nr = -1;
    forall(iter,name){
      if (geowin_choose_string==iter) { 
	export_handler_nr = i; break; 
      }
      i++;
    }  
    if (export_handler_nr == -1) return;
  
    // get export object ...
  
    list_item it = (cur->export_objects).get_item(export_handler_nr);
    geowin_export* eo = (cur->export_objects)[it];
    eo->operator()(cur, filename);  
  }
    
  os_file.close();
}

void GeoWin::file_open()
{
  string str1("no name");
  string str2(".");
  
  GeoWp = this;
  file_panel P(get_window(),str1, str2);
  P.set_load_handler(::geo_read_scene);
  P.set_cancel_handler(::geo_dummy);
  P.set_pattern("Geo Files (*.geo)","*.geo");
  
  panel& pan = P.get_panel();
  P.init_panel();
  pan.bool_item("Read without Header : ", read_without_header);
  pan.bool_item("Add objects to scene: ", add_obj_to_scene);

  open_panel(P); 
  GeoWp = 0;
}

void GeoWin::save_screen()
{
  string fname,dname(".");
 
#if defined(__win32__)
                 fname = "gw_screenshot.wmf";
#else
                 fname = "gw_screenshot.ps";
#endif

  GeoWp = this;
  file_panel P(get_window(),fname, dname);
  P.set_save_handler(::geo_write_screen);
  P.set_cancel_handler(::geo_dummy);

#if defined(__win32__)
  P.set_save_string("Write Metafile Screenshot");
  P.set_pattern("WMF Files (*.wmf)","*.wmf");
#else
  P.set_save_string("Write Postscript Screenshot");
  P.set_pattern("PS Files (*.ps)","*.ps");
#endif
         
  panel& pan = P.get_panel();
  P.init_panel();
  pan.bool_item("full color",s_full_color);

  open_panel(P);  
  GeoWp = 0;
}


// ******************************************************************
// *****************  View Operations  ******************************

void GeoWin::zoom_up()   {  
 int gmode;
 double gdist;
 grid_style gstyle;
 gw_save_grid_paras(Wp,gdist,gmode,gstyle);
 Wp->zoom(2.0,gw_anim_steps); 
 gw_restore_grid_paras(Wp,gdist,gmode,gstyle); 
}

void GeoWin::zoom_down() {  
 int gmode;
 double gdist;
 grid_style gstyle;
 gw_save_grid_paras(Wp,gdist,gmode,gstyle);
 Wp->zoom(0.5,gw_anim_steps); 
 gw_restore_grid_paras(Wp,gdist,gmode,gstyle); 
}


void GeoWin::zoom_area()
{
  if (!cur) return;

  point p1,p2;
  if (Wp->read_zoom_rect(p1,p2)){
     Wp->zoom_area(p1.xcoord(),p1.ycoord(),p2.xcoord(),p2.ycoord(),gw_anim_steps);
     Wp->set_show_coord_handler(GeoWin::coord_handler);
     update_status_line();
     redraw();
  }
  else 
     Wp->redraw();
}

void GeoWin::reset_window()
{
     Wp->zoom_area(START_XMIN,START_YMIN,START_XMAX,START_YMAX,gw_anim_steps);
}

void GeoWin::fill_window()
{ 
  double x_min =  MAXINT;
  double y_min =  MAXINT;
  double x_max = -MAXINT;
  double y_max = -MAXINT;


  geo_scene sc;
  forall( sc, scenes ) 
  { geo_editor ed = sc->type_editor();
    if (ed && sc->is_visible()) 
    { double x_min_akt, x_max_akt, y_min_akt, y_max_akt;
      if (sc->fill_rect(x_min_akt, x_max_akt, y_min_akt, y_max_akt))
      { if (x_min_akt < x_min) x_min = x_min_akt;
        if (x_max_akt > x_max) x_max = x_max_akt;
        if (y_min_akt < y_min) y_min = y_min_akt;
        if (y_max_akt > y_max) y_max = y_max_akt;
       }
     }
   }

  if (x_min < MAXINT) 
    Wp->zoom_area(x_min, y_min, x_max, y_max, gw_anim_steps);
}


string GeoWin::set_bg_pixmap(string pix_name)
{
  string old = bg_pixmap;
  color bg_col = Wp->get_bg_color();
  color fg_col = Wp->get_fg_color();

  list<string> pm_list;

  int i;
  for(i=0; i < num_texture; i++) pm_list.append( name_texture[i] );
  pm_list.append("none");
  
  i=0;
  string str;
  forall(str, pm_list) 
    if( str == pix_name ) break;
    else i++;
  
  char* pr = 0;
  if( i >= 0 && i < num_texture)
    pr = Wp->create_pixrect(xpm_texture[i]);
  else return old;
  
  pr = Wp->set_bg_pixmap(pr);
  if( pr ) Wp->del_pixrect(pr);
  bg_pixmap = pix_name;
  
  Wp->set_fg_color(fg_col);
  Wp->set_bg_color(bg_col);
  
  return old;
}

color GeoWin::get_bg_color() const
{
  return get_window().get_bg_color();
}

color GeoWin::set_bg_color(const color& c)
{
  return get_window().set_bg_color(c);
}

bool  GeoWin::get_show_grid() const
{
  return show_grid;
}
  
bool  GeoWin::set_show_grid(bool sh)
{
  double new_scale = Wp->scale();
  bool prev = show_grid;
  show_grid = sh;
  
  double g = grid;
  if( show_grid && g*new_scale < 4 ) g = 6/new_scale;
  if( show_grid ) { if (g>0) g=-g;  Wp->set_grid_dist(g); }
  else Wp->set_grid_mode(0);
  
  redraw();
  return prev;
}

double GeoWin::get_grid_dist() const
{
  return grid;
}
  
double GeoWin::set_grid_dist(double g)
{
  double new_scale = Wp->scale();
  double prev = grid;
  grid = g;
  if( show_grid && g*new_scale < 4 ) g = 6/new_scale;
  if( show_grid ) { if (g>0) g=-g;  Wp->set_grid_dist(g); }
  else Wp->set_grid_mode(0);
  
  redraw();
  return prev;
}

grid_style  GeoWin::get_grid_style() const
{
  return Wp->get_grid_style();
}

  
grid_style  GeoWin::set_grid_style(grid_style gr_style)
{
  double new_scale = Wp->scale();
  double gh = grid;
  if( show_grid && gh*new_scale < 4 ) gh = 6/new_scale;
  if( show_grid ) { if (gh>0) gh=-gh;  Wp->set_grid_dist(gh); }
  else            Wp->set_grid_mode(0);
  grid_style prev = Wp->set_grid_style(gr_style);
  redraw();
  return prev;
}

bool  GeoWin::get_show_position() const
{
  return show_coords;
}
  
bool  GeoWin::set_show_position(bool sp)
{
  bool prev = show_coords;
  if( !show_coords && sp ) 
    Wp->set_show_coord_handler(GeoWin::coord_handler);
  if( show_coords && !sp )  Wp->set_show_coord_handler(0);  
  show_coords = sp;
  redraw();
  return prev;
}

bool  GeoWin::get_show_quick_access_window() const
{
  return show_quick_access_window;
}
  
bool  GeoWin::set_show_quick_access_window(bool b)
{
  bool prev = show_quick_access_window;
  show_quick_access_window = b; 
  return prev;
}

bool GeoWin::get_d3_elimination() const
{ return gw_elim; }
  
bool GeoWin::set_d3_elimination(bool b)
{ bool h = gw_elim; gw_elim = b; return h; }
  
bool GeoWin::get_d3_solid() const
{ return gw_solid; }
  
bool GeoWin::set_d3_solid(bool b)  
{ bool h = gw_solid; gw_solid = b; return h; }

int GeoWin::get_d3_animation_speed() const
{ return gw_anim_speed; }
  
int GeoWin::set_d3_animation_speed(int s)
{ int h = gw_anim_speed; gw_anim_speed = s; return h; }

// --- new : use the scene color in d3 output ?
bool GeoWin::get_d3_use_scene_color(geo_scene sc) const
{  return sc->d3_use_scene_color; }

bool GeoWin::set_d3_use_scene_color(geo_scene sc, bool b)
{  bool prev = sc->d3_use_scene_color;
   sc->d3_use_scene_color = b;
   return prev;
}

bool GeoWin::get_d3_show_edges() const
{ return gw_d3edges; }
  
bool GeoWin::set_d3_show_edges(bool b)
{ bool h = gw_d3edges; gw_d3edges = b; return h; }

void GeoWin::z_scenes()
{
  if (scenes.empty()) return;
   
  int *znr= new int[scenes.length()];
  int i=scenes.length() -1;
  geo_scene sc;
  
  panel P("Z Order Scenes");
  P.buttons_per_line(2);
  
  P.text_item("\\bf\\blue Scenes/z-value");
  P.text_item("");
  P.text_item("");
  
  list_item it;
  for(it=scenes.last();it!=NULL;it=scenes.pred(it))  
    { 
      sc= scenes[it];
      znr[i] = sc->z;
      if (sc != cur) P.int_item(sc->get_name(), znr[i], 0, scenes.length()-1);
      else {
	string strh = sc->get_name();
        P.text_item("\\red "+strh);
	P.int_item("(active scene)", znr[i], 0, scenes.length()-1);
      }
      i--;
    }
  
  P.fbutton("Apply", APPLY_BUTTON);
  P.button("Cancel", CANCEL_BUTTON);

  window& W = get_window();
  W.disable_panel();
  int but = P.open();
  W.enable_panel();
  redraw_inp_panel();
  
  if (but) {
      i = 0;
      forall(sc , scenes) sc->z = znr[i++]; 
      redraw();
  }
  
  delete znr;
}

void GeoWin::visible_scenes()
{
  if( !cur ) return;
  
  bool *vis= new bool[scenes.length()];

  int i=0;
  geo_scene sc;
  
  panel P("Visible Scenes");
  P.buttons_per_line(2);
  
  P.text_item("\\bf\\blue Visible Scenes");
  P.text_item("");
  P.text_item("");
  
  forall(sc , scenes)  
    { 
      vis[i] = sc->get_visible();  
      
      if (sc != cur) P.bool_item( sc->get_name(), vis[i] );
      else {
        string strh = sc->get_name();
        P.text_item("\\red "+strh);
	P.bool_item( "(active scene)", vis[i] );
      }
      
      i++;
    }
  P.text_item("");
  P.fbutton("Apply", APPLY_BUTTON);
  P.button("Cancel", CANCEL_BUTTON);

  window& W = get_window();
  W.disable_panel();
  int but = P.open();
  W.enable_panel();
  redraw_inp_panel();
  
  if (but){ change_visibility(vis); init_quick_access_window(); }

  delete vis;
}

void GeoWin::change_visibility(bool* vis)
{
      geo_scene sc;
      int i = 0;
      
      forall(sc , scenes)  { sc->set_visible(vis[i++]); }  

      if (vis_buttons){     
       //reset buttons in vis_buttons ...
       geo_scene it;
    
       VIS_STATE = 0;
       i = 0;
    
       forall(it,VAR_1) {
         if (get_visible(it)) { VIS_STATE = VIS_STATE + (1 << i);  } 
         i++;
       }
       i=0;
       forall(it,VAR_2) {
         if (get_visible(it)) { VIS_STATE = VIS_STATE + (1 << i);  } 
         i++;
       }       
       
       if (! VAR_1.empty()) Wp->redraw_panel(pnit1);
       if (! VAR_2.empty()) Wp->redraw_panel(pnit2);
              
      }
      
      redraw();        
}

void GeoWin::scenegroup_dialog(geo_scenegroup gsg)
{
  panel P("Scene group");
  P.buttons_per_line(2);

  string st("\\bf\\blue Group name:");
  st = st + gsg->name;
  int i=0;
  bool *vis= new bool[gsg->group_members.length()];
  
  P.text_item(st);  
  P.text_item("");
  P.text_item("\\bf\\blue Visibility of group members:");
  P.text_item("");
  P.text_item("");  
  
  geo_scene iter;
  forall(iter,gsg->group_members)
  {
      vis[i] = iter->get_visible();  
      P.bool_item( iter->get_name(), vis[i] );
      i++;    
  }
  
  P.button("All visible",1);
  P.button("All invisible",2);
  P.button("Apply",3);
  P.button("Cancel", CANCEL_BUTTON); //0

  int but = open_panel(P);
  
  if (but != 0) {
   switch (but){
    case 1:
    { set_visible(gsg,true); break;}
    case 2:
    { set_visible(gsg,false); break; }
    case 3:
    {
      i = 0;
      forall(iter , gsg->group_members)  { iter->set_visible(vis[i++]);
        if (iter->is_visible()) iter->update(); 
      } 
      break;     
    }
   }
   
   if (vis_buttons){     
       geo_scene it;
       VIS_STATE = 0;
       i = 0;   
       forall(it,VAR_1) {
         if (get_visible(it)) { VIS_STATE = VIS_STATE + (1 << i);  } 
         i++;
       }
       i=0;
       forall(it,VAR_2) {
         if (get_visible(it)) { VIS_STATE = VIS_STATE + (1 << i);  } 
         i++;
       }      
       if (!VAR_1.empty()) Wp->redraw_panel(pnit1);     
       if (!VAR_2.empty()) Wp->redraw_panel(pnit2);               
      }
      redraw();   
      
     init_quick_access_window();
  }
 
  delete vis;
}


// ******************************************************************
// *****************  protected members  ****************************

void GeoWin::scene_defaults()
{
  int own_but;
 
  panel P("Scene Defaults");
  P.set_inf(&own_but);
   
  color c0 = DEFAULT_fill_color; 
  color c1 = DEFAULT_color; 
  color c2 = DEFAULT_color2;
  int aw = DEFAULT_active_line_width;
  int bw = DEFAULT_line_width;
  line_style lst = DEFAULT_line_style;
  point_style pst= DEFAULT_point_style;
  
  P.text_item("\\bf\\blue Colors");
  P.color_item(string("default interior"), c0);
  P.color_item(string("default boundary"), c1);
  P.color_item(string("default selected"), c2);
  P.text_item("");

  P.text_item("\\bf\\blue Points"); 
  P.pstyle_item("default point style", pst);
  
  P.text_item("\\bf\\blue Lines");
  P.lstyle_item("default line style", lst);
  P.lwidth_item("default line width", bw);
  P.lwidth_item("default active width", aw);
  
  P.fbutton("apply",   APPLY_BUTTON);
  P.button("cancel",   CANCEL_BUTTON);
  
  int& but = *((int*)(P.get_inf()));
  but = open_panel(P);

  if( but == CANCEL_BUTTON ) return;
  
  DEFAULT_fill_color = c0;
  DEFAULT_color = c1;
  DEFAULT_color2 = c2;
  DEFAULT_active_line_width = aw;
  DEFAULT_line_width = bw;
  DEFAULT_line_style = lst;
  DEFAULT_point_style = pst;
}

void GeoWin::advanced_global_options()
{
  bool scribble = scribble_mode;
  bool cscn = algo_create_scenes;
  bool tme= time_meas;
  bool mpc= mouse_pos_check;
  bool up = update_while_object_dragging;
  bool qa = show_quick_access_window;
  
  bool allow_mfg = allow_moving_from_grid;
  
  string rep_chars = replace_by_blanks;
  
  int  prec = rat_point::default_precision;
  bool  fi   = (rat_point::use_filter != 0);
  bool  check_existance = file_existence_test;

  // new -error handling
  //bool win_error = window::get_show_error_panel();

  panel P;
  P.set_inf(this);    
  P.text_item("\\bf\\blue Rational Kernel Options"); 
  P.int_item("precision", prec, 0, 53);
  P.bool_item("use filter",fi);
  P.text_item("");  
  
  P.text_item("\\bf\\blue File Options"); 
  P.bool_item("saving - file existance check", check_existance);  
  P.text_item("");
  P.text_item("Export - replace these characters by blanks");
  P.string_item("Characters:", rep_chars);
  P.text_item("");

  P.text_item("\\bf\\blue Miscellaneous");   
  P.bool_item("update while moving/rotating objects", up);
  P.bool_item("measure time", tme);
  P.bool_item("scribble mode", scribble);
  P.bool_item("buildin algorithms - create scenes",cscn);
  P.bool_item("current scene- check mouse position", mpc);
  P.bool_item("show quick access window", qa);
  P.bool_item("allow moving from grid in gridmode", allow_mfg);
  //P.bool_item("error handler - window panel", win_error);
  
  P.text_item("");
  
  P.fbutton("apply",   APPLY_BUTTON);
  P.button("defaults", DEFAULTS_BUTTON);
  P.button("cancel",   CANCEL_BUTTON);

  int but;
  
  while ((but=open_panel(P)) == DEFAULTS_BUTTON )
    { 
      up        = true;
      mpc       = true;
      prec      = 10;
      fi        = true;
      tme       = DEFAULT_time_meas;
      qa        = true;
      allow_mfg = false;
      rep_chars = ""; 
      check_existance  = true;
      //win_error = true;
    }
  
  if( but == CANCEL_BUTTON ) return;
  
  Wp->set_flush(false);
  
  // new
  //window::set_show_error_panel(win_error);
    
  scribble_mode=scribble;
  if (scribble_mode)  scribble_first=true;
  
  algo_create_scenes = cscn;
  
  update_while_object_dragging = up;

  rat_point::default_precision=prec;
  rat_point::use_filter=fi;
  time_meas = tme;
  mouse_pos_check = mpc;
  show_quick_access_window = qa;
  allow_moving_from_grid = allow_mfg;
  
  replace_by_blanks = rep_chars;
  file_existence_test = check_existance;
  
  init_quick_access_window();

  Wp->set_flush(true);
  
  redraw();      
  update_status_line();
}


void GeoWin::global_options()
{
  bool sh_grid     = show_grid;
  double new_scale = Wp->scale();
  int  grid_dist = (int)grid;
  int  max_dist  = (int)((Wp->xmax()-Wp->xmin())/4);
  int  min_dist  = (int)(new_scale/4);
  int  gr_style     = (int)Wp->get_grid_style();
  
  if (grid_dist > max_dist) max_dist = grid_dist;
  if (grid_dist < min_dist) min_dist = grid_dist;
  
  color bg_col = Wp->get_bg_color();
  
  string pm_name   = bg_pixmap;
  list<string> pm_list;
  
  int i;
  for(i=0; i < num_texture; i++) pm_list.append( name_texture[i] );
  pm_list.append("none");

  bool coords = show_coords;
  
  // change - get nrat from window ...
  bool nrat = Wp->get_normalize_rational_objects();
  
  int  mfct = (int) mouse_fct;
  
  int an_steps = gw_anim_steps;
  
  panel P;
  P.set_inf(this);  
  P.text_item("\\bf\\blue Grid Settings");
  P.bool_item("show grid", sh_grid);
  P.int_item("grid dist", grid_dist, min_dist, max_dist); 
  P.choice_item("grid style",gr_style,"invisible", "points", "lines");
  P.text_item("");
  
  P.text_item("\\bf\\blue Window Options");
  P.string_item("bg texture", pm_name, pm_list, 5);
  P.color_item("bg color",bg_col);
  P.int_item("animation", an_steps,1,30);
  P.bool_item("normalize before output", nrat);
  P.text_item("");
  P.bool_item("show position", coords);
  P.text_item("");


  P.text_item("\\bf\\blue Mouse Options"); 
  P.int_item("capture width", mfct, 1, 15);
  P.text_item("");  
  
  P.fbutton("apply",   APPLY_BUTTON);
  P.button("defaults", DEFAULTS_BUTTON);
  P.button("cancel",   CANCEL_BUTTON);
  
  int but;
  
  while ((but=open_panel(P)) == DEFAULTS_BUTTON )
    { 
      grid_dist = (int)(10/Wp->scale());
      sh_grid   = false;
      gr_style  = point_grid;
      
      bg_col    = white;
      pm_name   = "none";
      
      coords    = DEFAULT_show_coords;

      mfct      = 1;
      nrat      = false;
      
      an_steps  = DEFAULT_anim_steps;
    }
  
  if( but == CANCEL_BUTTON ) return;
  
  Wp->set_flush(false);
  
  if( sh_grid ) { grid = grid_dist; show_grid = true; }
  else         { grid = grid_dist; show_grid = false; }
  double g = grid;
  if( show_grid && g*new_scale < 4 ) g = 6/new_scale;
  if( show_grid ) { if (g>0) g=-g;  Wp->set_grid_dist(g); }
  else            Wp->set_grid_mode(0);
  Wp->set_grid_style((grid_style)gr_style);
  
  Wp->set_bg_color(bg_col);
  i=0;
  string str;
  forall(str, pm_list) 
    if( str == pm_name ) break;
    else i++;
  
  char* pr = 0;
  if( i >= 0 && i < num_texture)
    pr = Wp->create_pixrect(xpm_texture[i]);
  
  pr = Wp->set_bg_pixmap(pr);
  if( pr ) Wp->del_pixrect(pr);
  bg_pixmap = pm_name;
  
  if( !show_coords && coords ) 
    Wp->set_show_coord_handler(GeoWin::coord_handler);
  if( show_coords && !coords )  Wp->set_show_coord_handler(0);
  
  gw_anim_steps = an_steps;
  
  // handle nrat ...
  Wp->set_normalize_rational_objects(nrat);
  
  show_coords = coords;

  mouse_fct = mfct;
  
  Wp->set_flush(true);
  
  redraw();      
  update_status_line();
}


void GeoWin::scaling_options()
{
  panel P;  

  P.text_item("\\bf\\blue Window Coordinates");

  double  x0  = Wp->xmin();
  double  x1  = Wp->xmax();
  double  y0  = Wp->ymin();
  double ox0=x0, ox1=x1, oy0=y0;
  
  P.double_item("minimum xcoord", x0);
  P.double_item("maximum xcoord", x1);
  P.double_item("minimum ycoord", y0);

  P.fbutton("apply",   APPLY_BUTTON);
  P.button("defaults", DEFAULTS_BUTTON);
  P.button("cancel",   CANCEL_BUTTON);
  
  int but;
  
  while ((but=open_panel(P)) == DEFAULTS_BUTTON  )
    {
      x0 =   0;
      x1 = 100;
      y0 =   0;
    }
  
  if( but == CANCEL_BUTTON ) return;
  
  if( but && ( x0!=ox0 || x1!=ox1 || y0!=oy0) )
    {
      init(x0, x1, y0);
      redraw();
    }
}

void GeoWin::d3_options()
{
  panel P;  
  bool elim= gw_elim, solid= gw_solid, show_cs = gw_d3_cs;
  int sp= gw_anim_speed;
  string xstr;
  list<string> projection_parameters;
  int zv = (int) gw_d3_zvalue;
  
  projection_parameters.append(string("xy")); projection_parameters.append(string("xz")); 
  projection_parameters.append(string("yz")); 
  
  switch(GeoWin::get_projection_mode()){
   case 0:  { xstr = string("xy"); break; }
   case 1:  { xstr = string("xz"); break; }
   default:{ xstr = string("yz"); break; }
  }

  P.text_item("\\bf\\blue Options for 3d Output");
  P.text_item("");
  P.text_item("\\bf\\blue Elimination"); 
  P.bool_item("Eliminate ?", elim);
  P.text_item("");

  P.text_item("\\bf\\blue Solid Output"); 
  P.bool_item("Solid ?", solid);
  P.text_item("");

  P.text_item("\\bf\\blue Animation Speed"); 
  P.int_item("Speed", sp, 1, 32);
  P.text_item("");

  P.text_item("\\bf\\blue Projection");
  P.string_item("into plane:", xstr, projection_parameters, 8);
  P.text_item(""); 
  
  P.text_item("\\bf\\blue Z-value for input of new objects");
  P.int_item("z:",zv,-200,200);
  
  P.text_item("\\bf\\blue Coordinate system");
  P.bool_item("show?", show_cs);
  P.text_item(""); 

  P.fbutton("apply",   APPLY_BUTTON);
  P.button("defaults", DEFAULTS_BUTTON);
  P.button("cancel",   CANCEL_BUTTON);
  
  int but;
  
  while ((but=open_panel(P)) == DEFAULTS_BUTTON  ) {
    elim= DEFAULT_gw_elim; solid= DEFAULT_gw_solid;
    sp= DEFAULT_gw_anim_speed;
    show_cs=true;
    zv = 0;
  }
  if( but == CANCEL_BUTTON ) return;

  gw_elim= elim; gw_solid= solid; 
  gw_anim_speed= sp;
  gw_d3_cs = show_cs;
  gw_d3_zvalue = (double) zv;
  
  if (xstr == string("xy")) GeoWin::projection_mode=0;
  else
   if (xstr == string("xz")) GeoWin::projection_mode=1;
   else
    if (xstr == string("yz")) GeoWin::projection_mode=2;
    
  redraw();
}

void GeoWin::add_text(geo_scene sc, const geowin_text& gt)
{
 (sc->additional_objects).append(gt);
}

void GeoWin::remove_texts(geo_scene sc)
{
 (sc->additional_objects).clear();
}

void GeoWin::msg_open(string str)
{
  if( !msg_win->is_open() )  msg_win->open(0, 0);
  msg_win->message(str);
}

void GeoWin::message(string str)
{
  if (str==string("")){
   double x0, x1, y0, y1;
   cur_info = string("");
  
   x0 = Wp->xmin();
   x1 = Wp->xmax(); 
   y0 = Wp->ymin();
   y1 = Wp->ymax();
  
   redraw(x0, y0, x1, y1);  
  }
  else {
   if( !(cur_info == string("")) ) message("");
   cur_info = str;
   draw_info();
  }
}


void GeoWin::draw_info()
{
  double x0, x1, y;
  x0 = Wp->xmin() + Wp->pix_to_real(12);
  x1 = Wp->xmax() - Wp->pix_to_real(12);
  y  = Wp->ymax() - Wp->pix_to_real(2);
  
  Wp->text_box(x0, x1, y, cur_info );    
}

bool GeoWin::get_pin_point( point& p ) const
{
  if (pin_point) 
    {
      p = *pin_point;
      return true;
    }
  return false;
}

void GeoWin::set_pin_point( point p )
{
  if (pin_point)  delete pin_point;
  pin_point = new point(p); 
}

void GeoWin::del_pin_point( )
{
  if (pin_point)  delete pin_point;
  pin_point = 0; 
}

void GeoWin::draw_pin()
{
  if( pin_point )
    {
      int oldw = Wp->set_line_width(1);

      double r = Wp->pix_to_real(4);
      double d = Wp->pix_to_real(3);
      double x = pin_point->xcoord();
      double y = pin_point->ycoord();
      Wp->draw_segment(x-d, y-d, x+d, y+d, black);
      Wp->draw_segment(x-d, y+d, x+d, y-d, black);
      Wp->draw_circle(x, y, r, black);
      
      Wp->set_line_width(oldw);
    }
}


void GeoWin::rotate()
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;

  point p1;
  if( get_pin_point(p1) ) rotate(p1);
  else
    {
      string str("\\bf \\blue ROTATE : \\black click the rotation center.");
      message(str);

      Wp->read_mouse(p1);
      set_pin_point(p1);
      draw_pin();
      rotate(p1);
      del_pin_point();
      redraw();
    }
}

void GeoWin::rotate(const point& p1)
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;
  
  string str("\\bf \\blue ROTATE : \\black select a second point and drag.");
  message(str);
  
  point p2, p3;
  Wp->read_mouse(p2);
  rotate(p1, p2);

  message("");
}

void GeoWin::rotate(const point& p1, const point& p)
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;
  if( !cur->start_changing() )
    error_handler(0, "Changing selected objects is not allowed !");

  point p2 = p;
  point p3;
  segment s1, s2;
  double x1,y1;
  unsigned long t = Wp->button_press_time();
  int event;
  
  do 
    { 
      int val;
      unsigned long t1;
      event = Wp->read_event(val,x1,y1,t1);
      if (event == motion_event) if (t1 - t < 50) continue;
      p3 = point(x1, y1);
      s1 = segment(p1, p2);
      s2 = segment(p1, p3);
      ed->move_or_rotate_sel(p1.xcoord(), p1.ycoord(), s1.angle(s2), false, update_while_object_dragging);
      
      if (! update_while_object_dragging) redraw_and_update_status_line();  
      
      t = t1;      
      p2 = p3;
    }
  while( event != button_release_event);
  
  p3 = point(x1, y1);
  s1 = segment(p1, p2);
  s2 = segment(p1, p3);
  ed->move_or_rotate_sel(p1.xcoord(), p1.ycoord(), s1.angle(s2), false, true);
  
  cur->stop_changing();
}

void GeoWin::move()
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;
  
  string str("\\bf \\blue MOVE : \\black determine a point and drag.");
  message(str);

  point p1;
  Wp->read_mouse(p1);
  move(p1);
  
  message("");
}

void GeoWin::show_sel()
{
   if (cur) {
     cur->call_buffer_fcn(0,false);
   }
}

void GeoWin::paste_from()
{
   if (cur) {
    if (edit_scene) { 
     panel P("Paste from scene");
     P.buttons_per_line(1);
     
     // select one of these scenes, get the typename
     P.text_item("\\bf\\blue Paste from ");
     
     int i=0;
     geo_scene sc;
     
     list<geo_scene> scl;
     
     forall(sc , scenes)  
     { 
      if( ((sc->owner).empty()) && (sc->get_type_name() == edit_scene->get_type_name()) ) 
        scl.append(sc);
     }
     
     if (! scl.empty()) {
      forall(sc,scl) P.button(sc->get_name(),i++);

      P.button("cancel",   i);

      window& W = get_window();
      W.disable_panel();
      int but = P.open();
      W.enable_panel();
      redraw_inp_panel();
      
      if (but==i) return;
 
      geo_scene from_scene = scl[scl.get_item(but)];
      edit_scene->copy_clip_buffer(from_scene);
      
      // compare typename of selected scene with typename of edit_scene ...
     }
    }
   }
}


void GeoWin::move(const point& p)
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;
  if( !cur->start_changing() )
    error_handler(0, "Changing selected objects is not allowed !");

  double x = p.xcoord();
  double y = p.ycoord();

  double dx, dy, x1,y1;
  unsigned long t = Wp->button_press_time();
  int event;
     
  do 
  { 
    int val;
    unsigned long t1;
    event = Wp->read_event(val,x1,y1,t1);
    if (event == motion_event) if (t1 - t < 50) continue;
    dx = x1-x; 
    dy = y1-y;
    ed->move_or_rotate_sel(dx, dy, 0, true, update_while_object_dragging);
    
    if (! update_while_object_dragging) redraw_and_update_status_line(); 
    
    t = t1;
    x = x1;
    y = y1;
  }
  while( event != button_release_event);

  // get_window().get_mouse(x1,y1);

  dx = x1-x; 
  dy = y1-y;
  ed->move_or_rotate_sel(dx, dy, 0, true, true);
  
  cur->stop_changing();
}

void GeoWin::select()
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;
  
  string str("\\bf \\blue SELECT : \\black determine a corner and drag.");
  message(str);
  
  point p;
  if ( Wp->read_mouse(p) != MOUSE_BUTTON(1) ) return;
  select(p);

  message("");
}


void GeoWin::select_all()
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return; 

  cur->select_in_rect(false,0,0,0,0,true);
}

void GeoWin::select(const point& p0)
{
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( !ed ) return;
  
  double x0 = p0.xcoord();
  double y0 = p0.ycoord();
  //unsigned long t = Wp->button_press_time();
  int event;
  double x1,y1;
  
  char* buf = Wp->get_window_pixrect();
  Wp->start_buffering();
  Wp->set_line_width(1);
  
  do
    { 
      
      int val;
      unsigned long t1;
      event = Wp->read_event(val,x1,y1,t1);
      
      Wp->put_pixrect(buf);
      Wp->draw_rectangle(x0, y0, x1, y1); 
      Wp->flush_buffer();
    }
  while(event != button_release_event);
  
  Wp->put_pixrect(buf);  
  Wp->flush_buffer();
  
  Wp->del_pixrect(buf); 
  Wp->stop_buffering();
  
  double tmp;
  if( x1 < x0 ) { tmp = x0; x0 = x1; x1 = tmp; }
  if( y1 < y0 ) { tmp = y0; y0 = y1; y1 = tmp; }
  
  cur->select_in_rect(false, x0, y0, x1, y1,false); 
}

void GeoWin::unselect()
{ 
  geo_editor ed = cur ? cur->type_editor() : 0;
  if( ed ) cur->select_in_rect(true,0,0,0,0,true);
}

// d3 - support in GeoWin using d3_window ...

void GeoWin::show_d3()
{
  window& W = get_window();
  int nw = W.get_node_width();
  
  GRAPH<d3_point,int> H;
  // create d3 window ...
  d3_win = new d3_window(W,H);
  
  //new - init node and edge color map ...
  node_color_map.clear();
  edge_color_map.clear();
  
  d3_win->set_node_width(4);
  d3_win->set_elim(gw_elim);
  d3_win->set_solid(gw_solid);
  d3_win->set_speed(gw_anim_speed);
  d3_win->set_message("");
  
  //new 
  d3_win->set_draw_edges(gw_d3edges);
  
  geo_scene scn;
  int but; 

  forall(scn,scenes) {
    if (get_visible(scn) || is_active(scn)){
       GRAPH<d3_point,int> H2;

       color cact = get_color(scn);
       
       scn->compute_d3_output(scn, (*d3_win), H2);
       int nn = H2.number_of_nodes();
       int ne = H2.number_of_edges();       
       H.join(H2);
       
       if (get_d3_use_scene_color(scn)){
         // get node and edge lists of H ...
         list<node>& Lnode = (list<node>&) H.all_nodes();
         list<edge>& Ledge = (list<edge>&) H.all_edges();
        
         // set node and edge colors ...
         list_item it;
         it = Lnode.last();	
         for(;nn>0; it = Lnode.cyclic_pred(it),nn--) {  
           if (! node_color_map.defined(Lnode[it])) node_color_map[Lnode[it]] = cact; 
         }
       
         it = Ledge.last();
         for(;ne>0; it = Ledge.cyclic_pred(it),ne--) {  
           if (! edge_color_map.defined(Ledge[it])) edge_color_map[Ledge[it]] = cact; 
         }  
       }    
    }   
  }
  
  // coordinate system ...
  edge e1=NULL,e2,e3=NULL,e4,e5=NULL,e6;
  if (gw_d3_cs) {
  
  double v1 = get_xmax() - get_xmin();
  double v2 = get_ymax() - get_ymin();
  
  double w = (v1 > v2) ? v2/1.5 : v1/1.5 ;
  double val2 = w/25;
  
  node nv1 = H.new_node(d3_point(0,0,0));
  node nv2 = H.new_node(d3_point(w,0,0));
  node nv3 = H.new_node(d3_point(0,w,0));
  node nv4 = H.new_node(d3_point(0,0,w));  
  
  e1 = H.new_edge(nv1,nv2), e2 = H.new_edge(nv2,nv1);
  H.set_reversal(e1,e2);
  e3 = H.new_edge(nv1,nv3), e4 = H.new_edge(nv3,nv1);
  H.set_reversal(e3,e4);
  e5 = H.new_edge(nv1,nv4), e6 = H.new_edge(nv4,nv1);  
  H.set_reversal(e5,e6);
  

  // x, y, z ...
  node xn1 = H.new_node(d3_point(w-val2,-1.7*val2,0));
  node xn2 = H.new_node(d3_point(w,-1.7*val2,0));
  node xn3 = H.new_node(d3_point(w-val2,-0.5*val2,0));  
  node xn4 = H.new_node(d3_point(w,-0.5*val2,0));  
  edge xe1 = H.new_edge(xn1,xn4), xe2 = H.new_edge(xn4,xn1);
  edge xe3 = H.new_edge(xn2,xn3), xe4 = H.new_edge(xn3,xn2);  
  H.set_reversal(xe1,xe2);
  H.set_reversal(xe3,xe4);
  node yn1 = H.new_node(d3_point(-1.5*val2,w - 1.5*val2,0));
  node yn2 = H.new_node(d3_point(-0.5*val2,w,0));
  node yn3 = H.new_node(d3_point(-1.0*val2,w -0.75*val2,0));
  node yn4 = H.new_node(d3_point(-1.5*val2,w,0));
  edge ye1 = H.new_edge(yn1,yn2), ye2 = H.new_edge(yn2,yn1);
  edge ye3 = H.new_edge(yn3,yn4), ye4 = H.new_edge(yn4,yn3);
  H.set_reversal(ye1,ye2);
  H.set_reversal(ye3,ye4);
  node zn1 = H.new_node(d3_point(val2/2,-1.7*val2,w-val2/2));
  node zn2 = H.new_node(d3_point(-val2/2,-1.7*val2,w-val2/2));
  node zn3 = H.new_node(d3_point(val2/2,-0.5*val2,w-val2/2)); 
  node zn4 = H.new_node(d3_point(-val2/2,-0.5*val2,w-val2/2));
  edge ze1 = H.new_edge(zn1,zn2), ze2 = H.new_edge(zn2,zn1);
  edge ze3 = H.new_edge(zn3,zn2), ze4 = H.new_edge(zn2,zn3);
  edge ze5 = H.new_edge(zn3,zn4), ze6 = H.new_edge(zn4,zn3);
  H.set_reversal(ze1,ze2);
  H.set_reversal(ze3,ze4);
  H.set_reversal(ze5,ze6);
  }
 
  
  W.disable_panel();
  string msg("\\bf\\blue ~~~~~left:  \\black zoom up \
              \\bf\\blue ~~~~~middle:\\black zoom down \
              \\bf\\blue ~~~~~right: \\black return to geowin");

  W.set_show_coord_handler(nil);
  W.set_status_string(msg);
 
  node_array<vector> pos(H);
  node v;
  forall_nodes(v,H) pos[v] = H[v].to_vector();
  d3_win->init(pos);

  // new ...
  edge e;
  forall_defined(e,edge_color_map) {  d3_win->set_color(e, edge_color_map[e]); }
  forall_defined(v,node_color_map) {  d3_win->set_color(v, node_color_map[v]); } 

  //arrows 
  if (gw_d3_cs) {
  d3_win->set_arrow(e1,true);
  d3_win->set_arrow(e3,true);
  d3_win->set_arrow(e5,true);
  }
  
  
  d3_win->draw();
  
  do {
     but = d3_win->move();
  }
  while (but != MOUSE_BUTTON(3));
  
  delete d3_win;
  d3_win = NULL;

  W.set_show_coord_handler(GeoWin::coord_handler);
  W.enable_panel();
  redraw_inp_panel();
  W.set_node_width(nw);
  
  
  update_status_line();
  redraw();
}


void GeoWin::hide_d3()
{
  if (d3_win != NULL) 
  { delete d3_win; 
    d3_win=NULL; 
  }
}


// war version.c
double GeoWin::version()             { return GeoWinVersion; }
double GeoWin::fileformatversion()   { return GeoWinFileFormat; }

// setting contents function ...
void GeoWin::set_contents_fcn(geo_scene sc, bool (*f)(geo_scene, list<string>&))
{
  sc->contents_fcn = f;
}

// Scene group operations ...
geo_scenegroup GeoWin::new_scenegroup(string name)
{
  geo_scenegroup gsg = new GeoSceneGroup(name);
  scene_groups.append(gsg);
  return gsg;
}

geo_scenegroup GeoWin::new_scenegroup(string name, const list<geo_scene>& LS)
{
  geo_scenegroup gsg = new GeoSceneGroup(name, LS);
  scene_groups.append(gsg);
  return gsg;
}

void GeoWin::insert(geo_scenegroup GS, geo_scene sc)
{ 
  GS->group_members.append(sc);
  geo_scene iter;
  forall(iter,GS->dependent_scenes){
    add_dependence(sc,iter);
  }
}

bool GeoWin::del(geo_scenegroup GS, geo_scene sc)
{
  list_item rmv=NULL,it;
  forall_items(it,GS->group_members)
  {
   geo_scene act = (GS->group_members)[it];
   if (act == sc) rmv=it;
  }
  
  if (rmv == NULL) return false; // scene is not in the group ...
  
  // remove dependences ...
  GS->group_members.del_item(rmv);
  geo_scene iter;
  forall(iter,GS->dependent_scenes){
    del_dependence(sc,iter);
  }  
  
  return true;
}

void GeoWin::call_update_fcn(void (*fcn)(...),void* input, void* output,const list<void*>& infl)
{
  //cout << "call update function !\n";
  int sz = infl.size();
  
  switch(sz){
   case 0: { fcn(input,output); break; }
   case 1: {
     void* c1 = infl[infl.get_item(0)];
     fcn(input,c1,output);
     break;
   }
   case 2: {
     void* c1 = infl[infl.get_item(0)];
     void* c2 = infl[infl.get_item(1)];
     fcn(input,c1,c2,output);
     break;
   }
   case 3: {
     void* c1 = infl[infl.get_item(0)];
     void* c2 = infl[infl.get_item(1)];
     void* c3 = infl[infl.get_item(2)];   
     fcn(input,c1,c2,c3,output);   
     break;
   }
   case 4: {
     void* c1 = infl[infl.get_item(0)];
     void* c2 = infl[infl.get_item(1)];
     void* c3 = infl[infl.get_item(2)];
     void* c4 = infl[infl.get_item(3)];    
     fcn(input,c1,c2,c3,c4,output);
     break;
   }
   case 5: {
     void* c1 = infl[infl.get_item(0)];
     void* c2 = infl[infl.get_item(1)];
     void* c3 = infl[infl.get_item(2)];
     void* c4 = infl[infl.get_item(3)];
     void* c5 = infl[infl.get_item(4)];    
     fcn(input,c1,c2,c3,c4,c5,output);   
     break;
   }  
   case 6: {
     void* c1 = infl[infl.get_item(0)];
     void* c2 = infl[infl.get_item(1)];
     void* c3 = infl[infl.get_item(2)];
     void* c4 = infl[infl.get_item(3)];
     void* c5 = infl[infl.get_item(4)];    
     void* c6 = infl[infl.get_item(5)];
     fcn(input,c1,c2,c3,c4,c5,c6,output);   
     break;
   }     
   case 7: {
     void* c1 = infl[infl.get_item(0)];
     void* c2 = infl[infl.get_item(1)];
     void* c3 = infl[infl.get_item(2)];
     void* c4 = infl[infl.get_item(3)];
     void* c5 = infl[infl.get_item(4)];    
     void* c6 = infl[infl.get_item(5)];
     void* c7 = infl[infl.get_item(6)];
     fcn(input,c1,c2,c3,c4,c5,c6,c7,output);   
     break;
   }       
    
   default:{ cout << "Too many influencer scenes !\n"; }
  }
}


void GeoWin::add_import_object(geo_scene sc, geowin_import& io, string name, string desc)
{
  io.name=name; io.description=desc;
  sc->import_objects.append(&io);
}
  
void GeoWin::add_export_object(geo_scene sc, geowin_export& eo, string name, string desc)
{
  eo.name=name; eo.description=desc;
  sc->export_objects.append(&eo);
}


void GeoWin::set_font(geowin_font_type font_t, double size, string fn_user)
{ 
  if (! Wp) return;
  
  int sz = Wp->real_to_pix(size);
  string font_name;
  switch (font_t) {
    case roman_font:  font_name = string("T%d",sz);
                      break;
    case bold_font:   font_name = string("B%d",sz);
                      break;
    case italic_font: font_name = string("I%d",sz);
                      break;
    case fixed_font:  font_name = string("F%d",sz);
                      break;
    case user_font:   font_name = fn_user;
                      break;
  }
  Wp->set_font(font_name);
}

// former geo_rem (text objects)

void geowin_text::construct()
{
  x_offs = 0;  y_offs = 0;
  size = 10; 
  remark = "";
  ft_type = roman_font; 
  user_font = "";
  clr = black;  
  position_fcn = NULL;
}

geowin_text::geowin_text(string t, double ox, double oy, geowin_font_type ft, double sz, string uf, color c)
{
  construct();
  remark = t;
  x_offs = ox; y_offs = oy;
  ft_type = ft;
  size = sz;  
  user_font = uf;
  clr = c;
}

geowin_text::geowin_text(string t, geowin_font_type ft, double sz)
{
  construct();
  remark = t;
  ft_type = ft;
  size = sz;
}

string geowin_text::get_text() const 
{ return remark; }

string geowin_text::set_text(const string& n)
{ string prev = remark;
  remark = n;
  return prev;
}  

double geowin_text::get_size() const 
{ return size; }

double geowin_text::set_size(double n)
{ double prev = size;
  size = n;
  return prev;
}  

geowin_font_type geowin_text::get_font_type() const
{ return ft_type; }

geowin_font_type geowin_text::set_font_type(geowin_font_type n)
{ geowin_font_type prev = ft_type;
  ft_type = n;
  return prev;
}  

string geowin_text::get_user_font() const
{ return user_font; }

string geowin_text::set_user_font(const string& n)
{ string prev = user_font;
  user_font = n;
  return prev;
}  

color geowin_text::get_color() const
{ return clr; }

color geowin_text::set_color(color n)
{ color prev = clr;
  clr = n;
  return prev;
}  

double geowin_text::get_x_offset() const 
{ return x_offs; }

double geowin_text::set_x_offset(double n)
{ double prev = x_offs;
  x_offs = n;
  return prev;
}  

double geowin_text::get_y_offset() const 
{ return y_offs; }

double geowin_text::set_y_offset(double n)
{ double prev = y_offs;
  y_offs = n;
  return prev;
} 

GEOWIN_END_NAMESPACE








