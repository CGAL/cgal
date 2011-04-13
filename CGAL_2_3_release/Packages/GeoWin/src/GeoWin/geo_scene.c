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
// file          : src/GeoWin/geo_scene.c
// package       : GeoWin (1.2.2)
// revision      : 1.2.2
// revision_date : 30 January 2001 
// author(s)     : Matthias Baesken, Ulrike Bartuschka, Stefan Naeher
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ============================================================================

#include<LEDA/geowin.h>
#include "geo_localdefs.h"

GEOWIN_BEGIN_NAMESPACE

// ***************************************************************
// ************* GeoScenePrototypes ******************************

int GeoScenePrototypes::SceneCount = 0;

bool GeoScenePrototypes::st_alloc = false; 

static GeoScenePrototypes* tname_map_ptr = 0;

GeoScenePrototypes::GeoScenePrototypes() 
  : allscenes(MAXIMUM_SCENETYPES)
{
  st_alloc=false;
  allscenes.init((GeoScene*)0);
}

GeoScenePrototypes::~GeoScenePrototypes()
{
  st_alloc=false; 
  for( int i = 0; i < SceneCount; i++ )
    delete allscenes[i];
}


int GeoScenePrototypes::RegisterScene(geo_scene sc)
{
  if (tname_map_ptr == 0) { 
                            tname_map_ptr = new GeoScenePrototypes; 
                            st_alloc=true; }

  if( SceneCount < MAXIMUM_SCENETYPES )
    {
      tname_map_ptr->allscenes[SceneCount] = sc;
      string str("Type%dScene", SceneCount+1);
      sc->set_base_name(str);

      return SceneCount++;
    }
  else error_handler(0, "Maximum scene count reached");
  
  return -1;
}

void GeoScenePrototypes::FreeMem()
{
  if (st_alloc) delete tname_map_ptr;
}

geo_scene GeoScenePrototypes::get_prototype(int id) 
{ 
  if (tname_map_ptr && id >= 0 && id < SceneCount) 
     return tname_map_ptr->allscenes[id];
  return 0;
}

geo_scene GeoScenePrototypes::get_prototype(string str) 
{ 
  if (tname_map_ptr == 0) return 0;

  list<string> choice;
  string tmp;
  
  for( int i = 0; i < SceneCount; i++ )
    {
      tmp = GeoScenePrototypes::get_type_name(i);
      if( tmp == str ) return  tname_map_ptr->allscenes[i];
      choice.append(tmp);
    }
  
  cout << "There is no prototype for a scene of type " << str << " known. ";
  cout << endl << "Possible types are : " << endl;
  forall(tmp, choice) cout << tmp << endl;
  
  return 0;
}

void GeoScenePrototypes::set_type_name(int id, string str)
{ 
  if(tname_map_ptr && id>=0 && id<SceneCount) 
      tname_map_ptr->allscenes[id]->set_base_name(str);
}

string GeoScenePrototypes::get_type_name(int id) 
{ 
  if(tname_map_ptr && id>=0 && id<SceneCount) 
    return tname_map_ptr->allscenes[id]->get_base_name();
  return string("");
}

list<string> GeoScenePrototypes::get_all_editables()
{
  list<string> L;
  if (tname_map_ptr == 0) return L;
  for( int i = 0; i < SceneCount; i++ )
    { geo_scene  sc = tname_map_ptr->allscenes[i];
      geo_editor ed = sc->type_editor();
      if( ed ) L.append(sc->get_base_name());
    }
  return L;
}

// ***************************************************************
// ****************** GeoScene ***********************************

  
void GeoScene::init_default_values()
{
 col1=gw->DEFAULT_color;
 col2=gw->DEFAULT_color2; 
 filling_color=gw->DEFAULT_fill_color; 
 text_color=gw->DEFAULT_color;     
 back_line_width=gw->DEFAULT_line_width;       
 active_line_width=gw->DEFAULT_active_line_width;  
 l_style=gw->DEFAULT_line_style;        
 p_style=gw->DEFAULT_point_style;   
 
 //selection parameters ...
 sel_filling_color = gw->DEFAULT_fill_color; 
 sel_l_style = gw->DEFAULT_line_style;
 sel_line_width = gw->DEFAULT_line_width;      
}


GeoScene::GeoScene(GeoWin* win) : gw(win) 
{
  set_name(string(""));
  description       = string("");
  
  pos               = 0;
  edit_mode         = 0;
  z                 = 0;
  fcn_state         = 0;
  draw_mode_number  = 0;
  col1              = black;        
  col2              = red;   
  filling_color     = ivory;
  cyclic_colors     = false;         
  back_line_width   = 1;     
  active_line_width = 1;      
  l_style           = solid;      
  p_style           = cross_point;  
  copy_attr         = false;    
  visible           = false;
  active            = false;
  edit_menu_type    = -1;
  while_dragging    = false;
  mouse_changing    = false;
  redraw_pt         = NULL;
  debug_mode        = false;
  debug_file_name   = string("debug.geo");
  d3_use_scene_color= true;
  
  //selection parameters ...
  sel_filling_color = ivory;
  sel_l_style = solid;
  sel_line_width = 1;
  
  
  contents_fcn      = NULL;
  
  limit = -1;
  
  int i;
  for(i = 0; i<geo_max_ids; i++) ids[i] = -1;
  
  // client data ...
  for(i=0; i<16; i++) cl_data[i] = NULL;
  
  labels[geo_col0_label] = "interior";
  labels[geo_col1_label] = "boundary";
  labels[geo_col2_label] = "selected";

  init_d3_window = NULL;
  
  // removed from GeoBaseScene / GeoEditScene
   myobjs        = true;
   base_mode     = true;
   change_mode   = -1;			   
   high_light = false;
   mouse_obj_exist = false;
   base_mode = false;  
   
   short_name = "";
   
   draw_fcn = 0;
    
   color1_map = 0;
   color2_map = 0;
   lst_map    = 0;
   lw_map     = 0;
   label_map  = 0;
   text_map   = 0;
   original   = 0;    
   
   box_intersection_fcn = 0;
   get_bounding_box_fcn = 0;
   move_fcn = 0;
   rotate_fcn = 0;  
   info_fcn = 0;
   objects_init_d3_window = 0;
   move_point_fcn = 0; 
   edit_obj_fcn = 0;
   start_change_fcn = 0;
   end_change_fcn = 0;
   
   defining_points_fcn = NULL;
   handle_defining_points = geowin_highlight;
}

GeoScene::~GeoScene() 
{  
  additional_objects.clear();
  
  geo_scene scn;
  forall(scn,owner) {
     scn->del_dependence(this);
  }

  if( gw )    gw->remove(this);
  while( !results.empty() ) 
    {
      geo_scene sc = results.pop();
      // sc already deleted ?
      //(sc->owner) = 0;
      list_item del=0,lit;
      forall_items(lit,sc->owner)
        if ((sc->owner)[lit] == this) { del = lit; break; }
      if (del) (sc->owner).del(del);
      if ((sc->owner).empty()) delete sc;
    }
}

geo_scene GeoScene::new_scene(string str)
{
  geo_scene proto = GeoScenePrototypes::get_prototype(str);
  if (proto == NULL) return NULL;
  return proto->clone();
}

geo_scene GeoScene::new_scene(int id)
{
  geo_scene proto = GeoScenePrototypes::get_prototype(id);
  if (proto == NULL) return NULL;
  return proto->clone();
}

void GeoScene::set_name(string nm)
{
 name=nm;
}

int GeoScene::get_generator_number() const
{ return generator_number; }

void GeoScene::set_generator_number(int n)
{ generator_number = n; }

list<string> GeoScene::get_generator_names() const 
{ return generator_names; }

string GeoScene::get_name() const { return name; }

void GeoScene::set_base_name(string nm) { base_name = nm; }

string GeoScene::get_base_name() const { return base_name; }

void GeoScene::set_type_name(string str)
{
  int id = GeoSceneId();
  GeoScenePrototypes::set_type_name(id, str);
}

string GeoScene::get_type_name() const
{
  int id = GeoSceneId();
  //cout << "get_type_name - id:" << id << "\n";
  return GeoScenePrototypes::get_type_name(id);
}

color GeoScene::get_default_color1()          
{ 
 return col1; 
}

color GeoScene::get_sel_color()               
{ 
  return col2;
}

color GeoScene::get_default_color2()          
{
 return filling_color; 
}

point_style GeoScene::get_default_point_style() 
{ return p_style; }

color GeoScene::get_default_text_color()        
{ return text_color; }

line_style GeoScene::get_default_line_style() 
{ return l_style; }

int  GeoScene::get_default_line_width()       
{ return get_line_width(); }


bool GeoScene::IsMyKindOfScene(geo_scene sc)
{
  return sc->IsSuperOrEqual(this);
}

void GeoScene::init_from_prototype(geo_scene pr)
{
  set_base_name(pr->get_base_name());
  draw_object_parameters = pr->draw_object_parameters;
  info_fcn = pr->info_fcn;
  objects_init_d3_window = pr->objects_init_d3_window;
  box_intersection_fcn = pr->box_intersection_fcn;
  get_bounding_box_fcn = pr->get_bounding_box_fcn;    
  move_fcn = pr->move_fcn;
  rotate_fcn = pr->rotate_fcn;
  move_point_fcn = pr->move_point_fcn;  
  edit_obj_fcn = pr->edit_obj_fcn;
  start_change_fcn = pr->start_change_fcn;
  end_change_fcn = pr->end_change_fcn;  
  defining_points_fcn = pr->defining_points_fcn;
}

void GeoScene::init_data()
{
  update();
  if (gw){
    gw->redraw();
    gw->update_status_line();
  }
}

void GeoScene::edit()
{ if ( gw) gw->redraw(); }


// client data operations ...

void* GeoScene::set_client_data(void* cld,int i)
{
 if ((i>-1) && (i<16)){
   void* prev = cl_data[i]; 
   cl_data[i] = cld;
   return prev;
 }
 // i not in the expected range ...
 else return NULL; 
}

void* GeoScene::get_client_data(int i)
{
 if ((i>-1) && (i<16)) return cl_data[i];
 // i not in the expected range ...
 else return NULL;
}

void GeoScene::add_owner(GeoScene* sc)
{
 owner.append(sc);
}

void GeoScene::add_dependence(GeoScene* sc)
{
  results.append(sc);
  (sc->owner).append(this);
}

void GeoScene::del_dependence(GeoScene* sc)
{
  //cout << "remove " << sc->get_name() << " from " << this->get_name() << "\n";
  list_item lit, del=0;
  forall_items( lit, results )
    if( results[lit] == sc ) { del = lit; break; }
  if( del )
    {
      results.del_item(del);
      del=0;
      forall_items(lit,sc->owner)
        if ((sc->owner)[lit] == this) { del = lit; break; }
      if (del) (sc->owner).del(del);
      if ((sc->owner).empty()) delete sc;
    }
  else 
    {
      //cout << "error in remove " << sc->get_name() << " from " << this->get_name() << "\n";
    }
}

void GeoScene::redraw(window* w, double x1, double y1, double x2, double y2)
{
  // ----- moved code from base scenes ---------------
  //set calling GeoWin and scene 
  GeoWin::call_win = gw;
  GeoWin::call_scene = this;
  
  oldwl  = w->set_line_width(get_line_width());
  oldstl = w->set_line_style(l_style);
  oldpstl = w->set_point_style(p_style);
  bool b=true;

  call_redraw_objects(b,w,x1,y1,x2,y2);  
  
  // iterate on the scene container ...
  if (b) { // iteration on the container storing the objects ...
      init_iteration();
      void* hadr;
   
      while (no_end_reached()) {
        hadr = get_next_obj_adr();
        draw_object_with_properties(hadr);
        iterate();
      }   
  }    
  
  // -------------------------------------------------

  geo_scene sc;
  
  //draw all text objects ...
  geowin_text gt;
  forall(gt, additional_objects) {
    geowin_draw_text_object(*w, gw->get_xmin(), gw->get_xmax(), gw->get_ymin(), gw->get_ymax(),(void*) this, gt, false);
  }  

  forall( sc, results )
    {
      if( sc->is_visible() && sc->gw != gw )   sc->gw->redraw();
      if( sc->get_active() )  sc->gw->update_status_line();
    }

  // reset point style , line style, line width
  w->set_point_style(oldpstl);
  w->set_line_style(oldstl);
  w->set_line_width(oldwl);  
}

void GeoScene::show_points(const list<point>& L)
{
  window& w = gw->get_window();
  point_style pold=w.set_point_style(rect_point);
  color cold = w.set_fill_color(ivory);
    
  point piter;
  forall(piter,L) w << piter;
     
  w.set_point_style(pold);
  w.set_fill_color(cold);  
}

void GeoScene::write_postscript(ps_file& f)
{
  window& w = gw->get_window();
  geowin_text gt;
  forall(gt, additional_objects) {
    geowin_draw_text_object(w, f, gw->get_xmin(), gw->get_xmax(), gw->get_ymin(), gw->get_ymax(),(void*) this, gt, false);
  }    
}

void GeoScene::update()
{
  geo_scene sc;
  forall( sc, results ) {
    //cout << "  " << sc->get_name() << "\n";
    if( !while_dragging || sc->is_visible() || (! (sc->results).empty()))  sc->update();
  }
    
  // clean up lists of added/ deleted objects ...
  clear_op_objs();
}

void GeoScene::update_and_redraw()
{
  update();
  gw->redraw_and_update_status_line();
}

string GeoScene::information()
{
  string str = string("\\blue \\bf ~~%s~: ~~~\\black \\tt ", get_name());
  void* objs = get_untyped_objptr();
  
  if( info_fcn && objs!=NULL) str += info_fcn(objs);
  return str;  
}

void GeoScene::mouse_at_pos(double, double, long& MASK)
{
  gw->set_cursor(-1);
}

void GeoScene::split_and_append(string s, list<string>& LS, int cnt)
{
  int posact=0, posmax= s.length()-1, pos2=0;
  
  for(;posact< posmax+1;posact=posact+cnt){
    pos2 = posact+cnt-1;
    if (pos2 > posmax) pos2=posmax;
    LS.append(s(posact,pos2)); 
  }  
  LS.append(string(" "));
}

void GeoScene::contents()
{
  list<string> Ls;
  bool b = false;

  if (contents_fcn != NULL) {
    b = contents_fcn(this, Ls);
  }
  
  if (b || (contents_fcn==NULL))
  {
    string_ostream str_bf;
    write(str_bf,"\n");
    str_bf << ends;
    string h = str_bf.str();
    string sh("\n");
  
    int posakt=0, posmax= h.length()-1,act;
    act=h.pos(sh,posakt);
  
    while(act != -1){
     split_and_append(h(posakt,act-1),Ls,77);
     posakt=act+1;
     act=h.pos(sh,posakt);
    }
  
    split_and_append(h(posakt,posmax),Ls,77);
  }
  LedaFileViewer(Ls);
}

void GeoScene::scene_description()
{
  if (description == "") return;
 
  if (gw) {
    panel P("LEDA PANEL");
    int start =0, posi;
    do {
     posi = description.pos(string("\n"),start);
     if (posi != -1) { P.text_item(description(start,posi)); P.text_item(""); }
     else P.text_item(description(start,description.length()-1));
     
     start = posi+1;
    }
    while (posi != -1);
    
    P.button("OK");
    gw->open_panel(P);
  }
}


bool GeoScene::options(panel* p)
{ 
  panel* P = p;
  int own_but;
  
  if( !p )
    {
      P = new panel("Scene Options");
      P->set_inf(&own_but);
    }
   
  color c0 = filling_color; 
  color c1 = col1; 
  color c2 = col2;
  color c3 = sel_filling_color;
  
  bool enum_col = cyclic_colors;
  bool dbg = debug_mode;
  string dbs = debug_file_name;
  string xstr("");

  int aw = active_line_width, bw = back_line_width;
  line_style lst = l_style;
  point_style pst= p_style;
  
  int sel_lw = sel_line_width;
  line_style sel_lst = sel_l_style;
  
  P->text_item("\\bf\\blue Colors");
  P->color_item(labels[geo_col0_label], c0);
  P->color_item(labels[geo_col1_label], c1);
  P->color_item("selected interior", c3);
  P->color_item("selected boundary", c2);
  P->bool_item("enum colors in edit scenes",enum_col);
  P->text_item("");

  P->text_item("\\bf\\blue Points"); 
  P->pstyle_item("point style", pst);
  
  P->text_item("\\bf\\blue Lines");
  P->lstyle_item("line style", lst);
  P->lstyle_item("selected line style", sel_lst);
  P->lwidth_item("line width", bw);
  P->lwidth_item("selected line width", sel_lw);
  P->lwidth_item("active width", aw);
  P->text_item("");
  
  if (! draw_object_parameters.empty()) {
   xstr = draw_object_parameters[draw_object_parameters.get_item(draw_mode_number)];
   P->text_item("\\bf\\blue Draw function");
   P->string_item("Drawing mode:", xstr, draw_object_parameters, 8);
   P->text_item("");
  }
  
  P->text_item("\\bf\\blue Debug");      
  P->bool_item("Output input scene", dbg);
  P->string_item("File name", dbs);
  P->text_item("");  
  
  P->fbutton("apply",   APPLY_BUTTON);
  P->button("cancel",   CANCEL_BUTTON);
  
  int& butt = *((int*)(P->get_inf()));
 
  butt = gw->open_panel(*P);
  
  if( butt == CANCEL_BUTTON ) return false;
  
  filling_color = c0;
  sel_filling_color = c3;
  col1 = c1;
  col2 = c2;
  active_line_width = aw;
  back_line_width = bw;
  sel_line_width = sel_lw;
  l_style = lst;
  sel_l_style = sel_lst;
  p_style = pst;
  debug_mode = dbg;
  cyclic_colors = enum_col;
  debug_file_name = dbs;
  
  if (! draw_object_parameters.empty()) {
   if (xstr.length() > 0){ // a new draw mode was choosen
    string iter; 
    int nb=0, counter=0;
    forall(iter, draw_object_parameters){
      if (iter == xstr) nb=counter; counter++;
    }
    draw_mode_number = nb;
   }
  }

  if(!p) gw->redraw();  
  return true;
}

  
void GeoScene::get_export_object_names_and_descriptions(list<string>& N, list<string>& D)
{  geowin_export* iter;
   forall(iter, export_objects) { N.append(iter->get_name()); D.append(iter->get_description()); }
}
  
void GeoScene::get_import_object_names_and_descriptions(list<string>& N, list<string>& D)
{  geowin_import* iter;
   forall(iter, import_objects) { N.append(iter->get_name()); D.append(iter->get_description()); }
}   

// was GeoProp ------------------------------------------------------------------

color GeoScene::get_obj_color(void* o) 
{
  color c = get_default_color1();
  void* obj = (original && original->defined(o)) ? (*original)[o] : o;
  if( color1_map && color1_map->defined(obj)) c = (*color1_map)[obj];
  return c; 
}

color GeoScene::set_obj_color(void* obj, color c) 
{
  color cold = get_obj_color(obj);
  if( !color1_map ) color1_map = new map<void*, color>;
  (*color1_map)[obj] = c;
  return cold;
}

color GeoScene::get_obj_fill_color(void* o) 
{
  color c = get_default_color2();
  void* obj = (original && original->defined(o)) ? (*original)[o] : o;
  if( color2_map && color2_map->defined(obj)) { c = (*color2_map)[obj]; }

  return c; 
}

color GeoScene::set_obj_fill_color(void* obj, color c) 
{
  color cold = get_obj_fill_color(obj);
  if( !color2_map ) color2_map = new map<void*, color>;
  (*color2_map)[obj] = c;
  return cold;
}

  
line_style GeoScene::get_obj_line_style(void* o)   
{
  line_style lst = get_default_line_style();
  void* obj = (original && original->defined(o)) ? (*original)[o] : o;
  if( lst_map && lst_map->defined(obj)) lst = (*lst_map)[obj];
  return lst;
}

line_style GeoScene::set_obj_line_style(void* obj, line_style lst)
{
  line_style lold = get_obj_line_style(obj);
  if(!lst_map) lst_map = new map<void*,line_style>;
  (*lst_map)[obj] = lst;
  return lold;
}

int GeoScene::get_obj_line_width(void* o)
{ 
  int lw = get_default_line_width();
  void* obj = (original && original->defined(o)) ? (*original)[o] : o;
  if ( lw_map && lw_map->defined(obj) ) lw = (*lw_map)[obj];
  return  lw;
}

int GeoScene::set_obj_line_width(void* obj, int w)
{
  int lold = get_obj_line_width(obj);
  if(!lw_map) lw_map = new map<void*, int >;
  (*lw_map)[obj] = w;
  return lold;
}

bool GeoScene::get_obj_label(void* o,string& label)
{ 
  if (! label_map) return false;
  void* obj = (original && original->defined(o)) ? (*original)[o] : o;
  if ( label_map && label_map->defined(obj) )  { label=  (*label_map)[obj]; return true; }
  return false;
}

void GeoScene::set_obj_label(void* obj, string label)
{
  if(!label_map) label_map = new map<void*, string >;
  (*label_map)[obj] = label;
}  

bool GeoScene::get_obj_text(void* o, geowin_text& t)
{
  if (! text_map) return false; 
  void* obj = (original && original->defined(o)) ? (*original)[o] : o;
  if ( text_map && text_map->defined(obj) )  { t =  (*text_map)[obj]; return true; }
  return false;
}

void GeoScene::set_obj_text(void* obj, geowin_text t)  
{
  if(!text_map) text_map = new map<void*, geowin_text >;
  (*text_map)[obj] = t;
}

void GeoScene::set_original_properties(void* o, void* orig)
{
    if( color1_map && color1_map->defined(orig) )
      (*color1_map)[o] = (*color1_map)[orig];
    if( color2_map && color2_map->defined(orig) )
      (*color2_map)[o] = (*color2_map)[orig];
    if( lst_map && lst_map->defined(orig) )
      (*lst_map)[o] = (*lst_map)[orig];
    if( lw_map && lw_map->defined(orig) )
      (*lw_map)[o] = (*lw_map)[orig];
    if( label_map && label_map->defined(orig) )
      (*label_map)[o] = (*label_map)[orig];
    if( text_map && text_map->defined(orig) )
      (*text_map)[o] = (*text_map)[orig];      
}

void GeoScene::set_def_properties(void* o)
{
    if( color1_map && color1_map->defined(o) )
      (*color1_map)[o] = get_default_color1();
    if( color2_map && color2_map->defined(o) )
      (*color2_map)[o] = get_default_color2();
    if( lst_map && lst_map->defined(o) )
      (*lst_map)[o] = get_default_line_style();
    if( lw_map && lw_map->defined(o) )
      (*lw_map)[o] = get_default_line_width();
    if( label_map && label_map->defined(o) )
      (*label_map)[o] = string("");
    if( text_map && text_map->defined(o) )
      (*text_map)[o] = geowin_text();      
}

void GeoScene::original_properties(void* o)
{
    void* obj = original ? (*original)[o] : o;
    if(o==obj )  return; //was identical...

    set_original_properties(o,obj);
}
  
void GeoScene::set_default_properties(void* o)
{
    set_def_properties(o);
    if( original && original->defined(o))  (*original)[o] = o;
}

void GeoScene::undefine_colors1() { if(color1_map) delete color1_map; color1_map = 0;}

void GeoScene::undefine_colors2() { if(color2_map) delete color2_map; color2_map = 0; }

void GeoScene::undefine_line_styles()  { if(lst_map) delete lst_map; lst_map = 0; }

void GeoScene::undefine_line_width()  { if(lw_map) delete lw_map; lw_map = 0; }

void GeoScene::undefine_labels()  { 
   if(label_map) delete label_map; 
   label_map = 0; 
}

void GeoScene::undefine_texts()  { 
   if(text_map) delete text_map; 
   text_map = 0; 
} 


void GeoScene::undefine_all()
{
    undefine_colors1();
    undefine_colors2();
    undefine_line_styles();
    undefine_line_width();
    undefine_labels();
    undefine_texts();
}

int GeoScene::keyboard_panel(window& w, string& nvl)
{ 
    panel p;
    p.text_item("\\bf\\blue Input new object:");
    p.text_item("");
    p.string_item("New value:",nvl );
	
    p.fbutton("apply",   APPLY_BUTTON); 
    p.button("cancel",   CANCEL_BUTTON);
    
    int bt = p.open(w);
    return bt;
}

void GeoScene::setup_focus_dialog(window& w,void* obj)
{
   panel p("Object Properties");
   
   color c1      = get_obj_color(obj);
   color c2      = get_obj_fill_color(obj);
   int   lw     = get_obj_line_width(obj);
   line_style ls = get_obj_line_style(obj);

   string label,old_label;
   if (! get_obj_label(obj,label)) label="";
	
    p.text_item("\\bf\\blue Colors");
    p.color_item("interior", c2);   
    p.color_item("boundary", c1);   
    p.text_item("");
	
    p.text_item("\\bf\\blue Line Style");
    p.lstyle_item("line style", ls);
    p.lwidth_item("line width", lw);
    p.text_item("");  

    p.string_item("Object Label",label );
	
    p.fbutton("apply",   APPLY_BUTTON);
    p.button("default",  DEFAULTS_BUTTON);
    p.button("cancel",   CANCEL_BUTTON);

    int bt;
    while( (bt=p.open(w)) == DEFAULTS_BUTTON )
    {
	    c1 = get_default_color1();
	    c2 = get_default_color2();
	    lw = get_default_line_width();
	    ls = get_default_line_style();
    }

   if( bt == CANCEL_BUTTON ) return;
	
   if( c1 != get_obj_color(obj) ) set_obj_color(obj, c1);
   if( c2 != get_obj_fill_color(obj) ) set_obj_fill_color(obj, c2);
   if( lw != get_obj_line_width(obj)) set_obj_line_width(obj, lw);
   if( ls != get_obj_line_style(obj)) set_obj_line_style(obj, ls);

   if (! get_obj_label(obj,old_label)) old_label="";        
   if( label != old_label) { set_obj_label(obj, label); }
}
  
int GeoScene::focus_dialog(string& nvl, window& w)
{
    panel p;
    p.set_item_width(240);
    p.text_item("\\bf\\blue Object");
    p.text_item("");
    p.text_item(nvl); p.text_item(""); 
    p.string_item("new value",nvl );
	
    p.fbutton("apply",   APPLY_BUTTON); 
    p.button("cancel",   CANCEL_BUTTON);
    int bt=p.open(w);	
    return bt;
}

void GeoScene::set_window_params(GeoWin* gwin,window& w, void* adr, PresentationType pt)
{
    color c = get_obj_color(adr);
    if (pt != geowin_normal) c =  get_sel_color();

    color fc = get_obj_fill_color(adr);
    if (pt != geowin_normal) fc = sel_filling_color;

    oldcl   = w.set_color(c);
    oldfl   = w.set_fill_color(fc);
    
    line_style l = get_obj_line_style(adr);
    if (pt != geowin_normal) l = sel_l_style;
    
    int wd = get_obj_line_width(adr);
    if (pt != geowin_normal) wd = sel_line_width;
    
    old_ls  = w.set_line_style(l);
    oldw    = w.set_line_width(wd);  
    
    text_clr = get_default_text_color(); 
    old_ps = w.set_point_style(get_default_point_style());
    
    gwin->set_font(fixed_font, 10, "");  
    //setup_font(w,3);
}

void GeoScene::restore_window_params(window& w)
{
    w.set_color(oldcl);
    w.set_fill_color(oldfl);
    w.set_point_style(old_ps);
    w.set_line_style(old_ls);
    w.set_line_width(oldw);
}  
  
void GeoScene::set_ps_params(ps_file& F, void* adr)
{
    F.set_line_width((int)get_obj_line_width(adr));
    F.set_line_style(get_obj_line_style(adr));
    F.set_fill_color(get_obj_fill_color(adr));
    F.set_color(get_obj_color(adr));
}
 
void GeoScene::draw_object(void* adr, window& w, PresentationType pt)
{ 
  set_window_params(gw,w,adr,pt);

  if (draw_fcn) { // if draw_fcn is set, use it ...
   draw_fcn((void*)(&w), adr, draw_mode_number);
  }
  else window_output(w, adr); 

  string label;
  geowin_text gt;
  double xw1,xw2,yw1,yw2;
    
  if (get_obj_text(adr, gt)){
    get_bounding_box_fcn(adr,(void*)(&xw1),(void*)(&xw2),(void*)(&yw1),(void*)(&yw2));
    geowin_draw_text_object(w, xw1,xw2,yw1,yw2, adr, gt, true);      
  }
  
  // set font
  int wsz = w.real_to_pix(12);
  if (wsz < 7) wsz=7;
  string font_name = string("F%d",wsz);
  w.set_font(font_name); 
      
    
  if (get_obj_label(adr,label))
  { // get bounding box and draw label
    get_bounding_box_fcn(adr,(void*)(&xw1),(void*)(&xw2),(void*)(&yw1),(void*)(&yw2));     
    w.draw_ctext((xw1+xw2)/2,(yw1+yw2)/2,label, text_clr);
  }     
  
  restore_window_params(w); 
}      

void GeoScene::ps_draw_text(ps_file& F, void* adr)
{
 string label;
 geowin_text gt; 
 window& w = gw->get_window();
 double xw1, xw2, yw1, yw2;
 
 double wsz = (double) w.real_to_pix(12);  
 if (wsz < 7) wsz=7; 

 if (get_obj_text(adr, gt)){
    get_bounding_box_fcn(adr,(void*)(&xw1),(void*)(&xw2),(void*)(&yw1),(void*)(&yw2));
    geowin_draw_text_object(w, F, xw1,xw2,yw1,yw2, adr, gt, true);      
 }	
	
 if (get_obj_label(adr,label))
 { // get bounding box and draw label
    get_bounding_box_fcn(adr,(void*)(&xw1),(void*)(&xw2),(void*)(&yw1),(void*)(&yw2));
    double fold = F.set_font_size(wsz/45);
    F.set_text_font(string("Helvetica"));
    F.draw_ctext((xw1+xw2)/2,(yw1+yw2)/2,label, text_clr);
    F.set_font_size(fold);
 }      	
}


void GeoScene::move_point(void* obj_adr, double dx, double dy)
{ if (move_point_fcn) move_point_fcn(obj_adr,dx,dy,obj_pnr); } 


bool GeoScene::fill_rect(double& x1, double& x2, double& y1, double& y2)
{
  init_iteration();
  void* hadr;
   
  if (no_end_reached()) {
      hadr = get_next_obj_adr();
      get_bounding_box_fcn(hadr,(void*)(& x1),(void*)(& x2),(void*)(& y1),(void*)(& y2) );
      iterate();
  }  
  else return false;

  double x1akt,x2akt,y1akt,y2akt;

  while(no_end_reached()) {
       hadr = get_next_obj_adr();
       get_bounding_box_fcn(hadr,(void*)(& x1akt),(void*)(& x2akt),(void*)(& y1akt),(void*)(& y2akt) );
       iterate();
       if (x1akt<x1) x1=x1akt;
       if (x2akt>x2) x2=x2akt;
       if (y1akt<y1) y1=y1akt;
       if (y2akt>y2) y2=y2akt;
  }
  return true;
}


void GeoScene::set_select(void* adr, bool b){}
  
void GeoScene::set_select(const list<void*>& L, bool b)
{ void* adr;
  forall(adr, L) set_select(adr, b); 
}


bool GeoScene::get_select(void* obj_adr){ return false; }

void GeoScene::select_in_rect(bool uflag, double x1, double y1, double x2, double y2,bool b)
// uflag false : select in rect, true: unselect
{
  bool sel = false;
  list<void*> tmp;
  void* act;
    
  if (! uflag){
   init_iteration(); 
    
   while(no_end_reached() ) 
   {
       act = get_next_obj_adr();
       if(b || box_intersection_fcn(act, x1, y1, x2, y2,true) ) tmp.push(act); 
       iterate();
   }
   forall(act, tmp)    if( !get_select(act))   { sel = true; break; }
  }
  else  tmp = sel_objs; 
    
  set_select(tmp, sel);
  if( !tmp.empty() && gw ) gw->redraw();
}
  

void GeoScene::compute_d3_output(GeoScene* sc, d3_window& W, GRAPH<d3_point,int>& H)
{ 
  if (init_d3_window!=NULL) init_d3_window(sc,W,H); 
  else if (objects_init_d3_window!=NULL) {
    void* objs = get_untyped_objptr();
    objects_init_d3_window(objs,(void*)(&W),(void*)(&H)); 
  }   
}

void GeoScene::obj_edit(){
  void* mouse_obj_ptr = get_untyped_mouseobjptr();
  if (mouse_obj_exist && edit_obj_fcn) edit_obj_fcn((void*)gw, mouse_obj_ptr,edit_mode);
}

bool GeoScene::start_changing()
{
  original = new map<void*, void*>;
  void* act_adr;
    
  void* m_adr = get_untyped_mouseobjptr();
    
  if( mouse_obj_exist  && !get_select(m_adr) )
      {
	if( !start_change_fcn || start_change_fcn((void*)(&gw), m_adr) )
	  mouse_changing = true;
	else return false;
	(*original)[m_adr] = m_adr;
      }
  else
      {
	mouse_changing = false;
	
	forall(act_adr, sel_objs) 
	  {
	    if( start_change_fcn && !start_change_fcn((void*)(&gw), act_adr)) return false;
	    (*original)[act_adr] = act_adr;
	  }
      }
  while_dragging = true;
  return true;
}

void GeoScene::stop_changing()
{
  void* act_adr;
  void* m_adr = get_untyped_mouseobjptr();

  if( mouse_changing )
      {
	if ( end_change_fcn ) end_change_fcn((void*)(&gw), m_adr);
	mouse_changing = false;
	original_properties(m_adr); 
      }
  else
      forall(act_adr, sel_objs)
      {
	if ( end_change_fcn ) end_change_fcn((void*)(&gw), act_adr);
	original_properties(act_adr);
      }
    
  while_dragging = false;
  delete original;
  original = 0;
  update_and_redraw();
}

bool GeoScene::del_mouse_object(bool b)
{ return false; }

void GeoScene::add_mouse_object_copy(void*& new_adr)
{ new_adr = NULL; }

void GeoScene::del_focus(geo_editor ed)
{  
  void* m_adr = get_untyped_mouseobjptr();

  if (mouse_obj_exist)
  { if (get_select(m_adr))  ed->del_sel(); 
    else {
        if (del_mouse_object(true))
        { 
	  change_mode=1;
	  update_and_redraw();
         }
    }
  }
}

void GeoScene::setup_focus_or_raise(bool flag)
{
  if( mouse_obj_exist )
      {
	void* m_adr = get_untyped_mouseobjptr();
	
	if (flag) { // setup focus ...
	 window& w = gw->get_window();
         setup_focus_dialog(w, m_adr);
	}
	else {	
         void* new_adr = NULL;

         add_mouse_object_copy(new_adr);
	 
	 if (new_adr != NULL){ 
          set_original_properties(new_adr ,m_adr);
          del_mouse_object(false);
         } 
         update();
	}
	
	gw->redraw();
      }
}

void GeoScene::toggle_selection() 
{
    if( mouse_obj_exist )
      {
        void* m_adr = get_untyped_mouseobjptr();      
	set_select(m_adr, !get_select(m_adr));
	if( gw ) draw_object_with_properties(m_adr);
      }
}

void GeoScene::draw_object_with_properties(void* obj_adr)  
{
    window& w = gw->get_window();
    bool draw_flag=false;
    
    PresentationType pt = geowin_normal;
    if( get_select(obj_adr) )  pt = geowin_selected;   
    
    if( mouse_obj_exist ) {  
      void* mptr = get_untyped_mouseobjptr();
    
      if ( mptr == obj_adr ) 
      {  if (high_light) pt = geowin_focus; 
         if (handle_defining_points==geowin_highlight) draw_flag=true; 
      } 
    }

    draw_object(obj_adr, w, pt);
    list<point> hcont; 
    
    if (defining_points_fcn && (handle_defining_points==geowin_show || draw_flag)) 
       defining_points_fcn(obj_adr,(void*)(&hcont));
    show_points(hcont);    
} 
 
 
 
bool GeoScene::find_object(double x, double y)
{
    if ( !gw ) return false;
    window& w = gw->get_window();
    double d  = w.pix_to_real(1)* gw->get_mfct();
    double x1=x-d, x2=x+d, y1=y-d, y2=y+d;
    
    init_iteration();
    void* hadr;
    
    list<point> hcont;  
    
    void* old_mouse_objptr = get_untyped_mouseobjptr();
    bool old_mouse_obj_exist = mouse_obj_exist;

    mouse_obj_exist = false;
    obj_defining_points_change = false;
    
    if (defining_points_fcn && (handle_defining_points!=geowin_hide)){
    
      while (no_end_reached()){
        hadr = get_next_obj_adr();
      
        defining_points_fcn(hadr,(void*)(&hcont));
        point piter;
	int cnt=0;
        forall(piter,hcont){
	  double xw = piter.xcoord(), yw = piter.ycoord();
	  if ((xw>=x1) && (xw<=x2) && (yw>=y1) && (yw<=y2)) {
	     obj_defining_points_change=true;
	     obj_pnr= cnt;
	  }
	  cnt++;
	}
        hcont.clear();	
	iterate();
      }
      
    }

    init_iteration();
    
    while(no_end_reached())
      {
        hadr = get_next_obj_adr();
	if( box_intersection_fcn(hadr, x1, y1, x2, y2,true) )
	  {
	    mouse_obj_exist = true;
	    //set the mouse object ...
	    set_mouse_object(hadr);
	  }
	iterate();
      }

    if ((high_light || (defining_points_fcn && (handle_defining_points == geowin_highlight))) && 
        (old_mouse_objptr != get_untyped_mouseobjptr() || old_mouse_obj_exist != mouse_obj_exist))
    { gw->redraw();
      if (mouse_obj_exist) draw_object_with_properties(get_untyped_mouseobjptr()); 
     }
       
    return mouse_obj_exist;
} 
 
 
GEOWIN_END_NAMESPACE







