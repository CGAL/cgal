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
// file          : src/GeoWin/geo_menu.c
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

// ****************** menu handling ************************************

void GeoWin::call_handler(int n) 
{
   GeoFunctions* mf = menu_fcn[n];
   mf->call(*this); 
}

void geo_call_handler(int n) 
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->call_handler(n);
}

void GeoWin::activate(int id)
{
  geo_scene sc;
  forall(sc, scenes)
    if( sc->ids[geo_activ_id] == id ) { activate(sc); return; }
}

void geo_activate(int n)
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->activate(n);
}

void GeoWin::scene_options(int id)
{
  geo_scene sc;
  forall(sc, scenes)
    if( sc->ids[geo_opt_id] == id )  { sc->options();  return; }
}

void GeoWin::scene_description(int id)
{
  geo_scene sc;
  forall(sc, scenes)
    if( sc->ids[geo_desc_id] == id )  { sc->scene_description();  return; }
}


void GeoWin::scene_group_menu(int id)
{
  //cout << id << "\n";
  geo_scenegroup gsg;
  forall(gsg, get_scenegroups()){
    if (gsg->id == id) {
       //cout << "Found\n";
       scenegroup_dialog(gsg);
       return;
    }
  }
}

void GeoWin::scene_contents(int id)
{
  geo_scene sc;
  forall(sc, scenes)
    if( sc->ids[geo_cont_id] == id )  { sc->contents();  return; }
}

void GeoWin::active_scene_options()
{
  if(cur) cur->options();
}

void geo_scene_option(int n)
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->scene_options(n);
}

void geo_scene_description(int n)
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->scene_description(n); 
}

void geo_scene_algo(int n)
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  geo_scene cur = gw->get_active_scene();
  
  if (cur) cur->call_buffer_fcn(n,true);
}

void geo_scene_groups(int n)
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->scene_group_menu(n);
}

void geo_scene_contents(int n)
{
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->scene_contents(n);
}


void GeoWin::new_scene(int n, double dummy)
{
  dummy=dummy+1;
  string type_name = editables[editables[n]];
  //cout << type_name << "\n";
  geo_scene sc = GeoScene::new_scene(type_name);
  insert_scene(sc);
  sc->init_default_values();
  if( !sc->active) activate(sc);
}

void geo_new_scene_handler(int n)
{ 
  window*    wp = window::get_call_window();
  GeoWin* gw = (GeoWin*)wp->get_inf();
  gw->new_scene(n, 0.5); 
}

void geo_call(GeoWin& gw, void (*fkt)(GeoWin&), string)
{
  fkt(gw);
}

void geo_call(GeoWin&, void (*fkt)(), string)
{
  fkt();
}

int GeoWin::add_call(GeoFunctions* f, string label, window* m)
{
  if (f)
    {
      menu_fcn[menu_fcn_count] = f;
      if (m) m->button(label, menu_fcn_count++, geo_call_handler);
      else  
	{
	  menus[run_menu]->button(label, menu_fcn_count++, geo_call_handler);
	  user_fcn_count++;
	}
      return menu_fcn_count-1;
    }
  return -1;
}

int GeoWin::add_button_call(GeoFunctions* f, string label,window* m, char* pmap)
{
  if (f)
    {
      menu_fcn[menu_fcn_count] = f;
      m->button(pmap,pmap,label, menu_fcn_count++, geo_call_handler);

      return menu_fcn_count-1;
    }
  return -1;
}



void GeoWin::init_menu()
{
  if (menu_init) return;
  menu_init = true;

  Wp->buttons_per_line(8);  
  int but,but2, but3,but4, but5;
  
  // File menu
  add_call( new GeowinMember(&GeoWin::file_open), "&Load Scene",
	    menus[io_menu]);
  but = add_call( new GeowinMember(&GeoWin::file_save_active), "&Save active Scene",menus[io_menu]);
  
  add_call( new GeowinMember(&GeoWin::file_save), "Save Scene",menus[io_menu]);
 
  but5= add_call( new GeowinMember(&GeoWin::file_import), "Import Data",menus[io_menu]);
 
  but2= add_call( new GeowinMember(&GeoWin::file_export), "E&xport Data",menus[io_menu]);
  menus[io_menu]->new_line();

  but3= add_call( new GeowinMember(&GeoWin::save_ps), "PS &active scene",menus[io_menu]);
  but4= add_call( new GeowinMember(&GeoWin::save_ps_vis), "PS &visible scenes",menus[io_menu]);

  add_call( new GeowinMember(&GeoWin::save_screen), "&Screenshot",menus[io_menu]);

  menus[io_menu]->new_line();
   
  buttons[GEO_CM_IMPORT]=but5;
  button_info[GEO_CM_IMPORT]=1;

  buttons[GEO_CM_EXPORT]=but2;
  button_info[GEO_CM_EXPORT]=1;

  buttons[GEO_CM_WRITE]     = but;
  button_info[GEO_CM_WRITE] = 1;

  buttons[GEO_CM_POST]      = but3;
  button_info[GEO_CM_POST]  = 1;

  buttons[GEO_CM_POSTVIS]   = but4;
  button_info[GEO_CM_POSTVIS]  = 1;

  but = add_call( new GeowinMember(&GeoWin::quit_button_pressed), 
		  "&Quit", menus[io_menu]);  
  buttons[GEO_CM_QUIT]     = but;
  button_info[GEO_CM_QUIT] = 1;

  // Edit menu1
  add_call(new GeoSceneMember(&GeoScene::edit), "&Edit",  menus[edit_menu1]);
  
  // Edit menu2
  add_call( new GeowinMember(&GeoWin::select_menu),   "Select", menus[mark_menu]);
  // change
  add_call( new GeowinMember(&GeoWin::select_all),   "Select all", menus[mark_menu]); 
 
  add_call( new GeowinMember(&GeoWin::unselect), "Unselect",  menus[mark_menu]);

  menus[edit_menu2]->button("Selection", menu_const + mark_menu, 
			    *(menus[mark_menu]));
 
  // operate submenu
  add_call( new GeoEditorMember(&GeoEditor::copy), "Copy", 
	    menus[operate_menu]);
  add_call( new GeoEditorMember(&GeoEditor::paste), "Paste",
	    menus[operate_menu]);
  add_call( new GeowinMember(&GeoWin::paste_from), "Paste from",
            menus[operate_menu]);
  add_call( new GeowinMember(&GeoWin::move_menu), "Move", menus[operate_menu]);
  add_call( new GeowinMember(&GeoWin::rotate_menu), "Rotate",  
	    menus[operate_menu]); 
  add_call( new GeowinMember(&GeoWin::show_sel), "Show selected",  
	    menus[operate_menu]);     
  add_call( new GeoEditorMember(&GeoEditor::del_sel), "Delete", 
	    menus[operate_menu]);
	
  menus[edit_menu2]->button("Edit Selection", menu_const + edit_menu2,
			    *(menus[operate_menu]));

  add_call(new GeowinMember(&GeoWin::call_mouse_read_object), 
	   "Mouse Read Object",  menus[readobject_menu]);
  add_call(new GeowinMember(&GeoWin::call_keyboard_read_object), 
	   "Keyboard Read Object",  menus[readobject_menu]);
	   
  menus[edit_menu2]->button("Read Object", menu_const + edit_menu2,
			    *(menus[readobject_menu]));
  menus[edit_menu2]->new_line();  
  add_call(new GeoEditorMember(&GeoEditor::generate), "Generate", 
	   menus[edit_menu2]);

  void (GeoWin::*za)() = &GeoWin::zoom_area;

  // Window menu
  add_call( new GeowinMember(&GeoWin::redraw_menu) , "Redraw",
	    menus[view_menu]);
  add_call(new GeowinMember(&GeoWin::zoom_up) ,   "Zoom &Up",
	   menus[view_menu]);
  add_call(new GeowinMember(&GeoWin::zoom_down) , "Zoom &Down", 
	   menus[view_menu]); 
  add_call(new GeowinMember(za) , "Zoom Area",
	   menus[view_menu]);
  add_call(new GeowinMember(&GeoWin::fill_window),"Fill &Window", 
	   menus[view_menu]);
  add_call(new GeowinMember(&GeoWin::reset_window),"Reset", 
	   menus[view_menu]);
  menus[view_menu]->new_line();
  add_call(new GeowinMember(&GeoWin::show_d3),"Show 3d Output", 
	   menus[view_menu]);  
  
  // Options menu
  add_call(new GeowinMember(&GeoWin::global_options), "General", 
	   menus[option_menu]);
  add_call(new GeowinMember(&GeoWin::advanced_global_options), "Advanced", 
	   menus[option_menu]);	   
  add_call(new GeowinMember(&GeoWin::d3_options), "3d Output", 
	   menus[option_menu]);
  add_call(new GeowinMember(&GeoWin::scaling_options), "Scaling", 
	   menus[option_menu]);
  but = add_call(new GeowinMember(&GeoWin::active_scene_options), 
		 "Active Scene", menus[option_menu]);
  buttons[GEO_CM_ACTIVE_OPT] = but;
  button_info[GEO_CM_ACTIVE_OPT] = 1;
  
  menus[option_menu]->new_line();
  add_call(new GeowinMember(&GeoWin::scene_defaults), "Scene defaults", 
	   menus[option_menu]);  
  
  // Scene menu  
  if( editables.empty() )
    {
      buttons[GEO_CM_NEW] = 
	add_call( new GeowinMember(&GeoWin::redraw_menu), "open",
		  menus[scene_menu]);
      menus[scene_menu]->disable_button(buttons[GEO_CM_NEW]);
      button_info[GEO_CM_NEW] = 2;
    }
  else
    {
      string str;
      forall(str, editables)
	menus[scene_type_menu]->button(str, geo_new_scene_handler);
      
      buttons[GEO_CM_NEW] = 
	menus[scene_menu]->button("open", menu_const + scene_type_menu,
				  *(menus[scene_type_menu]));
      button_info[GEO_CM_NEW] = 1;
    }
  
  buttons[GEO_CM_CLOSE] = 
    add_call( new GeowinMember(&GeoWin::destroy_menu), "close",
	      menus[scene_menu]);
  button_info[GEO_CM_CLOSE] = 1;

  buttons[GEO_CM_CLOSE_ALL] = 
    add_call( new GeowinMember(&GeoWin::destroy_all), "close all",
	      menus[scene_menu]);
  button_info[GEO_CM_CLOSE_ALL] = 1;


  buttons[GEO_CM_CLEAR] =
  add_call( new GeoEditorMember(&GeoEditor::clear), "clear",
	    menus[scene_menu]);
  button_info[GEO_CM_CLEAR] = 1;

  buttons[GEO_CM_VISIBLE]=
    add_call( new GeowinMember(&GeoWin::visible_scenes) , "visibility", 
	    menus[scene_menu]);
  button_info[GEO_CM_VISIBLE]=1;
  
  buttons[GEO_CM_Z]=
    add_call( new GeowinMember(&GeoWin::z_scenes) , "z-order", 
	    menus[scene_menu]);
  button_info[GEO_CM_Z]=1;
  
  menus[scene_menu]->button("desription", menu_const + scene_desc_menu,
			     *(menus[scene_desc_menu]));  

  menus[scene_menu]->button("options", menu_const + scene_opt_menu,
			     *(menus[scene_opt_menu]));

  menus[scene_menu]->button("activate", menu_const + scene_list_menu,
			    *(menus[scene_list_menu])); 

  menus[scene_menu]->button("groups", menu_const + scene_groups_menu,
			    *(menus[scene_groups_menu])); 
			    
  menus[scene_menu]->new_line();
  
  menus[scene_menu]->button("contents", menu_const + scene_cont_menu,
			     *(menus[scene_cont_menu]));  

  //Help menu
  add_call(new GeowinMember(&GeoWin::help_about), "About", 
	   menus[help_menu]);
  add_call(new GeowinMember(&GeoWin::help_news), "News", 
	   menus[help_menu]);  
  // help
  add_help_text("geowin_buttons");
  add_help_text("geowin_file");
  add_help_text("geowin_edit");
  add_help_text("geowin_scenes");
  add_help_text("geowin_scene_options");
  add_help_text("geowin_window");
  add_help_text("geowin_options1");
  add_help_text("geowin_options2");
	     
  if (!help_list.empty())
  { menus[help_menu]->new_line(); 
    string hlp;
    forall(hlp,help_list) 
        add_call(new GeowinMember(&GeoWin::help_user), hlp, menus[help_menu]); 
   }
   
  //Algorithm menu
  menus[algo_menu]->button("Algorithm", menu_const + algo_algorithm_menu,
			     *(menus[algo_algorithm_menu]));

  menus[algo_menu]->button("Algorithm options", menu_const + algo_options_menu,
			    *(menus[algo_options_menu])); 
  
  // Main Menu
  Wp->button("File", menu_const + io_menu, *(menus[io_menu]));
  Wp->button("Edit", menu_const + edit_menu1,*(menus[edit_menu1]));
  Wp->button("Scenes", menu_const + scene_menu, *(menus[scene_menu]));
  Wp->button("Window", menu_const + view_menu, *(menus[view_menu])); 
  if (user_fcn_count > 0) 
    Wp->button("User", menu_const + run_menu, *(menus[run_menu]));
  but = add_call( new GeowinMember(&GeoWin::done_button_pressed), 
		  "Done", Wp ); 
  Wp->button("Options", menu_const + option_menu, *(menus[option_menu]));
  Wp->button("Algorithms", menu_const + algo_menu, *(menus[algo_menu]));
  Wp->button("Help", menu_const + help_menu, *(menus[help_menu]));
  buttons[GEO_CM_OK]     = but;

  // set default state
  disable_done_button();
  
  add_call(new GeowinMember(&GeoWin::msg_clear) , "Clear", msg_win);
  add_call(new GeowinMember(&GeoWin::msg_close) , "Close", msg_win);
  
  Wp->make_menu_bar(2);
}

void   GeoWin::make_scene_menu()
{
  if( is_open )
    {
      menu* m1 = menus[scene_list_menu] = new menu;
      m1->set_inf(this);
      
      menu* m6 = menus[scene_desc_menu] = new menu;
      m6->set_inf(this);
      
      menu* m2 = menus[scene_opt_menu] = new menu;
      m2->set_inf(this);
      menu* m3 = menus[scene_cont_menu] = new menu;
      m3->set_inf(this);     
      menu* m4 = menus[scene_groups_menu] = new menu;
      m4->set_inf(this);
      menu* m5 = menus[algo_algorithm_menu] = new menu;
      m5->set_inf(this);
      
      int count = 0;
      
      if (cur){
        int dummy; 
        list<string> LS = cur->get_fcn_names(dummy, false);
	//cout << LS << "\n";
	string iter;
	forall(iter,LS) m5->button(iter,geo_scene_algo);
      }
      
      if( cur )
	{
	  geo_scene sc;
	  forall(sc, scenes)
	    {
	      if( (sc->owner).empty() )
		{
		  sc->ids[geo_activ_id] =
		    m1->button(sc->get_name(),geo_activate);
		  count++;
		}
	      sc->ids[geo_opt_id] = 
		m2->button(sc->get_name(),geo_scene_option);
	      sc->ids[geo_cont_id] = 
		m3->button(sc->get_name(),geo_scene_contents);	
	      sc->ids[geo_desc_id] = 
		m6->button(sc->get_name(),geo_scene_description);			
		
	    }
	  // scene groups ...
	  geo_scenegroup gsg;
	  
	  forall(gsg,scene_groups)
	  {
	     if (gsg->id == -1) {
	       gsg->id = m4->button(gsg->name,geo_scene_groups);
	       //cout << gsg->name << "\n";
	     }
	  }  
	    
	}
	
      if( count )
	{
	  delete menus[scene_menu]->set_window( menu_const+scene_list_menu,m1);
	  menus[scene_menu]->enable_button(menu_const + scene_list_menu);
	}
      else menus[scene_menu]->disable_button(menu_const + scene_list_menu);
      delete menus[scene_menu]->set_window(menu_const + scene_opt_menu,  m2);
      delete menus[scene_menu]->set_window(menu_const + scene_desc_menu,  m6);
      delete menus[scene_menu]->set_window(menu_const + scene_cont_menu,  m3);
      delete menus[scene_menu]->set_window(menu_const + scene_groups_menu, m4);
      delete menus[algo_menu]->set_window(menu_const + algo_algorithm_menu, m5);
    }
}

void GeoWin::make_algo_menu(geo_scene sc)
{
  if (is_open){
      menu* m5 = menus[algo_algorithm_menu] = new menu;
      m5->set_inf(this);
      
      int dummy;
      list<string> LS = sc->get_fcn_names(dummy,false);
      string iter;
      forall(iter,LS) m5->button(iter,geo_scene_algo);
      delete menus[algo_menu]->set_window(menu_const + algo_algorithm_menu, m5);
  }
}

void GeoWin::enable_done_button()
{
  Wp->enable_button(buttons[GEO_CM_OK]);
  button_info[GEO_CM_OK] = 1;
}

void GeoWin::disable_done_button()
{
  Wp->disable_button(buttons[GEO_CM_OK]);
  button_info[GEO_CM_OK] = 0;
}

void GeoWin::enable_close_button()
{
  menus[scene_menu]->enable_button(buttons[GEO_CM_CLOSE]);
  menus[scene_menu]->enable_button(buttons[GEO_CM_CLOSE_ALL]);
  button_info[GEO_CM_CLOSE] = 1;
  button_info[GEO_CM_CLOSE_ALL] = 1;  
}

void GeoWin::disable_close_button()
{
  menus[scene_menu]->disable_button(buttons[GEO_CM_CLOSE]);
  menus[scene_menu]->disable_button(buttons[GEO_CM_CLOSE_ALL]);  
  button_info[GEO_CM_CLOSE] = 0;
  button_info[GEO_CM_CLOSE_ALL] = 0;
}

void GeoWin::enable_new_button()
{
  if( button_info[GEO_CM_NEW] < 2 )
    {
      menus[scene_menu]->enable_button(buttons[GEO_CM_NEW]);
      button_info[GEO_CM_NEW] = 1;
    }
}

void GeoWin::disable_new_button()
{
 if( button_info[GEO_CM_NEW] < 2 )
   {
     menus[scene_menu]->disable_button(buttons[GEO_CM_NEW]);
     button_info[GEO_CM_NEW] = 0;
   }
}

void GeoWin::disable_menus()
{
  Wp->disable_button(menu_const +  io_menu);          
  Wp->disable_button(menu_const +  edit_menu1);    
  Wp->disable_button(menu_const +  edit_menu2);        
  Wp->disable_button(menu_const +  scene_menu);     
  Wp->disable_button(menu_const +  view_menu);       
  Wp->disable_button(menu_const +  zoom_menu);       
  Wp->disable_button(menu_const +  option_menu);    
  Wp->disable_button(menu_const +  help_menu);      
}

void GeoWin::enable_menus()
{
  Wp->enable_button(menu_const +  io_menu);          
  Wp->enable_button(menu_const +  edit_menu1); 
  Wp->disable_button(menu_const +  edit_menu2);             
  Wp->enable_button(menu_const +  scene_menu);     
  Wp->enable_button(menu_const +  view_menu);       
  Wp->enable_button(menu_const +  zoom_menu);       
  Wp->enable_button(menu_const +  option_menu);    
  Wp->enable_button(menu_const +  help_menu);        
}

void GeoWin::no_scene_on()
{
  if( !is_open ) return;
  menus[io_menu]->disable_button(buttons[GEO_CM_WRITE]);
  menus[io_menu]->disable_button(buttons[GEO_CM_EXPORT]);
  menus[io_menu]->disable_button(buttons[GEO_CM_IMPORT]);
  menus[io_menu]->disable_button(buttons[GEO_CM_POST]);

  menus[scene_menu]->disable_button(menu_const + scene_list_menu);
  menus[scene_menu]->disable_button(buttons[GEO_CM_CLEAR]);
  menus[scene_menu]->disable_button(buttons[GEO_CM_Z]);  
  menus[scene_menu]->disable_button(buttons[GEO_CM_VISIBLE]);
  menus[scene_menu]->disable_button(menu_const + scene_desc_menu);  
  menus[scene_menu]->disable_button(menu_const + scene_opt_menu);
  menus[scene_menu]->disable_button(menu_const + scene_cont_menu); 
  menus[scene_menu]->disable_button(menu_const + scene_groups_menu);

  menus[option_menu]->disable_button(buttons[GEO_CM_ACTIVE_OPT]);
  menus[scene_menu]->disable_button(buttons[GEO_CM_CLOSE]);
  menus[scene_menu]->disable_button(buttons[GEO_CM_CLOSE_ALL]);
  Wp->disable_button(menu_const + edit_menu1);
}

void GeoWin::no_scene_off()
{
  if( !is_open ) return;
  menus[io_menu]->enable_button(buttons[GEO_CM_WRITE]);
  menus[io_menu]->enable_button(buttons[GEO_CM_EXPORT]);
  menus[io_menu]->enable_button(buttons[GEO_CM_IMPORT]);
  menus[io_menu]->enable_button(buttons[GEO_CM_POST]);

  menus[scene_menu]->enable_button(menu_const + scene_list_menu);
  menus[scene_menu]->enable_button(buttons[GEO_CM_CLEAR]);
  menus[scene_menu]->enable_button(buttons[GEO_CM_Z]);
  menus[scene_menu]->enable_button(buttons[GEO_CM_VISIBLE]);
  menus[scene_menu]->enable_button(menu_const + scene_desc_menu);  
  menus[scene_menu]->enable_button(menu_const + scene_opt_menu);
  menus[scene_menu]->enable_button(menu_const + scene_cont_menu);  
  menus[scene_menu]->enable_button(menu_const + scene_groups_menu);  

  menus[option_menu]->enable_button(buttons[GEO_CM_ACTIVE_OPT]);
  menus[scene_menu]->enable_button(buttons[GEO_CM_CLOSE]);  
  menus[scene_menu]->enable_button(buttons[GEO_CM_CLOSE_ALL]);  
  if( cur->edit_menu_type != -1)
    Wp->enable_button(menu_const + edit_menu1);
} 

// ************************** mouse handling *****************************

void GeoWin::clear_actions()
{ 
  mouse_actions.clear();
} 

void GeoWin::set_action(long mask, geo_action f)
{
  mouse_actions[mask] = f;
}

geo_action GeoWin::get_action(long mask)
{
  return mouse_actions[mask];
}

// ********************* Default Actions ***************************

void geo_mouse_input(GeoWin& gw, const point& p)
{
  window& W = gw.get_window();

  geo_scene  sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  
  if (ed)
    {
      //restore_events();
      W.set_point_buffer(p);
      ed->scene_read_object(0);
    }
}

void geo_scroll_scene(GeoWin& gw, const point& p)
{ window& W = gw.get_window();
  if (gw.show_grid) { W.set_grid_dist(-W.get_grid_dist()); }
  W.scroll_window(p,GeoWin::redraw_geowin);
  if (gw.show_grid) { W.set_grid_dist(-W.get_grid_dist()); }
}



void geo_select(GeoWin& gw, const point& p)
{
  geo_scene  sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  
  if (ed) gw.select(p);
}

void geo_toggle_selection(GeoWin& gw , const point&)
{
  geo_scene  sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  
  if (ed) sc->toggle_selection();
}

void geo_move(GeoWin& gw, const point& p)          
{
  geo_scene  sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  if (ed) gw.move(p);
}

void geo_rotate(GeoWin& gw, const point& p)        
{
  geo_scene  sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  
  if (ed) gw.rotate(p);
}

void geo_object_dragging(GeoWin& gw, const point& p)
{
  geo_scene sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  if( !ed ) return;
 
  point pp;
  if( gw.get_pin_point(pp) ) gw.rotate(pp, p);
  else                       gw.move(p);
}

void geo_local_object_menu(GeoWin& gw, const point& p)
{
  window& w = gw.get_window();
  geo_scene sc = gw.get_active_scene();
  geo_editor ed = sc ? sc->type_editor() : 0;
  if( !ed ) return;
  
  menu  M;
  M.button("setup",  1);
  M.button("(de)select", 2);
  M.button("delete", 3);
  M.button("object", 4);
  M.button("raise",  5);
  M.button("edit",   6);
  
  int x = w.xpix(p.xcoord()) - 1;
  int y = w.ypix(p.ycoord()) + 2;
  
  int ymax = w.ypix(gw.get_ymin());
  int xmax = w.xpix(gw.get_xmax());
  
  if (y > ymax-100) y=y-100;
  if (x > xmax-65) x=x-65;

  int but = M.open(w, x, y);
  
  switch (but) 
    {
    case 1 :  sc->setup_focus_or_raise(true); break;  // setup focus
    case 2 :  sc->toggle_selection(); break;
    case 3 :  sc->del_focus(ed); break;
    case 4 :  ed->obj_focus(); break;
    case 5 :  sc->setup_focus_or_raise(false); break; // raise
    case 6 :  sc->obj_edit(); break;
    default :  break;
    }
  
}

void geo_active_scene_options(GeoWin& gw, const point&)
{
  geo_scene sc = gw.get_active_scene();
  if(sc) sc->scene_options();
}

void geo_active_scene_input_options(GeoWin& gw, const point&)
{
  geo_scene sc = gw.get_active_scene();
  if(sc) sc->scene_input_options();
}


void geo_pin_point(GeoWin& gw, const point& p)
{
  point oldp;
  
  if( gw.get_pin_point(oldp) )
    {
      double d =  gw.get_window().pix_to_real(4);
      double xo = oldp.xcoord();
      double yo = oldp.ycoord();
      double x =  p.xcoord();
      double y =  p.ycoord();
      
      if( x-d <= xo && x+d >= xo && y-d <= yo && y+d >= yo)
	{
	  gw.del_pin_point();
	  gw.redraw();
	  return;
	}
    }
  gw.set_pin_point(p);
  gw.redraw();
}


void GeoWin::reset_actions()
{
  mouse_actions.clear();
  
  set_action( A_LEFT , geo_mouse_input);
  set_action( A_LEFT | A_OBJECT, geo_mouse_input);
  set_action( A_LEFT | A_DRAG, geo_scroll_scene);
  set_action( A_LEFT | A_DRAG | A_OBJECT, geo_object_dragging);
  
  set_action( A_LEFT | A_CTRL, geo_pin_point);
  set_action( A_LEFT | A_CTRL | A_OBJECT, geo_pin_point);
  

  set_action( A_MIDDLE |  A_OBJECT, geo_toggle_selection);
  set_action( A_MIDDLE |  A_DRAG, geo_select); 
  set_action( A_MIDDLE |  A_DRAG | A_OBJECT, geo_select); 
  
  set_action( A_RIGHT | A_CTRL, geo_active_scene_input_options);
  set_action( A_RIGHT | A_IMMEDIATE, geo_active_scene_options);
  set_action( A_RIGHT | A_IMMEDIATE | A_OBJECT, geo_local_object_menu);
  
}

GEOWIN_END_NAMESPACE
