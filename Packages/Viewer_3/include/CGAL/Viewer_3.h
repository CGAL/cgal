// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Viewer_3.h
// revision      : $Revision$
//
// author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#include <CGAL/GL_win.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
extern "C" {
#include <stdio.h>
}
#include <FL/x.H> // platform-indep. fttk wrapper for X11/Windows
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Radio_Light_Button.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Dial.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Roller.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Multi_Browser.H>
#include <FL/Fl_Repeat_Button.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Counter.H>

CGAL_BEGIN_NAMESPACE
class Viewer_3;


void custom_win(GL_win* , Viewer_3*);



void* mainloop(Viewer_3 * );

const int size=600;

class Viewer_3
{
public:
typedef void (*User_ctr_win)(GL_win *,Viewer_3 *);

private:
typedef Point_3<Cartesian<double> >   Point3;
  User_ctr_win default_impl;
  int *scale;


#ifdef USE_THREAD
  pthread_t thr1;
#endif
  int group ;

  int clip_plane;
  Color obj_color1;
  Color obj_color2;
  Size  obj_size;
  Style obj_style;
  Precision obj_precision;
  
  bool boucle ;

  Fl_Window        *form;
  GL_win           *canvas ;
  Fl_Button        *close_but;
  Fl_Button        *exit_but;
  Fl_Button        *reset_but;
  Fl_Dial          *angle_sld;
  Fl_Dial          *deep_sld;
  Fl_Round_Button  *ortho_but;
  Fl_Round_Button  *persp_but;
  Fl_Choice        *mode_but;
  Fl_Roller        *mvplan_sld;
  Fl_Multi_Browser *scene;
  Fl_Menu_Button   *insert_but;
  Fl_Repeat_Button *mvobjup_sld; 
  Fl_Repeat_Button *mvobjdn_sld; 
  Fl_Choice        *group_but;
  Fl_Menu_Button   *group_men;
  Fl_Button        *light_but;
  Fl_Button        *ps_but;
  Fl_Button        *user_but;
  Fl_Window        *light_win;
  Fl_Window        *para_win;
  Fl_Button        *para_but;
  Fl_Window        *clip_win;
  Fl_Button        *clip_but;
  Fl_Menu_Button   *style_but;



  // Private functions ------------------


void perspective_view();

void orthogonal_view();

  //removes line in scene graph
void remove_line(int);

  // inserts a new line for drawable object 
void insert_in_group(int, Drawable_object*);

void rebuild_graph();

void update_color();

  // PUBLIC PART
public:

void reset();

  // sets the mouse events for the user mode
void set_mouse_push_handler(GL_win::Mouse_click);

void set_mouse_grab_handler(GL_win::Mouse_grab);

void set_mouse_release_handler(GL_win::Mouse_click);

  // for all groups, sets all local transformations to identity
void reset_groups();

int get_line_number(int, int);

int get_first_selected_obj();

  // moves object indiced by i to next group (down if second
  // parameter is true, up if not)
void move_to_next_group(int, int);

  // adds a new group in the scene graph
void add_group();

  // adds drawable object in group (default is 1)
void add_drawable(Drawable_object* , int);

  // removes dravable from group for indice
void remove_drawable(int, int);

  // idem, but destroy definitively the drawable (actually not working)
void delete_drawable(int, int);

  // deletes selected drawables
void delete_selection();

void delete_group(int); //a commenter?????

pthread_t get_window_thread();

  // Inits the window thread
void init_window_thread() ;

  // Displays the scene on the screen
void display();

  // Initializes the viewer (should it be private??)
void init_window();

  // set and get default size
void set_size(Size );
Size get_size();

  // for the style
void set_style(Style);
Style get_style();

  // for the precision
void set_precision(Precision );
Precision get_precision();

  //for the color
void set_color(Color,int);
  Color get_color(int); // 1 for the first color, anything for the second.

  // the number of groups in the scene graph
int get_group();

  // set the user defined widgets panel
void set_custom_panel(User_ctr_win);

  // return the GL_win member
GL_win* get_window();

  // The main loop 
void main_loop();

void set_style_in_selection(Style);

Viewer_3() : group(1){

  scale = new int[6];
  scale[0] = scale[2] = scale[4] = 0; 
  scale[1] = scale[3] = scale[5] = 500; 
  init_window();obj_color1=RED; obj_color2=BLACK;obj_size=10 ;
  obj_precision = 20 ; obj_style = FILL;
}

~Viewer_3(){
  delete scale;
  for (int i = scene->size(); i > 0; i--)
    canvas->reset_group(((char*) scene->data(i))[0]);
}

Viewer_3(GLsizei s) : group(1)  {
  scale = new int[6];
  scale[0] = scale[2] = scale[4] = 0;
  scale[1] = scale[3] = scale[5] = s; 
  init_window(); obj_color1=RED;obj_color2=BLACK; obj_size=10 ;
  obj_precision = 20 ; obj_style = FILL; default_impl=custom_win ;
}

Viewer_3(GLsizei x_max, GLsizei y_max, GLsizei z_max) : group(1)  {
  scale = new int[6];
  scale[0] = scale[2] = scale[4] = 0;
  scale[1] = x_max; 
  scale[3] = y_max; 
  scale[5] = z_max; 
  init_window(); obj_color1=RED;obj_color2=BLACK; obj_size=10 ;
  obj_precision = 20 ; obj_style = FILL; default_impl=custom_win ;
}

Viewer_3(GLsizei x_min, GLsizei x_max, GLsizei y_min, 
	 GLsizei y_max, GLsizei z_min, GLsizei z_max) : group(1)  {
  scale = new int[6];
  scale[0] = x_min; scale[1] = x_max; 
  scale[2] = y_min; scale[3] = y_max; 
  scale[4] = z_min; scale[5] = z_max; 
  init_window(); obj_color1=RED;obj_color2=BLACK; obj_size=10 ;
  obj_precision = 20 ; obj_style = FILL; default_impl=custom_win ;
}



// #########   All the callbacks #####################################



  // Callback for clip_win

static void closec_cb(Fl_Widget* w, void* v)
{ 
  Viewer_3* W= (Viewer_3*) v;
  W->clip_win->hide();
}

static void orient_cb(Fl_Widget* w, void* v)
{ 
  Viewer_3* W= (Viewer_3*) v;
  Fl_Choice* c = (Fl_Choice*) w;
  W->clip_plane=c->value();
  W->canvas->set_clip_planes(c->value());
  W->display();
}

static void wid_cb(Fl_Widget* w, void* v)
{ 
  Viewer_3* W= (Viewer_3*) v;
  Fl_Counter* c = (Fl_Counter*) w;
  W->canvas->set_clip_width(c->value(),W->clip_plane);
  W->display();
}

static void move_cb(Fl_Widget* w, void* v)
{ 
  Viewer_3* W= (Viewer_3*) v;
  Fl_Roller* c = (Fl_Roller*) w;
  W->canvas->set_clip_move(c->value(),W->clip_plane);
  W->display();
}



// Calback for parameters ##########################################


static void bg_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Choice* c = (Fl_Choice*) w;
    switch (c->value()) {
    case 0:
      W->canvas->set_bgcolor(BLACK);
      break;
    case 1:
      W->canvas->set_bgcolor(WHITE);
      break;
    case 2:
      W->canvas->set_bgcolor(GRAY);
      break;
    case 3:
      W->canvas->set_bgcolor(RED);
      break;
    case 4:
      W->canvas->set_bgcolor(YELLOW);
      break;
    case 5:
      W->canvas->set_bgcolor(ORANGE);
      break;
    case 6:
      W->canvas->set_bgcolor(VIOLET);
      break;
    case 7:
      W->canvas->set_bgcolor(PURPLE);
      break;
    case 8:
      W->canvas->set_bgcolor(DEEPBLUE);
      break;
    case 9:
      W->canvas->set_bgcolor(BLUE);
      break;
    case 10:
      W->canvas->set_bgcolor(GREEN);
      break;
    }
    W->display();
  }

static void col_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Choice* c = (Fl_Choice*) w;
    switch (c->value()) {
    case 0:
      W->set_color(BLACK,1);
      break;
    case 1:
      W->set_color(WHITE,1);
      break;
    case 2:
      W->set_color(GRAY,1);
      break;
    case 3:
      W->set_color(RED,1);
      break;
    case 4:
      W->set_color(YELLOW,1);
      break;
    case 5:
      W->set_color(ORANGE,1);
      break;
    case 6:
      W->set_color(VIOLET,1);
      break;
    case 7:
      W->set_color(PURPLE,1);
      break;
    case 8:
      W->set_color(DEEPBLUE,1);
      break;
    case 9:
      W->set_color(BLUE,1);
      break;
    case 10:
      W->set_color(GREEN,1);
      break;
    }
    W->display();
  }

static void col2_cb(Fl_Widget* w, void* v)
{
  Viewer_3* W= (Viewer_3*) v;
  Fl_Choice* c =(Fl_Choice*) w;
  switch(c->value()) {
    case 0:
      W->set_color(BLACK,2);
      break;
   case 1:
       W->set_color(WHITE,2);
      break;
   case 2:
       W->set_color(GRAY,2);
      break;
   case 3:
       W->set_color(RED,2);
      break;
   case 4:
       W->set_color(YELLOW,2);
      break;
   case 5:
       W->set_color(ORANGE,2);
      break;
   case 6:
       W->set_color(VIOLET,2);
      break;
   case 7:
       W->set_color(PURPLE,2);
      break;
   case 8:
       W->set_color(DEEPBLUE,2);
      break;
   case 9:
       W->set_color(BLUE,2);
      break;
   case 10:
       W->set_color(GREEN,2);
      break;
  }
  W->display();

}

static void size_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Counter* c = (Fl_Counter*) w;
    W->set_size(c->value());
  }
static void prec_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Counter* c = (Fl_Counter*) w;
    W->set_precision(c->value());
  }

static void sty_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Choice* c = (Fl_Choice*) w;
    switch(c->value()) {
    case 0:
      W->set_style(FILL);
      break;
    case 1:
      W->set_style(WIRE);
      break;
    case 2:
      W->set_style(RAW);
      break;
    }
  }





// Callback for light ##########################################

static void on_off_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->canvas->switch_light();
     W->display();
  }




static void close_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->light_win->hide();
    W->para_win->hide();
  }

static void XL_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Slider* s =(Fl_Slider*) w;
    W->canvas->set_X_light_pos(s->value());
    W->display();
  }

static void YL_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Slider* s =(Fl_Slider*) w;
    W->canvas->set_Y_light_pos(s->value());
    W->display();
  }

static void ZL_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Slider* s =(Fl_Slider*) w;
    W->canvas->set_Z_light_pos(s->value());
    W->display();
  }

static void white_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->canvas->set_light_diff(1,1,1);
    W->display();
  }

static void yell_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->canvas->set_light_diff(1,1,0);
    W->display();
  }
static void red_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->canvas->set_light_diff(1,0,0);
    W->display();
  }
static void green_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->canvas->set_light_diff(0,1,0);
    W->display();
  }
static void blue_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->canvas->set_light_diff(0,0,1);
    W->display();
  }

static void amb_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Slider* s= (Fl_Slider*) w;
    W->canvas->set_light_variation(s->value());
    W->display();
  
  }

static void shy_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Slider* s= (Fl_Slider*) w;

    W->canvas->set_light_shy(s->value());
    W->display();
  
  }

//###################################################################3

static void style_cb(Fl_Widget* w, void* v)
{
  static int old_val = 0;
  Viewer_3* W= (Viewer_3*) v;
  Fl_Menu_Button* b = (Fl_Menu_Button*) w;
  switch (b->value()) {
  case 0: 
    W->set_style_in_selection(FOS1);
    old_val = 0;
    break;
  case 1: 
    W->set_style_in_selection(FOS2);
    old_val = 1;
    break;
  case 2: 
    W->set_style_in_selection(FOS3);
    old_val = 2;
    break;
  case 3: 
    W->set_style_in_selection(FOS4);
    old_val = 3;
    break;
  case 4: 
    W->set_style_in_selection(FOS5);
    old_val = 4;
    break;
  case 5:
    W->update_color();
    b->value(old_val);
    break;
  }
  W->display();
}

static void user_cb(Fl_Widget* w, void* v)
{
  Viewer_3* W= (Viewer_3*) v;
  W->default_impl(W->canvas,W);
}


static void ps_cb(Fl_Widget* w, void* v)
{
  Viewer_3* W= (Viewer_3*) v;
  W->canvas->draw_ps();
}

static void light_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->light_win->show();
  }

static void para_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->para_win->show();
  }

  
    
static void clip_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    W->clip_win->show();
  }




static void graph_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Menu_Button* b = (Fl_Menu_Button*) w;
    switch (b->value()) {
    case 0: 
      W->add_group();
      break;
    case 1:
      W->delete_selection();
      break;
    case 2:
      W->reset_groups();
      break;
    }
    W->display();
  }

static void group_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    Fl_Choice* b = (Fl_Choice*) w;
    W->canvas->change_group(b->value());
    if (b->value() == 3)
      W->boucle=true;
    else
      W->boucle=false;
  }

static void mvobjdn_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    int obj = W->get_first_selected_obj();
    if (obj)
      W->move_to_next_group(obj,0);
  }

static void mvobjup_cb(Fl_Widget* w, void* v)
  { 
    Viewer_3* W= (Viewer_3*) v;
    int obj = W->get_first_selected_obj();
    if (obj)
      W->move_to_next_group(obj,1);
  }


static void scene_cb(Fl_Widget* w, void * v)
{  
  Viewer_3* W= (Viewer_3*) v;
  Fl_Multi_Browser* s= (Fl_Multi_Browser*) w;
  char* va= new char[2];
  int i;
  for (i=1 ; i<=s->size() ; i++) {
    va = (char*) s->data(i);
    int g_val = va[0];
    if ((va[1] == 0) && (!s->selected(i)))
      W->canvas->group_visible(false,g_val);
    if ((va[1] == 0) && (s->selected(i)) && (i<s->size())) {
      
      int j=  i+1;
      va = (char*) s->data(j);
      W->canvas->make_visible(g_val);
      while ((va[0] == g_val) && (j <= s->size())) {
	s->select(j);
        j++;
	if (j <= s->size())
	  va = (char*) s->data(j);
      }
    }
    else if (va[1] != 0) { 
	   if (s->selected(i))
	     W->canvas->change_visibility(g_val,va[1],true);   
           else 
	     W->canvas->change_visibility(g_val,va[1],false);   
    }
  }
  W->display();
}   

static void set_mode_cb(Fl_Widget* w, void* v)
{
  Viewer_3* W= (Viewer_3*) v;
  W->canvas->set_mode(W->mode_but->value());
  if (W->mode_but->value() != 1) {
    W->mvplan_sld->value(0);
    W->canvas->add_zplan(0);
  }
}

static void set_plan_cb(Fl_Widget* w, void* v)
{
  Viewer_3* W= (Viewer_3*) v;
  W->canvas->add_zplan(W->mvplan_sld->value());
}

static void insert_cb(Fl_Widget* w, void* v)
{
  int i=1;
  Viewer_3* W= (Viewer_3*) v;
  Fl_Menu_Button* b = (Fl_Menu_Button*) w;
  std::vector<double> p = W->canvas->get_point();
  if (W->mode_but->value() == 1) {
    switch (b->value()) {
    case 0: 
      for (i=1 ; i <= W->scene->size() ; i++)
	if ((W->scene->selected(i)) && (((char*) W->scene->data(i))[1] != 0))
	  W->canvas->add_point_to_object(((char*) W->scene->data(i))[0],((char*) W->scene->data(i))[1],p);
      break;
    case 1:
      while ((i<=W->scene->size()) && (!((W->scene->selected(i)) && (((char*)
					     W->scene->data(i))[1]== 0))))
	i++;
      if (i<=W->scene->size()) {
          std::vector<double> v = W->canvas->get_real_point(((char*)
					     W->scene->data(i))[0],p);
          Drawable_point_3<Point3>* dp=new
	    Drawable_point_3<Point3>(v[0],v[1],v[2],W->obj_color1,
				   W->obj_style,W->obj_size,W->obj_precision);
	  W->add_drawable(dp,((char*) W->scene->data(i))[0]);
      }
	break;
    case 2:
      std::vector<double> v = W->canvas->get_real_point(W->group+1,p);
      Drawable_point_3<Point3>* dp=new
	Drawable_point_3<Point3>(v[0],v[1],v[2],W->obj_color1,
				   W->obj_style,W->obj_size,W->obj_precision);
      W->add_drawable(dp,W->group +1);
    }
    W->canvas->redraw();
  }
}

};



// ################ Implementations ############################//
void Viewer_3::reset()
{
 canvas->reset();
 canvas->redraw();
}


void Viewer_3::perspective_view()
{
  ortho_but->value(0);
  canvas->set_projection(true);
  canvas->reshape();
  canvas->redraw();
}


void Viewer_3::orthogonal_view()
{
  persp_but->value(0);
  canvas->set_projection(false);
  canvas->reshape();
  canvas->redraw();
}

void Viewer_3::set_mouse_push_handler(GL_win::Mouse_click hand)
{      
  canvas->set_mouse_push_handler(hand);
}

void Viewer_3::set_mouse_release_handler(GL_win::Mouse_click hand)
{      
  canvas->set_mouse_release_handler(hand);
}

void Viewer_3::set_mouse_grab_handler(GL_win::Mouse_grab hand)
{      
  canvas->set_mouse_grab_handler(hand);
}

void Viewer_3::reset_groups()
    {
      for (int i=1; i <= scene->size() ; i++)
	if ((scene->selected(i)) && (((char*)
				      scene->data(i))[1]==0))
	  canvas->reset_group(((char*) scene->data(i))[0]);
    }

int Viewer_3::get_first_selected_obj()
    {
      for (int i=1; i <= scene->size() ; i++)
	if ((scene->selected(i)) && (((char*)
				      scene->data(i))[1]!=0)) 
	 return i; 
	  
      return 0;
    }

void Viewer_3::move_to_next_group(int i, int w)
    {
      int gr =  ((char*) scene->data(i))[0];
      int ind = ((char*) scene->data(i))[1];
  
      Drawable_object* obj = canvas->get_drawable(gr,ind);
      if (group !=1) {
	remove_drawable(gr,ind);
        remove_line(i);
	if (w) {
	  if (gr == 1) 
	    add_drawable(obj,group);
	  else 
	    add_drawable(obj,gr-1);
	}
	else {
	  if (gr == group) 
	    add_drawable(obj,1);
	  else 
	    add_drawable(obj,gr+1);
	}
      }
    }
	
void Viewer_3::remove_line(int i)
    {
      char nb[30];
      char* v;
      scene->remove(i);
      while ((i<=scene->size()) && (((char*) scene->data(i))[1] != 0)){

	v=(char*) scene->data(i);
	v[1]=v[1]-1;   
	Drawable_object* obj= canvas->get_drawable(v[0],v[1]);
	scene->data(i,v);
	sprintf(nb,"@s  %d.%d %s",v[0],v[1],obj->type);
	scene->text(i,nb);
	i++;
      }
    }

int Viewer_3::get_line_number(int gr, int i)
{
  if ((gr>group) || (gr<=0) || (group==0)) return 0;
  int ind=1;
  for (; ind<=scene->size() ; ind++)
    if ((((char*) scene->data(ind))[0] == gr) && (((char*)
						  scene->data(ind))[1]==i))
      return ind;

      
  return 0;
}

void Viewer_3::insert_in_group(int i, Drawable_object* obj)
    {
      char nb[30];
      char* v = new char(2);
      int j=1;
      v =  (char*)  scene->data(1);
      while ((i != v[0] - 1) && (j<=scene->size())) {
	v = (char*) scene->data(j);
	if (i != v[0] - 1)
         j++;
      }
        v =  (char*) scene->data(j-1);

	sprintf(nb,"@s  %d.%d %s",v[0],v[1] + 1,obj->type);
        char* v1= new char(2);
        v1[0]=v[0]; v1[1] = v[1] + 1;
	scene->insert(j,nb,v1);
	scene->select(j);
    }	   
	  
void Viewer_3::add_group()
    {
      group=group+1;
      char nb[30];
      sprintf(nb,"%s %d", "@bGroup" , group);
      char* v= new char(2);
      v[0]=group; v[1] = 0;
      scene->add(nb,v);
#ifdef USE_THREAD
      pthread_mutex_lock(&Synchronizer::sgMutex);
      canvas->add_new_group();
      pthread_mutex_unlock(&Synchronizer::sgMutex);
#else
      canvas->add_new_group();
#endif
    }

void Viewer_3::delete_group(int gr)
{
  if ((gr <= group) && (gr>0)) {
     int l=get_line_number(gr,0);
     scene->deselect();
     while ( (l<=scene->size()) && (((char*) scene->data(l))[0] == gr)){
        scene->select(l);
	l++;
     }
     delete_selection();
  }
}
void Viewer_3::add_drawable(Drawable_object* obj, int i=1)
    { 
#ifdef USE_THREAD
      pthread_mutex_lock(&Synchronizer::sgMutex);
      canvas->add_drawable(obj,i);
      pthread_mutex_unlock(&Synchronizer::sgMutex);
#else
      canvas->add_drawable(obj,i);
#endif
      if (i > group) {
        group++;
  	char nb[30];
        sprintf(nb,"%s %d", "@bGroup" , group);
	char* v= new char(2);
        v[0]=group; v[1] = 0;
        scene->add(nb, v);
	scene->select(scene->size());
	sprintf(nb,"@s  %d.%d %s",group,1,obj->type);
        v= new char(2);
	v[0]=group; 
	v[1] = 1;
        scene->add(nb,v);
	scene->select(scene->size());
      }
      else {
            insert_in_group(i,obj);
      }
      
    }

  
void Viewer_3::remove_drawable(int gr, int i)
   { 
#ifdef USE_THREAD
     pthread_mutex_lock(&Synchronizer::sgMutex);
     canvas->remove_drawable(gr,i);
     pthread_mutex_unlock(&Synchronizer::sgMutex);
#else
     canvas->remove_drawable(gr,i);
#endif
     int l=get_line_number(gr,i);
     if (l!=0) remove_line(l);
   }

void Viewer_3::delete_drawable(int gr, int i)
   {
#ifdef USE_THREAD     
     pthread_mutex_lock(&Synchronizer::sgMutex);
     canvas->free_drawable(gr,i);
     pthread_mutex_unlock(&Synchronizer::sgMutex);
#else
     canvas->free_drawable(gr,i);
#endif
   }


void Viewer_3::delete_selection()
    {
      int j=0;
      int ref=0;
      for (int i=1 ; i<=scene->size() ; i++) {
	if (scene->selected(i)) {
	  int gr =  ((char*) scene->data(i))[0];
	  if (gr!=ref) {j=0;ref=gr;}
	  int ind = ((char*) scene->data(i))[1]-j;
	  if (ind != 0) {
	    delete_drawable(gr,ind);
	    ++j;
	  }
	}
      }
      rebuild_graph();
    }

void Viewer_3::set_style_in_selection(Style s)
{
  char* v;
  for (int i=1 ; i<=scene->size() ; i++) {
    if ((scene->selected(i)) && (((char*) scene->data(i))[1] !=0)) {
      v=(char*) scene->data(i);   
      Drawable_object* obj= canvas->get_drawable(v[0],v[1]);
      obj->set_style(s);
    }
  }
}

void Viewer_3::update_color()
{
  char* v;
  for (int i=1 ; i<=scene->size() ; i++) {
    if ((scene->selected(i)) && (((char*) scene->data(i))[1] !=0)) {
      v=(char*) scene->data(i);
      std::cerr << v[0] << " " << v[1] << std::endl;   
      Drawable_object* obj= canvas->get_drawable(v[0],v[1]);
      obj->set_colors(obj_color1,obj_color2);
    }
  }

}

void Viewer_3::rebuild_graph()
{
  scene->clear();
  group=0;
  int o;
  char nb[30];
  char* v;
  canvas->clean_graph();
  Scene_graph* scg=canvas->get_scene_graph();
  Scene_graph::iterator it;
  Scene_group::iterator git;
  for (it=scg->begin() ; it!=scg->end(); it++) {
    group++;
    sprintf(nb,"%s %d", "@bGroup" , group);
    v= new char(2);
    v[0]=group; v[1] = 0;
    scene->add(nb, v);
    o=1;
    for (git=it->begin() ; git!=it->end(); git++) {
      sprintf(nb,"@s  %d.%d %s",group,o,(*git)->type);
      v= new char(2);
      v[0]=group; v[1] = o;
      scene->add(nb, v);
      o++;
    }
  }
  if (group==0) {
    group=1;    
    sprintf(nb,"%s %d", "@bGroup" , group);
    v=new char(2);
    v[0]=group; v[1] = 0;
    scene->add(nb,v);
  }
}


#ifdef USE_THREAD
pthread_t Viewer_3::get_window_thread() 
  { 
  return thr1; 
  } 

void Viewer_3::init_window_thread() 
{ 
  pthread_create(&thr1, NULL, reinterpret_cast<void *(*)(void *)>(&mainloop), this); 
} 
#endif

void Viewer_3::display()
{
  canvas->redraw();
}


void Viewer_3::init_window()
{
#ifdef USE_THREAD
#ifndef WIN32  // a temporary solution ?
  XInitThreads();
#endif
#endif
  form = new Fl_Window(size,size,"CGAL Viewer");
  
  new Fl_Box(FL_DOWN_FRAME,107,97,size-119,size-144,"");
  canvas = new GL_win(110,100,size-125,size-150,scale,0);
  //    form->resizable(canvas);

  close_but = new Fl_Button(size-85,size-35,80,25,"Quit"); 

  exit_but = new Fl_Button(size-255,size-35,80,25,"Exit thread");
  reset_but = new Fl_Button(size-170,size-35,80,25,"Reset");


  ortho_but = new Fl_Round_Button(size-210, 10, 80, 25, "Orthogonal");

  persp_but = new Fl_Round_Button(size-210, 35, 80, 25, "Perspective");
  persp_but->value(1);

  angle_sld = new Fl_Dial(size-55, 5, 50,50, "Angle");
  angle_sld->range(10,100);
  angle_sld->step(1);
  angle_sld->value(40);
  angle_sld->color(40);
  angle_sld->color2(FL_RED);

  deep_sld = new Fl_Dial(size-110, 5, 50,50, "Deep");
  deep_sld->range((scale[5] - scale[4])/10,
		  4*(scale[5] - scale[4]));
  deep_sld->step(1);
  deep_sld->value(2*(scale[5] - scale[4]));
  deep_sld->color(40);
  deep_sld->color2(FL_RED);

  mode_but= new Fl_Choice(size-105,70, 100, 25,0);
  mode_but->down_box(FL_DOWN_BOX);
  mode_but->add("View");
  mode_but->add("Insert Point");
  mode_but->add("Slice");
  mode_but->add("User mode");
  mode_but->callback(set_mode_cb,(void*)this);
  mode_but->value(0);

  insert_but = new Fl_Menu_Button(size-210, 70, 100, 25,"Insert");
  insert_but->add("in object");
  insert_but->add("in group");
  insert_but->add("in new group");
  insert_but->callback(insert_cb,(void*)this);

  style_but = new Fl_Menu_Button(size-305, 70, 120, 25,"Facet object");
  style_but->add("Wire");
  style_but->add("Wire & Vertices");
  style_but->add("Wire hidden lines");
  style_but->add("Facets");
  style_but->add("Tubes");
  style_but->add("updates colors");
  style_but->callback(style_cb,(void*)this);
  
  mvplan_sld = new Fl_Roller(size - 360, size -35, 100, 25, 0);
  mvplan_sld->type(FL_HORIZONTAL);
  mvplan_sld->callback(set_plan_cb, (void*)this);
  mvplan_sld->value(0);
  mvplan_sld->range(-1000,1000);
  mvplan_sld->step(1);

  group_but = new Fl_Choice(5,5,100,25,0);
  group_but->callback(group_cb, (void*)this);
  group_but->add("View all groups");
  group_but->add("View selected groups");
  group_but->add("Move selected groups");
  group_but->add("Highlight selection");
  group_but->value(0);

  scene = new Fl_Multi_Browser(20,97,85,size-120, "Scene Graph");
  scene->color(45);
  scene->color2(50);
  char* v =new char(2);
  v[0]=1; v[1] = 0;
  scene->add("@bGroup 1", v);
  scene->select(1);
  scene->callback(scene_cb, (void*) this);

  mvobjup_sld = new Fl_Repeat_Button (2,97,14,50,"^");
  mvobjup_sld->color(40);
  mvobjup_sld->callback(mvobjup_cb, (void*) this);

  mvobjdn_sld = new Fl_Repeat_Button (2,147,14,50,"V");
  mvobjdn_sld->color(40);
  mvobjdn_sld->callback(mvobjdn_cb, (void*) this);

  group_men = new Fl_Menu_Button(5, 35, 100, 25,"Graph");
  group_men->add("Add new group");
  group_men->add("Delete selection");
  group_men->add("Reset groups");
  group_men->callback(graph_cb,(void*)this);

  light_but = new Fl_Button(110,5,80,25,"Set light");
  light_but->callback(light_cb,(void*)this);

  ps_but = new Fl_Button(110,35,80,25,"PostScript");
  ps_but->callback(ps_cb,(void*)this);

  para_but = new Fl_Button(195,5,80,25,"Parameters");
  para_but->callback(para_cb,(void*)this);

  clip_but = new Fl_Button(280,5,80,25,"Cliping");
  clip_but->callback(clip_cb,(void*)this);

  user_but = new Fl_Button(size - 465, size -35, 100, 25, "User Panel");
  user_but->callback(user_cb,(void*)this);
  // clip plane box ________________________________________-

  clip_win = new Fl_Window(150,220,"Cliping Planes");
  clip_win->hide();

  Fl_Button* closec_but = new Fl_Button(25,185,100,25,"close");
  closec_but->callback(closec_cb,(void *) this);

  Fl_Choice* orient_but = new Fl_Choice(25, 10, 100, 25,0);
  orient_but->add("X-Y plane");
  orient_but->add("X-Z plane");
  orient_but->add("Y-Z plane");
  orient_but->callback(orient_cb,(void*)this);
  orient_but->value(0);


  Fl_Counter* wid_ct = new Fl_Counter(25,100,100,25,"Width");
  wid_ct->callback(wid_cb,(void*)this);
  wid_ct->value((float) scale[1]/5);
  wid_ct->step((float) scale[1]/250);
  wid_ct->lstep((float) scale[1]/100);
  wid_ct->range(0,scale[1]);
  wid_ct->align(FL_ALIGN_TOP);

  Fl_Roller* move_sld = new Fl_Roller(25,150,100,25,"Move");
  move_sld->type(FL_HORIZONTAL);
  move_sld->callback(move_cb, (void*)this);
  //  move_sld->value(0);
  move_sld->value((scale[0]+scale[1])/2);
  //
  move_sld->range(2*scale[0] - scale[1], 2*scale[1] - scale[0]);
  move_sld->step((float) (scale[1]-scale[0])/500);
  move_sld->align(FL_ALIGN_TOP);

  clip_win->end();
  // Objects parameters


  para_win = new Fl_Window(380,140,"Parameters");
  para_win->hide();

  Fl_Button* closep_but = new Fl_Button(150,105,80,30,"close");
  closep_but->callback(close_cb,(void*)this);


  Fl_Choice* bgcol_but = new Fl_Choice(130,20, 120 ,25,"Background");
  bgcol_but->down_box(FL_DOWN_BOX);
  bgcol_but->callback(bg_cb,(void*)this);
  bgcol_but->add("Black");
  bgcol_but->add("White");
  bgcol_but->add("Gray");
  bgcol_but->add("Red");
  bgcol_but->add("Yellow");
  bgcol_but->add("Orange");
  bgcol_but->add("Violet");
  bgcol_but->add("Purple");
  bgcol_but->add("Deep Blue");
  bgcol_but->add("Blue");
  bgcol_but->add("Green");
  bgcol_but->value(0);
  bgcol_but->align(FL_ALIGN_TOP);
 
  Fl_Choice* col_but = new Fl_Choice(5,20, 120 ,25,"First Color");
  col_but->down_box(FL_DOWN_BOX);
  col_but->callback(col_cb,(void*)this);
  col_but->add("Black");
  col_but->add("White");
  col_but->add("Gray");
  col_but->add("Red");
  col_but->add("Yellow");
  col_but->add("Orange");
  col_but->add("Violet");
  col_but->add("Purple");
  col_but->add("Deep Blue");
  col_but->add("Blue");
  col_but->add("Green");
  col_but->value(0);
  col_but->align(FL_ALIGN_TOP);

  Fl_Choice* col2_but = new Fl_Choice(255,20, 120 ,25,"Second Color");
  col2_but->down_box(FL_DOWN_BOX);
  col2_but->callback(col2_cb,(void*)this);
  col2_but->add("Black");
  col2_but->add("White");
  col2_but->add("Gray");
  col2_but->add("Red");
  col2_but->add("Yellow");
  col2_but->add("Orange");
  col2_but->add("Violet");
  col2_but->add("Purple");
  col2_but->add("Deep Blue");
  col2_but->add("Blue");
  col2_but->add("Green");
  col2_but->value(0);
  col2_but->align(FL_ALIGN_TOP);

  Fl_Counter* size_ct = new Fl_Counter(5,70,120,30,"Size");
  size_ct->callback(size_cb,(void*)this);
  size_ct->value(10);
  size_ct->step(1);
  size_ct->range(0,500);
  size_ct->align(FL_ALIGN_TOP);


  Fl_Counter* prec_ct = new Fl_Counter(130,70,120,30,"Precision");
  prec_ct->callback(prec_cb,(void*)this);
  prec_ct->value(20);
  prec_ct->step(1);
  prec_ct->range(1,100);
  prec_ct->align(FL_ALIGN_TOP);

  Fl_Choice* sty_but= new Fl_Choice(255,70, 120 ,30,0);
  sty_but->down_box(FL_DOWN_BOX);
  sty_but->add("FILL");
  sty_but->add("WIRE");
  sty_but->add("RAW");
  sty_but->callback(sty_cb,(void*)this);
  sty_but->value(0);



  para_win->end();



//------------------------Light window ------------------------
  light_win = new Fl_Window(350,200,"Light");
  light_win->hide();
  Fl_Button* closel_but = new Fl_Button(200,175,80,20,"close");
  closel_but->callback(close_cb,(void*)this);

  Fl_Round_Button* on_off_but = new Fl_Round_Button(5,175,80,20,"ON/OFF");
  on_off_but->callback(on_off_cb,(void*)this);
  on_off_but->value(1);

 Fl_Round_Button* white_but = new Fl_Round_Button(5,5,80,20,"White");
  white_but->callback(white_cb,(void*)this);
  white_but->labelcolor(FL_WHITE);
  white_but->type(FL_RADIO_BUTTON);
  white_but->value(1);

 Fl_Round_Button* yell_but = new Fl_Round_Button(5,30,80,20,"Yellow");
  yell_but->callback(yell_cb,(void*)this);
  yell_but->labelcolor(FL_YELLOW);
  yell_but->type(FL_RADIO_BUTTON);

 Fl_Round_Button* red_but = new Fl_Round_Button(5,55,80,20,"Red");
  red_but->callback(red_cb,(void*)this);
  red_but->labelcolor(FL_RED);
  red_but->type(FL_RADIO_BUTTON);

 Fl_Round_Button* green_but = new Fl_Round_Button(5,80,80,20,"Green");
  green_but->callback(green_cb,(void*)this);
  green_but->labelcolor(FL_GREEN);
  green_but->type(FL_RADIO_BUTTON);

 Fl_Round_Button* blue_but = new Fl_Round_Button(5,105,80,20,"Blue");
  blue_but->callback(blue_cb,(void*)this);
  blue_but->labelcolor(FL_BLUE);
  blue_but->type(FL_RADIO_BUTTON);

  Fl_Slider* amb = new Fl_Slider(90,20,230,15,"Intensity");
  amb->callback(amb_cb,(void*)this);
  amb->type(FL_HORIZONTAL);
  amb->value(1);
  amb->range(10,0.3);
  amb->align(FL_ALIGN_TOP);


  Fl_Slider* shy = new Fl_Slider(90,55,230,15,"Shyniness");
  shy->callback(shy_cb,(void*)this);
  shy->type(FL_HORIZONTAL);
  shy->value(50);
  shy->range(128,0);
  shy->step(1);
  shy->align(FL_ALIGN_TOP);

  Fl_Slider* XL = new Fl_Slider(90,100,230,15,"X : ");
  XL->callback(XL_cb,(void*)this);
  XL->type(FL_HORIZONTAL);
  XL->value(0);
  XL->range(-10*size,10*size);
  XL->step(1);
  XL->align(FL_ALIGN_LEFT);

  Fl_Slider* YL = new Fl_Slider(90,120,230,15,"Y : ");
  YL->callback(YL_cb,(void*)this);
  YL->type(FL_HORIZONTAL);
  YL->value(0);
  YL->range(-10*size,10*size);
  YL->step(1);
  YL->align(FL_ALIGN_LEFT);

  Fl_Slider* ZL = new Fl_Slider(90,140,230,15,"Z : ");
  ZL->callback(ZL_cb,(void*)this);
  ZL->type(FL_HORIZONTAL);
  ZL->value(2*size);
  ZL->range(-10*size,10*size);
  ZL->step(1);
  ZL->align(FL_ALIGN_LEFT);


  light_win->end();



  form->end();

}

void Viewer_3::set_size(Size s)
    {
      obj_size = s;
    }

Size Viewer_3::get_size()
    {
      return obj_size;
    }


void Viewer_3::set_style(Style s)
    {
      obj_style = s;
    }
Style Viewer_3::get_style()
    {
      return obj_style;
    }

void Viewer_3::set_precision(Precision p)
    {
      obj_precision = p;
    }
Precision Viewer_3::get_precision()
    {
      return obj_precision;
    }

void Viewer_3::set_color(Color c, int i)
    {
      if (i==1)
	obj_color1=c;
      else
	obj_color2=c;
    }
Color Viewer_3::get_color(int i)
    {
      if (i==1)
	return obj_color1;
      return obj_color1;
    }

int Viewer_3::get_group()
{
  return group;
}


void Viewer_3::set_custom_panel(User_ctr_win fct)
{
default_impl=fct;
}


GL_win* Viewer_3::get_window()
{
  return canvas;
}


void Viewer_3::main_loop()
{
  form->show();
  canvas->show();
  while(Fl::wait()) {
    Fl_Widget *obj = Fl::readqueue();
    if (obj == close_but) break;
    if (obj == exit_but) {
#ifdef USE_THREAD
      sendSignal();
#endif
    }
    if (obj == reset_but) reset();
    if (obj == ortho_but) orthogonal_view();
    if (obj == persp_but) perspective_view();
    if (obj == angle_sld) {canvas->set_angle(angle_sld->value());
    canvas->reshape();canvas->redraw();}
    if (obj == deep_sld) {canvas->set_deep(deep_sld->value());
    canvas->reshape();canvas->redraw();}
    if (group_but->value()==3) 
      while( (Fl::wait(1.0)) && (boucle)) 
	canvas->redraw();
	
     
   }
   exit(1);
}


// ----------------------------------------------------------------


void* mainloop(Viewer_3* W)
{
  W->main_loop();
  return NULL;
}


//------------------------------------------------------------------

void close_cb(Fl_Widget* w, void* v)
{

  Fl_Window* win = (Fl_Window*) v;
  delete win;
  
}


void custom_win(GL_win* win, Viewer_3* view)
{


 Fl_Window* flwin = new Fl_Window(150,100,"Empty win");


 Fl_Button* close = new Fl_Button(35,25,80,25,"Done");
 close->callback(close_cb,(void *) flwin);
    flwin->end();
    flwin->show();

}
CGAL_END_NAMESPACE
