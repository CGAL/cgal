// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : src/CGALWin/_base_window.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================



// defines the BASE_WINDOW operations declared in <LEDA/base_window.h>
// using the basic graphics routines from <LEDA/impl/x_basic.h>

#include <CGAL/LEDA/base_window.h>
#include <CGAL/LEDA/impl/x_basic.h>


#include <cstdio>
#include <cstdarg>
#include <string>
#include <cctype>
#include <ctime>
#include <cassert>

#if defined(__BORLANDC__)
using std::tolower;
using std::sprintf;
using std::isdigit;
using std::isspace;
using std::time_t;
#endif

#if defined(__unix__)
#include <csignal>
#include <sys/time.h>
#endif


const char* event_name[] = { 
 "key_press_event", 
 "key_release_event", 
 "button_press_event", 
 "button_release_event",
 "panel_press_event",
 "configure_event", 
 "exposure_event", 
 "motion_event", 
 "destroy_event", 
 "timer_event",
 "no_event" 
};


namespace CGAL {


inline double HyPot(double x, double y) { return sqrt(x*x + y*y); }


const char* LEDA::copyright_string= "powered by CGAL and LEDA"; 
const char* LEDA::copyright_window_string= "powered by CGAL and LEDA";



int  BASE_WINDOW::real_to_pix(double d) const
{ if (d > 0)
    return int(d*scaling + 0.5);
  else
    return int(d*scaling - 0.5);
}

double BASE_WINDOW::pix_to_real(int p) const { return p*one_over_scaling; }


int BASE_WINDOW::xpix(double coord) const
{ return int(xorigin + xoffset + coord*scaling + 0.5); }

int BASE_WINDOW::ypix(double coord) const
{ return int(yorigin + yoffset - coord*scaling + 0.5); }


double BASE_WINDOW::xreal(int pix) const    
{ pix -= xoffset;
  return double(pix-xorigin)*one_over_scaling; 
}

double BASE_WINDOW::yreal(int pix) const
{ pix -= yoffset;
  return double(yorigin-pix)*one_over_scaling; 
}



void BASE_WINDOW::REDRAW_FUNC(void* p, int x, int y, int wi, int he)
{ 
  BASE_WINDOW* w = (BASE_WINDOW*)p;

  // call configure() if window size has changed
  if (w->window_width  != x_window_width(w->draw_win)  ||  
      w->window_height != x_window_height(w->draw_win)) w->configure();  

  if (y <= w->panel_height) 
  { w->draw_frame(x,y,x+wi+1,y+he+1);
    x_flush_display();
    y = w->panel_height+1;
   }



  x_set_clip_rectangle(w->draw_win,x-1,y-1,wi+2,he+2);

  double x0 = w->xreal(x-1); 
  double y0 = w->yreal(y+he+1); 
  double x1 = w->xreal(x+wi+1); 
  double y1 = w->yreal(y-1); 

  if (w->redraw0 || w->redraw1 || w->redraw2 || w->redraw2_ptr || w->redraw1_ptr ) 
  { 
    drawing_mode save_mode = x_set_mode(w->draw_win,src_mode);

    if (w->grid_mode != 0) w->draw_grid(x0,y0,x1,y1);

    if (w->redraw2 || w->redraw2_ptr) {
      if (w->redraw2_ptr) {
         w->redraw2_ptr->set_params((window*)w,NULL,0,x0,y0,x1,y1);
         w->redraw2_ptr->operator()();
      }
      else w->redraw2(w,x0,y0,x1,y1);
    }
    else 
      if (w->redraw1 || w->redraw1_ptr) {
        if (w->redraw1_ptr) {
	   w->redraw1_ptr->set_params((window*)w, NULL);
	   w->redraw1_ptr->operator()();
	}
        else w->redraw1(w);
      }
      else 
        w->redraw0();
     x_set_mode(w->draw_win,save_mode);
    }
   else 
      w->clear(x0,y0,x1,y1);

  if (w->mesg_count > 0) w->draw_messages();

  w->clipping(2);
}


void BASE_WINDOW::redraw()
{
    // new - background redraw added
    if (clear_redraw_ptr) {
      clear_redraw_ptr->set_params((window*)this,NULL,0,xmin(),ymin(),xmax(),ymax());
      clear_redraw_ptr->operator()();
    }
    else if (clear_redraw) clear_redraw(this,xmin(),ymin(),xmax(),ymax());

    if (redraw2 || redraw2_ptr) {
       if (redraw2_ptr) {
         redraw2_ptr->set_params((window*)this,NULL,0,xmin(),ymin(),xmax(),ymax());
         redraw2_ptr->operator()();
       }
       else redraw2(this,xmin(),ymin(),xmax(),ymax());
    }
    else 
      if (redraw1 || redraw1_ptr){ 
         if (redraw1_ptr) {
	   redraw1_ptr->set_params((window*)this,NULL);
	   redraw1_ptr->operator()();
	 }
         else redraw1(this);
      }	 
      else 
         if (redraw0) redraw0(); // new - check added
 }



char* BASE_WINDOW::access_str(void* p)
{ //fprintf(stderr, "Use of BASE_WINDOW::access_str() is illegal.");
  return (char*)p;
 }

void BASE_WINDOW::assign_str(void* p, const char* str)
{ //fprintf(stderr, "Use of BASE_WINDOW::assign_str() is illegal.");
  strcpy((char*)p,str);
}


void BASE_WINDOW::set_border_width(int w) { x_set_border_width(draw_win,w); }
void BASE_WINDOW::set_border_color(int c) { x_set_border_color(draw_win,c); }

char* BASE_WINDOW::set_bg_pixrect(char* pr, double xorig, double yorig)
{ char* save = bg_pixrect; 
  bg_pixrect = pr; 
  x_set_bg_pixmap(draw_win,pr);
  bg_xoff = xorig;
  bg_yoff = yorig;
  return save;
}

char* BASE_WINDOW::set_bg_pixrect(char* pr)
{ char* save = bg_pixrect; 
  bg_pixrect = pr; 
  x_set_bg_pixmap(draw_win,pr);
  return save;
}


int BASE_WINDOW::set_bg_color(int c) 
{ int save = bg_color; 
  bg_color = c; 
  x_set_bg_color(draw_win,c);
  return save;
 }


int BASE_WINDOW::set_fg_color(int c)
{ int save = fg_color; 
  fg_color = c; 
  return save;
 }

int BASE_WINDOW::set_fill_color(int c)
{ int save = fill_color; 
  fill_color = c; 
  return save;
 }


int  BASE_WINDOW::mono() const { return x_display_depth() == 1; }


void BASE_WINDOW::set_buttons(int b0, int b1, int b2, int b3,
                              int s0, int s1, int s2, int s3)
{ button_table[0] = b0;
  button_table[1] = b1;
  button_table[2] = b2;
  button_table[3] = b3;
  shift_table[0]  = s0;
  shift_table[1]  = s1;
  shift_table[2]  = s2;
  shift_table[3]  = s3;
}


void BASE_WINDOW::set_buttons(int* new_value, int* save_values)
{ for(int i=0; i<4; i++) 
  { if (save_values)
    { save_values[i]   = button_table[i];
      save_values[i+4] = shift_table[i];
     }
    button_table[i] = new_value[i];
    shift_table[i]  = new_value[i+4];
  }
}

void BASE_WINDOW::std_buttons(int* save_values)
{ if (save_values)
    for(int i=0; i<4; i++) 
    { save_values[i]   = button_table[i];
      save_values[i+4] = shift_table[i];
     }
   set_buttons(NO_BUTTON,MOUSE_BUTTON(1),MOUSE_BUTTON(2),MOUSE_BUTTON(3),
               NO_BUTTON,MOUSE_BUTTON(1),MOUSE_BUTTON(2),MOUSE_BUTTON(3));
}

void BASE_WINDOW::old_buttons(int* save_values)
{ if (save_values)
    for(int i=0; i<4; i++) 
    { save_values[i]   = button_table[i];
      save_values[i+4] = shift_table[i];
     }
   set_buttons(0,1,2,3,0,-1,-2,-3);
}



void BASE_WINDOW::flush() { x_flush_display(); }


void BASE_WINDOW::clipping(int mode)
{ 
   switch (mode) {

   case 0: // entire window
     x_set_clip_rectangle(draw_win,xoffset,yoffset,window_width,window_height);
     break;

   case 1: // panel area
     //x_set_clip_rectangle(draw_win,xoffset,yoffset,panel_width,panel_height);
     x_set_clip_rectangle(draw_win,xoffset,yoffset,panel_width,panel_height+1);
     break;

   case 2: // drawing area
     x_set_clip_rectangle(draw_win,xoffset,yoffset+panel_height,window_width,
                                                  window_height-panel_height);
     break;
   }
}



void BASE_WINDOW::set_color(int c) 
{ if (!is_open() && buf_level == 0) 
{   std::cerr << "Error: Window has to be displayed before drawing.\n";
    quit_action(1);
   }
  if (c == DEF_COLOR) c = fg_color;
  x_set_color(draw_win,c); 
 }


void BASE_WINDOW::set_stipple(char* bits, int c) 
{ x_set_stipple(draw_win,bits,c); }


int BASE_WINDOW::is_open() const
{ return x_window_opened(draw_win); }

int BASE_WINDOW::is_closed() const 
{ return x_window_opened(draw_win)==0; }

int BASE_WINDOW::set_cursor(int c) 
{ return x_set_cursor(draw_win,c); }

void BASE_WINDOW::set_frame_label(const char* label) 
{ strncpy(default_frame_label,label,128);
  default_frame_label[127] = '\0';
  x_set_label(draw_win,default_frame_label); 
  x_set_icon_label(draw_win,default_frame_label); 
}

void BASE_WINDOW::set_tmp_label(const char* label) 
{ x_set_label(draw_win,label); }


void BASE_WINDOW::set_icon_label(const char* label) 
{ x_set_icon_label(draw_win,label); }


int BASE_WINDOW::load_text_font(const char* fname) 
{ return x_load_text_font(draw_win,fname); }

int BASE_WINDOW::load_italic_font(const char* fname) 
{ return x_load_italic_font(draw_win,fname); }

int BASE_WINDOW::load_bold_font(const char* fname) 
{ return x_load_bold_font(draw_win,fname); }

int BASE_WINDOW::load_fixed_font(const char* fname) 
{ return x_load_fixed_font(draw_win,fname); }

int BASE_WINDOW::load_button_font(const char* fname) 
{ return x_load_button_font(draw_win,fname); }



void BASE_WINDOW::set_text_font()  { x_set_text_font(draw_win); }
void BASE_WINDOW::set_italic_font(){ x_set_italic_font(draw_win); }
void BASE_WINDOW::set_bold_font()  { x_set_bold_font(draw_win); }
void BASE_WINDOW::set_fixed_font() { x_set_fixed_font(draw_win); }

int BASE_WINDOW::set_font(const char* fn) { return x_set_font(draw_win,fn); }


void BASE_WINDOW::reset_frame_label() 
{ set_frame_label(default_frame_label); }

double BASE_WINDOW::set_grid_dist(double d) 
{ double d0 = grid_mode;
  grid_mode = d;
  return d0;
 }

grid_style BASE_WINDOW::set_grid_style(grid_style s) 
{ grid_style s0 = g_style;
  g_style = s;
  return s0;
 }

int BASE_WINDOW::set_grid_mode(int i) 
{ int save = (int)grid_mode;
  if (grid_mode != i) init(xmin(),xmax(),ymin(),i); 
  return save;
 }

int BASE_WINDOW::set_node_width(int w)
{ if (w < 1) w = 1;
  int save = node_width;
  node_width = w;
  return save;
 }


drawing_mode BASE_WINDOW::set_mode(drawing_mode m) 
{ return x_set_mode(draw_win,m); }

int BASE_WINDOW::set_line_width(int w)
{ if (w < 1) w = 1;
  return x_set_line_width(draw_win,w);
 }

point_style BASE_WINDOW::set_point_style(point_style ps)
{ point_style ps_old = pt_style;
  pt_style = ps;
  return ps_old;
}

line_style BASE_WINDOW::set_line_style(line_style ls)
{ return x_set_line_style(draw_win,ls); }


int BASE_WINDOW::set_join_style(int js)
{ return x_set_join_style(draw_win,js); }

text_mode BASE_WINDOW::set_text_mode(text_mode m)
{ return x_set_text_mode(draw_win,m); }



coord_handler_func BASE_WINDOW::set_coord_handler(coord_handler_func f) 
{ coord_handler_func old = coord_handler;
  coord_handler= f;  
  return old;
}

const window_handler* BASE_WINDOW::set_coord_object(const window_handler& obj)
{
  const window_handler* ptr = coord_handler_ptr;
  coord_handler_ptr = & (window_handler&) obj;
  return ptr;
}


win_delete_handler_func BASE_WINDOW::set_window_delete_handler(win_delete_handler_func f) 
{ win_delete_handler_func old = win_delete_handler;
  win_delete_handler= f;  
  return old;
}

const window_handler*   BASE_WINDOW::set_window_delete_object(const window_handler& obj)
{
  const window_handler* prev = win_delete_ptr;
  win_delete_ptr = & (window_handler&) obj;
  return prev;
}

event_handler_func BASE_WINDOW::set_event_handler(event_handler_func f) 
{ event_handler_func old = user_event_handler;
  user_event_handler= f;  
  return old;
}

void BASE_WINDOW::set_redraw(win_redraw_func2 f) 
{ redraw2 = f;  
  redraw0 = 0;
  redraw1 = 0;
}

void BASE_WINDOW::set_redraw2(const window_handler& obj)
{
  redraw2_ptr = & (window_handler&) obj;
  redraw0 = 0;
  redraw1 = 0;  
}

void BASE_WINDOW::set_redraw(win_redraw_func1 f) 
{ redraw1 = f; 
  redraw0 = 0;
  redraw2 = 0;
  redraw2_ptr = 0;
}

void BASE_WINDOW::set_redraw(const window_handler& obj) 
{ 
  redraw1_ptr = & (window_handler&) obj;
  redraw0 = 0;
  redraw2 = 0;
  redraw2_ptr = 0;  
}

void BASE_WINDOW::set_redraw(win_redraw_func0 f) 
{ redraw0 = f;
  redraw1 = 0;
  redraw2 = 0;
  redraw2_ptr = 0;
}


int BASE_WINDOW::get_line_width() const
{ return x_get_line_width(draw_win); }

line_style BASE_WINDOW::get_line_style() const
{ return x_get_line_style(draw_win); }

text_mode BASE_WINDOW::get_text_mode() const
{ return x_get_text_mode(draw_win);  }

drawing_mode BASE_WINDOW::get_mode() const 
{ return x_get_mode(draw_win);       }

int BASE_WINDOW::get_cursor() const
{ return x_get_cursor(draw_win); }


double  BASE_WINDOW::text_width(const char* s) const
{ int tw_max = 0;
  const char* p = s;
  for(;;)
  { const char* q = p;
    while (*q && *q != '\n') q++;
    int tw1 = x_text_width(draw_win,p,q-p);
    if (tw1 > tw_max) tw_max = tw1;
    if (*q == '\0') break;
    p = q+1;
   }
  return pix_to_real(tw_max);
}
 
double  BASE_WINDOW::text_height(const char* s) const
{ int lc = 1;
  const char* p = s;
  while (*p)
    if (*p++ == '\n') lc++;
  return pix_to_real(lc*x_text_height(draw_win,s));
}

int BASE_WINDOW::get_border_width() const
{ return x_get_border_width(draw_win); }

int BASE_WINDOW::get_border_color() const
{ return x_get_border_color(draw_win); }



// frame window (window manager)

void BASE_WINDOW::frame_box(int& x0, int& y0, int& x1, int& y1) const
{ x_window_frame(draw_win,&x0,&y0,&x1,&y1); }

int BASE_WINDOW::xpos() const
{ int x0,y0,x1,y1;
  x_window_frame(draw_win,&x0,&y0,&x1,&y1); 
  return x0;
}

int BASE_WINDOW::ypos() const
{ int x0,y0,x1,y1;
  x_window_frame(draw_win,&x0,&y0,&x1,&y1); 
  return y0;
}

int BASE_WINDOW::frame_width() const
{ int x0,y0,x1,y1;
  x_window_frame(draw_win,&x0,&y0,&x1,&y1); 
  return x1-x0+1;
}

int BASE_WINDOW::frame_height() const
{ int x0,y0,x1,y1;
  x_window_frame(draw_win,&x0,&y0,&x1,&y1); 
  return y1-y0+1;
}




int BASE_WINDOW::read_event(int& k, double& x, double& y)
{ BASE_WINDOW* w=0;
  int e = event_handler(w,1);
  while (w != this || e == no_event) e = event_handler(w,1);
  x = mouse_xreal;
  y = mouse_yreal;
  k = mouse_key;
  return e;
 }

int BASE_WINDOW::read_event(int& k, double& x, double& y, unsigned long& t)
{ BASE_WINDOW* w=0;
  int e = event_handler(w,1);
  while (w != this || e == no_event) e = event_handler(w,1);
  x = mouse_xreal;
  y = mouse_yreal;
  k = mouse_key;
  t = event_time;
  return e;
 }


int BASE_WINDOW::read_event(int& k, double& x, double& y, unsigned long& t,
                                                          int timeout)
{ BASE_WINDOW* w=this;
  int e = event_handler(w,timeout);
  if (w == this) 
  { x = mouse_xreal;
    y = mouse_yreal;
    k = mouse_key;
    t = event_time;
   }
  //else e = no_event;
  return e;
 }


int BASE_WINDOW::get_event(int& k, double& x, double& y)
{ BASE_WINDOW* w=0;
  int e = event_handler(w,0);
  while (e != no_event && w != this) e = event_handler(w,0);
  x = mouse_xreal;
  y = mouse_yreal;
  k = mouse_key;
  return e;
 }


static BASE_WINDOW* help_win = 0;
static void help_win_redraw()
{ help_win->draw_ctext((char*)help_win->get_inf()); }

void BASE_WINDOW::open_help_window()
{ 
  if (help_win) return;

  panel_item   it = active_button;
  BASE_WINDOW* wp = active_button_win;
  const char* str = it->help_str;
  if (str == 0) return;

  int w = x_text_width(wp->draw_win,str)  + 6;
  int h = x_text_height(wp->draw_win,str) + 4;

  if (w < it->width) w = it->width;

  help_win = new BASE_WINDOW(w,1);
  help_win->set_bg_color(ivory);
  help_win->init(0,w,0);
  help_win->set_inf((void*)str);

  int x = it->xcoord;
  int y = it->ycoord + it->height + 10;

  if (y + h > wp->window_height)
    { x_window_to_screen(wp->draw_win,&x,&y);
      help_win->display(x,y,help_win);
     }
  else
     help_win->display(x,y,wp);

  for(int j=1; j<h; j++) 
  { help_win->resize(x,y,w,j);
    leda_wait(0.005);
   }

  help_win->set_redraw(help_win_redraw);

  help_win_redraw();
}

void BASE_WINDOW::close_help_window()
{ if (help_win)
  { delete help_win; 
    help_win = 0;
  }
}



int BASE_WINDOW::event_handler(BASE_WINDOW*& w, int blocking)
{
  int win = 0;
  int val,x,y,e;
  unsigned long t;

  if (blocking == 1)
   { e = no_event;
     win = 1;

     if (active_button_win) 
     { if (active_button != last_active_button)
       { active_button_win->draw_frame(); 
         last_active_button = active_button;
        }
       if (active_button == 0) active_button_win = 0;
       if (active_button != help_last_item)
       { if (help_win) close_help_window();
         if (active_button == 0)
            help_last_item = 0;
         else
            { e = x_get_next_event(&win, &val, &x, &y, &t,800);
              if (e == no_event && win == 0) //timeout
              { help_last_item = active_button;
                open_help_window();
                e = x_get_next_event(&win, &val, &x, &y, &t,2400);
                if (e == no_event && win == 0) close_help_window();
               }
             }
        }
       else
       { if (help_win)
         { e = x_get_next_event(&win, &val, &x, &y, &t,2400);
           if (e == no_event) close_help_window();
          }
        }
      }
      
     if (e == no_event) e = x_get_next_event(&win, &val, &x, &y, &t);
    }
  else
     if (blocking > 1)
       { win = w->draw_win;
         e = x_get_next_event(&win, &val, &x, &y, &t,blocking);
        }
     else
        e = x_check_next_event(&win, &val, &x, &y, &t);

  //printf("w = %d k = %d  x = %d y = %d v = %d t = %d \n",win,e,x,y,val,t);

  if (win == 0 || e == no_event) 
  { w = 0;
    return no_event;
   }

  w = (BASE_WINDOW*)x_window_inf(win);

  if (e == timer_event) 
  { if (w->timer_action || w->timer_action_ptr) 
    { drawing_mode save = x_set_mode(w->draw_win,src_mode);
    
      if (w->timer_action_ptr) {
        w->timer_action_ptr->set_params((window*)w,NULL);
        w->timer_action_ptr->operator()();
      }
      else w->timer_action(w);
      
      x_set_mode(w->draw_win,save);
     }
    return e;
   }

  int alt_down = (val >> 8) == 4;


  if (alt_down && (e == key_press_event || e == key_release_event))
  { int key_char = (val & 0xFF);
    if (e == key_press_event)
    { // wait for key release event and put it back
      int w,v,x,y;
      unsigned long t;
      while (x_get_next_event(&w,&v,&x,&y,&t) != key_release_event);
      x_put_back_event();
     }
    for(int i = 0; i < w->item_count; i++)
    { panel_item it = w->Item[i];
      if (it->kind != Button_Item || it->shortcut < 0) continue;
      char c = tolower(it->label_str[it->shortcut]);
      if ( c == key_char)
      { x = it->xcoord + w->button_w/2;
        y = it->ycoord + w->button_h/2;
        if (e == key_press_event)
          e = button_press_event;
        else
           e = button_release_event;
        break;
       }
     }
   }
    

  if (w->panel_enabled == 0 && y > 0 && y < w->panel_height) 
      return (blocking > 1) ? e : no_event;


  if (w->scroll_bar)
  { BASE_WINDOW* sw = w->scroll_bar;
    int xx = x;
    int yy = y;
    x_window_to_screen(w->draw_win,&xx,&yy);
    x_screen_to_window(sw->draw_win,&xx,&yy);
    if (xx > 0 && xx < sw->window_width && yy > 0 && yy < sw->window_height)
    { w = sw;
      x = xx;
      y = yy;
     }
   }


  if (e != configure_event && e != exposure_event &&  e != destroy_event &&
      (e != key_press_event || (val != KEY_PRINT && val != KEY_PRINT1)))
  { 
    int panel_closed = 0;
    BASE_WINDOW* wp = BASE_WINDOW::active_window;

    while (wp && wp->p_win)
    { 
      if (e == motion_event && wp->p_win->panel_menu_mode == 0) break;

      if (x > 0 && x < w->window_width && y > 0 && y < w->window_height) break;

      if (e != button_press_event &&
          x > w->hotx1 && x < w->hotx2 && y > w->hoty1 && y < w->hoty2) break;

      panel_closed++;
      x_window_to_screen(wp->draw_win,&x,&y);
      wp = wp->p_win;
      x_screen_to_window(wp->draw_win,&x,&y);
      if (wp->p_root != wp) x_grab_pointer(wp->draw_win);
      w = wp;
      BASE_WINDOW::active_window = wp;
    }

   int v = w->panel_event_handler(w->draw_win,e,val,x,y,t); 

   if (v >= 0)
   { w = w->p_root;
     w->mouse_key = v;
     w->mouse_press_time = t;
     return button_press_event;
    }

   if (panel_closed || (y > 0 && y < w->panel_height)) 
      return (blocking > 1) ? e : no_event;
  }

  // w may have been changed in panel event handler
  w = (BASE_WINDOW*)x_window_inf(win);



  switch (e) {

  case exposure_event:
      { int wi = val;
        int he = (int)t;
        //printf("EXPOSE x0 = %d  y0 = %d w = %d h = %d \n",x,y,wi,he);
        REDRAW_FUNC(w,x,y,wi,he);
        break; 
       }


 case configure_event: // window position changed
           w->window_xpos = x;
           w->window_ypos = y;
           break; 

  case key_press_event:
  case key_release_event:

           if (e == key_press_event && (val == KEY_PRINT || val == KEY_PRINT1)) 
           { char fname[256];
             bool full_color = 0;
#if defined(__win32__)
             strcpy(fname,"screenshot.wmf");
#else
             strcpy(fname,"screenshot.ps");
#endif

             BASE_WINDOW P(-1,-1,"Screenshot Panel");
             P.text_item("\\bf\\blue Save Window Screenshot");
             P.string_item("to file",fname,0,"");
             P.bool_item("full color",&full_color,panel_action_func(0),"");
             int b = P.button("write",0,panel_action_func(0));
             P.set_focus_button(b);
             P.button("cancel",1,panel_action_func(0));


             BASE_WINDOW* disp_w = w;
             if (w->window_height == w->panel_height) disp_w = 0;

             int but = 0;

             if (val != KEY_PRINT1)
             { char* pr = w->get_window_pixrect();
               but = P.panel_open(BASE_WINDOW::center,
                                  BASE_WINDOW::center,disp_w);
               w->put_pixrect(pr);
               w->del_pixrect(pr);
               w->draw_frame();
              }

             if (but == 0) w->screenshot(fname,full_color);
            }

           w->shift_key_state = val >> 8;
           w->mouse_key = val & 0xFF;
           w->mouse_xpix = x;
           w->mouse_ypix = w->ydots-y-1;
           break;


  case destroy_event:
         { if ((w->win_delete_handler == 0) && (w->win_delete_ptr ==0)) 
           { 
#if defined(__win32__)
             x_close_display(); 
#endif
             exit(0);
            }
           BASE_WINDOW::call_window = w; 
	   
	   if (w->win_delete_ptr) {
	     w->win_delete_ptr->set_params((window*)w, NULL);
	     w->win_delete_ptr->operator()();
	   }
           else w->win_delete_handler(w);
	   
           /* e = button_press_event;
              w->mouse_key = CLOSE_BUTTON;
              w->mouse_xpix = 0;
              w->mouse_ypix = 0;
              w->mouse_press_time = t;
            */
           break;
         }



  case button_press_event:
           if (val >> 8)
              w->mouse_key = w->shift_table[val & 0xFF];
           else
              w->mouse_key = w->button_table[val & 0xFF];

           w->shift_key_state = val >> 8;

           w->mouse_xpix = x;
           w->mouse_ypix = w->ydots-y-1;
           w->mouse_press_time = t;
           break;

  case button_release_event:
           if (val >>8 )
              w->mouse_key = w->shift_table[val & 0xFF];
           else
              w->mouse_key = w->button_table[val & 0xFF];
           w->shift_key_state = val >> 8;
           w->mouse_xpix = x;
           w->mouse_ypix = w->ydots-y-1;
           w->mouse_release_time = t;
           break;

  case motion_event:

         w->mouse_xpix = x;
         w->mouse_ypix = w->ydots-y-1;
         break;
    }


  if (e==motion_event || e==button_press_event || e==button_release_event)
  { if (w->grid_mode > 0) w->cursor();

    w->mouse_xreal =  w->min_xcoord + 
                      ((double)w->mouse_xpix)*w->one_over_scaling;

    w->mouse_yreal =  w->min_ycoord + 
                      ((double)w->mouse_ypix)*w->one_over_scaling;
          
    if (blocking)
    { 
      if (w->grid_mode > 0)
      { int g = (int)w->grid_mode;
        w->mouse_xreal=g*(int)(w->mouse_xreal/g+((w->mouse_xreal>0)?0.5:-0.5));
        w->mouse_yreal=g*(int)(w->mouse_yreal/g+((w->mouse_yreal>0)?0.5:-0.5));
        w->cursor();
      }
       
      // coordinate handling ...
      if (w->coord_handler_ptr){ 
       w->coord_handler_ptr->set_params((window*)w,NULL,0,w->mouse_xreal, w->mouse_yreal);
       w->coord_handler_ptr->operator()();
      }
      else { 
       if (w->coord_handler) 
        (w->coord_handler)(w,w->mouse_xreal, w->mouse_yreal);
      }
            
      if (w == read_window && w->mouse_action) /* user defined action */
      { w->mouse_action(w->mouse_last_xreal,w->mouse_last_yreal);
        w->mouse_action(w->mouse_xreal,w->mouse_yreal);
      }
    }
    
    w->mouse_last_xreal = w->mouse_xreal;
    w->mouse_last_yreal = w->mouse_yreal;
  }

  w->event_time = t;

  if (w->user_event_handler)
    e = w->user_event_handler(w,e,w->mouse_key,w->mouse_xreal,w->mouse_yreal,t);

  return e;
}


int BASE_WINDOW::read_mouse(int kind, double xstart, double ystart, double &x, 
                                                                    double &y)
{ int result = 0;

  int buf = x_test_buffer(draw_win);
  if (buf) x_stop_buffering(draw_win);

  drawing_mode dm = x_set_mode(draw_win,xor_mode);

  switch(kind) {

  case  0: // point
           result = read_mouse_action(mouse_default_action,xstart,ystart,x,y);
           break;

  case  1: // segment
           result = read_mouse_action(mouse_segment_action,xstart,ystart,x,y);
           break;

  case  2: // ray
           result = read_mouse_action(mouse_ray_action,xstart,ystart,x,y);
           break;

  case  3: // line
           result = read_mouse_action(mouse_line_action,xstart,ystart,x,y);
           break;

  case  4: // circle
           result = read_mouse_action(mouse_circle_action,xstart,ystart,x,y);
           break;

  case  5: // rectangle
           result = read_mouse_action(mouse_rect_action,xstart,ystart,x,y);
           break;

  default: result = read_mouse_action(mouse_default_action,xstart,ystart,x,y);
           break;
  }

  x_set_mode(draw_win,dm);
  if (buf) x_start_buffering(draw_win);

  return result;
}


int BASE_WINDOW::get_mouse(double& x, double& y)
{ 
  if (!is_open() && buf_level == 0) 
  { std::cerr << "Error: Window has to be displayed before mouse input.\n";
    quit_action(1);
   }

  int but = button_table[0]; // no button

  BASE_WINDOW* w;

  int e = event_handler(w,0);

  while ((w != this || e != button_press_event) && e != no_event) 
   e = event_handler(w,0);

  if (e != no_event) but = mouse_key;

  x = mouse_xreal;
  y = mouse_yreal;

  return but;
 }


unsigned long BASE_WINDOW::button_press_time()  { return mouse_press_time; }
unsigned long BASE_WINDOW::button_release_time(){ return mouse_release_time; }

int BASE_WINDOW::shift_key_down() { return shift_key_state == 1; }
int BASE_WINDOW::ctrl_key_down()  { return shift_key_state == 2; }
int BASE_WINDOW::alt_key_down()   { return shift_key_state == 4; }




#if defined(__unix__)
static BASE_WINDOW* timer_win;
 
#if defined(SIG_PF) || defined(sgi) || defined(__SUNPRO_CC)
#define leda_sig_pf SIG_PF
#else
typedef void (*leda_sig_pf)(int);
#endif

void BASE_WINDOW::timer_handler(int)
{ drawing_mode save_mode = timer_win->set_mode(src_mode);

  if (timer_win->timer_action_ptr) {
    timer_win->timer_action_ptr->set_params((window*)timer_win,NULL);
    timer_win->timer_action_ptr->operator()();
  }
  else timer_win->timer_action(timer_win);
  
  timer_win->set_mode(save_mode);
  signal(SIGALRM,(leda_sig_pf)timer_handler);
 }
#endif


void BASE_WINDOW::start_timer(int msec, win_redraw_func1 F)
{ timer_action = F;
  x_start_timer(draw_win,msec); 
#if defined(__unix__)
  timer_win = this;
  signal(SIGALRM,(leda_sig_pf)timer_handler);

  static itimerval it;
  float fsec  = msec/1000.0;
  float usec  = 1000000*(fsec - int(fsec));
  it.it_interval.tv_sec  = int(fsec);
  it.it_interval.tv_usec = int(usec);
  it.it_value.tv_sec  = int(fsec);
  it.it_value.tv_usec = int(usec);
  setitimer(ITIMER_REAL, &it, NULL);
#endif
}

void BASE_WINDOW::start_timer(int msec, const window_handler& obj)
{ timer_action_ptr = & (window_handler&) obj;
  x_start_timer(draw_win,msec); 
#if defined(__unix__)
  timer_win = this;
  signal(SIGALRM,(leda_sig_pf)timer_handler);

  static itimerval it;
  float fsec  = msec/1000.0;
  float usec  = 1000000*(fsec - int(fsec));
  it.it_interval.tv_sec  = int(fsec);
  it.it_interval.tv_usec = int(usec);
  it.it_value.tv_sec  = int(fsec);
  it.it_value.tv_usec = int(usec);
  setitimer(ITIMER_REAL, &it, NULL);
#endif
}

void BASE_WINDOW::stop_timer() 
{ timer_action = 0;
  timer_action_ptr = 0;
  x_stop_timer(draw_win); 
#if defined(__unix__)
  static itimerval it;
  setitimer(ITIMER_REAL, &it, NULL);
#endif
}


int BASE_WINDOW::read_mouse_action(mouse_action_func action, double xstart, 
                                                             double ystart, 
                                                             double &x, 
                                                             double &y)
{
  if (!is_open() && buf_level == 0)
  { std::cerr << "Error: Window has to be displayed before mouse input.\n";
    quit_action(1);
   }


  read_window = this;

  mouse_action = action;

  mouse_key = 0;

  mouse_start_xreal = xstart;
  mouse_start_yreal = ystart;

  //mouse_last_xreal = xstart;
  //mouse_last_yreal = ystart;

  if (grid_mode > 0) cursor();
  
  if (coord_handler_ptr) {
    coord_handler_ptr->set_params((window*)this,NULL,0,mouse_xreal,mouse_yreal);
    coord_handler_ptr->operator()();
  }
  else { if (coord_handler) coord_handler(this,mouse_xreal,mouse_yreal); }

  if (mouse_action) mouse_action(mouse_last_xreal,mouse_last_yreal);

  BASE_WINDOW* w;
  while (event_handler(w,1) != button_press_event || w != this);

  if (mouse_action) mouse_action(mouse_xreal,mouse_yreal);

  if (grid_mode > 0) cursor();

  x = mouse_xreal;
  y = mouse_yreal;

  mouse_action = mouse_default_action;

  return mouse_key;
}

void BASE_WINDOW::iconify() { x_iconify_window(draw_win); }


void BASE_WINDOW::close() 
{ 
  if (draw_win == 0) return;

  if (active_button_win == this) active_button_win = 0;

  x_delete_buffer(draw_win);

  if (scroll_bar)
  { delete scroll_bar;
    scroll_bar = 0;
   }

  destroy_status_window();

  if (this == active_window) active_window = 0;
  // close sub-windows
  for(int i=0; i<item_count; i++)
  { panel_item it=Item[i];
    if (it->kind == Button_Item && it->ref)
    { BASE_WINDOW* wp = (BASE_WINDOW*)it->ref;
      if (wp != this && wp->is_open()) wp->close();
     }
   }
  if (is_open()) x_close_window(draw_win); 
}


BASE_WINDOW::BASE_WINDOW(int width, int height, const char* frame_label)
{ create(width,height,frame_label); }


BASE_WINDOW::BASE_WINDOW(const char* frame_label)
{ int l  =  screen_height() - 30;
  if (l > 600) l = 600;
  create(l,l,frame_label); 
 }

BASE_WINDOW::BASE_WINDOW()
{ int l  =  screen_height() - 30;
  if (l > 600) l = 600;
  create(l,l,""); 
 }


BASE_WINDOW::~BASE_WINDOW() 
{ close();
  if (status_str) delete[] status_str;
  if (owner_item) 
  { owner_item->menu_win = 0;
    owner_item->ref = 0;
   }
  for(int i = 0; i<item_count; i++) delete Item[i];
  x_destroy_window(draw_win);
  if (--win_count==0) x_close_display(); 
 }



void BASE_WINDOW::create(int w_width, int w_height, const char* label)
{
  if (win_count==0) x_open_display();
  win_count++;

  grawin_ptr = 0;
  geowin_ptr = 0;

  for(int i=0; i<16; i++) data[i] = 0;

/*
  float fy = float(x_display_height()-28)/w_height;
  if (0 < fy && fy < 1)
  { w_height *= fy;
    w_width  *= fy;
   }

  //float fx = float(x_display_width()-16)/w_width;
  float fx = float(x_display_width())/w_width;
  if (0 < fx && fx < 1)
  { w_height *= fx;
    w_width  *= fx;
   }
*/

  active_window = this;
  p_root = this;
  p_win = 0;

  win_parent = 0;

  state = 1;
  window_width  = int(w_width);
  window_height = int(w_height);
  bg_pixrect = 0;
  bg_xoff = 0;
  bg_yoff = 0;
  bg_color = white;
  fg_color = black;
  fill_color = invisible;

  strncpy(default_frame_label,label,128);
  default_frame_label[127] = '\0';

  panel_init();

  draw_win = x_create_window(this, window_width, window_height, bg_color, 
                             default_frame_label,default_frame_label,0,
                             REDRAW_FUNC);

  panel_win = 0;

  
  if (draw_win)
  { x_set_text_font(draw_win);
    x_set_color(draw_win,black);
    x_set_mode(draw_win,src_mode);
    x_set_line_style(draw_win,solid);
    x_set_line_width(draw_win,1);
    x_set_text_mode(draw_win,transparent);
   }

  pt_style = cross_point;
  
  show_grid_cursor = 1;

  clear_on_resize = 1;

  node_width = 11;

  redraw0 = 0;
  
  redraw1 = 0;
  redraw1_ptr = 0;
  
  redraw2 = 0;
  redraw2_ptr = 0;

  clear_redraw = 0;  
  clear_redraw_ptr = 0;

  coord_handler = 0;
  coord_handler_ptr = 0;
  
  win_delete_handler = 0;
  win_delete_ptr = 0;
  
  user_event_handler = 0;

  timer_action = 0;
  timer_action_ptr = 0;

  hotx1 = 0;
  hotx2 = 0;
  hoty1 = 0;
  hoty2 = 0;

  mesg_count = 0;
  win_flush = 1;

  scaling = 1;
  one_over_scaling = 1;
  scaling_prec = 12;

  xoffset = 0;
  yoffset = 0;
  xdots = 0;
  ydots = 0;
  min_xcoord = 0;
  min_ycoord = 0;
  max_xcoord = 100;
  max_ycoord = 100;
  mouse_xpix = 0;
  mouse_ypix = 0;
  mouse_xreal = 0;
  mouse_yreal = 0;
  mouse_last_xreal = 0;
  mouse_last_yreal = 0;
  shift_key_state = 0;
  grid_mode  = 0;
  g_style  = point_grid;
  owner_item = 0;
  panel_enabled = true;
  buf_level = 0;

d3_view_pos_x = 0;
d3_view_pos_y = 0;
d3_view_pos_z = 100;
d3_view_norm_x = 0;
d3_view_norm_y = 0;
d3_view_norm_z = 100;

  status_str = 0;
  status_win = 0;

  scroll_bar = 0;

  set_buttons(NO_BUTTON,MOUSE_BUTTON(1),MOUSE_BUTTON(2),MOUSE_BUTTON(3),
              NO_BUTTON,MOUSE_BUTTON(1),MOUSE_BUTTON(2),MOUSE_BUTTON(3));

  if (coord_handler_ptr) {
    coord_handler_ptr->set_params((window*)this ,NULL,0,mouse_xreal, mouse_yreal);
    coord_handler_ptr->operator()();
  }
  else {
   if (coord_handler) 
    (coord_handler)(this,mouse_xreal,mouse_yreal);
  }
}


void BASE_WINDOW::resize(int xpos, int ypos,int width, int height)
{ window_xpos = xpos;
  window_ypos = ypos;
  x_resize_window(draw_win,xpos,ypos,width,height, (win_parent != 0));
  configure();
/*
  int win,val,x,y;
  unsigned long t;
  int e = x_check_next_event(&win, &val, &x, &y, &t);
  while (e != no_event && e != button_press_event && win == draw_win)
     e = x_check_next_event(&win, &val, &x, &y, &t);
  if (e != no_event) x_put_back_event();
*/
  clipping(0);
 }


void BASE_WINDOW::display(int xpos, int ypos, BASE_WINDOW* w)
{
  if (x_window_opened(draw_win)) 
  { clipping(0);
    redraw_panel(0);
    clipping(2);
    return;
   }

  x_open_display();

  win_parent = w;

  int pw = 0;

  if (w)
  { if (w != this)  
      pw = w->draw_win;
    else
     { pw = -1;
       w = 0;
      }
   }


  // open parent if closed
  if (w && w->is_closed()) 
      w->display(BASE_WINDOW::center,BASE_WINDOW::center,0);


  if (item_count > 0)
  { panel_width = window_width;
    place_panel_items();
    window_width  = panel_width;
    if (window_height < panel_height) window_height = panel_height;
   }


  int parent_width = (w==0) ? screen_width() : w->width();
  int parent_height= (w==0) ? screen_height() : w->height();

  int width  = window_width;
  int height = window_height;


  if (xpos > 0xFFFF)
  {  switch (xpos) {

     case BASE_WINDOW::min:    
                  xpos = 0;
                  break;

     case BASE_WINDOW::center: 
                  xpos = -(parent_width+window_width)/2;
                  if (pw == 0) xpos -= 5;
                  break;

     case BASE_WINDOW::max:    
                  xpos = -(parent_width-1);
                  break;
     }
   }
   


  if (ypos > 0xFFFF)
  {  switch (ypos) {

     case BASE_WINDOW::min:   
                 ypos = 0;
                 break;

     case BASE_WINDOW::center: { 
                 int dy = 0;
                 // parent with panel: center on drawing area 
                 if (w && w->window_height > w->panel_height) 
                   dy = w->panel_height;
                 ypos = -(parent_height+window_height+dy)/2;
                 if (pw == 0) ypos -= 16; 
                 break;
                }

     case BASE_WINDOW::max:    
                 ypos = -(parent_height-1);
                 break;
     }
  }


  if (xpos < 0)
  { xpos  = -xpos;
    width = -width;
   }
    
  if (ypos < 0)
  { ypos   = -ypos;
    height = -height;
   }
    
  window_xpos = xpos;
  window_ypos = ypos;


  if (window_height == panel_height) // panel
     x_set_bg_color(draw_win,panel_bg_color);


 x_open_window(draw_win,window_xpos,window_ypos,width,height,pw);


#if !defined(__unix__)
  x_set_focus(draw_win);
#endif

  if (pw == 0)  // frame
  { window_xpos = 0;
    window_ypos = 0;
    x_window_to_screen(draw_win,&window_xpos,&window_ypos);
   }

  configure();
  draw_frame();
  clear();

/*
  if (focus_button) 
  { int x = focus_button->xcoord + button_w/2;
    int y = focus_button->ycoord + yskip/2 + button_h/2;
    x_move_pointer(draw_win,x,y);
   }
*/

  if (getenv("LEDA_OPEN_ICONIFIED") && window_height > panel_height) 
    iconify();
}


void BASE_WINDOW::move_pointer(double x, double y)
{ x_move_pointer(draw_win,xpix(x),ypix(y)); }


void BASE_WINDOW::quit_action(int flag)
{
#if defined(__win32__)
  x_close_display();
#endif
  if (flag == 0)
    exit(0);
  else
    abort();
 }


int BASE_WINDOW::panel_open(int x, int y, BASE_WINDOW* w)
{ if (but_count==0) // panel without buttons
  { button("continue",1,panel_action_func(0));
    button("quit",0,quit_action);
   }
  display(x,y,w);
  set_focus();
  int b = read_mouse();
  close();
  return b;
 }


int BASE_WINDOW::menu_open(int xpos, int ypos, BASE_WINDOW* w)
{ 
  x_set_border_width(draw_win,0);
  p_win = w;

  place_panel_items();
  if (xpos + panel_width > w->width()) xpos = w->width() - panel_width;
  if (ypos + panel_height > w->height()) ypos = w->height() - panel_height;

  display(xpos,ypos,w);
  p_win = 0;

  BASE_WINDOW* wp;
  double x,y;
  w->panel_menu_mode = 1;
  int b = read_mouse(wp,x,y);

  if (wp != this)  b = NO_BUTTON;

    
{
  // read all remaining events 
  int e,ww,v,x,y;
  unsigned long t;
  while (x_check_next_event(&ww, &v, &x, &y, &t) != no_event);

  close();

  if (!x_display_bits_saved() && win_parent != this)
    { while (x_get_next_event(&ww, &v, &x, &y, &t) != exposure_event);
      x_put_back_event();
      event_handler(wp,1);
     }
  else
    while ((e=x_check_next_event(&ww, &v, &x, &y, &t)) != no_event)
    { if (e == exposure_event)
      { x_put_back_event();
        event_handler(wp,1);
       }
     }
 }

  return b;
}



void BASE_WINDOW::draw_frame(int x0, int y0, int x1, int y1)
{ 
  if (item_count == 0) return;

  clipping(0);
  start_buffering();
  clipping(0);
  draw_panel_items();
  if (!clear_on_resize && y1 > panel_height) y1 = panel_height;
  if (buf_level == 1) x_flush_buffer(draw_win,x0,y0,x1,y1);
  stop_buffering();
  clipping(2);
}


void BASE_WINDOW::draw_frame()
{ draw_frame(0,0,panel_width,panel_height); }


void BASE_WINDOW::set_offset(double dx, double dy)
{  xoffset = real_to_pix(dx);
   yoffset = real_to_pix(dy);
}

int BASE_WINDOW::set_precision(int prec)
{ int old_prec = scaling_prec;
  scaling_prec = prec;

  scaling = double(xdots)/(max_xcoord-min_xcoord);
  scaling = truncate(scaling,scaling_prec);

  one_over_scaling = double(max_xcoord-min_xcoord)/xdots;
  one_over_scaling = truncate(one_over_scaling,scaling_prec);

  return old_prec;
}


void BASE_WINDOW::configure()
{
  int panel_only = (panel_height == window_height); 

  window_width  = x_window_width(draw_win);
  window_height = x_window_height(draw_win);


  if (item_count > 0) 
  { panel_width = window_width;
    place_panel_items();
   }

  if (panel_only) panel_height = window_height;

  bool size_changed = (xdots != window_width || ydots != window_height);

  xdots = window_width;
  ydots = window_height;

  scaling = double(xdots)/(max_xcoord-min_xcoord);
  one_over_scaling = double(max_xcoord-min_xcoord)/xdots;

  scaling = truncate(scaling,scaling_prec);
  one_over_scaling = truncate(one_over_scaling,scaling_prec);

  if ((grid_mode > 0) && (grid_mode*scaling < 4))
  { // at least grid distance of 4 pixels
    grid_mode=0; 
	std::cerr << "warning: grid distance to small.\n";
   }

  if (grid_mode > 0)
    { max_xcoord = min_xcoord+int(xdots*one_over_scaling);
      max_ycoord = min_ycoord+int((ydots-panel_height)*one_over_scaling);
     }
  else
    { max_xcoord = min_xcoord+xdots*one_over_scaling;
      max_ycoord = min_ycoord+(ydots-panel_height)*one_over_scaling;
     }

  xorigin = (int)(-min_xcoord*scaling);
  yorigin = (int)(ydots+min_ycoord*scaling);

  mouse_xreal = 0;
  mouse_yreal = 0;

  if (!is_open() || !size_changed) return;

  int buffering = x_test_buffer(draw_win);
  x_stop_buffering(draw_win);
  x_delete_buffer(draw_win);


  if (status_win)
  { int w = window_width+3;
    int h = status_win->window_height;
    int x = -1; 
    int y = window_height-h;
    status_win->resize(x,y,w,h);
   }

  if (scroll_bar)
  { panel_item it0 = scroll_bar->Item[0];
    panel_item it1 = scroll_bar->Item[1];
    panel_item it2 = scroll_bar->Item[2];
    panel_action_func f0 = it0->action;
    panel_action_func f1 = it1->action;
    panel_action_func f2 = it2->action;
    double h = scroll_bar->panel_height - it0->height - it1->height;
    double sz = it2->height/h;
    double y = it2->ycoord - (it2->height - scroll_bar->yskip)/2;
    double pos = (y - it0->height)/(h - it2->height);
    close_scrollbar();
    open_scrollbar(f0,f1,f2,sz,pos);
   }

  if (buffering) 
    x_start_buffering(draw_win);

}



void BASE_WINDOW::init(double x0, double x1, double y0, int g_mode, int erase)
{ 
  if (x0 >= x1)
  { std::cerr << "illegal arguments in W.init: x0  >= x1\n" << x0 << " " << x1 << "\n";
    quit_action(1);
   }

  min_xcoord = x0;
  max_xcoord = x1;
  min_ycoord = y0;
  grid_mode  = g_mode;

  configure();

  if (erase)
  { //configure();
    if (is_open()) clear();
   }
}



void BASE_WINDOW::init0(double x0, double x1, double y0, int bg_off_restore)
{ 
  int bgx = xpix(bg_xoff);
  int bgy = ypix(bg_yoff);

  min_xcoord = x0;
  max_xcoord = x1;
  min_ycoord = y0;

  scaling = double(xdots)/(max_xcoord-min_xcoord);
  one_over_scaling = double(max_xcoord-min_xcoord)/xdots;

  scaling = truncate(scaling,scaling_prec);
  one_over_scaling = truncate(one_over_scaling,scaling_prec);

  if (grid_mode > 0)
    { max_xcoord = min_xcoord+int(xdots*one_over_scaling);
      max_ycoord = min_ycoord+int((ydots-panel_height)*one_over_scaling);
     }
  else
    { max_xcoord = min_xcoord+xdots*one_over_scaling;
      max_ycoord = min_ycoord+(ydots-panel_height)*one_over_scaling;
     }

  xorigin = (int)(-min_xcoord*scaling);
  yorigin = (int)(ydots+min_ycoord*scaling);

  //if (bg_pixrect)
  if (bg_off_restore)
  { bg_xoff = xreal(bgx);
    bg_yoff = yreal(bgy);
   }
}





int BASE_WINDOW::query_pix(double x, double y)
{ return x_get_pixel(draw_win,xpix(x),ypix(y)); }


void BASE_WINDOW::draw_pix(double x, double y, int col)
{ if (col == invisible) return;
  set_color(col);
  x_pixel(draw_win,xpix(x),ypix(y));
  if (win_flush) flush();
}

void BASE_WINDOW::draw_pixels(int n, double* xcoord, double* ycoord, int col)
{ if (col == invisible) return;
  set_color(col);
  int* x = new int[n];
  int* y = new int[n];
  int i;
  for(i=0;i<n;i++)
  { x[i] = xpix(xcoord[i]);
    y[i] = ypix(ycoord[i]);
   }
  x_pixels(draw_win,n,x,y); 
  delete[] x;
  delete[] y;
  if (win_flush) flush();
}



void BASE_WINDOW::draw_point(double x, double y, int col)
{ if (col == invisible) return;
  int X = xpix(x);
  int Y = ypix(y);
  int ws = x_set_line_width(draw_win,1);
  line_style ls = x_set_line_style(draw_win,solid);
  set_color(col);
  switch (pt_style) {

   case pixel_point:  x_pixel(draw_win,X,Y);
                      break;

   case cross_point:  x_point(draw_win,X,Y);
                      break;

   case plus_point:   x_plus(draw_win,X,Y);
                      break;

   case circle_point: set_color(fill_color);
                      x_fill_circle(draw_win,X,Y,2);
                      set_color(col);
                      x_circle(draw_win,X,Y,2);
                      break;

   case rect_point:   set_color(fill_color);
                      x_box(draw_win,X-2,Y-2,X+2,Y+2);
                      set_color(col);
                      x_rect(draw_win,X-2,Y-2,X+2,Y+2);
                      break;

   case disc_point:   x_fill_circle(draw_win,X,Y,2);
                      x_circle(draw_win,X,Y,2);
                      break;

   case box_point:    x_box(draw_win,X-2,Y-2,X+2,Y+2);
                      break;

  }

  x_set_line_width(draw_win,ws);
  x_set_line_style(draw_win,ls);

  if (win_flush) flush();
}


void BASE_WINDOW::draw_segment(double x1, double y1, double x2, double y2, int col)
{ if (col == invisible) return;
  set_color(col);
  x_line(draw_win, xpix(x1), ypix(y1), xpix(x2), ypix(y2));
  if (win_flush) flush();
}


void BASE_WINDOW::draw_segments(int n, double* xcoord1, double* ycoord1, 
                                       double* xcoord2, double* ycoord2, int col)
{ if (col == invisible) return;
  set_color(col);

  int* x1 = new int[n];
  int* y1 = new int[n];
  int* x2 = new int[n];
  int* y2 = new int[n];
  int i;

  for(i=0;i<n;i++)
  { x1[i] = xpix(xcoord1[i]);
    y1[i] = ypix(ycoord1[i]);
    x2[i] = xpix(xcoord2[i]);
    y2[i] = ypix(ycoord2[i]);
   }

  x_lines(draw_win,n,x1,y1,x2,y2); 

  delete[] x1;
  delete[] y1;
  delete[] x2;
  delete[] y2;
  if (win_flush) flush();
}


void BASE_WINDOW::draw_line(double x1, double y1, double x2, double y2, int col)
{
  double dx = x2 - x1;
  double dy = y2 - y1;

  if (dx == 0 && dy == 0)
  { draw_pix(x1,y1,col);
    return;
   }

/*
  double xl = xmin();
  double yl = ymin();
  double xr = xmax();
  double yr = ymax();
*/

  double xl = xreal(0);
  double xr = xreal(window_width);
  double yl = yreal(window_height);
  double yr = yreal(0);


  if (fabs(dy) < fabs(dx))
  { yl = y1 + (xl-x1)*dy/dx;
    yr = y1 + (xr-x1)*dy/dx;
   }
  else
  { xl = x1 + (yl-y1)*dx/dy;
    xr = x1 + (yr-y1)*dx/dy;
   }

  BASE_WINDOW::draw_segment(xl,yl,xr,yr,col);

}



void BASE_WINDOW::draw_ray(double x1, double y1, double x2, double y2, int col)
{
  double dx = x2 - x1;
  double dy = y2 - y1;

  if (dx == 0 && dy == 0)
  { draw_pix(x1,y1,col);
    return;
   }

  double xmin = xreal(0);
  double xmax = xreal(window_width);
  double ymin = yreal(window_height);
  double ymax = yreal(0);

  double x,y;

  if (fabs(dy) < fabs(dx))
    { x = (x1 < x2) ? xmax : xmin;
      y = y1 + (x-x1)*dy/dx;
     }
  else
    { y = (y1 < y2) ? ymax : ymin;
      x = x1 + (y-y1)*dx/dy;
     }

  BASE_WINDOW::draw_segment(x1,y1,x,y,col);

}



void BASE_WINDOW::draw_arc(double x0, double y0, double r1, double r2, 
                                                            double start, 
                                                            double angle, 
                                                            int col)
{ if (col == invisible) return;
  set_color(col);
  int R1 = real_to_pix(r1);
  int R2 = real_to_pix(r2);
  x_arc(draw_win,xpix(x0),ypix(y0),R1,R2,start,angle);
  if (win_flush) flush();
}

void BASE_WINDOW::draw_filled_arc(double x0, double y0, double r1, double r2, 
                                                                   double start,
                                                                   double angle,
                                                                   int col)
{ if (col == invisible) return;
  set_color(col);
  int R1 = real_to_pix(r1);
  int R2 = real_to_pix(r2);
  x_fill_arc(draw_win,xpix(x0),ypix(y0),R1,R2,start,angle);
  if (win_flush) flush();
}


void BASE_WINDOW::draw_node(double x0, double y0, int col)
{ if (col == invisible) return;
  int save = x_set_line_width(draw_win,1);
  double R = pix_to_real(node_width);
  draw_circle(x0,y0,R,col);
  x_set_line_width(draw_win,save);
 }

void BASE_WINDOW::draw_filled_node(double x0, double y0, int col)
{ if (col == invisible) return;
  set_color(col);
  int X = xpix(x0);
  int Y = ypix(y0);
  x_fill_circle(draw_win,X,Y,node_width);
  int save = x_set_line_width(draw_win,1);
  x_circle(draw_win,X,Y,node_width);
  set_color(black);
  x_circle(draw_win,X,Y,node_width);
  x_set_line_width(draw_win,save);
  if (win_flush) flush();
 }


void BASE_WINDOW::draw_text_node(double x0, double y0, const char *s, int col)
{ text_mode t_save = x_set_text_mode(draw_win,transparent);

  if (col == DEF_COLOR) col = bg_color;

  if (mono() && col!=black) col = white;
  draw_filled_node(x0,y0,col);
  
  draw_ctext(x0,y0,s,col);

  if (col == black || col == blue || col == violet || col == brown) 
     draw_ctext(x0,y0,s,white);
  else
     draw_ctext(x0,y0,s,black);

  x_set_text_mode(draw_win,t_save);
}

void BASE_WINDOW::draw_int_node(double x0, double y0, int i, int col)
{ char buf[16];
  CGAL_CLIB_STD::sprintf(buf,"%d",i);
  draw_text_node(x0,y0,buf,col);
 }


void BASE_WINDOW::draw_edge(double x1, double y1, double x2, double y2, int col)
{ 
  if (col == invisible) return;

  double dx = x2-x1;
  double dy = y2-y1;
  double L  = scaling*HyPot(dx,dy);

  if (L > 2*node_width)
  { set_color(col);
    double l  = double(node_width+1)/L;
    x1 += l*dx;
    x2 -= l*dx;
    y1 += l*dy;
    y2 -= l*dy;
    x_line(draw_win,xpix(x1),ypix(y1),xpix(x2),ypix(y2));
    if (win_flush) flush();
   }
}
 

void BASE_WINDOW::draw_circle(double x0, double y0, double r, int col)
{ if (col == invisible) return;
  set_color(col);
  x_circle(draw_win,xpix(x0),ypix(y0),real_to_pix(r));
  if (win_flush) flush();
 }


void BASE_WINDOW::draw_filled_circle(double x0, double y0, double r, int col)
{ if (col == invisible) return;
  set_color(col);
  int R = real_to_pix(r);
  if (R > 0)
     x_fill_circle(draw_win,xpix(x0),ypix(y0),R);
  else
     x_pixel(draw_win,xpix(x0),ypix(y0));
  if (win_flush) flush();
 }


void BASE_WINDOW::draw_ellipse(double x0, double y0, double a, double b, int col)
{ if (col == invisible) return;
  set_color(col);
  int R1 = real_to_pix(a);
  int R2 = real_to_pix(b);
  x_ellipse(draw_win,xpix(x0),ypix(y0),R1,R2);
  if (win_flush) flush();
 }


void BASE_WINDOW::draw_filled_ellipse(double x0, double y0, double a, double b, int col)
{ if (col == invisible) return;
  set_color(col);
  int R1 = real_to_pix(a);
  int R2 = real_to_pix(b);
  x_fill_ellipse(draw_win,xpix(x0),ypix(y0),R1,R2);
  if (win_flush) flush();
 }



void BASE_WINDOW::plot_xy(double x0, double x1, win_draw_func f, int col)
{
  if (col == invisible) return;
  set_color(col);

  int *xcoord;
  int *ycoord;

  int x = xpix(x0);
  int y_old = ypix((*f)(x0));
  int i,y_new;
  int size = 0;
  int n = 0;


  for(x = xpix(x0)+1; x <= xpix(x1); x++)
  { y_new = ypix((*f)(xreal(x)));
    if (y_new > y_old)
       size += (y_new-y_old+1);
    else
       size += (y_old-y_new+1);
    y_old = y_new;
   }

  xcoord = new int[size];
  ycoord = new int[size];

  y_old = ypix((*f)(x0));

  for(x = xpix(x0)+1; x <= xpix(x1); x++)
  { y_new = ypix((*f)(xreal(x)));
    if (y_new > y_old)
      for(i=y_old; i<=y_new; i++) 
      { xcoord[n] = x;
        ycoord[n] = i;
        n++;
       }
    else
      for(i=y_old; i>=y_new; i--) 
      { xcoord[n] = x;
        ycoord[n] = i;
        n++;
       }
    y_old = y_new;
  }

 x_pixels(draw_win,size,xcoord,ycoord);
 
 delete[] xcoord;
 delete[] ycoord;

  if (win_flush) flush();
}


void BASE_WINDOW::plot_yx(double y0, double y1, win_draw_func f, int col)
{
  if (col == invisible) return;
  set_color(col);

  int *xcoord;
  int *ycoord;

  int y;
  int i,x_new;
  int x_old = xpix((*f)(y0));
  int size = 0;
  int n = 0;


  for(y = ypix(y0)-1; y >= ypix(y1); y--)
  { x_new = xpix((*f)(yreal(y)));
    if (x_new > x_old)
       size += (x_new-x_old+1);
    else
       size += (x_old-x_new+1);
    x_old = x_new;
   }

  xcoord = new int[size];
  ycoord = new int[size];

  x_old = xpix((*f)(y0));

  for(y = ypix(y0)-1; y >= ypix(y1); y--)
  {
    x_new = xpix((*f)(yreal(y)));
    if (x_new > x_old)
      for(i=x_old; i<=x_new; i++) 
      { xcoord[n] = i;
        ycoord[n] = y;
        n++;
       }
    else
      for(i=x_old; i>=x_new; i--) 
      { xcoord[n] = i;
        ycoord[n] = y;
        n++;
       }
    x_old = x_new;
  }

 x_pixels(draw_win,size,xcoord,ycoord);
 
 delete[] xcoord;
 delete[] ycoord;

 if (win_flush) flush();
}

void BASE_WINDOW::adjust_polyline(int n, double *xcoord, double *ycoord)
{
 int* x = new int[n];
 int* y = new int[n];
 int i;
 for(i=0;i<n;i++)
 { x[i] = xpix(xcoord[i]);
   y[i] = ypix(ycoord[i]);
  }
 x_polyline(draw_win,n,x,y,1);
 xcoord[0] = xreal(x[0]);
 ycoord[0] = yreal(y[0]);
 xcoord[n-1] = xreal(x[n-1]);
 ycoord[n-1] = yreal(y[n-1]);
 delete[] x;
 delete[] y;
}


void BASE_WINDOW::draw_polyline(int n, double *xcoord, double *ycoord, int col)
{
 if (col == invisible) return;
 set_color(col);

 int* x = new int[n];
 int* y = new int[n];
 int i;

 for(i=0;i<n;i++)
 { x[i] = xpix(xcoord[i]);
   y[i] = ypix(ycoord[i]);
  }

 x_polyline(draw_win,n,x,y);

 delete[] x;
 delete[] y;

 if (win_flush) flush();
}


void BASE_WINDOW::draw_polygon(int n, double *xcoord, double *ycoord, int col)
{
 if (col == invisible) return;
 set_color(col);

 int* x = new int[n+1];
 int* y = new int[n+1];
 int i;

 for(i=0;i<n;i++)
 { x[i] = xpix(xcoord[i]);
   y[i] = ypix(ycoord[i]);
  }

 x[n] = x[0];
 y[n] = y[0];

 x_polyline(draw_win,n+1,x,y);

 delete[] x;
 delete[] y;

  if (win_flush) flush();
}




void BASE_WINDOW::draw_filled_polygon(int n, double *xcoord, double *ycoord, int col)
{
 if (col == invisible) return;
 set_color(col);

 int* x = new int[n];
 int* y = new int[n];
 int i;

 for(i=0;i<n;i++)
 { x[i] = xpix(xcoord[i]);
   y[i] = ypix(ycoord[i]);
  }

 x_fill_polygon(draw_win,n,x,y);

 delete[] x;
 delete[] y;

  if (win_flush) flush();
}



void BASE_WINDOW::draw_rectangle(double x1, double y1, double x2, double y2, int col)
{ if (col == invisible) return;
  set_color(col);
  x_rect(draw_win,xpix(x1),ypix(y1),xpix(x2),ypix(y2));
  if (win_flush) flush();
 }


void BASE_WINDOW::draw_filled_rectangle(double x1, double y1, double x2, double y2, int col)
{ if (col == invisible) return;
  set_color(col);
  x_box(draw_win,xpix(x1),ypix(y1),xpix(x2),ypix(y2));
  if (win_flush) flush();
 }



inline void swap_coords(int* X, int* Y, int i, int j)
{ int  tmp;
  tmp = X[i]; X[i] = X[j]; X[j] = tmp;
  tmp = Y[i]; Y[i] = Y[j]; Y[j] = tmp;
}

int circle_quadrant(int* X, int* Y, int x0, int y0, int r, int f1, int f2)
{
  int y = r;
  int x = 0;
  int e = 3-2*y;
  int n = 0;

  while (x < y)
  { X[n] = x;
    Y[n] = y;
    n++;
    x++;
    if (e >= 0) { y--; e -= 4*y; }
    e += 4*x+2;
   }

  int j = n-1;

  while(j >= 0)
  { X[n] = Y[j];
    Y[n] = X[j];
    n++;
    j--;
   }

  for(j=0; j<n; j++)
  { X[j] = x0 + f1*X[j];
    Y[j] = y0 + f2*Y[j];
   }

  if (f1 != f2)
    for(j=0; j<n/2; j++) swap_coords(X,Y,j,n-j-1);

  return n;
}



int BASE_WINDOW::construct_roundrect(int*& xc, int*& yc, double x0, double y0, 
                                                         double x1, double y1, 
                                                         double rndness)
{ if (x0 > x1) { double t = x0; x0 = x1; x1 = t; }
  if (y0 > y1) { double t = y0; y0 = y1; y1 = t; }
  double w = x1 - x0;
  double h = y1 - y0;
  double r;
  if (w > h)
     r = rndness*h/2;
  else
     r = rndness*w/2;

  int R = real_to_pix(r);

  int X0 = xpix(x0);
  int Y0 = ypix(y0);
  int X1 = xpix(x1);
  int Y1 = ypix(y1);

  xc = new int[8*(R+1)+1];
  yc = new int[8*(R+1)+1];

  int n = 0;

  n += circle_quadrant(xc+n,yc+n,X1-R,Y0-R,R,+1,+1);
  n += circle_quadrant(xc+n,yc+n,X1-R,Y1+R,R,+1,-1);
  n += circle_quadrant(xc+n,yc+n,X0+R,Y1+R,R,-1,-1);
  n += circle_quadrant(xc+n,yc+n,X0+R,Y0-R,R,-1,+1);

  if (n > 8*(R+1))
  { std::cerr << "construct round rect: internal error\n";
    x_close_display();
    exit(1);
  }

  xc[n] = xc[0];
  yc[n] = yc[0];

  return n+1;
}


void BASE_WINDOW::draw_roundbox(double x0, double y0, 
                                double x1, double y1, double rndness, int col)
{ if (col == invisible) return;
  int* xc;
  int* yc;
  int n = construct_roundrect(xc,yc,x0,y0,x1,y1,rndness);
  set_color(col);
  x_fill_polygon(draw_win,n,xc,yc);
  if (win_flush) flush();
  delete[] xc;
  delete[] yc;
}


void BASE_WINDOW::draw_roundrect(double x0, double y0, 
                                 double x1, double y1, double rndness, int col)
{ if (col == invisible) return;
  int* xc;
  int* yc;
  int n = construct_roundrect(xc,yc,x0,y0,x1,y1,rndness);
  set_color(col);
  x_polyline(draw_win,n,xc,yc);
  if (win_flush) flush();
  delete[] xc;
  delete[] yc;
}



char* BASE_WINDOW::create_bitmap(int w, int h, unsigned char* pr_data)
{ return x_create_bitmap(draw_win,w,h,pr_data); }

char* BASE_WINDOW::get_bitmap(double x1, double y1, double x2, double y2)
{ clipping(0);
  char* pr = x_create_bitmap(draw_win, xpix(x1), ypix(y2), xpix(x2), ypix(y1));
  clipping(2);
  return pr;
 }


char* BASE_WINDOW::create_pixrect(int w,int h, unsigned char* pr_data,
                                               int fcol, int bcol)
{ return x_create_pixrect(draw_win,w,h,pr_data,fcol,bcol); }

char* BASE_WINDOW::create_pixrect(const char** xpm)
{ return x_create_pixrect(draw_win,xpm); }


char* BASE_WINDOW::get_pixrect(double x1, double y1, double x2, double y2)
{ clipping(0);
  char* pr = x_create_pixrect(draw_win, xpix(x1), ypix(y2), xpix(x2), ypix(y1));
  clipping(2);
  return pr;
 }


char* BASE_WINDOW::get_window_pixrect()
{ clipping(0);
  char* pr = x_create_pixrect(draw_win,0,0, window_width-1, window_height-1);
  clipping(2);
  return pr;
 }

static BASE_WINDOW* progress_win;
static int progress_val = 0;

static void open_progress_win(void* win, int max)
{ BASE_WINDOW* P = new BASE_WINDOW(-1,-1,"pixmap2ps");
  progress_val = 0;
  P->text_item("\\bf\\blue Converting X11 Pixmap to Postscript");
  P->slider_item("line",&progress_val,0,max,0,"");
  P->display(BASE_WINDOW::center,BASE_WINDOW::center,(BASE_WINDOW*)win);
  progress_win = P;
}

static void close_progress_win()
{ leda_wait(1.5);
  delete progress_win; 
}

static void show_progress(int i)
{ if (progress_val != i)
  { progress_val = i;
    progress_win->draw_frame();
   }
 }


void BASE_WINDOW::pixrect_to_matrix(char* pmap, int* matrix)
{ x_pixrect_to_matrix(draw_win,pmap,matrix); }


void BASE_WINDOW::pixrect_to_ps(char* pmap, std::ostream& psout, bool full_color,
                                                            bool animate)
{ if (pmap == 0) return;

  int w = get_width(pmap);
  int h = get_height(pmap);
  int count = 0;

  int   buf_sz = (full_color) ? (6*w*h) : (2*w*h);

  char* buf = new char[buf_sz+1];

  if (animate) 
   { open_progress_win(this,h-1);
     x_pixrect_to_ps(draw_win,pmap,buf,full_color,show_progress);
     close_progress_win();
    }
  else
    x_pixrect_to_ps(draw_win,pmap,buf,full_color,0);

  for(int i=0; i<buf_sz; i++)
  { psout << buf[i];
    if (++count == 72) 
    { count = 0;
	  psout << std::endl;
     }
   }
  psout << std::endl;
  delete[] buf;
}

  
void BASE_WINDOW::screenshot(const char* fn, bool full_color)
{
  int fnn = strlen(fn);

  char fname[256];
  strcpy(fname,fn);

  char* p = fname + fnn;

  while (p != fname && *p != '.') p--;

  if (strcmp(p,".wmf") != 0 && strcmp(p,".ps") != 0)
  { p =  fname + fnn;
#if defined(__win32__)
  CGAL_CLIB_STD::sprintf(p,".wmf");
#else
  CGAL_CLIB_STD::sprintf(p,".ps");
#endif
   }

  if (strcmp(p,".wmf") == 0) screenshot_wmf(fname,full_color);
  if (strcmp(p,".ps") == 0) screenshot_ps(fname,full_color);
}
  



void BASE_WINDOW::screenshot_ps(const char* fname, bool full_color)
{ 
  std::ofstream o(fname);

  if (o.fail()) 
  { acknowledge("Cannot write file",fname);
    return;
   }

  int x0,y0,x1,y1;

  if (win_parent == 0)
    x_window_frame(draw_win,&x0,&y0,&x1,&y1);
  else
  { x0 = 0; 
    y0 = 0;
    x_window_to_screen(draw_win,&x0,&y0);
    x1 = window_width;
    y1 = window_height;
    x_window_to_screen(draw_win,&x1,&y1);
   }

  char* pmap = x_root_pixrect(x0,y0,x1,y1);

  int w = get_width(pmap);
  int h = get_height(pmap);

  int dw = 20;
  int dh = 20;

  time_t t;
  time(&t);

  o << "%!PS-Adobe-2.0" << std::endl;
  o << "%%Creator: LEDA Window" << std::endl;
  o << "%%CreationDate: " << CGAL_CLIB_STD::ctime(&t);
  o << "%%Pages:  1" << std::endl;
  o << "%%BoundingBox: " << " 0 0 " << (w+2*dw) << " " << (h+2*dh) << std::endl; 
  o << "%%EndComments" << std::endl;
  o << std::endl;

  o << "/draw_pixmap {" << std::endl;
  o << " 3 dict begin" << std::endl;
  o << "   /h exch def /w exch def" << std::endl;
  o << "   /pix w 3 mul string def" << std::endl;
  o << "   w h scale" << std::endl;
  o << "   w h 8 [w 0 0  0 h sub  0 h]" << std::endl;
  o << "   {currentfile pix readhexstring pop}" << std::endl;
  o << "   false 3 colorimage" << std::endl;
  o << "} def" << std::endl;
  o << std::endl;

  o << "/draw_grey_pixmap {" << std::endl;
  o << " 3 dict begin" << std::endl;
  o << "   /h exch def /w exch def" << std::endl;
  o << "   /pix w string def" << std::endl;
  o << "   w h scale" << std::endl;
  o << "   w h 8 [w 0 0  0 h sub  0 h]" << std::endl;
  o << "   {currentfile pix readhexstring pop}" << std::endl;
  o << "   image" << std::endl;
  o << "} def" << std::endl;
  o << std::endl;


  o << "gsave" << std::endl;
  o << std::endl;
  o << dw << " " << dh << " translate" << std::endl;
  o << std::endl;

  if (full_color)
	  o << w << " " << h << " draw_pixmap" << std::endl;
  else
	  o << w << " " << h << " draw_grey_pixmap" << std::endl;

  pixrect_to_ps(pmap,o,full_color,true);


  o << std::endl;
  o << "grestore" << std::endl;
  o << "showpage" << std::endl;
  o << std::endl;

  x_delete_pixrect(pmap);
}

void BASE_WINDOW::screenshot_wmf(const char* fname, bool full_color)
{ 
#if !defined(__win32__)
  std::ofstream out(fname);
  out << "WMF-format available only for win32 platform" << std::endl;
#else

  int x0,y0,x1,y1;

  if (win_parent == 0)
    x_window_frame(draw_win,&x0,&y0,&x1,&y1);
  else
  { x0 = 0; 
    y0 = 0;
    x_window_to_screen(draw_win,&x0,&y0);
    x1 = window_width;
    y1 = window_height;
    x_window_to_screen(draw_win,&x1,&y1);
   }

  char* pmap = x_root_pixrect(x0,y0,x1,y1);

  open_metafile(fname);
  put_pixrect(0,0,pmap);
  close_metafile();
  x_delete_pixrect(pmap);

#endif
}


void BASE_WINDOW::pixrect_to_clipboard(char* pmap)
{ x_pixrect_to_clipboard(draw_win,pmap); }


char* BASE_WINDOW::pixrect_from_clipboard()
{ return x_pixrect_from_clipboard(draw_win); }


void BASE_WINDOW::open_metafile(const char* fname)
{ if (strlen(fname) > 0) 
     x_open_metafile(draw_win,fname); 
  else
     x_open_metafile(draw_win,NULL); 
}

void BASE_WINDOW::close_metafile()
{ x_close_metafile(draw_win); }

void BASE_WINDOW::metafile_to_clipboard()
{ x_metafile_to_clipboard(draw_win); }


void BASE_WINDOW::load_metafile(double x0, double y0, double x1, 
                                                      double y1,
                                                      const char* fname)
{ x_load_metafile(draw_win,xpix(x0),ypix(y1),xpix(x1),ypix(y0),fname); }






void BASE_WINDOW::put_bitmap(double x, double y, char* bm, int col)
{ x_set_color(draw_win,col);
  x_insert_bitmap(draw_win, xpix(x), ypix(y), bm); 
  if (win_flush) flush();
 }


void  BASE_WINDOW::put_pixrect(double x, double y, char* prect)
{ x_insert_pixrect(draw_win, xpix(x), ypix(y), prect); 
  if (win_flush) flush();
 }

void  BASE_WINDOW::center_pixrect(double x, double y, char* prect)
{ int w,h;
  x_pixrect_dimensions(prect,&w,&h);
  int X = xpix(x - pix_to_real(w)/2);
  int Y = ypix(y - pix_to_real(h)/2);
  x_insert_pixrect(draw_win,X,Y,prect); 
  if (win_flush) flush();
 }


void  BASE_WINDOW::put_pixrect(double x, double y, char* rect, int x0, 
                                                               int y0, 
                                                               int width, 
                                                               int height)
{ x_insert_pixrect(draw_win, xpix(x), ypix(y), rect, x0, y0, width,height);
  if (win_flush) flush();
 }


void  BASE_WINDOW::put_pixrect(char* rect)
{ x_insert_pixrect(draw_win,rect); 
  if (win_flush) flush();
 }

void  BASE_WINDOW::del_bitmap(char* rect) { x_delete_bitmap(rect); }
void  BASE_WINDOW::del_pixrect(char* rect) { x_delete_pixrect(rect); }

int BASE_WINDOW::get_width(char* rect) const
{ int w,h;
  x_pixrect_dimensions(rect,&w,&h);
  return w;
}

int BASE_WINDOW::get_height(char* rect) const
{ int w,h;
  x_pixrect_dimensions(rect,&w,&h);
  return h;
}


void BASE_WINDOW::copy_rect(double x1, double y1, double x2, double y2, double x, double y)
{ x_copy_pixrect(draw_win,xpix(x1),ypix(y2),xpix(x2),ypix(y1),xpix(x),ypix(y));
  if (win_flush) flush();
 }



void BASE_WINDOW::draw_grid_lines(double xmin,double ymin, 
                                  double xmax, double ymax, 
                                  double xorig, double yorig, 
                                  double d, bool axis)
{ 
  if (xmin == xmax || ymin == ymax) return;

  if (d < 0) d = -d;

  if (xmin > xmax)
  { double t = xmin;
    xmin = xmax;
    xmax = t;
   }

  if (ymin > ymax)
  { double t = ymin;
    ymin = ymax;
    ymax = t;
   }

  int nx = int((xmax - xmin)/d) + 2;
  int ny = int((ymax - ymin)/d) + 2;

  double x0 = xorig;
  while (x0 >  xmin) x0 -= d;
  while (x0 <  xmin) x0 += d;

  double y0 = yorig;
  while (y0 >  ymin) y0 -= d;
  while (y0 <  ymin) y0 += d;

  int N = nx + ny;

  int* x1 = new int[N];
  int* y1 = new int[N];
  int* x2 = new int[N];
  int* y2 = new int[N];

  int xorig_pix = xpix(xorig);
  int yorig_pix = ypix(yorig);

  int n = 0;

  for(double i = x0; i < xmax; i += d) 
  { x1[n] = xpix(i); y1[n] = ypix(ymin);
    x2[n] = xpix(i); y2[n] = ypix(ymax);
    n++;
   }
  for(double j = y0; j < ymax; j += d)
  { x1[n] = xpix(xmin); y1[n] = ypix(j);
    x2[n] = xpix(xmax); y2[n] = ypix(j);
    n++;
   }

  int r,g,b;
  x_get_rgb(grey1,&r,&g,&b);

  int c1 = x_new_color(r+12,g+12,b+12);
  if (c1 == -1) c1 = grey1;

  int c2 = white;

  if (bg_color == white)
  { c2 = x_new_color(r-12,g-12,b-12);
    if (c2 == -1) c2 = grey1;
   }

  x_set_color(draw_win,c1);
  x_lines(draw_win,n,x1,y1,x2,y2); 

  if (axis) 
  { x_set_color(draw_win,c2);

    if (xorig_pix > xpix(xmin) && xorig_pix < xpix(xmax))
       x_line(draw_win,xorig_pix,ypix(ymin),xorig_pix,ypix(ymax)); 

    if (yorig_pix > ypix(ymax) && yorig_pix < ypix(ymin))
       x_line(draw_win,xpix(xmin),yorig_pix,xpix(xmax),yorig_pix); 


    x_set_color(draw_win,grey3);

    if (xorig_pix > xpix(xmin) && xorig_pix < xpix(xmax)  &&
        yorig_pix > ypix(ymax) && yorig_pix < ypix(ymin))
    { x_pixel(draw_win,xorig_pix,yorig_pix);
      x_pixel(draw_win,xorig_pix-1,yorig_pix);
      x_pixel(draw_win,xorig_pix+1,yorig_pix);
      x_pixel(draw_win,xorig_pix,yorig_pix-1);
      x_pixel(draw_win,xorig_pix,yorig_pix+1);
    }
   }

  x_flush_display();

  delete[] x1;
  delete[] y1;
  delete[] x2;
  delete[] y2;
}



void BASE_WINDOW::draw_grid_points(double xmin, double ymin, 
                                   double xmax, double ymax, 
                                   double xorig, double yorig, 
                                   double d, bool axis)
{

  if (xmin == xmax || ymin == ymax) return;

  if (d < 0) d = -d;

  if (xmin > xmax)
  { double t = xmin;
    xmin = xmax;
    xmax = t;
   }

  if (ymin > ymax)
  { double t = ymin;
    ymin = ymax;
    ymax = t;
   }

  int nx = int((xmax - xmin)/d) + 1;
  int ny = int((ymax - ymin)/d) + 1;

  double x0 = xorig;
  while (x0 <= xmin) x0 += d;
  while (x0 >= xmin) x0 -= d;

  double y0 = yorig;
  while (y0 <= ymin) y0 += d;
  while (y0 >= ymin) y0 -= d;

  int xorig_pix = xpix(xorig);
  int yorig_pix = ypix(yorig);

  int N = (nx+2)*(ny+2);

  int* xc = new int[N];
  int* yc = new int[N];

  int n = 0;
  double i,j;

  for(i = x0; i < xmax; i += d)
    for(j = y0; j < ymax; j += d)
    { xc[n] = xpix(i); 
      yc[n] = ypix(j); 
      n++;
     }

  set_color(text_color(bg_color));
  x_pixels(draw_win,n,xc,yc);

  if (axis) 
  { n = 0;
    for(i = x0; i < xmax; i += d)
    { int xp = xpix(i); 
      if (xp > xorig_pix) { xc[n] = xp+1; yc[n] = yorig_pix; n++; }
      if (xp < xorig_pix) { xc[n] = xp-1; yc[n] = yorig_pix; n++; }
     }
  
    for(j = y0; j < ymax; j += d)
    { int yp = ypix(j); 
      if (yp > yorig_pix) { xc[n] = xorig_pix; yc[n] = yp+1; n++; }
      if (yp < yorig_pix) { xc[n] = xorig_pix; yc[n] = yp-1; n++; }
     }

    x_pixels(draw_win,n,xc,yc);
    x_point(draw_win,xorig_pix,yorig_pix);
  }

  delete[] xc;
  delete[] yc;

  x_flush_display();
}



void BASE_WINDOW::draw_grid(double xmin, double ymin, double xmax, double ymax, 
                            double xorig, double yorig, double d)
{ 
  bool axis = true;
  line_style ls = x_set_line_style(draw_win,solid);
  x_set_color(draw_win,black);
  if (g_style == point_grid)
     draw_grid_points(xmin,ymin,xmax,ymax,xorig,yorig,d,axis);
  if (g_style == line_grid)
     draw_grid_lines(xmin,ymin,xmax,ymax,xorig,yorig,d,axis);
  x_set_line_style(draw_win,ls);
 }


void BASE_WINDOW::draw_grid(double xmin, double ymin, double xmax, double ymax)
{ if (grid_mode) draw_grid(xmin,ymin,xmax,ymax,bg_xoff,bg_yoff,grid_mode); }


void BASE_WINDOW::draw_grid()
{ draw_grid(xmin(),ymin(),xmax(),ymax()); }


void BASE_WINDOW::draw_copy_right()
{

/* 
   This function generates a copyright  message in every LEDA window 
   (and the classes derived from it). It is part of the LEDA license 
   conditions that this function is neither modified nor disabled. 
*/
  if (window_width < 150 || window_height < 100) return;

  const char*  copy_right = LEDA::copyright_window_string;

  int col = text_color(bg_color);

  if (bg_color == white) col = grey3;

  int save_col = x_set_color(draw_win,col);
  x_set_font(draw_win,"T12");

  int dx = 4;
  int dy = 3;

  if (status_win && status_win->is_open()) 
      dy += status_win->window_height;
      
  if (scroll_bar) 
      dx += (scroll_bar->window_width + 1);

  x_text(draw_win,window_width-x_text_width(draw_win,copy_right) - dx,
                  window_height-x_text_height(draw_win,copy_right) - dy,
                  copy_right);
  x_set_text_font(draw_win);
  x_set_color(draw_win,save_col);
}





void BASE_WINDOW::clear(double x0, double y0, double x1, double y1)
{ int X0 = xpix(x0);
  int Y0 = ypix(y1);
  int X1 = xpix(x1);
  int Y1 = ypix(y0);

  if (X0 > X1) { int tmp = X0; X0 = X1; X1 = tmp; }
  if (Y0 > Y1) { int tmp = Y0; Y0 = Y1; Y1 = tmp; }

  x_clear_window(draw_win,X0,Y0,X1,Y1,xpix(bg_xoff),ypix(bg_yoff));

  if (grid_mode != 0) 
  { x_set_clip_rectangle(draw_win,X0,Y0,X1-X0+1,Y1-Y0+2);
    draw_grid(x0,y0,x1,y1);
    x_set_clip_rectangle(draw_win,0,panel_height,window_width, 
                                            window_height-panel_height);
   }

  if (clear_redraw_ptr) {
     clear_redraw_ptr->set_params((window*)this,NULL,0,x0,y0,x1,y1);
     clear_redraw_ptr->operator()();
  }
  else if (clear_redraw) clear_redraw(this,x0,y0,x1,y1);
  
  if (win_parent == 0) draw_copy_right();
  
  flush();
}


void BASE_WINDOW::clear(double x0, double y0, double x1, double y1, double xo,
                                                                    double yo)
{ bg_xoff = xo;
  bg_yoff = yo;
  clear(x0,y0,x1,y1);
}


void BASE_WINDOW::clear()
{ del_messages();

  int X0 = xoffset;
  int Y0 = yoffset+panel_height;
  int X1 = window_width-1;
  int Y1 = window_height-1;

  int h = window_height;

  if (x_test_buffer(draw_win)) 
  { Y0 = 0;
    h += (panel_height+1);
  }

  x_clear_window(draw_win,X0,Y0,X1,Y1,xpix(bg_xoff),ypix(bg_yoff));

  if (grid_mode != 0) 
  { x_set_clip_rectangle(draw_win,X0,Y0,window_width, h);
    draw_grid(xreal(X0),yreal(Y1),xreal(X1),yreal(Y0));
    x_set_clip_rectangle(draw_win,0,panel_height,window_width, 
                                          window_height-panel_height);
  }

  if (clear_redraw_ptr) {
    clear_redraw_ptr->set_params((window*)this,NULL,0,xmin(),ymin(),xmax(),ymax());
    clear_redraw_ptr->operator()();
  }
  else if (clear_redraw) clear_redraw(this,xmin(),ymin(),xmax(),ymax());
  
  if (win_parent == 0) draw_copy_right();
  
  flush();
}


void BASE_WINDOW::clear(double xo, double yo)
{ bg_xoff = xo;
  bg_yoff = yo;
  clear();
}


void BASE_WINDOW::clear(int col)
{ color bg_col   = set_bg_color(col); 
  char* bg_prect = set_bg_pixrect(0); 
  clear();
  set_bg_color(bg_col);
  set_bg_pixrect(bg_prect);
}



#define DRAW_MESSAGE(i,th) { x_text(draw_win,10,panel_height+2+th*i,mesg_list[i]); }


void BASE_WINDOW::message(const char *s)
{ 
  if (mesg_count == 64) return;

  drawing_mode save_dm = x_set_mode(draw_win,xor_mode);
  text_mode    save_tm = x_set_text_mode(draw_win,transparent);

  mesg_list[mesg_count] = new char[strlen(s)+1];
  strcpy(mesg_list[mesg_count],s);

  set_color(black);
  int h = int(1.2 * x_text_height(draw_win,"H"));

  DRAW_MESSAGE(mesg_count,h);
  mesg_count++;

  x_set_text_font(draw_win);
  x_set_mode(draw_win,save_dm);
  x_set_text_mode(draw_win,save_tm);
  flush();
}

void BASE_WINDOW::draw_messages(void)
{ drawing_mode save_dm = x_set_mode(draw_win,xor_mode);
  text_mode    save_tm = x_set_text_mode(draw_win,transparent);
  set_color(black);
  int h = int(1.2 * x_text_height(draw_win,"H"));
  for(int i=0; i<mesg_count; i++) DRAW_MESSAGE(i,h);
  x_set_text_font(draw_win);
  x_set_mode(draw_win,save_dm);
  x_set_text_mode(draw_win,save_tm);
  flush();
}


void BASE_WINDOW::del_messages(void)
{ if (mesg_count > 0) draw_messages();
  while(mesg_count > 0) delete[] mesg_list[--mesg_count];
 }



/*
void BASE_WINDOW::draw_text(double x, double y, const char *s, int col)
{ if (col == invisible) return;
  set_color(col);
  x_text(draw_win,xpix(x),ypix(y),s);
  if (win_flush) flush();
}
void BASE_WINDOW::draw_ctext(double x, double y, const char *s, int col)
{ if (col == invisible) return;
  set_color(col);
  x_ctext(draw_win,xpix(x),ypix(y),s);
  if (win_flush) flush();
}
*/


void BASE_WINDOW::draw_text(double x, double y, const char *s, int col)
{ if (col == invisible) return;
  set_color(col);
  double tht = pix_to_real(x_text_height(draw_win,s));
  const char* p = s;
  for(;;)
  { const char* q = p;
    while (*q && *q != '\n') q++;
    x_text(draw_win,xpix(x),ypix(y),p,q-p);
    y -= tht; 
    if (*q == '\0') break;
    p = q+1;
   }
  if (win_flush) flush();
}



void BASE_WINDOW::draw_ctext(double x, double y, const char *s, int col)
{ double tw1 = text_width(s);
  double th1 = text_height(s);
  draw_text(x-tw1/2,y+th1/2,s,col);
}

void BASE_WINDOW::draw_ctext(const char* s, int col)
{ double x = (xmax() - xmin())/2;
  double y = (ymax() - ymin())/2 + pix_to_real(1);
  draw_ctext(x,y,s,col); 
}



static void set_attrib(const char* s, int, int& attrib)
{ int col = attrib & 0x0F;
  int att = attrib & 0xF0;

  if (strncmp(s,"\\tf",    3) == 0) { attrib = col;          return; }
  if (strncmp(s,"\\rm",    3) == 0) { attrib = col;          return; }
  if (strncmp(s,"\\bf",    3) == 0) { attrib = col | 0x10;   return; }
  if (strncmp(s,"\\tt",    3) == 0) { attrib = col | 0x20;   return; }
  if (strncmp(s,"\\it",    3) == 0) { attrib = col | 0x30;   return; }

  if (strncmp(s,"\\black", 6) == 0) { attrib = att | black;  return; }
  if (strncmp(s,"\\white", 6) == 0) { attrib = att | white;  return; }
  if (strncmp(s,"\\blue2", 6) == 0) { attrib = att | blue;   return; }
  if (strncmp(s,"\\blue",  5) == 0) { attrib = att | blue;   return; }
  if (strncmp(s,"\\red",   4) == 0) { attrib = att | red;    return; }
  if (strncmp(s,"\\green2",7) == 0) { attrib = att | green2; return; }
  if (strncmp(s,"\\green", 6) == 0) { attrib = att | green;  return; }
  if (strncmp(s,"\\yellow",7) == 0) { attrib = att | yellow; return; }
  if (strncmp(s,"\\violet",7) == 0) { attrib = att | violet; return; }
  if (strncmp(s,"\\orange",7) == 0) { attrib = att | orange; return; }
  if (strncmp(s,"\\cyan",  5) == 0) { attrib = att | cyan;   return; }
  if (strncmp(s,"\\brown", 6) == 0) { attrib = att | brown;  return; }
  if (strncmp(s,"\\grey1", 6) == 0) { attrib = att | grey1;  return; }
  if (strncmp(s,"\\grey2", 6) == 0) { attrib = att | grey2;  return; }
  if (strncmp(s,"\\grey3", 6) == 0) { attrib = att | grey3;  return; }
}



int BASE_WINDOW::split_text(const char* s, int& argc, char**& argv)
{
  int wc = 1;

  //count words
  const char* x = s;
  while (*x)
  { 
    while (*x && isspace(*x)) x++;

    const char* p=x;

    if (*p == '$' || *p == '.' /*|| *p == ','*/) 
       p++;
    else
      { if (*p == '\\') p++;
        while (*p)
        { if (isspace(*p)) break;
          if (*p == '$' || (*p == '.' && (*p+1) == ' ') /*|| *p == ','*/) break;
          if (*p == '\\'  && *(p+1) !='~') break;
          p++;
         }
       }

    if (p > x) 
    { wc++;
      x = p;
     }
   }

  // split into words
  argc = 0;
  argv = new char*[wc];

  int attrib = black;

  while (*s)
  { 
    if (argc >= wc) 
    { std::cerr << "BASE_WINDOW::error in split_text\n";
      quit_action(1);
     }

    while (*s && isspace(*s)) s++;

    const char* p = s;

    if (*p == '$' || *p == '.' /* || *p == ',' */) 
       p++;
    else
      { if (*p == '\\') p++;
        while (*p)
        { if (isspace(*p)) break;
          if (*p == '$' || (*p == '.'&& (*p+1) == ' ') /*|| *p == ','*/) break;
          if (*p == '\\'  && *(p+1) !='~') break;
          p++;
         }
       }

    if (p > s) 
    { int n = p-s;

      if (s[0] == '\\' && s[1] != '~')
      { if ( s[1] == 'n' || 
             isdigit(s[1]) || 
             (s[1] == 'c' && (s[2] == ' ' || isdigit(s[2]) || s[2] == '\0')))
           { 
             char* w = new char[2];
             w[0] = 0;
             int j = 1;
             if (s[1] == 'c')
             { w[0] = char(0x80);
               j = 2;
              }
             if ( isdigit(s[j]) ) w[0] |= ((char)atoi(s+j) & 0x7F);
             w[1] = 0;
             argv[argc++] = w;
            }
         else
             set_attrib(s,n,attrib);
         s = p;
         continue;
       }

      char* w = new char[n+2];
      int j = 0;
      w[j++] = (char)attrib;
      for(int i=0; i<n; i++) 
      { char c = s[i];
        if (c == '\\' && s[i+1] == '~' ) continue;
        if (c == '~'  && (i == 0 || s[i-1] != '\\')) c = ' ';
        w[j++] = c;
       }
      w[j++] = '\0';
      argv[argc++] = w;
      s = p;
     }
   }

  return argc;
}



int BASE_WINDOW::format_text(int argc, char** argv,
                             int xmin, int xmax, int ymin, int dy, int win) 
{
  int   y = ymin;
  int   i = 0;

  bool math = false;

  while (i < argc && (win==0 || y < window_height))
  { int j = i;
    int x = xmin;
    int x_off = 0;
    int w = x_text_width(draw_win,argv[i]+1);

    if (w == 0)  
       w = xmax-x;
    else
      if (i < argc-1) 
      { char c = argv[i+1][1];
        if (c != ',' && c != '.') w += x_text_width(draw_win," ");
       }

    if (x+w > xmax) i++;

    while(x+w <= xmax)
    { x += w;
      if (++i >= argc) break;
      char* s = argv[i];

      if (!math)
      { int att = s[0] & 0xF0;
        switch (att) {
          case 0x10: x_set_bold_font(draw_win);
                     break;
          case 0x20: x_set_fixed_font(draw_win);
                     break;
          case 0x30: x_set_italic_font(draw_win);
                     break;
          default:   x_set_text_font(draw_win);
                     break;
         }
       }

      if (s[1] == '$') 
       { math = !math;
         if (math)
            x_set_italic_font(draw_win);
         else
            x_set_text_font(draw_win);
         w = 0;
        }
      else
       { w = x_text_width(draw_win,s+1);
         if (w==0)  // newline
         { w = xmax-x;
           if (s[0] & 0x80) // center
             x_off = w/2;
          }
        else 
          if (i < argc-1) 
          { char c = argv[i+1][1];
            if (c != ',' && c != '.') w += x_text_width(draw_win," ");
           }
       }

      if (math)
      { char* p=s+1;
        while (*p && *p != '_') p++;
        if (*p == '_') w -= x_text_width(draw_win,"_"); 
       }

     }

    int d = 0;
    int r = 0;

    if (i < argc && i-j > 1)
    { d = (xmax-x)/(i-j-1);
      r = (xmax-x)%(i-j-1);
     }

 
    x = xmin + x_off;

    int dy1 = dy;

    for(int k=j; k < i; k++)
    { char* s = argv[k];
      if (s[1] == 0) // special text item 
      { dy1 += (s[0] & 0x7F);
        continue;
       }

      int dyfix = 0;

      if (!math)
      { int att = s[0] & 0xF0;
        switch (att) {
          case 0x10: x_set_bold_font(draw_win);
                     dyfix = 0;
                     break;
          case 0x20: x_set_fixed_font(draw_win);
                     dyfix = 1;
                     break;
          case 0x30: x_set_italic_font(draw_win);
                     dyfix = 0;
                     break;
          default:   x_set_text_font(draw_win);
                     dyfix = 0;
                     break;
         }
       }

      if (s[1] == '$') 
      { math = !math;
        if (math)
          x_set_italic_font(draw_win);
        else
          x_set_text_font(draw_win);
        continue;
       }

      int w = x_text_width(draw_win,s+1);

      if (win) x_set_color(win,color(s[0]&0x0F));

      char* p=s+1;
      while (*p && *p != '_') p++;

      if (math && *p == '_') 
        { if (win) 
          { *p = '\0';
            x_text(win,x,y,s+1); 
            int y_sub = y + x_text_height(draw_win,s+1)/3;
            x_text(win,x+x_text_width(win,s+1),y_sub,p+1); 
            *p = '_';
           }
          w -= x_text_width(draw_win,"_"); 
         }
      else
        if (win) x_text(win,x,y+dyfix,s+1); 

      x += w;

      if (k < i-1) 
      { char c = argv[k+1][1];
        if (c != ',' && c != '.')
        { x += x_text_width(draw_win," ");
          x += d;
          if (r-- > 0) x++;
         }
       }
     }
    y += dy1;
   }

  x_set_text_font(draw_win);
  x_set_color(draw_win,black);

  return y;
}



double BASE_WINDOW::text_box(double x0, double x1, double y1, const char* s, bool draw)
{ char** argv;
  int    argc;
  int    win = draw ? draw_win : 0;
  x_set_text_font(draw_win);
  split_text(s,argc,argv);
  int dy = x_text_height(draw_win,s);
  int y0 = format_text(argc, argv, xpix(x0),xpix(x1),ypix(y1),dy,win);
  for(int i=0; i<argc; i++) delete[] argv[i];
  delete[] argv;
  return yreal(y0);
}
  




void BASE_WINDOW::cursor(void)
{ if (show_grid_cursor && this == read_window)
  { x_set_read_gc(draw_win);
    int X = xpix(mouse_xreal);
    int Y = ypix(mouse_yreal);
    x_line(draw_win, X,Y,X+8,Y);
    x_line(draw_win, X,Y,X-8,Y);
    x_line(draw_win, X,Y,X,Y+8);
    x_line(draw_win, X,Y,X,Y-8);
    x_reset_gc(draw_win);
   }
}

// define static members

void BASE_WINDOW::do_not_open_display(int x)
{ x_do_not_open_display(x); }

char* BASE_WINDOW::root_pixrect(int x0, int y0, int x1, int y1)
{ return x_root_pixrect(x0,y0,x1,y1); }



BASE_WINDOW* BASE_WINDOW::active_window = 0;
BASE_WINDOW* BASE_WINDOW::read_window = 0;
int          BASE_WINDOW::win_count = 0;

panel_item   BASE_WINDOW::active_button = 0;
panel_item   BASE_WINDOW::last_active_button = 0;
BASE_WINDOW* BASE_WINDOW::active_button_win = 0;
panel_item   BASE_WINDOW::help_last_item = 0;

int  BASE_WINDOW::screen_width(void)  { return x_display_width(); }
int  BASE_WINDOW::screen_height(void) { return x_display_height(); }

void BASE_WINDOW::mouse_default_action(double,double) 
{ /* do nothing */}

void BASE_WINDOW::mouse_segment_action(double x, double y) 
{ double x0 = read_window->mouse_start_xreal;
  double y0 = read_window->mouse_start_yreal;
  x_set_read_gc(read_window->draw_win);
  read_window->draw_segment(x0,y0,x,y,black); 
  x_reset_gc(read_window->draw_win);
 }

void BASE_WINDOW::mouse_line_action(double x, double y) 
{ double x0 = read_window->mouse_start_xreal;
  double y0 = read_window->mouse_start_yreal;
  int dw = read_window->draw_win;
  int X = read_window->xpix(x0);
  int Y = read_window->ypix(y0);
  x_set_read_gc(read_window->draw_win);
  read_window->draw_line(x0,y0,x,y,black); 
  int csave = x_set_color(dw,black);
  x_circle(dw,X,Y,2);
  x_set_color(dw,csave);
  x_reset_gc(read_window->draw_win);
}

void BASE_WINDOW::mouse_ray_action(double x, double y) 
{ double x0 = read_window->mouse_start_xreal;
  double y0 = read_window->mouse_start_yreal;
  x_set_read_gc(read_window->draw_win);
  read_window->draw_ray(x0,y0,x,y,black); 
  x_reset_gc(read_window->draw_win);
}


void BASE_WINDOW::mouse_rect_action(double x, double y)
{ double x0 = read_window->mouse_start_xreal;
  double y0 = read_window->mouse_start_yreal;
  x_set_read_gc(read_window->draw_win);
  read_window->draw_rectangle(x0,y0,x,y,black);
  x_reset_gc(read_window->draw_win);
 }

void BASE_WINDOW::mouse_circle_action(double x, double y)
{ double x0 = read_window->mouse_start_xreal;
  double y0 = read_window->mouse_start_yreal;
  int dw = read_window->draw_win;
  int X = read_window->xpix(x0);
  int Y = read_window->ypix(y0);
  x_set_read_gc(read_window->draw_win);
  read_window->draw_circle(x0,y0,HyPot(x-x0,y-y0),black);
  int csave = x_set_color(dw,black);
  x_point(dw,X,Y);
  x_set_color(dw,csave);
  x_reset_gc(read_window->draw_win);
 }


int BASE_WINDOW::text_color(int col)
{ 

  if (col==black || col==red || col==blue || col==violet || col==brown 
      || col==pink || col==blue2 || col == grey2 || col==grey3) 
   return white;

  int r,g,b;
  x_get_rgb(col,&r,&g,&b);

  if (r+g+b < 256)
     return white;
  else
     return black;
}



void BASE_WINDOW::start_buffering() 
{ if (buf_level++ > 0) return; 
  if (x_start_buffering(draw_win)) clear();
 }

void BASE_WINDOW::start_buffering(int w, int h) 
{ buf_level = 1;
  if (window_height > panel_height) h += panel_height;
  x_start_buffering(draw_win,w,h);
  window_width  = w;
  window_height = h;
  clipping(0);
 }


void BASE_WINDOW::set_buffer(char* pr) { x_set_buffer(draw_win,pr); }

bool BASE_WINDOW::is_buffering() { return x_test_buffer(draw_win) != 0; }


void BASE_WINDOW::stop_buffering()  
{ if (buf_level > 0) buf_level--;
  if (buf_level == 0)
  { x_stop_buffering(draw_win); 
    window_width  = x_window_width(draw_win);
    window_height = x_window_height(draw_win);
    clipping(2);
   }
}

void BASE_WINDOW::stop_buffering(char*& pr)  
{ if (buf_level > 0) buf_level--;
  if (buf_level == 0)
  { x_stop_buffering(draw_win,&pr); 
    window_width  = x_window_width(draw_win);
    window_height = x_window_height(draw_win);
    clipping(2);
   }
}


void BASE_WINDOW::delete_buffer() { x_delete_buffer(draw_win); }  


void BASE_WINDOW::flush_buffer(double x0, double y0, double x1, double y1)
{ //if (buf_level <= 1)
    x_flush_buffer(draw_win,xpix(x0),ypix(y0),xpix(x1),ypix(y1)); 
  }

void BASE_WINDOW::flush_buffer() 
{ if (buf_level > 1) return;
  int x0 = 0;
  int y0 = panel_height;
  int x1 = window_width-1;
  int y1 = window_height-1;
  x_flush_buffer(draw_win,x0,y0,x1,y1,xoffset,yoffset);
}

void BASE_WINDOW::flush_buffer(double dx, double dy) 
{ if (buf_level > 1) return;
  int x0 = 0;
  int y0 = panel_height;
  int x1 = window_width-1;
  int y1 = window_height-1;
  x_flush_buffer(draw_win, x0, y0, x1, y1, xoffset+real_to_pix(dx),
                                           yoffset+real_to_pix(dy)); 
}

void BASE_WINDOW::flush_buffer(double dx, double dy, 
                               double x0, double y0, double x1, double y1)
{ if (buf_level <= 1)
    x_flush_buffer(draw_win,xpix(x0),ypix(y0),xpix(x1),ypix(y1),
                                              real_to_pix(dx),real_to_pix(dy)); 
 }



// static members


void BASE_WINDOW::default_coord_handler(BASE_WINDOW* win, double x, double y)
{ char s[128];
  CGAL_CLIB_STD::sprintf(s,"%8.2f %8.2f", x,y);
  text_mode save_tm = x_set_text_mode(win->draw_win,opaque);
  drawing_mode save_dm = x_set_mode(win->draw_win,src_mode);
  int save_col = x_set_color(win->draw_win, black);
  x_text(win->draw_win,win->window_width-x_text_width(win->draw_win,s)-2, 
                                                 win->panel_height+2,s);
  x_set_text_mode(win->draw_win,save_tm);
  x_set_mode(win->draw_win,save_dm);
  x_set_color(win->draw_win, save_col);
}

int BASE_WINDOW::read_event(BASE_WINDOW*& wp, int& k, double& x, double& y)
{ read_window = 0;
  int e = BASE_WINDOW::event_handler(wp,1);
  if (wp)
  { x = wp->mouse_xreal;
    y = wp->mouse_yreal;
    k = wp->mouse_key;
   }
  return e;
 }


int BASE_WINDOW::get_event(BASE_WINDOW*& wp, int& k, double& x, double& y)
{ read_window = 0;
  int e = BASE_WINDOW::event_handler(wp,0);
  if (wp)
  { x = wp->mouse_xreal;
    y = wp->mouse_yreal;
    k = wp->mouse_key;
   }
  return e;
 }


void BASE_WINDOW::put_back_event() { x_put_back_event(); }


int BASE_WINDOW::read_mouse(BASE_WINDOW*& w, double& x, double& y)
{
  read_window = 0;

  int e = BASE_WINDOW::event_handler(w,1);
 
  while (e != button_press_event && e != key_press_event) 
      e = BASE_WINDOW::event_handler(w,1);

  x = w->mouse_xreal;
  y = w->mouse_yreal;
  return w->mouse_key;
}


int BASE_WINDOW::get_mouse(BASE_WINDOW*& w, double& x, double& y)
{ 
  int but;

  if (BASE_WINDOW::event_handler(w,0) == button_press_event) 
      but = w->mouse_key;
  else
      but = w->button_table[0];  // no button

  x = w->mouse_xreal;
  y = w->mouse_yreal;

  return but;
}

void BASE_WINDOW::grab_mouse()   { x_grab_pointer(draw_win); }

void BASE_WINDOW::ungrab_mouse() { x_ungrab_pointer(); }



void* BASE_WINDOW::get_inf(int i) const { return data[i]; }

void* BASE_WINDOW::set_inf(void* x, int i) 
{ void* tmp = data[i]; 
  data[i] = x; 
  return tmp; 
}


void BASE_WINDOW::set_icon_pixrect(char* pr) { x_set_icon_pixmap(draw_win,pr); }

void BASE_WINDOW::set_icon_window(BASE_WINDOW& icon_win) 
{ x_set_icon_window(draw_win,icon_win.draw_win); }




void BASE_WINDOW::set_clip_rectangle(double x0, double y0, double x1, double y1)
{ int X0 = xpix(x0);
  int X1 = xpix(x1);
  int Y0 = ypix(y0);
  int Y1 = ypix(y1);
  if (Y1 < panel_height) Y1 = panel_height;
  x_set_clip_rectangle(draw_win,X0,Y1,X1-X0+1,Y0-Y1+1);
}


void BASE_WINDOW::set_clip_ellipse(double x0, double y0, double r1, double r2)
{ int R1 = real_to_pix(r1);
  int R2 = real_to_pix(r2);
  x_set_clip_ellipse(draw_win,xpix(x0),ypix(y0),R1,R2);
}

void BASE_WINDOW::set_clip_polygon(int n, double *xc, double *yc)
{
 int* x = new int[n+1];
 int* y = new int[n+1];
 int i;

 for(i=0;i<n;i++)
 { x[i] = xpix(xc[i]);
   y[i] = ypix(yc[i]);
  }

 x[n] = x[0];
 y[n] = y[0];

 x_set_clip_polygon(draw_win,n+1,x,y);

 delete[] x;
 delete[] y;
}


void BASE_WINDOW::reset_clipping() 
{ x_set_clip_rectangle(draw_win,0,panel_height,window_width,
                                          window_height-panel_height);
}


// status window

//static
void BASE_WINDOW::status_redraw(BASE_WINDOW* wp) 
{ BASE_WINDOW* win = (BASE_WINDOW*)wp->get_inf();
  char* str = win->status_str;
  if (str == 0) return;

  bool fstring = false;
  for(unsigned int i=0; i<strlen(str); i++)
       if (str[i] == '\\') fstring = true;


  double x0 = wp->xmin();
  double x1 = wp->xmax();
  double y0 = wp->ymin();
  double y1 = wp->ymax();

  wp->start_buffering();
  wp->clear();
  if (fstring)
     wp->text_box(x0,x1,y1,str);
  else
   { double d = (y1 - y0 - wp->text_height(str))/4;
     wp->draw_text(x0,y1-d,str);
    }
  wp->flush_buffer();
  wp->stop_buffering();
}



BASE_WINDOW* BASE_WINDOW::create_status_window(color col, int h)
{ if (status_win == 0)
  { status_win = new BASE_WINDOW(window_width+3,h);
    status_win->set_inf(this);
    status_win->set_redraw(status_redraw);
    status_win->set_bg_color(col);
   }
  return status_win;
}


BASE_WINDOW* BASE_WINDOW::open_status_window(color col, int h)
{ create_status_window(col,h);
  if (!status_win->is_open())
  { status_win->display(BASE_WINDOW::center,window_height-h,this);
    status_redraw(status_win);
   }
  return status_win;
}


void BASE_WINDOW::close_status_window()
{ if (status_win && status_win->is_open()) status_win->close(); }


void BASE_WINDOW::destroy_status_window()
{ if (status_win) 
  { delete status_win;
    status_win = 0;
   }
}

void BASE_WINDOW::set_status_string(const char* str)
{ if (status_str) delete[] status_str;
  if (str == 0)  str = "";
  status_str = new char[strlen(str) +1]; 
  strcpy(status_str,str);
  if (status_win) status_redraw(status_win);
}



void BASE_WINDOW::draw_text_cursor(double x, double y, int col)
{ if (col == invisible) return;
  int xc = xpix(x);
  int yc = ypix(y) + 1;
  int X[3];
  int Y[3];
  int d = x_text_height(draw_win,"H")/4 + 2;
  X[0] = xc-d+1;
  Y[0] = yc;
  X[1] = xc+d-1;
  Y[1] = yc;
  X[2] = xc;
  Y[2] = yc-d;
  set_color(col);
  x_fill_polygon(draw_win,3,X,Y);
  if (win_flush) flush();
}



void BASE_WINDOW::set_focus()
{ x_set_focus(draw_win); }


void BASE_WINDOW::string_edit(double x , double y, void* s)
{ 
  int string_h0 = string_h;
  string_h = string_h+6;

  int string_w0 = string_w;
  string_w = 3*string_w/4;

  double x1 = x - pix_to_real(string_w/2+2);
  double y1 = y - pix_to_real(string_h/2);
  double x2 = x + pix_to_real(string_w/2+2);
  double y2 = y + pix_to_real(string_h/2);
  char* pm =  get_pixrect(x1,y1,x2,y2);

  panel_item it = string_item0(x1,y2,s);

  act_str_item = it;

  draw_string_item(it,0);
  while (panel_text_edit(it))
  { assign_str(it->ref,it->data_str);
    it->data_str = access_str(it->ref);
   }
  put_pixrect(x1,y1,pm);
  del_pixrect(pm);
  act_str_item = 0;
  delete it;
  string_h = string_h0;
  string_w = string_w0;
}
  


void BASE_WINDOW::acknowledge(const char* s1, const char* s2)
{ BASE_WINDOW P(-1,-1,"s1");
  P.text_item("\\bf\\blue");
  P.text_item(s1);
  P.text_item("");
  if (strlen(s2) > 0)
  { P.text_item("\\tt");
    P.text_item(s2);
   }
  P.button("ok",0,panel_action_func(0));
  P.panel_open(BASE_WINDOW::center,BASE_WINDOW::center,this);
}




void BASE_WINDOW::reset_clip_mask()
{ x_clip_mask_polygon(draw_win,0,0,0,0); }



void BASE_WINDOW::clip_mask_polygon(int n, double* xc, double* yc, int c)
{ 
 int* x = new int[n];
 int* y = new int[n];

 for(int i=0;i<n;i++)
 { x[i] = xpix(xc[i]);
   y[i] = ypix(yc[i]);
  }

 x_clip_mask_polygon(draw_win,n,x,y,c);

 delete[] x;
 delete[] y;
}



void BASE_WINDOW::clip_mask_rectangle(double x0, double y0, double x1, 
                                                            double y1, int c)
{ int X0 = xpix(x0);
  int Y0 = ypix(y0);
  int X1 = xpix(x1);
  int Y1 = ypix(y1);
  if (X0 > X1) { int tmp = X0; X0 = X1; X1 = tmp; }
  if (Y0 > Y1) { int tmp = Y0; Y0 = Y1; Y1 = tmp; }
  int x[4], y[4];
  x[0] = X0; y[0] = Y0;
  x[1] = X1; y[1] = Y0;
  x[2] = X1; y[2] = Y1;
  x[3] = X0; y[3] = Y1;
  x_clip_mask_polygon(draw_win,4,x,y,c);
}


void BASE_WINDOW::clip_mask_ellipse(double x,double y,double r1,double r2,int c)
{ int R1 = real_to_pix(r1);
  int R2 = real_to_pix(r2);
  x_clip_mask_ellipse(draw_win,xpix(x),ypix(y),R1,R2,c); 
}


void BASE_WINDOW::set_drop_handler(void (*func)(void*,const char*,int,int))
{ x_set_drop_handler(draw_win,func); }



void BASE_WINDOW::set_d3_view_point(double px, double py, double pz, 
                                    double nx, double ny, double nz)
{ d3_view_pos_x = px;
  d3_view_pos_y = py;
  d3_view_pos_z = pz;
  d3_view_norm_x = nx;
  d3_view_norm_y = ny;
  d3_view_norm_z = nz;
}

void BASE_WINDOW::get_d3_view_point(double& px, double& py, double& pz, 
                                    double& nx, double& ny, double& nz)
{ px = d3_view_pos_x;
  py = d3_view_pos_y;
  pz = d3_view_pos_z;
  nx = d3_view_norm_x;
  ny = d3_view_norm_y;
  nz = d3_view_norm_z;
}


void BASE_WINDOW::project_d3_point(double& x, double& y, double& z)
{
  //q = normal_project(p)

  double px = d3_view_pos_x;
  double py = d3_view_pos_y;
  double pz = d3_view_pos_z;
  double nx = d3_view_norm_x;
  double ny = d3_view_norm_y;
  double nz = d3_view_norm_z;

  // v = p - point1
  double vx = x - px;
  double vy = y - py;
  double vz = z - pz;

  double W = nx*nx + ny*ny + nz*nz;
  double A = (nx*vx + ny*vy + nz*vz)/W;

  double qx = x-A*nx;
  double qy = y-A*ny;
  double qz = z-A*nz;

  x = qx;
  y = qy;
  z = qz;
}


int BASE_WINDOW::get_open_cmd(const char* suffix, char* buf, unsigned long sz)
{ return x_get_open_cmd(suffix,buf,sz); }


}
