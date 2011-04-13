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
// file          : include/CGAL/LEDA/base_window.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_BASE_WINDOW_H
#define CGAL_WINDOW_BASE_WINDOW_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/system.h>
#include <CGAL/LEDA/impl/x_window.h>
#include <CGAL/LEDA/color.h>


#define MOUSE_BUTTON(i) (-i)       

#define CLOSE_BUTTON  -8
#define NO_BUTTON    -16


namespace CGAL {


class __exportC BASE_WINDOW;

class __exportC window;

typedef double (*win_draw_func) (double);
typedef void   (*win_redraw_func0) ();
typedef void   (*win_redraw_func1) (BASE_WINDOW*);
typedef void   (*win_redraw_func2) (BASE_WINDOW*,double,double,double,double);
typedef void   (*mouse_action_func) (double,double);
typedef void   (*panel_action_func)(int);
typedef void   (*panel_str_action_func)(char*);

typedef void   (*coord_handler_func)(BASE_WINDOW*,double,double);

typedef void   (*win_delete_handler_func)(BASE_WINDOW*);

typedef int    (*event_handler_func)(BASE_WINDOW*,int,int,double,double,
                                                                 unsigned long);

// ----------------------------------------------------------------
// alternative for window handlers ...
// ----------------------------------------------------------------


class __exportC window_handler {
protected:
  double  d[4];
  int     i;
  window* win_ptr;
  char*   char_ptr;
  
public:
  window_handler() : win_ptr(NULL),char_ptr(NULL) {}

  virtual void operator()() { }
  
  double  get_double(int nr) const { return d[nr]; }
  int     get_int() const { return i; } 
  window* get_window_ptr() const { return win_ptr; }
  char*   get_char_ptr() const { return char_ptr; }
  
  void   set_params(window* wp, char* cp, int w = 0, double v1 =0, double v2=0, double v3=0, double v4=0)
  { win_ptr=wp; char_ptr=cp; i=w; d[0]=v1; d[1]=v2; d[2]=v3; d[3]=v4; }
};


// ----------------------------------------------------------------




//------------------------------------------------------------------------------
// panel items
//------------------------------------------------------------------------------

enum { Text_Item, 
       String_Item, 
       String_Menu_Item, 
       Int_Item, 
       Slider_Item, 
       Float_Item, 
       Button_Item, 
       Choice_Item,
       Choice_Mult_Item,
       Bool_Item,
       Color_Item,
       Bitmap_Choice_Item,
       Bitmap_Choice_Mult_Item,
       Newline_Item,
       Space_Item,
       Fill_Item };


class __exportC Panel_Item {

bool   enabled;
int    index;
int    kind;
char*  label_str;
char*  data_str;
char*  help_str;
void*  ref;
int    xcoord;
int    ycoord;
int    width;
int    height;
int    dat1;     // min (slider)
int    dat2;     // max (slider)  step (choice)
int    offset;   // choice item
int    argc;     // choice item, string_menu_item
char** argv;  

int   shortcut;

BASE_WINDOW*  menu_win;

Panel_Item*   menu_but;

panel_action_func action;
window_handler* action_ptr;

panel_str_action_func str_action;
window_handler* str_action_ptr;

Panel_Item(int,char*,void*,int,int);
~Panel_Item();

friend class __exportC BASE_WINDOW;

friend void str_menu_selection(int sel);

};

typedef Panel_Item* panel_item;

#define MAX_ITEM_NUM 256 


#if defined(GEOWIN_USE_NAMESPACE)
#define GeoWinTypeName CGAL::GeoWin
namespace CGAL {
  class _exportC GeoWin;
}
#else
#define GeoWinTypeName GeoWin
#endif


class __exportC BASE_WINDOW {

friend class __exportC window;
friend class __exportC panel;
friend class __exportC menu;
friend class __exportC Panel_Item;
friend class __exportC file_panel;

friend class __exportC GraphWin;
friend class _exportC GeoWinTypeName;


friend void str_menu_selection(int sel);

private: 

// LEDA windows cannot be copied, use reference parameters! 

BASE_WINDOW(const BASE_WINDOW&) { }
BASE_WINDOW& operator=(const BASE_WINDOW&) { return *this; }


// functions telling me how to access strings

virtual char* access_str(void*);
virtual void  assign_str(void*, const char*);

BASE_WINDOW* win_parent;   // window parent
BASE_WINDOW* p_win;        // window parent (panel/menu)
BASE_WINDOW* p_root;       // window root   (panel/menu)

int draw_win;

int panel_win;

char default_frame_label[128];

mouse_action_func       mouse_action;

win_redraw_func0        redraw0;
win_redraw_func1        redraw1;
window_handler          *redraw1_ptr;

win_redraw_func2        redraw2;
window_handler          *redraw2_ptr;

win_redraw_func1        timer_action;
window_handler          *timer_action_ptr;

coord_handler_func      coord_handler;
window_handler          *coord_handler_ptr;

event_handler_func      user_event_handler;

win_delete_handler_func win_delete_handler;
window_handler          *win_delete_ptr;

win_redraw_func2        clear_redraw;
window_handler          *clear_redraw_ptr;

int clear_on_resize;

int show_grid_cursor;

// pixel coordinates
int xdots;
int ydots;
int xorigin;
int yorigin;
int xoffset;
int yoffset;

// mouse button return values

int button_table[4];  
int shift_table[4];

// mouse data

int mouse_key;
int mouse_xpix;
int mouse_ypix;

unsigned long mouse_press_time;
unsigned long mouse_release_time;
unsigned long event_time;

double mouse_xreal;
double mouse_yreal;
double mouse_last_xreal;
double mouse_last_yreal;
double mouse_start_xreal;
double mouse_start_yreal;

// shift key state of last event
int shift_key_state;  


char* mesg_list[64];
int   mesg_count;

//window coordinates
double max_xcoord;
double min_xcoord;
double max_ycoord;
double min_ycoord;

double scaling;
double one_over_scaling;

int scaling_prec;


// d3 view
double d3_view_pos_x;
double d3_view_pos_y;
double d3_view_pos_z;
double d3_view_norm_x;
double d3_view_norm_y;
double d3_view_norm_z;


//window geometry
int window_xpos;
int window_ypos;
int window_width;
int window_height;

// hot area
int hotx1;
int hotx2;
int hoty1;
int hoty2;


double grid_mode;

grid_style  g_style;
point_style pt_style;

int node_width;
int win_flush;

int panel_enabled;


// status window

char*   status_str;
BASE_WINDOW* status_win;

// scroll bars

BASE_WINDOW* scroll_bar;


// user defined data slots

void* data[16];

GraphWin*       grawin_ptr;
GeoWinTypeName* geowin_ptr;


// static data members

static void REDRAW_FUNC(void* p,int,int,int,int);
static BASE_WINDOW*  read_window; 
static BASE_WINDOW*  active_window; 
static int win_count;

static panel_item   active_button;
static panel_item   last_active_button;
static BASE_WINDOW* active_button_win;
static panel_item help_last_item;

static void scroll_down_action(int);
static void scroll_up_action(int);
static void scroll_drag_action(int);

protected:

int buf_level;

static BASE_WINDOW*  call_window; 
static panel_item    call_item; 

panel_item  owner_item;


private:

void cursor(void);

static void open_help_window();
static void close_help_window();

static int event_handler(BASE_WINDOW*&,int blocking = 1);

#if defined(__unix__)
static void timer_handler(int);
#endif

static void mouse_default_action(double,double);
static void mouse_line_action(double x, double y) ;
static void mouse_ray_action(double x, double y) ;
static void mouse_segment_action(double x, double y) ;
static void mouse_rect_action(double x, double y);
static void mouse_circle_action(double x, double y);

void set_color(int c);
void draw_messages(void);
void clipping(int mode);

void draw_grid_points(double,double,double,double,double,double,double,bool);
void draw_grid_lines(double,double,double,double,double,double,double,bool);

void draw_copy_right();

public:

enum { min = 0x10000, center = 0x20000, max = 0x30000 };


void draw_grid(double,double,double,double,double,double,double);
void draw_grid(double,double,double,double);
void draw_grid();


void set_stipple(char* bist, int c=black);

void* get_inf(int=0) const;
void* set_inf(void*,int=0);

GraphWin*      get_graphwin() const { return grawin_ptr; }
GeoWinTypeName* get_geowin() const  { return geowin_ptr; }

int   state;
int   fg_color;
int   bg_color;
int   fill_color;

char* bg_pixrect;

double  bg_xoff;
double  bg_yoff;

static void do_not_open_display(int);

static int screen_width(void);
static int screen_height(void);

static char* root_pixrect(int,int,int,int);

static void quit_action(int);

// set and read mouse button return values

void set_buttons(int b0, int b1, int b2, int b3, 
                 int s0=0, int s1=-1, int s2=-2, int s3=-3);
void set_buttons(int* new_value, int* save_values=0);
void std_buttons(int* save_values=0);
void old_buttons(int* save_values=0);

int  mouse_button_val(int i) { return button_table[i]; }
int  mouse_shift_val(int i) { return shift_table[i]; }


// setting parameters
int set_precision(int prec);

win_delete_handler_func set_window_delete_handler(win_delete_handler_func f);
const window_handler*    set_window_delete_object(const window_handler& obj);

event_handler_func set_event_handler(event_handler_func f);
event_handler_func get_event_handler() const { return user_event_handler; }

coord_handler_func set_coord_handler(coord_handler_func f);
coord_handler_func get_coord_handler() const { return coord_handler; }


const window_handler* set_coord_object(const window_handler& obj);
const window_handler* get_coord_object_ptr() const { return coord_handler_ptr; }


void set_show_coordinates(int x) 
{ coord_handler = (x) ? default_coord_handler : 0; }

void set_show_cursor(int x)      { show_grid_cursor = x; }
void set_clear_on_resize(int x)  { clear_on_resize = x; }

void set_focus();

void set_redraw(win_redraw_func0);
void set_redraw(win_redraw_func1);
void set_redraw(win_redraw_func2);

void set_redraw(const window_handler&);
void set_redraw2(const window_handler&);

void set_bg_redraw(win_redraw_func2 func) { clear_redraw = func; }

void set_bg_redraw(const window_handler& obj) 
{ clear_redraw_ptr = & (window_handler&) obj; }

void set_drop_handler(void (*func)(void*,const char*,int,int));

int load_text_font(const char* fname);
int load_italic_font(const char* fname);
int load_bold_font(const char* fname);
int load_fixed_font(const char* fname);
int load_button_font(const char* fname);

int   set_bg_color(int);
int   set_fg_color(int);
int   set_fill_color(int);

char* set_bg_pixrect(char*);
char* set_bg_pixrect(char*,double,double);


void         set_border_width(int);
void         set_border_color(int);

void         set_text_font();
void         set_bold_font();
void         set_italic_font();
void         set_fixed_font();

int          set_font(const char* s);
int          set_grid_mode(int d);
grid_style   set_grid_style(grid_style s);
double       set_grid_dist(double d);
point_style  set_point_style(point_style s);
line_style   set_line_style(line_style s);
int          set_join_style(int s);
int          set_line_width(int w);
drawing_mode set_mode(drawing_mode m);
int          set_node_width(int w);
text_mode    set_text_mode(text_mode m);
int          set_cursor(int c);

void         set_frame_label(const char* s);
void         set_tmp_label(const char* s);
void         set_icon_label(const char* s);
void         reset_frame_label();
void         set_flush(int b) { win_flush = b; }


void         set_d3_view_point(double px, double py, double pz, 
                               double nx, double ny, double nz);

void         get_d3_view_point(double& px, double& py, double& pz, 
                               double& nx, double& ny, double& nz);

void         project_d3_point(double& x, double& y, double& z);


void set_icon_pixrect(char*);
void set_icon_window(BASE_WINDOW&);

int          get_line_width() const;
point_style  get_point_style() const { return pt_style; }
line_style   get_line_style() const;
drawing_mode get_mode() const;
text_mode    get_text_mode() const;
int          get_cursor() const;
grid_style   get_grid_style() const { return g_style; }

int    get_node_width() const { return node_width; }
int    get_grid_mode()  const { return (int)grid_mode; }
double get_grid_dist()  const { return grid_mode; }
int    get_fg_color()   const { return fg_color; }
int    get_fill_color() const { return fill_color; }
int    get_bg_color()   const { return bg_color; }
char*  get_bg_pixrect() const { return bg_pixrect; }
int    get_border_width() const;
int    get_border_color() const;

int get_panel_height() const { return panel_height; }
int get_panel_width()  const { return panel_width;  }

double  text_width(const char*) const;
double  text_height(const char*) const;


int    mono() const;

double xmin()  const { return min_xcoord; }
double xmax()  const { return max_xcoord; }
double ymin()  const { return min_ycoord; }
double ymax()  const { return max_ycoord; }
double scale() const { return scaling; }

// frame window
void frame_box(int& x0, int& y0, int& x1, int& y1) const;
int xpos() const;
int ypos() const;
int frame_width() const;
int frame_height() const;

// client window
int width()  const { return window_width; }
int height() const { return window_height; }


int xpix(double xcoord) const;
int ypix(double ycoord) const;

double xreal(int xpix) const;
double yreal(int ypix) const;

double pix_to_real(int p) const;
int    real_to_pix(double d) const;

void create(int w_width, int w_height, const char* label);
void display(int w_xpos, int w_ypos, BASE_WINDOW* parent = 0);

void resize(int xpos, int ypos, int width, int height);

void iconify();

void close();

void init0(double x0, double x1, double y0, int bgoff_restore);

void init(double x0, double x1, double y0, int g_mode=0, int erase=1);

void redraw();

void set_offset(double dx, double dy);
void configure();
void flush();


int panel_open(int,int,BASE_WINDOW*);
int menu_open(int,int,BASE_WINDOW*);

int  is_open() const;
int  is_closed() const;


 BASE_WINDOW(int width,int height, const char *frame_label = "");
 BASE_WINDOW(const char* frame_label);
 BASE_WINDOW();

virtual ~BASE_WINDOW();



// events & mouse input

void grab_mouse();
void ungrab_mouse();

void move_pointer(double,double);

int  read_event(int& val, double& x, double& y);
int  read_event(int& val, double& x, double& y, unsigned long& t);
int  read_event(int& val, double& x, double& y, unsigned long& t, int timeout);
int  get_event(int& val, double& x, double& y);

int read_mouse_action(mouse_action_func action, double xstart, 
                      double ystart, double& x, double& y);

int read_mouse(int kind, double xstart, double ystart, double& x, double& y);

int read_mouse(double& x, double& y) 
{ return read_mouse(0,0,0,x,y); }

int read_mouse() 
{ double x,y; 
  return read_mouse(x,y); 
 }

int get_mouse(double& x, double& y);

int get_mouse()
{ double x,y; 
  return get_mouse(x,y); 
 }



static void default_coord_handler(BASE_WINDOW*,double,double);

static int read_mouse(BASE_WINDOW*&,double&,double&);
static int get_mouse(BASE_WINDOW*&,double&,double&);

static int  read_event(BASE_WINDOW*&, int& val, double& x, double& y);
static int  get_event(BASE_WINDOW*&, int& val, double& x, double& y);

static void put_back_event();


static int get_open_cmd(const char*, char*, unsigned long);


unsigned long button_press_time();
unsigned long button_release_time();

// wait for panel button

int read_menu();


// get shift key info for last handled mouse button event

int shift_key_down();
int ctrl_key_down();
int alt_key_down();

// timer

void start_timer(int msec, win_redraw_func1 = 0);
void start_timer(int msec, const window_handler& obj);

void stop_timer();


int query_pix(double x, double y);

// drawing

void draw_pix(double x, double y, int col=black);
void draw_pixels(int, double*, double*, int col=black);
void draw_point(double x, double y, int col=black);
void draw_segment(double x1, double y1, double x2, double y2, int col=black);
void draw_segments(int, double*, double*, double*, double*, int col=black);
void draw_line(double x1, double y1, double x2, double y2, int col=black);
void draw_ray(double x1, double y1, double x2, double y2, int col=black);

void draw_node(double x0, double y0, int col=black);
void draw_filled_node(double x0, double y0, int col=black);
void draw_text_node(double x0, double y0, const char *s, int col=white);
void draw_int_node(double x0, double y0, int i, int col=white);

void draw_edge(double x1, double y1, double x2, double y2, int col=black);

void draw_circle(double x0, double y0, double r, int col=black);
void draw_filled_circle(double x0, double y0, double r, int col=black);
void draw_ellipse(double x0, double y0, double a, double b, int col=black);
void draw_filled_ellipse(double x0, double y0, double a, double b, int col=black);
void draw_arc(double x0, double y0, double r1, double r2, double start, double angle, int col=black);
void draw_filled_arc(double x0, double y0, double r1, double r2, double start, double angle, int col=black);

void adjust_polyline(int n, double *xcoord, double *ycoord);

void draw_polyline(int n, double *xcoord, double *ycoord, int col=black);
void draw_polygon(int n, double *xcoord, double *ycoord, int col=black);
void draw_filled_polygon(int n, double *xcoord, double *ycoord, int col=black);
void draw_rectangle(double x1, double y1, double x2, double y2, int col=black);
void draw_filled_rectangle(double x1, double y1, double x2, double y2, int col=black);

int construct_roundrect(int*&,int*&,double, double, double, double, double);

void draw_roundrect(double x1, double y1, double x2, double y2, double rndness,
                                                               int col=black);

void draw_roundbox(double x1, double y1, double x2, double y2, double rndness,
                                                               int col=black);


void plot_xy(double x0, double x1, win_draw_func f, int col=black);
void plot_yx(double y0, double y1, win_draw_func f, int col=black);

// text

void draw_text(double x, double y, const char* s, int col=black);
void draw_ctext(double x, double y, const char* s, int col=black);
void draw_ctext(const char* s, int col=black);

void draw_text_cursor(double x, double y, int col=black);

int split_text(const char* s, int& argc, char**& argv);

int format_text(int argc, char** argv,
                int xmin, int xmax, int ymin, int dy, int win);

double text_box(double x0, double x1, double y1, const char* s, bool draw=true);


// misc

void draw_frame(int,int,int,int);
void draw_frame();

void clear(double,double,double,double);
void clear(double,double,double,double,double xo,double yo);
void clear(double xo,double yo);


void get_bg_origin(double& x_orig, double& y_orig) 
{ x_orig = bg_xoff;
  y_orig = bg_yoff;
 }
 

void clear();
void clear(int col);

void message(const char *s);
void del_messages(void);

char* create_bitmap(int w, int h, unsigned char* data);
void  put_bitmap(double x, double y, char* bmap, int col=black);
void  del_bitmap(char*);

char* create_pixrect(int w, int h, unsigned char* data, int fc=black, 
                                                        int bc=white);
char* create_pixrect(const char**);

char* get_pixrect(double x1, double y1, double x2, double y2);
char* get_bitmap(double x1, double y1, double x2, double y2);

char* get_window_pixrect();
char* get_buffer_pixrect();

void  put_pixrect(double x, double y, char*);
void  center_pixrect(double x, double y, char*);
void  put_pixrect(double x, double y, char*,int x0,int y0,int w,int h);
void  put_pixrect(char*);
void  pixrect_to_ps(char* pmap, std::ostream& psout, bool full_color=true,
                                                bool animate=false);

void  pixrect_to_matrix(char* pmap, int* matrix);


void  screenshot_ps(const char* fname, bool full_color=true);
void  screenshot_wmf(const char* fname, bool full_color=true);

void  screenshot(std::ostream& ostr, bool full_color=true);
void  screenshot(const char* fname,  bool full_color=true);

void  pixrect_to_clipboard(char* pmap);
char* pixrect_from_clipboard();

void  open_metafile(const char* fname);
void  close_metafile();
void  metafile_to_clipboard();
void  load_metafile(double x0, double y0, double x1, double y1,
                                                     const char* fname);

void  del_pixrect(char*);
void  copy_rect(double x1, double y1, double x2, double y2, double x, double y);

int   get_width(char*) const;
int   get_height(char*) const;

void start_buffering();
void start_buffering(int,int);
void set_buffer(char*);

bool is_buffering();

void flush_buffer();
void flush_buffer(double,double,double,double);
void flush_buffer(double dx, double dy);
void flush_buffer(double dx, double dy, double,double,double,double);

void stop_buffering();
void stop_buffering(char*&);
void delete_buffer();

void set_clip_rectangle(double x0, double y0, double x1, double y1);
void set_clip_ellipse(double x0, double y0, double r1, double r2);
void set_clip_polygon(int n, double *xcoord, double *ycoord);
void reset_clipping();

void reset_clip_mask();

void clip_mask_polygon(int n,double*,double*,int);
void clip_mask_ellipse(double,double,double,double,int);
void clip_mask_rectangle(double,double,double,double,int);

void acknowledge(const char* s1, const char* s2="");

//------------------------------------------------------------------------------
// panel section
//------------------------------------------------------------------------------

private:

panel_item  Item[MAX_ITEM_NUM];

int but_count;
int orig_but_count;
int item_count;


int panel_bg_color;
int item_bg_color;
int d3_box_bg_color;
int shadow_color;
int press_color;
int bitmap_color0;
int bitmap_color1;
int disable_color;

int panel_width;
int panel_height;

int   th;               // text height
int   tw;               // text width
int   xoff;             // left and right boundary space
int   yoff;             // top and bottom boundary space
int   label_w;          // label width
int   sl_num_w;         // width of numerical fields in sliders
int   slider_w;         // slider item length
int   slider_h;         // height of slider items  
int   string_w;         // string item length
int   string_h;         // string item height
int   number_w;         // int/float item length
int   number_h;         // int/float item height
int   choice_h;         // height of choice items
int   color_h;          // height of color items
int   choice_w;         // choice field width
int   button_h;         // button height
int   button_w;         // button width
int   button_d;         // horizontal space between buttons
int   ytskip;           // height of text items
int   yskip;            // height of other items
int   buts_per_line;    // number of buttons in line
int   menu_button_style;     
int   menu_size;

int   center_button_label; // if true center text on buttons

panel_item act_str_item;
panel_item last_sel_button;
panel_item focus_button;


int panel_menu_mode;

int has_menu_bar;
int menu_bar_height;

panel_item new_panel_item(int,const char*,void*,int,const char*);

void draw_label(panel_item);
void draw_box_with_shadow(int,int,int,int,int,int);
void draw_d3_box(int x1,int y1,int x2,int y2, int pressed, int enabled=1);

void draw_string_item(panel_item i, const char* s=0);
void activate_string_item(panel_item,int);
bool panel_text_edit(panel_item i);

void put_text_item(int x, int y, const char* s, int t_len);
int  draw_text_item(panel_item i, int);
void draw_choice_item(panel_item i,int x);
void draw_bool_item(panel_item i);
void draw_color_button(int xt, int yt, int w, int col, int pressed, int enabled=1);
void draw_color_item(panel_item i);
void change_color_item(panel_item i, int j);
void draw_slider_item(panel_item i,int x);
void draw_bitmap_choice_item(panel_item i);
void change_bitmap_choice_item(panel_item i, int j);
void draw_button(panel_item i,int pressed);
void draw_menu_button(panel_item i,int pressed);
void draw_menu_item(panel_item i,int pressed);
void draw_down_menu_button(int x, int y,int pressed, int enabled);
void draw_up_menu_button(int x, int y,int pressed, int enabled);
void draw_right_menu_button(int x, int y,int pressed, int enabled);
void draw_left_menu_button(int x, int y,int pressed, int enabled);

void draw_separator(panel_item i);

void item_error();

void open_sub_panel(panel_item,bool=false);
void close_sub_panel(BASE_WINDOW*);
void scroll(int);

int  panel_event_handler(int w, int k, int b, int x, int y, unsigned long t);

void panel_init();
void place_panel_items();
void draw_panel_items();

public:

typedef panel_item item;


void redraw_panel(panel_item);

void make_menu_bar(int kind=1);

void buttons_per_line(int n) 
{ buts_per_line = n; }

void set_menu_style(int n)   
{ menu_button_style = n; }

void set_bitmap_colors(int c0, int c1)
{ bitmap_color0 = c0; bitmap_color1 = c1; }

void set_panel_bg_color(int c) { panel_bg_color = c; }

void set_item_width(int);
void set_item_height(int);
void set_item_space(int);
void set_button_space(int);


panel_item text_item(const char*);

void set_text(panel_item, const char*);

panel_item string_item0(double x, double y, void* s);

panel_item string_item(const char*, void* x, panel_str_action_func,
                                             const char* hlp);

void set_menu(panel_item,int,const char**, int sz);


panel_item string_menu_item(const char*, void* x, const char* menu_label, 
                            int argc, const char** argv, int size,
                            panel_str_action_func, const char* hlp);

panel_item int_item(const char*, int* x,const char* hlp);

panel_item slider_item(const char*,int* x,int low,int high, panel_action_func,
                                                            const char* hlp);

panel_item float_item(const char*, double* x, const char* hlp);

panel_item bool_item(const char* text, bool* x, panel_action_func,
                                                const char* hlp);


 
panel_item color_item(const char* text, int* x, panel_action_func, 
                                                const char* hlp);

panel_item pstyle_item(const char* text, point_style* x, panel_action_func, 
                                                         const char* hlp);

panel_item lstyle_item(const char* text, line_style* x, panel_action_func, 
                                                        const char* hlp);

panel_item lwidth_item(const char* text, int* x, panel_action_func, 
                                                 const char* hlp);



//choice items
panel_item choice_item(const char* text, int* x, int argc, 
                       const char** argv, int step, int offset, 
                       panel_action_func, const char* hlp);

panel_item choice_mult_item(const char* text,int* x,int argc,
                            const char** argv, panel_action_func, 
                            const char* hlp);

panel_item bitmap_choice_item(const char* label, int *x, int argc,
                              int width, int height, unsigned char **argv, 
                              panel_action_func, const char* hlp);

panel_item bitmap_choice_mult_item(const char* label, int* x, int argc,
                                  int width, int height, unsigned char **argv,
                                  panel_action_func, const char* hlp);



//buttons
int button(const char*, int val, panel_action_func F, const char* hlp=0);
int button(const char*, int val, BASE_WINDOW* p,const char* hlp=0);

// bitmap
int button(int w, int h, unsigned char* bm, const char* s, int val, 
                                                           panel_action_func F,
                                                           const char* hlp=0);

int button(int w, int h, unsigned char* bm, const char* s, int val, 
                                                           BASE_WINDOW* p,
                                                           const char* hlp=0);

// pixmaps
int button(char* pmap1, char* pmap2, const char* s, int val,panel_action_func F,
                                                            const char* hlp=0);

int button(char* pmap1, char* pmap2, const char* s, int val,BASE_WINDOW* p,
                                                            const char* hlp=0);


int menu_button(const char*, int val, BASE_WINDOW* p, const char* hlp=0);

void set_focus_button(int but);

void new_line();
void fill_line();

void separator();

void hspace(int);
void vspace(int);

// string edit

void string_edit(double,double,void*);


static BASE_WINDOW* get_call_window();
static panel_item   get_call_item();
static int          get_call_button(); 

// manipulation of panel items

panel_item get_item(const char* s);
panel_item get_button_item(int b);
int        get_button(const char* s);

BASE_WINDOW* get_window(panel_item it);
BASE_WINDOW* get_window(int but);

BASE_WINDOW*      set_window(int but, BASE_WINDOW* p);

panel_action_func set_action(int but, panel_action_func F);
const window_handler* set_action(int but, const window_handler& obj);


int               set_value(int but, int val);


void set_item_menu_func(panel_item p, panel_action_func F) 
{ p->menu_but->action = F; }

void set_item_func(panel_item p, panel_action_func F) 
{ p->action = F; }


void set_item_menu_object(panel_item p, const window_handler& obj) 
{ p->menu_but->action_ptr = & (window_handler&)obj; }

void set_item_object(panel_item p, const window_handler& obj) 
{ p->action_ptr = & (window_handler&) obj; }

void set_item_object_str(panel_item p, const window_handler& obj) 
{ p->str_action_ptr = & (window_handler&) obj; }


panel_item  get_owner_item() { return owner_item; }


void enable_panel();
void disable_panel(bool disable_items=true);

void enable_item(panel_item);
void disable_item(panel_item);

bool is_enabled(panel_item it) { return it->enabled; }

void        set_item_label(panel_item,const char*);
const char* get_item_label(panel_item);

void        set_button_repeat(int but, int msec);
void        set_button_label(int but, const char* str);
void        set_button_help_str(int but, const char* str);
void        set_button_pixrects(int but,char*pmap0,char*pmap1);

const char* get_button_label(int but);


int          get_kind(panel_item it) { return it->kind; }
const char*  get_label(panel_item it) { return it->label_str; }

panel_item  active_string_item() { return act_str_item; }


int   get_item_xpos(panel_item it) { return it->xcoord; }
int   get_item_ypos(panel_item it) { return it->ycoord; }

void   set_item_xpos(panel_item it, int x) { it->xcoord = x; }
void   set_item_ypos(panel_item it, int y) { it->ycoord = y; }



// iteration
panel_item first_item() const;
panel_item next_item(panel_item) const;




// scroll bar

void open_scrollbar(void (*scroll_up)(int),
                    void (*scroll_down)(int),
                    void (*scroll_drag)(int), double sz , double pos=0);
void close_scrollbar();

void set_scrollbar_pos(double pos);


static int text_color(int);



// status window

static void  status_redraw(BASE_WINDOW*);
BASE_WINDOW* create_status_window(color bc=ivory, int h=17);
BASE_WINDOW* open_status_window(color bc=ivory, int h=17);
void         close_status_window();
void         destroy_status_window();
void         set_status_string(const char*);
const char*  get_status_string() { return status_str ? status_str : ""; }
BASE_WINDOW* get_status_window() { return status_win; }

};

}

#endif
