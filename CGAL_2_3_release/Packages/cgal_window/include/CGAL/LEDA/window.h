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
// file          : include/CGAL/LEDA/window.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================



#ifndef CGAL_WINDOW_H
#define CGAL_WINDOW_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/window_point.h>
#include <CGAL/LEDA/base_window.h>
#include <list>
#include <string>

/*{\Moptions usesubscripts=yes}*/


namespace CGAL {

using std::string;


typedef void (*show_coord_handler_func)(window*,double,double);


/*{\Manpage {window} {} {Windows}}*/

class __exportC window : public BASE_WINDOW {

//friend class __exportC GraphWin;
//friend class _exportC  GeoWinTypeName;


/*{\Mdefinition
The data type $window$ provides an interface for graphical input and output 
of basic two-dimensional geometric objects.
Application programs using this data type have to be linked with {\em libCGALWin.a} 
and (on UNIX systems) with the X11 base library {\em libX11.a}
(cf.~section~1.6):

CC $prog.c$ -lCGALWin -lX11 -lm

An instance $W$ of type $window$ is an iso-oriented rectangular window in the 
two-dimensional plane. The default representation of $W$ on the screen is a 
square of maximal possible edge length positioned in the upper right corner 
of the display. 

In general, a window consists of two rectangular sections,
a {\em panel section} in the upper part and 
a {\em drawing section} in the rest of the window. 
The panel section  contains {\em panel items} such as sliders, choice
fields, string items and buttons. They have to be created before
the window is opened by special panel operations described in section
\ref{panel-operations}.

The drawing section can be used for the output of geometric objects such
as points, lines, segments, arrows, circles, polygons, graphs, \dots and 
for the input of all these objects using the mouse input device.  
All drawing and input operations in the drawing section use a
coordinate system that is defined by three parameters of type $double$: 
$xmin$, the minimal x-coordinate, $xmax$, the maximal x-coordinate, 
and $ymin$, the minimal y-coordinate. The two parameters $xmin$ and $xmax$
define the scaling factor $scaling$  as $w/(xmax-xmin)$, where $w$ is the width 
of the drawing section in pixels. The maximal y-coordinate $ymax$ of the
drawing section is equal to $ymin + h \cdot  scaling$ and depends on the 
actual shape of the window on the screen. Here, $h$ is the height of the
drawing section in pixels.


A list of all window parameters:

\begin{enumerate}

\item
The {\em foreground color} parameter (default $black$) defines 
the default color to be used in all drawing operations.
There are 18 predefined colors (enumeration in \<CGAL/LEDA/incl/impl/x\_window.h\>):
$black$, $white$, $red$, $green$, $blue$, $yellow$, $violet$, $orange$,
$cyan$, $brown$, $pink$, $green2$, $blue2$, $grey1$, $grey2$, $grey3$
$ivory$, and $invisible$. 
Note that all drawing operations have an optional color argument that can
be used to override the default foreground color.
On monochrome systems all colors different from $white$ are displayed as 
$black$. 
The color $invisible$ can be used for invisible (transparent) objects.


\item
The {\em background color} parameter (default $white$) defines
the default background color (e.g. used by $W$.clear()).

\item
The {\em text font} parameter defines the name of the font to be used
in all text drawing operations.

\item
Minimal and maximal coordinates of the drawing area
$xmin$ (default 0), $xmax$ (default 100), $ymin$ (default 0).

\item
The {\em grid width} parameter (default 0) defines the width of the
grid that is used in the drawing area. A grid width of 0 indicates
that no grid is to be used.

\item
The {\em frame label} parameter defines a string to be displayed
in the frame of the window.

\item
The {\em show coordinates} flag  (default $true$) determines whether
the current coordinates of the mouse cursor in the drawing section
are displayed in the upper right corner.

\item
The {\em flush output} flag  (default $true$) determines whether
the graphics output stream is flushed after each draw action.

\item
The {\em line width} parameter (default value 1 pixel) defines the width of all
kinds of lines (segments, arrows, edges, circles, polygons).

\item
The {\em line style} parameter defines the style of lines. Possible line
styles are |solid| (default), |dashed|, and |dotted|.

\item
The {\em point style} parameter defines the style points are drawn by 
the |draw_point| operation. Possible point styles are |pixel_point|, 
|cross_point| (default), |plus_point|, |circle_point|, |disc_point|, 
|rect_point|, and |box_point|.

\item
The {\em node width} parameter (default value 8 pixels) defines the diameter
of nodes created by the draw\_node and draw\_filled\_node operations.

\item
The {\em text mode} parameter defines how text is inserted into the window.
Possible values are |transparent| (default)  and |opaque|.

\item
The {\em drawing mode} parameter defines the logical operation that is used
for setting pixels in all drawing operations. Possible values are
|src_mode| (default) and |xor_mode|. In |src_mode| pixels are set to
the respective color value, in |xor_mode| the value is bitwise added to the
current pixel value.

\item
The {\em redraw function} parameter is a pointer to a function of type 
void ($*$F)(window$*$). It is called with a pointer to the corresponding window
as argument to redraw (parts of) the window whenever a redrawing 
is necessary, e.g., if the shape of the window is changed or previously
hidden parts of it become visible. 

\item
The {\em window delete handler} parameter is a pointer to a function of type 
void ($*$F)(window$*$). It is called with a pointer to the corresponding window
as argument when the window is to be closed by the window manager (e.g. by
pressing the $\times$-button on Windows-NT systems). The default
window delete handler closes the window and terminates the program.

\item
The {\em buttons per line} parameter (default $\infty$) defines
the maximal number of buttons in one line of the panel section.


\item
The {\em precision} parameter (default $16$) defines the precision
that is used for representing window coordinates, more precisely,
all $x$ and $y$ coordinates generated by window input operations 
are doubles whose mantissa are truncated after $precision-1$ bits
after the binary point.


\end{enumerate} 

In addition to call-back (handler) functions windows now also
support the usage of function objects. Function object classes have
to be derived from the |window_handler| base class.
\begin{verbatim}
class window_handler {
  ...
  virtual void operator()() { }
  
  // parameter access functions
  double  get_double(int nr) const;
  int     get_int() const; 
  window* get_window_ptr() const;
  char*   get_char_ptr() const;
};
\end{verbatim}

Derived classes have to implement the handling function in the definition
of the |operator()| method. The different |get_| methods can be
called to retrieve parameters.

If both, a handler function and an object for the same action is supplied
the object has higher priority.\\
}*/


private: 

/*
#if defined(__KCC)
// instantiate lists for this compiler
std::list<CGAL::window_point> __dummy_list__;
std::list<std::string> __dummy_list_2__;
std::list<CGAL::window_point>::const_iterator __dummy_iter__;
#endif
*/

// windows cannot be copied, use reference parameters! 

window(const window&) { }
window& operator=(const window&) { return *this; }

char* access_str(GenPtr p) { return (char*)(((string*)p)->c_str()); }
void  assign_str(GenPtr p, const char* s) { *((string*)p) = s; }

window_point* point_buffer;

public:

//enum { min = 0x10000, center = 0x20000, max = 0x30000 };


static color fgcol;
static color bgcol;


bool normalize_rat;


/*{\Mcreation W }*/ 

window();
/*{\Mcreate creates a squared window with maximal possible edge length
            (minimum of width and height of the display).}*/


window(const char* label);
/*{\Mcreate creates a maximal squared window with frame label $label$.}*/

window(int w, int h);
/*{\Mcreate creates a window $W$ of physical size $w$ pixels 
            $\times$ $h$ pixels .}*/


window(int w, int h, const char* label);
/*{\Mcreate creates a window $W$ of physical size $w$ pixels 
            $\times$ $h$ pixels and frame label $label$.}*/



/*{\Mtext

\medskip
All four variants initialize the coordinates of $W$ to $xmin = 0$,
$xmax = 100$ and $ymin = 0$. The $init$ operation (see  below) can later
be used to change the window coordinates and scaling. Please note, that a 
window is not displayed before the function $display$ is called for it.}*/


~window() {}


/*{\Moperations 3.1 4.4}*/

/*{\Mtext
{\bf 3.1 Initialization} }*/

void init(double x0, double x1, double y0) { BASE_WINDOW::init(x0,x1,y0); }
/*{\Mopl     sets $xmin$ to $x_0$, $xmax$ to $x_1$, and $ymin$ to $y_0$, 
             the scaling factor $scaling$ to $w / (xmax-xmin)$, and 
             $ymax$ to $ymin + h \cdot scaling$. Here $w$ and $h$ are the
             width and height of the drawing section in pixels. }*/

double set_grid_dist(double d) 
{ return BASE_WINDOW::set_grid_dist(d); }
/*{\Mopl     sets the grid distance of $W$ to $d$. }*/

grid_style set_grid_style(grid_style s) 
{ return BASE_WINDOW::set_grid_style(s); }
/*{\Mopl     sets the grid style of $W$ to $s$. }*/

int set_grid_mode(int d) 
{ return BASE_WINDOW::set_grid_mode(d); }
/*{\Mopl     sets the grid distance of $W$ to $d$ pixels. }*/

bool set_normalize_rational_objects(bool b)
{ bool prev = normalize_rat; 
  normalize_rat = b;
  return prev;
}

int set_precision(int prec) 
{ return BASE_WINDOW::set_precision(prec); }
/*{\Mopl     sets the precision of $W$ to $prec$. }*/


void init(double x0, double x1, double y0, int d, bool erase=true) 
{ BASE_WINDOW::init(x0,x1,y0,d,erase); }
/*{\Mopl     same as $W$.init($x_0,x_1,y_0$) followed by 
             $W$.set\_grid\_mode($d$). If the optional flag |erase| is
             set to |false| the window will not be erased. }*/


void display() { BASE_WINDOW::display(window::max,window::min,0); }
/*{\Mopl      opens $W$ and displays it with its right upper corner in the 
              upper right corner of the 
              screen. Note that $W$.display() has to be called before
              all drawing operations and that all operations adding panel 
              items to $W$ (cf. \ref{panel-operations}) have to be called 
              before the first call of $W$.display(). }*/

void display(int x, int y) { BASE_WINDOW::display(x,y,0); }
/*{\Mopl      opens $W$ and displays it with its left upper corner at position
              $(x,y)$. Special 
              values for $x$ and $y$ are  $window::min$, $window::center$, 
              and $window::max$ for positioning $W$ at the minimal or maximal 
              $x$ or $y$ coordinate or centering it in the $x$ or $y$ 
              dimension. }*/

/*{\Moptions nextwarning=no}*/
void display(window& w, int x, int y) { BASE_WINDOW::display(x,y,&w); }
void display(int x, int y, window& W0) { BASE_WINDOW::display(x,y,&W0); }
/*{\Mopl      opens $W$ and displays it with its left upper corner at position 
              $(x,y)$ relative to the upper left corner of  window $W_0$.   }*/ 
/*{\Mtext
$W.open$\dots can be used as a synonym for $W.display$\dots\ Note,
that the $open$ operation for panels (cf. \ref{Panels}) is defined 
slightly different.}*/

void open()                        { display(); }
void open(int x, int y)            { display(x,y); }
void open(window& w, int x, int y) { display(w,x,y); }
void open(int x, int y, window& w) { display(x,y,w); }

void close() { BASE_WINDOW::close(); }
/*{\Mopl      closes $W$ by removing it from the display. }*/


void clear() { BASE_WINDOW::clear(); }
/*{\Mopl      clears $W$ using the current background color or pixmap,
              i.e., if $W$ has a background pixmap defined it is tiled with
              $P$ such that the upper left corner is the tiling origin.
              Otherwise, it is filled with background color of $W$. }*/

void clear(double x0,double y0,double x1,double y1)
{ BASE_WINDOW::clear(x0,y0,x1,y1); }
/*{\Mopl      only clears the rectangular area $(x0,y0,x1,y1)$ of window $W$
              using the current background color or pixmap.}*/

void clear(double x0,double y0,double x1,double y1, double xorig, double yorig)
{ BASE_WINDOW::clear(x0,y0,x1,y1,xorig,yorig); }


void clear(color c) { BASE_WINDOW::clear(c); }
/*{\Mopl      clears $W$ with color $c$ and sets the background color of $W$
              to $c$. }*/

void clear(double xorig, double yorig) { BASE_WINDOW::clear(xorig,yorig); }
/*{\Mopl      clears $W$. If a background pixmap is defined the point
              $(xorig,yorig)$ is used as the origin of tiling. }*/

void redraw() { BASE_WINDOW::redraw(); }
/*{\Mop      repaints the drawing area if $W$ has a redraw function.}*/ 



/*{\Mtext
{\bf 3.2 Setting parameters}\label{parameters} }*/

color set_fg_color(color c) { return BASE_WINDOW::set_fg_color(c); }

color set_color(color c) { return BASE_WINDOW::set_fg_color(c); }
/*{\Mopl     sets the foreground color parameter to $c$ and
	     returns its previous value.}*/


color set_fill_color(color c) { return BASE_WINDOW::set_fill_color(c); }
/*{\Mopl     sets the fill color parameter (used by |<<| operators) 
             to $c$ and returns its previous value.}*/



color set_bg_color(color c) { return BASE_WINDOW::set_bg_color(c); }
/*{\Mopl     sets the background color parameter to $c$ and
	     returns its previous value.}*/

char* set_bg_pixmap(char* pr) { return BASE_WINDOW::set_bg_pixrect(pr); }
/*{\Mopl     sets the background pixmap to $pr$ and
	     returns its previous value.}*/


int set_line_width(int pix){return BASE_WINDOW::set_line_width(pix);}
/*{\Mopl     sets the line width parameter to $pix$ pixels and
	     returns its previous value.}*/

line_style set_line_style(line_style s)
{return BASE_WINDOW::set_line_style(s);}
/*{\Mopl     sets the line style parameter to $s$ and returns its
	     previous value.}*/

int set_node_width(int pix) {return BASE_WINDOW::set_node_width(pix);}
/*{\Mopl     sets the node width parameter to $pix$ pixels and
	     returns its previous value.}*/

text_mode set_text_mode(text_mode m)
{return BASE_WINDOW::set_text_mode(m);}
/*{\Mopl     sets the text mode parameter to $m$ and returns
	     its previous value.}*/

drawing_mode set_mode(drawing_mode m)
{return BASE_WINDOW::set_mode(m);}
/*{\Mopl     sets the drawing mode parameter to $m$ and returns 
	     its previous value.}*/

int          set_cursor(int cursor_id = -1)
{return BASE_WINDOW::set_cursor(cursor_id);}
/*{\Mopl     sets the mouse cursor of $W$ to |cursor_id|. Here |cursor_id| 
             must be one of the constants  predefined in |<X11/cursorfont.h>|
             or $-1$ for the system default cursor. Returns the previous 
             cursor. }*/

void set_show_coordinates(bool b) { BASE_WINDOW::set_show_coordinates(b); }
/*{\Mopl     sets the show coordinates flag to $b$. }*/

void set_frame_label(string s){ BASE_WINDOW::set_frame_label(s.c_str()); }
/*{\Mopl     makes $s$ the window frame label. }*/

void set_icon_label(string s){ BASE_WINDOW::set_icon_label(s.c_str()); }
/*{\Mopl     makes $s$ the window icon label. }*/

void reset_frame_label() { BASE_WINDOW::reset_frame_label(); }
/*{\Mop      restores the standard frame label. }*/


void set_window_delete_handler(void (*F)(window*)) 
{ BASE_WINDOW::set_window_delete_handler((win_delete_handler_func)F);}
/*{\Mopl     sets the window delete handler function parameter to $F$.}*/

void set_window_delete_object(const window_handler& obj)
{ BASE_WINDOW::set_window_delete_object(obj); }
/*{\Mopl     sets the window delete object parameter to $obj$.}*/

void set_event_handler(int (*F)(window*,int,int,double,double,unsigned long)) 
{ BASE_WINDOW::set_event_handler((event_handler_func)F);}

show_coord_handler_func get_show_coord_handler() const
{ return (show_coord_handler_func)BASE_WINDOW::get_coord_handler(); }

void set_show_coord_handler(void (*F)(window*,double,double)) 
{ BASE_WINDOW::set_coord_handler((coord_handler_func)F);}
/*{\Mopl     sets the show coordinate handler function parameter to $F$.}*/

void set_show_coord_object(const window_handler& obj)
{ BASE_WINDOW::set_coord_object(obj); }
/*{\Mopl     sets the show coordinate object parameter to $obj$.}*/

void set_redraw(void (*F)(window*)) 
{ BASE_WINDOW::set_redraw((win_redraw_func1)F);}
/*{\Mopl     sets the redraw function parameter to $F$.}*/

void set_redraw(const window_handler& obj)
{ BASE_WINDOW::set_redraw(obj); }
/*{\Mopl     sets the redraw object parameter to $obj$.}*/

void set_redraw(void (*F)(window*,double,double,double,double)=0) 
{ BASE_WINDOW::set_redraw((win_redraw_func2)F);}
/*{\Mopl     sets the redraw function parameter to $F$.}*/

void set_redraw2(const window_handler& obj)
{ BASE_WINDOW::set_redraw2(obj); }
/*{\Mopl     sets the redraw object parameter to $obj$.}*/


void set_redraw(void (*F)()) { BASE_WINDOW::set_redraw((win_redraw_func0)F);}
/*{\Xopl     for backward compatibility. }*/

void set_bg_redraw(void (*F)(window*,double,double,double,double)=0) 
{ BASE_WINDOW::set_bg_redraw((win_redraw_func2)F);}
/*{\Mopl     sets the background redraw function parameter to $F$.}*/

void set_bg_redraw(const window_handler& obj) 
{ BASE_WINDOW::set_bg_redraw(obj); }
/*{\Mopl     sets the background redraw object parameter to $obj$.}*/


void start_timer(int msec, void (*F)(window*)) 
{ BASE_WINDOW::start_timer(msec,(win_redraw_func1)F);}
/*{\Mopl   starts a timer that runs $F$ every |msec| milliseconds
           with a pointer to |\Mvar|.}*/
	   
void start_timer(int msec, const window_handler& obj)
{ BASE_WINDOW::start_timer(msec,obj); }
/*{\Mopl   starts a timer that runs the $operator()$ of $obj$ every |msec| milliseconds.}*/


void stop_timer() { BASE_WINDOW::stop_timer(); }
/*{\Mopl   stops the timer.}*/


void set_flush(bool b)     { BASE_WINDOW::set_flush(b); }
/*{\Mopl     sets the flush parameter to $b$.}*/

void set_icon_pixrect(char* pr)  { BASE_WINDOW::set_icon_pixrect(pr); }
/*{\Mopl     makes |pr| the new icon of |\Mvar|. }*/


void* set_client_data(void* p, int i=0) { return BASE_WINDOW::set_inf(p,i); }
/*{\Mopl    sets the $i$-th client data pointer of |\Mvar| to $p$ and 
            returns its previous value. \precond $i < 16$. }*/ 



bool load_text_font(string fn)    { return BASE_WINDOW::load_text_font(fn.c_str())!=0; }
/*{\Xopl      loads $X11$ font $fn$ as text font. Returns true on success and 
              false if the font is not available.}*/

bool load_bold_font(string fn)    { return BASE_WINDOW::load_bold_font(fn.c_str())!=0; }
/*{\Xopl      loads $X11$ font $fn$ as bold font. Returns true on success and 
              false if the font is not available.}*/

bool load_fixed_font(string fn) { return BASE_WINDOW::load_fixed_font(fn.c_str())!=0;  }
/*{\Xopl      load $X11$ font $fn$ as fixed font. Returns true on success 
              and false if the font is not available.}*/

bool load_button_font(string fn) { return BASE_WINDOW::load_button_font(fn.c_str())!=0;}
/*{\Xopl      load $X11$ font $fn$ as button font. Returns true on success 
              and false if the font is not available.}*/


/*{\Mtext
{\bf 3.3 Reading parameters} }*/

bool get_normalize_rational_objects() const { return normalize_rat; }


int get_line_width() const {return BASE_WINDOW::get_line_width();}
/*{\Mop      returns the current line width.}*/ 

line_style get_line_style() const {return BASE_WINDOW::get_line_style();}
/*{\Mop      returns the current line style.}*/ 

int get_node_width() const {return BASE_WINDOW::get_node_width();}
/*{\Mop      returns the current node width.}*/ 

text_mode get_text_mode() const {return BASE_WINDOW::get_text_mode();}
/*{\Mop      returns the current text mode.}*/ 

drawing_mode get_mode() const {return BASE_WINDOW::get_mode();}
/*{\Mop      returns the current drawing mode.}*/ 

int          get_cursor() const {return BASE_WINDOW::get_cursor();}
/*{\Mop      returns the id of the current cursor, i.e, one of the constants
             predefined in $<X11/cursorfont.h>$ or $-1$ for the default
             cursor. }*/ 

double xmin() const {return BASE_WINDOW::xmin();}
/*{\Mop      returns the minimal x-coordinate of the drawing area of $W$.}*/ 

double ymin() const {return BASE_WINDOW::ymin();}
/*{\Mop      returns the minimal y-coordinate of the drawing area of $W$.}*/ 

double xmax() const {return BASE_WINDOW::xmax();}
/*{\Mop      returns the maximal x-coordinate of the drawing area of $W$.}*/ 

double ymax() const {return BASE_WINDOW::ymax();}
/*{\Mop      returns the maximal y-coordinate of the drawing area of $W$.}*/ 

double scale() const {return BASE_WINDOW::scale();}
/*{\Mop      returns the scaling factor of the drawing area of $W$, i.e. the 
             number of pixels of a unit length line segment.}*/

double       get_grid_dist()  const { return BASE_WINDOW::get_grid_dist(); }
/*{\Mop      returns the width of the current grid
             (zero if no grid is used). }*/

grid_style   get_grid_style() const 
{ return (grid_style) BASE_WINDOW::get_grid_style(); }
/*{\Mop      returns the current grid style. }*/


int          get_grid_mode()  const { return BASE_WINDOW::get_grid_mode(); }
/*{\Mop      returns the width of the current grid in pixels (zero if
             no grid is used). }*/


void* get_client_data(int i=0) const { return BASE_WINDOW::get_inf(i); }
/*{\Mop   returns the $i$-th client data pointer of of |\Mvar|. 
          \precond $i < 16$. }*/ 


//GraphWin* get_graphwin() const { return BASE_WINDOW::get_graphwin(); }
/*{\Xop   returns a pointer to the |GraphWin| (see \ref{Graph Windows})
          that uses |\Mvar| as its display window or |NULL| if |\Mvar| 
          is not used by any |GraphWin|. }*/ 

//GeoWinTypeName*   get_geowin() const { return BASE_WINDOW::get_geowin(); }
/*{\Xop   returns a pointer to the |GeoWin| (see \ref{Geo Windows})
          that uses |\Mvar| as its display window or |NULL| if |\Mvar| 
          is not used by any |GeoWin|. }*/ 



int width() const { return BASE_WINDOW::width(); }
/*{\Mop      returns the width of $W$ in pixels.}*/ 

int height() const { return BASE_WINDOW::height(); }
/*{\Mop      returns the height of $W$ in pixels.}*/ 

int xpos() const { return BASE_WINDOW::xpos(); }
/*{\Mop      returns the x-coordinate of the upper left corner of 
             the frame of $W$. }*/ 

int ypos() const { return BASE_WINDOW::ypos(); }
/*{\Mop      returns the y-coordinate of the upper left corner of 
             the frame of $W$. }*/ 


operator void*() const { return (state==0) ? 0 : (void*)this; }

int  get_state() const { return state; }
/*{\Mop      returns the state of  $W$. }*/ 

void set_state(int stat) { state = stat; }
/*{\Mop      sets the state of $W$ to |stat|. }*/ 


bool contains(double x, double y) const;
/*{\Mop returns true if $(x,y)$ lies in the drawing area. }*/


/*{\Mtext
{\bf 3.4 Drawing Operations}

All drawing operations have an optional color argument at the end of the
parameter list. If this argument is omitted the current foreground color 
(cf. section \ref{parameters}) of $W$ is used.
}*/



// points

/*{\Mtext
{\bf 3.4.1 Drawing points} 
\setopdims{1.2cm}{4.4cm}
}*/

void draw_point(double x, double y, color c=window::fgcol);
/*{\Mopl     draws the point $(x,y)$ . }*/


/*{\Moptions nextwarning=no}*/
void draw_pix(double x, double y, color c=window::fgcol);
void draw_pixel(double x, double y, color c=window::fgcol) { draw_pix(x,y,c); }
/*{\Mopl     sets the color of the pixel at position $(x,y)$ to $c$. }*/



void draw_pixels(int n, double* xcoord, double* ycoord,
                                                color c=window::fgcol);
/*{\Mopl     draws all pixels $(xcoord[i],ycoord[i])$ for $0 \le i \le n-1$. }*/



// segments

/*{\Mtext
{\bf 3.4.2 Drawing line segments}
}*/

void draw_segment(double x1, double y1, double x2, double y2, color c=window::fgcol);
/*{\Mopl     draws a line segment from $(x_1,y_1)$ to $(x_2,y_2)$.}*/


// lines

/*{\Mtext
{\bf 3.4.3 Drawing lines}
}*/

void draw_line(double x1, double y1, double x2, double y2, color c=window::fgcol);
/*{\Mopl     draws a straight line passing through points $(x_1,y_1)$ and 
	     $(x_2,y_2)$.}*/

void draw_hline(double y, color c=window::fgcol);
/*{\Mopl     draws a horizontal line with y-coordinate $y$. }*/

void draw_vline(double x, color c=window::fgcol);
/*{\Mopl     draws a vertical line with x-coordinate $x$. }*/



/*{\Mtext
{\bf 3.4.4 Drawing Rays}
}*/

void draw_ray(double x1, double y1, double x2, double y2, color c=window::fgcol);
/*{\Mopl     draws a ray starting in $(x_1,y_1)$ and passing through 
	     $(x_2,y_2)$.}*/


/*{\Mtext
{\bf 3.4.5 Drawing Arcs and Curves}
}*/

static void compute_bezier(const std::list<window_point>&, int,double*,double*);
static void bezier_segments(const std::list<window_point>&, int,window_point&,window_point&,window_point&,window_point&);

void draw_bezier(const std::list<window_point>& L, int m, color c);
void draw_bezier(const std::list<window_point>&, int n, int arr, double d,color);



// arrows

/*{\Mtext
{\bf 3.4.6 Drawing arrows}
}*/

void draw_arrow(double x1, double y1, double x2, double y2, color c =window::fgcol);
/*{\Mopl     draws an arrow pointing from $(x_1,y_1)$ to $(x_2,y_2)$.}*/

void draw_arrow(const window_point& p, const window_point& q, color c=window::fgcol);
/*{\Mopl     draws an arrow pointing from point $p$ to point $q$.}*/


void draw_polyline_arrow(const std::list<window_point>& lp, color c=window::fgcol);
/*{\Mopl     draws a polyline arrow with vertex sequence $lp$.}*/



void draw_bezier_arrow(const std::list<window_point>& C, int n,  color c=window::fgcol);
/*{\Mopl     draws the bezier curve with control polygon $C$ by a polyline 
             with $n$ points, the last segment is drawn as an arrow. }*/



window_point arrow_head(const window_point&, double, double, double*, double*);

window_point draw_arrow_head(const window_point& p, double dir, double d, color c);

window_point draw_arrow_head(const window_point& p, double dir, color c=window::fgcol);
/*{\Mopl     draws an arrow head at position $p$ pointing to direction $dir$,
             where $dir$ is an angle from $[0,2\pi]$. }*/



//circles

/*{\Mtext
{\bf 3.4.7 Drawing circles}
}*/

void draw_circle(double x, double y, double r, color c=window::fgcol);
/*{\Mopl     draws the circle with center $(x,y)$ and radius $r$.}*/

void draw_circle(const window_point& p, double r, color c=window::fgcol);
/*{\Mopl     draws the circle with center $p$ and radius $r$.}*/


void draw_ellipse(double x, double y, double r1, double r2, color c=window::fgcol);
/*{\Mopl     draws the ellipse with center $(x,y)$ and radii $r_1$ and $r_2$.}*/

void draw_ellipse(const window_point& p, double r1, double r2, color c=window::fgcol);
/*{\Mopl     draws the ellipse with center $p$ and radii $r_1$ and $r_2$.}*/



/*{\Mtext
{\bf 3.4.8 Drawing discs}
}*/

void draw_disc(double x, double y, double r, color c=window::fgcol);
/*{\Mopl     draws a filled circle with center $(x,y)$ and radius $r$.}*/

void draw_disc(const window_point& p, double r, color c=window::fgcol);
/*{\Mopl     draws a filled circle with center $p$ and radius $r$.}*/


void draw_filled_ellipse(double x, double y, double r1, double r2, color c=window::fgcol);
/*{\Mopl  draws a filled ellipse with center $(x,y)$ and radii $r_1$ and $r_2$.}*/

void draw_filled_ellipse(const window_point& p, double r1, double r2, color c=window::fgcol);
/*{\Mopl  draws a filled ellipse with center $p$ and radii $r_1$ and $r_2$.}*/


//polygons 

/*{\Mtext
{\bf 3.4.9 Drawing polygons }
}*/

void draw_polyline(const std::list<window_point>& lp, int arrow, double d, color c);

void draw_polyline(const std::list<window_point>& lp, color c=window::fgcol);
/*{\Mopl     draws a polyline with vertex sequence $lp$.}*/

void draw_polyline(int n, double* xc, double* yc, color c=window::fgcol)
{ BASE_WINDOW::draw_polyline(n,xc,yc,c); }
/*{\Mopl     draws a polyline with vertex sequence 
             $(xc[0],yc[0]), \dots, (xc[n-1],yc[n-1])$.}*/


void draw_polygon(const std::list<window_point>& lp, color c=window::fgcol);
/*{\Mopl     draws the polygon with vertex sequence $lp$.}*/

void draw_oriented_polygon(const std::list<window_point>& lp, color c=window::fgcol);
/*{\Mopl     draws the polygon with vertex sequence |lp| and indicates the 
              orientation by an arrow.}*/


void draw_filled_polygon(const std::list<window_point>& lp, color c=window::fgcol);
/*{\Mopl     draws the filled polygon with vertex sequence $lp$.}*/



void draw_rectangle(double x0, double  y0, double x1, double y1, color=window::fgcol);
/*{\Mopl     draws a rectangle with lower left corner $(x_0,y_0)$ and upper
             right corner $(x_1,y_1)$.\\
             \precond $x_0 < x_1$ and $y_0 < y_1$. }*/

void draw_rectangle(window_point p, window_point q, color=window::fgcol);
/*{\Mopl     draws a rectangle with lower left corner $p$ and upper
             right corner $q$.\\
             \precond $p < q$. }*/


void draw_filled_rectangle(double,double,double,double,color=window::fgcol);

void draw_box(double x0, double y0, double x1, double y1,color c=window::fgcol)
{ draw_filled_rectangle(x0,y0,x1,y1,c); }
/*{\Mopl     draws a filled rectangle with lower left corner $(x_0,y_0)$ and upper
             right corner $(x_1,y_1)$.\\
             \precond $x_0 < x_1$ and $y_0 < y_1$. }*/

void draw_filled_rectangle(window_point p, window_point q, color=window::fgcol);
/*{\Mopl     draws a filled rectangle with lower left corner $p$ and upper
             right corner $q$.\\
             \precond $p < q$. }*/

void draw_box(window_point p, window_point q, color c=window::fgcol)
{ draw_filled_rectangle(p,q,c); }
/*{\Mopl     same as |draw_filled_rectangle(p,q,c)|. }*/


void draw_roundrect(double x0, double  y0, double x1, double y1, double rndness,
                                                      color col=window::fgcol)
{ BASE_WINDOW::draw_roundrect(x0,y0,x1,y1,rndness,col); }
/*{\Mopl     draws a rectangle $(x_0,y_0,x_1,y_1)$ with round corners.
             The |rndness| argument must be real number in the intervall
             $[0,1]$ and defined the ``roundness'' of the rectangle. }*/

void draw_rectangle(window_point p, window_point q, double rndness, color col=window::fgcol)
{ draw_roundrect(p.x(),p.y(),q.x(),q.y(),rndness,col); }
/*{\Mopl     draws a round rectangle with lower left corner $p$, upper
             right corner $q$, and roundness |rndness|. }*/


void draw_roundbox(double x0, double  y0, double x1, double y1, double rndness,
                                                      color col=window::fgcol)
{ BASE_WINDOW::draw_roundbox(x0,y0,x1,y1,rndness,col); }
/*{\Mopl     draws a filled rectangle $(x_0,y_0,x_1,y_1)$ with round corners.
             The |rndness| argument must be real number in the intervall
             $[0,1]$ and defined the ``roundness'' of the rectangle. }*/

void draw_roundbox(window_point p, window_point q, double rndness, color col=window::fgcol)
{ draw_roundbox(p.x(),p.y(),q.x(),q.y(),rndness,col); }
/*{\Mopl     draws a round filled rectangle with lower left corner $p$, upper
             right corner $q$, and roundness |rndness|. }*/



void draw_triangle(window_point a, window_point b, window_point c, color=window::fgcol);
/*{\Mopl     draws triangle $(a,b,c)$. }*/

void draw_filled_triangle(window_point a, window_point b, window_point c, color=window::fgcol);
/*{\Mopl     draws filled triangle $(a,b,c)$. }*/


// functions

/*{\Mtext
{\bf 3.4.10 Drawing functions}
}*/

void plot_xy(double x0, double x1, win_draw_func F, color c=window::fgcol);
/*{\Mopl     draws the graph of function $F$ in the x-range $[x_0,x_1]$, i.e., 
             all pixels $(x,y)$ with $y = F(x)$ and $x_0\le x\le x_1$.}*/

void plot_yx(double y0, double y1, win_draw_func F, color c=window::fgcol);
/*{\Mopl     draws the graph of function $F$ in the y-range $[y_0,y_1]$, i.e., 
             all pixels $(x,y)$ with $x = F(y)$ and $y_0\le y\le y_1$.}*/




/*{\Mtext
{\bf 3.4.11 Drawing text}
}*/

void draw_text_with_cursor(double x, double y, string s, int cursor, 
                                                         color c=window::fgcol);

void draw_text(double x, double y, string s, color c=window::fgcol);
/*{\Mopl     writes string $s$ starting at position $(x,y)$.}*/

void draw_text(const window_point& p, string s, color c=window::fgcol);
/*{\Mopl     writes string $s$ starting at position $p$.}*/

void draw_ctext(double x, double y, string s, color c=window::fgcol);
/*{\Mopl     writes string $s$ centered at position $(x,y)$.}*/

void draw_ctext(const window_point& p, string s, color c=window::fgcol);
/*{\Mopl     writes string $s$ centered at position $p$.}*/

void draw_ctext(string s, color c=window::fgcol);
/*{\Mopl     writes string $s$ centered in window $W$.}*/

double text_box(double x0, double x1, double y, string s, bool draw=true);
/*{\Mopl     formats and writes string $s$ into a box with its left border
             at x-coordinate $x0$, its right border at $x1$, and its upper 
             border at y-coordinate $y$. Some LaTeX-like formatting commands
             can be used: {\tt $\backslash$bf, $\backslash$tt, $\backslash$rm,
             $\backslash$n, $\backslash$c, $\backslash$<color>,} \dots 
             returns y-coordinate of lower border of box. If the optional
             last parameter |draw| is set to |false| no drawing takes place
             and only the lower y-coordinate of the box is computed. }*/

void text_box(string s);
/*{\Mopl     as above with |x0 = W.xmin()|, |x1 = W.xmax()|, 
             and |y = W.ymax()|. }*/



void message(string s) {BASE_WINDOW::message(s.c_str());};
/*{\Mop      displays the message $s$ (each call adds a new line).}*/

/*{\Moptions nextwarning=no}*/
void del_message()  { BASE_WINDOW::del_messages(); };
void del_messages() { BASE_WINDOW::del_messages(); };
/*{\Mop      deletes the text written by all previous message 
	     operations.}*/



// nodes

/*{\Mtext
{\bf 3.4.12 Drawing nodes}

Nodes are represented by circles of diameter $node\_width$.
}*/

void draw_node(double x0, double y0, color c=window::fgcol);
/*{\Mopl     draws a node at position $(x_0,y_0)$.}*/

void draw_node(const window_point& p, color c=window::fgcol);
/*{\Mopl     draws a node at position $p$.}*/

void draw_filled_node(double x0, double y0, color c=window::bgcol);
/*{\Mopl     draws a filled node at position $(x_0,y_0)$.}*/

void draw_filled_node(const window_point& p, color c=window::bgcol);
/*{\Mopl     draws a filled node at position $p$.}*/

void draw_text_node(double x, double y, string s, color c=window::bgcol);
/*{\Mopl     draws a node with label $s$ at position $(x,y)$. }*/

void draw_text_node(const window_point& p, string s, color c=window::bgcol);
/*{\Mopl     draws a node with label $s$ at position $p$. }*/

void draw_int_node(double x, double y, int i, color c=window::bgcol);
/*{\Mopl     draws a node with integer label $i$ at position 
	     $(x,y)$. }*/

void draw_int_node(const window_point& p, int i, color c=window::bgcol);
/*{\Mopl     draws a node with integer label $i$ at position  $p$. }*/



// edges

/*{\Mtext
{\bf 3.4.13 Drawing edges}

Edges are drawn as straigth line segments or arrows with a clearance
of $node\_width/2$ at each end.  }*/

void draw_edge(double x1, double y1, double x2, double y2, color c=window::fgcol);
/*{\Mopl     draws an edge from $(x_1,y_1)$ to $(x_2,y_2)$.}*/

void draw_edge(const window_point& p, const window_point& q, color c=window::fgcol);
/*{\Mopl     draws an edge from $p$ to $q$.}*/

void draw_edge_arrow(double x1, double y1, double x2, double y2, color c=window::fgcol);
/*{\Mopl     draws a directed edge from $(x_1,y_1)$ to $(x_2,y_2)$.}*/

void draw_edge_arrow(const window_point& p, const window_point& q, color c=window::fgcol);
/*{\Mopl     draws a directed edge from $p$ to $q$.}*/



/*{\Mtext
{\bf 3.4.14 Bitmaps and Pixrects}
}*/

char* create_bitmap(int w, int h, unsigned char* bm_data)
{ return BASE_WINDOW::create_bitmap(w,h,bm_data); }
/*{\Mopl   creates a bitmap (monochrome pixrect) of width $w$, height $h$, 
           from the bits in |data|. }*/

char* create_pixrect(const char** xpm_str)
{ return BASE_WINDOW::create_pixrect(xpm_str); }
/*{\Mopl   creates a pixrect from the {\bf xpm} data string |xpm_str|. }*/

char* create_pixrect(string xpm_file);
/*{\Mopl   creates a pixrect from the {\bf xpm} file |xpm_file|. }*/

char* create_pixrect(int w, int h, unsigned char* bm_data, int fg=window::fgcol,
                                                           int bg=window::bgcol)
{ return BASE_WINDOW::create_pixrect(w,h,bm_data,fg,bg); }
/*{\Mopl   creates a pixrect of width $w$, height $h$, foreground color 
           $fg$, and background color $bg$ from bitmap |data|. }*/


char* get_pixrect(double x1, double y1, double x2, double y2)
{ return BASE_WINDOW::get_pixrect(x1,y1,x2,y2); }
/*{\Mopl   creates a color pixrect of width $w=x_2-x_1$, height $h=y2-y1$,
           and  copies all pixels from the rectangular area $(x_1,x_2,y_1,y_2)$
           of $W$ into it. }*/

char* get_window_pixrect()
{ return BASE_WINDOW::get_window_pixrect(); }
/*{\Mopl   creates a pixrect copy of the current window contents. }*/


int get_width(char* pr) const { return BASE_WINDOW::get_width(pr); }
/*{\Mopl   returns the width (number of pixels in a row) 
           of pixrect |pr|. }*/

int get_height(char* pr) const { return BASE_WINDOW::get_height(pr); }
/*{\Mopl   returns the height (number of pixels in a column) 
           of pixrect |pr|. }*/


void  put_pixrect(double x, double y, char* pr)
{ BASE_WINDOW::put_pixrect(x,y,pr); }
/*{\Mopl   copies the contents of pixrect $pr$ with lower left corner at
           position $(x,y)$ into $W$. }*/

void  put_pixrect(window_point p, char* pr)
{ BASE_WINDOW::put_pixrect(p.x(),p.y(),pr); }
/*{\Mopl   copies the contents of pixrect $pr$ with lower left corner at
           position $p$ into $W$. }*/

void  center_pixrect(double x, double y, char* pr)
{ BASE_WINDOW::center_pixrect(x,y,pr); }
/*{\Mopl   copies the contents of pixrect $pr$ into $W$ such that its
           center lies on position $(x,y)$. }*/


void  put_pixrect(char* pr) { BASE_WINDOW::put_pixrect(pr); }
/*{\Mopl   copies the contents of pixrect $pr$ with lower left corner at
           position |(W.xmin(),W.ymin())| into $W$. }*/


void  put_bitmap(double x, double y, char* bm, color c=window::fgcol)
{ BASE_WINDOW::put_bitmap(x,y,bm,c); }
/*{\Mopl   draws all pixels corresponding to 1-bits in $bm$ with color $c$,
           here the lower left corner of $bm$ corresponds to the pixel at
           position $(x,y)$ in $W$. }*/

void  put_pixrect(double x, double y, char* pr,int x0, int y0, int w, int h)
{ BASE_WINDOW::put_pixrect(x,y,pr,x0,y0,w,h); }
/*{\Mopl   copies (pixel) rectangle $(x0,y0,x0+w,y0+h)$ of $pr$ with lower left 
           corner at position $(x,y)$ into $W$. }*/


void  del_bitmap(char* bm) { BASE_WINDOW::del_bitmap(bm); }
/*{\Mopl   destroys bitmap $bm$. }*/

void  del_pixrect(char* pr) { BASE_WINDOW::del_pixrect(pr); }
/*{\Mopl   destroys pixrect $pr$. }*/


void  copy_rect(double x0, double y0, double x1, double y1, double x, double y)
{ BASE_WINDOW::copy_rect(x0,y0,x1,y1,x,y); }
/*{\Mopl   copies all pixels of rectangle $(x_0,y_0,x_1,y_1)$ into
           the rectangle $(x,y,x+w,y+h)$, where $w = x_1-x_0$ and 
           $h = y_1-y_0$. }*/


void  screenshot(string fname, bool full_color=true);
/*{\Mopl   creates a screenshot of the current window. On unix
           systems suffix |.ps| is appended to |fname| and
           the output format is postscript. On windows systems
           the suffix |.wmf| is added and the format is windows
           metafile. If the flag |full_color| is set to |false|
           colors will be translated into grey scales. }*/  



/*{\Mtext
{\bf 3.4.15 Buffering}
}*/

void start_buffering() { BASE_WINDOW::start_buffering(); }
/*{\Mopl  starts buffering mode for $W$. If $W$ has no associated buffer
          a buffer pixrect $buf$ of the same size as the current drawing 
          area of $W$ is created. All subsequent drawing operations draw 
          into $buf$ instead of $W$ until buffering mode is ended by 
          calling |W.stop_buffering()|. }*/


void start_buffering(int w, int h) { BASE_WINDOW::start_buffering(w,h); }

void set_buffer(char* pr)          { BASE_WINDOW::set_buffer(pr); }



void flush_buffer() { BASE_WINDOW::flush_buffer(); }
/*{\Mopl copies the contents of the buffer pixrect into the drawing
         area of $W$. }*/

void flush_buffer(double dx, double dy) { BASE_WINDOW::flush_buffer(dx,dy); }
/*{\Mopl copies the contents of the buffer pixrect translated by vector
         $(dx,dy)$ into the drawing area of $W$. }*/


void flush_buffer(double x0, double y0, double x1, double y1) 
{ BASE_WINDOW::flush_buffer(x0,y0,x1,y1); }
/*{\Mopl  copies the contents of rectangle $(x0,y0,x1,y1)$ of the
          buffer pixrect into the corresponding rectangle of the
          drawin area. }*/

void flush_buffer(double dx,double dy,double x0,double y0,double x1,double y1) 
{ BASE_WINDOW::flush_buffer(dx,dy,x0,y0,x1,y1); }
/*{\Mopl  copies the contents of rectangle $(x0,y0,x1,y1)$ of the
          buffer pixrect into the corresponding rectangle of the
          drawin area translated by vector $(dx,dy)$. }*/


void stop_buffering() { BASE_WINDOW::stop_buffering(); }
/*{\Mopl ends buffering mode. }*/

void stop_buffering(char*& prect) { BASE_WINDOW::stop_buffering(prect); }
/*{\Mopl ends buffering mode and returns the current buffer pixrect
         in |prect|. }*/



/*{\Mtext
{\bf 3.4.16 Clipping}
}*/

void set_clip_rectangle(double x0, double y0, double x1, double y1)
{ BASE_WINDOW::set_clip_rectangle(x0,y0,x1,y1); }
/*{\Mopl sets the clipping region of $W$ to rectangle $(x0,y0,x1,y1)$. }*/


void set_clip_ellipse(double x0, double y0, double r1, double r2)
{ BASE_WINDOW::set_clip_ellipse(x0,y0,r1,r2); }
/*{\Mopl sets the clipping region of $W$ to ellipse with center $(x,y)$ and
         horizontal radius $r1$ and vertical radius $r2$.}*/

void set_clip_polygon(const std::list<window_point>&);



void reset_clipping() { BASE_WINDOW::reset_clipping(); }
/*{\Mopl restores the clipping region to the entire drawing area of $W$.}*/





// mouse input

/*{\Mtext
{\bf 3.5 Input}

The main input operation for reading positions, mouse clicks, and buttons
from a window $W$ is the operation $W$.read\_mouse(). This operation is 
blocking, i.e., waits for a button to be pressed which is either a ``real''
button on the mouse device pressed inside the drawing area of $W$ or a 
button in the panel section of $W$.
In both cases, the number of the selected button is returned. Mouse
buttons have pre-defined numbers MOUSE\_BUTTON(1) for the left button, 
MOUSE\_BUTTON(2) for the middle button, and MOUSE\_BUTTON(3) for the
right button. The numbers of the panel buttons can be defined by the
user. If the selected button has an associated action function or sub-window 
this function/window is executed/opened (cf. \ref{panel-operations} for
details).

There is also a non-blocking version $W$.get\_mouse() which returns 
the constant NO\_BUTTON if no button was pressed. 

The window data type also provides two more general input operations 
$W$.read\_event() and $W$.get\_event() for reading events. They return 
the event type (enumeration in \<CGAL/LEDA/impl/x\_window.h\>), the value of 
the event, the position of the event in the drawing section, and a time 
stamp of the event.

}*/

void set_point_buffer(const window_point p);

/*{\Mtext
{\bf 3.5.1 Read Mouse}
}*/

int read_mouse();
/*{\Mop      waits until a mouse button is pressed inside of the drawing 
             area or until a button of the panel section is selected.  
             In both cases, the number $n$ of the button is 
             returned which is one of the predefined constants 
             MOUSE\_BUTTON($i$) with $i\in\{1,2,3\}$ for mouse buttons and
             a user defined value (defined when adding the button with 
             $W$.button()) for panel buttons.
             If the button has an associated action function this function is 
             called with parameter $n$. If the button has an associated
             window $M$ it is opened and $M$.read\_mouse() is returned. }*/

int read() { return read_mouse(); }

int read_mouse(double& x, double& y);
/*{\Mopl     If a button is pressed inside the drawing area the 
             current position of the cursor is assigned to $(x,y)$. 
             The operation returns the number of the pressed button
             (see $W$.read\_mouse().)}*/ 

int read_mouse(window_point& p);
/*{\Mopl     as above, the current position is assigned to point |p|. }*/


int read_mouse_seg(double x0, double y0, double& x, double& y);
/*{\Mopl     displays a line segment from $(x_0,y_0)$ to the
	     current cursor position until a mouse button is
	     pressed inside the drawing section of $W$. When a 
             button is pressed the current position is assigned to $(x,y)$ 
             and the number of the pressed button is returned.}*/

int read_mouse_seg(const window_point& p, window_point& q);
/*{\Mopl     as above with |x0 = p.xcoord()| and |y0 = p.ycoord()|
             and the current position is assigned to |q|. }*/

int read_mouse_line(double x0, double y0, double& x, double& y);
/*{\Mopl     displays a line passing through $(x_0,y_0)$ and the
	     current cursor position until a mouse button is
	     pressed inside the drawing section of $W$. When a 
             button is pressed the current position is assigned to $(x,y)$ 
             and the number of the pressed button is returned.}*/

int read_mouse_line(const window_point& p, window_point& q);
/*{\Mopl     as above with |x0 = p.xcoord()| and |y0 = p.ycoord()|
             and the current position is assigned to |q|. }*/

int read_mouse_ray(double x0, double y0, double& x, double& y);
/*{\Mopl     displays a ray from $(x_0,y_0)$ passing through the
	     current cursor position until a mouse button is
	     pressed inside the drawing section of $W$. When a 
             button is pressed the current position is assigned to $(x,y)$ 
             and the number of the pressed button is returned.}*/

int read_mouse_ray(const window_point& p, window_point& q);
/*{\Mopl     as above with |x0 = p.xcoord()| and |y0 = p.ycoord()|
             and the current position is assigned to |q|. }*/

int read_mouse_rect(double x0, double y0, double& x, double& y);
/*{\Mopl     displays a rectangle with diagonal from $(x_0,y_0)$ 
     	     to the current cursor position until a mouse button 
	     is pressed inside the drawing section of $W$.
             When a button is pressed the current 
	     position is assigned to $(x,y)$ and the number of the 
             pressed button is returned.}*/

int read_mouse_rect(const window_point& p, window_point& q);
/*{\Mopl     as above with |x0 = p.xcoord()| and |y0 = p.ycoord()|
             and the current position is assigned to |q|. }*/

int read_mouse_circle(double x0, double y0, double& x, double& y);
/*{\Mopl     displays a circle with center $(x_0,y_0)$ passing 
	     through the current cursor position until a mouse 
	     button is pressed inside the drawing section of $W$.
             When a button is pressed the 
	     current position is assigned to $(x,y)$ and the 
	     number of the pressed button is returned.}*/

int read_mouse_circle(const window_point& p, window_point& q);
/*{\Mopl     as above with |x0 = p.xcoord()| and |y0 = p.ycoord()|
             and the current position is assigned to |q|. }*/


int read_mouse_action(mouse_action_func, double&, double&);
int read_mouse_action(mouse_action_func, window_point&);


int get_mouse();
/*{\Mop     non-blocking read operation, i.e., if a button was pressed 
            its number is returned, otherwise the constant NO\_BUTTON is 
            returned. }*/

int get_mouse(double& x, double& y);
/*{\Mopl    if a mouse button was pressed the corresponding position is
            assigned to $(x,y)$ and the button number is returned, 
            otherwise the constant NO\_BUTTON is returned. }*/

int get_mouse(window_point& p);
/*{\Mop     if a mouse button was pressed the corresponding position is
            assigned to $p$ and the button number is returned, 
            otherwise the constant NO\_BUTTON is returned. }*/

// for backward compatibility
int get_button() { return get_mouse(); }
int get_button(window_point& p) { return get_mouse(p); }
int get_button(double& x, double& y) { return get_mouse(x,y); }

int read_mouse(double x0, double y0, int timeout1, int timeout2, 
                       bool& double_click, bool& drag);
/*{\Mop  ... }*/

int read_mouse(const window_point& p, int timeout1, int timeout2, 
                       bool& double_click, bool& drag)
{ return read_mouse(p.x(),p.y(),timeout1,timeout2,double_click,drag);}
/*{\Mop  ... }*/


/*{\Mtext
\bigskip
{\bf 3.5.2 Events}
}*/

int  read_event(int& val, double& x, double& y, unsigned long& t)
{ return BASE_WINDOW::read_event(val,x,y,t); }
/*{\Mopl   waits for next event in window $W$ and returns it. 
           Assigns the button or key to $val$, the position in $W$ 
           to  $(x,y)$, and the time stamp of the event to $t$. 
           Possible events are (cf. \<CGAL/LEDA/impl/x\_window.h\>): 
           key\_press\_event, key\_release\_event,
           button\_press\_event, button\_release\_event, 
           configure\_event, motion\_event, destroy\_event. }*/

int  read_event(int& val, double& x, double& y, unsigned long& t, int timeout)
{ return BASE_WINDOW::read_event(val,x,y,t,timeout); }
/*{\Mopl   as above, but waits only $timeout$ milliseconds;
           if no event occured the special event $no\_event$ is 
           returned. }*/

int  read_event(int& val, double& x, double& y)
{ return BASE_WINDOW::read_event(val,x,y); }
/*{\Mopl   waits for next event in window $W$ and returns it. 
           Assigns the button or key to $val$ and the position
           in $W$ to  $(x,y)$. }*/

int  read_event()
{ int val; double x,y; return BASE_WINDOW::read_event(val,x,y); }
/*{\Mopl   waits for next event in window $W$ and returns it. }*/

int  get_event(int& val, double& x, double& y)
{ return BASE_WINDOW::get_event(val,x,y); }
/*{\Mopl   if there is an event for window $W$ in the event queue a 
           $W.read\_event$ operation is performed, otherwise the integer 
           constant $no\_event$ is returned. }*/


// get shift key info for last handled mouse button event

bool shift_key_down() { return BASE_WINDOW::shift_key_down()!=0; }
/*{\Mop   returns $true$ if a {\em shift} key was pressed during
          the last handled mouse button event. }*/ 

bool ctrl_key_down()  { return BASE_WINDOW::ctrl_key_down()!=0; }
/*{\Mop   returns $true$ if a {\em ctrl} key was pressed during
          the last handled mouse button event. }*/ 

bool alt_key_down()  { return BASE_WINDOW::alt_key_down()!=0; }
/*{\Mop   returns $true$ if an {\em alt} key was pressed during
          the last handled mouse button event. }*/ 


int button_press_time()
{ return (int)BASE_WINDOW::button_press_time(); }
/*{\Mop     returns the time-stamp (in msec) of the last button press event. }*/


int button_release_time()
{ return (int)BASE_WINDOW::button_release_time(); }
/*{\Mop     returns the time-stamp (in msec) of the last button release event. }*/
        
 

/*{\Mtext
\bigskip
{\bf 3.6 Panel Input}

The operations listed in this section are useful for simple input of
strings, numbers, and Boolean values.
}*/

bool    confirm(string s);
/*{\Mop      displays string $s$ and asks for confirmation. 
	     Returns true iff the answer was ``yes''.}*/

/*{\Moptions nextwarning=no}*/
void    notice(string s);
void    acknowledge(string s);
/*{\Mopl     displays string $s$ and asks for acknowledgement.}*/


int     read_panel(string h, int n, string* S);
/*{\Mopl     displays a panel with header $h$ and an array of $n$ 
             buttons with labels $S[0..n-1]$, returns the index of 
             the selected button.}*/

int     read_vpanel(string h, int n, string* S);
/*{\Mopl     like read\_panel with vertical button layout.}*/

string  read_string(string p);
/*{\Mop      displays a panel with prompt $p$ for string input, 
    	     returns the input.}*/

double  read_real(string p);
/*{\Mop      displays a panel with prompt $p$ for double input 
 	     returns the input.}*/

int     read_int(string p);
/*{\Mop      displays a panel with prompt $p$ for integer input, 
	     returns the input.}*/


void  panel_read(string p, int& x);
void  panel_read(string p, double& x);
void  panel_read(string p, string& x);




// I/O operators
/*{\Mtext
\bigskip
{\bf 3.7 Input and output operators}

}*/

/*{\Mtext
{\bf 3.7.1 Output}
\setopdims{2.5cm}{4.4cm}
}*/

// drawing points, segments, circles ...

// points ...
window& draw(const window_point& p,color c=window::fgcol);

// segments ...
window& draw(const window_point& ps, const window_point& pt, color c=window::fgcol);

// circles ...
window& draw(const window_point& m, double rad, color c=window::fgcol);



/*{\Mtext
\bigskip
{\bf 3.7.2 Input}
}*/

// reading points, segments, circles ...

// reading points ...
window& read(window_point&);

// reading segments ...
window& read(window_point&,window_point&);

// reading circles ...
window& read(window_point& m, double& rad);



std::list<window_point> read_polygon();
/*{\Mopl    as above, however, returns list of vertices. }*/



/*{\Mtext
\bigskip
{\bf 3.8 Non-Member Functions} 
}*/

friend int  read_event(window*& w, int& val, double& x, double& y)
{ return BASE_WINDOW::read_event((BASE_WINDOW*&)w,val,x,y); }
/*\{\Mfuncl waits for next event and returns it. Assigns the window to
            $w$, the button or key to $val$ and the position in $w$ to 
            $(x,y)$. Possible events are (cf. <CGAL/LEDA/impl/x_window.h>): 
            key_press_event, key_release_event,
            button_press_event, button_release_event, 
            configure_event,motion_event, destroy_event. }*/


friend int  get_event(window*& w, int& val, double& x, double& y)
{ return BASE_WINDOW::get_event((BASE_WINDOW*&)w,val,x,y); }
/*\{\Mfuncl if the event queue is not empty a $read_event$ operation
            is performed, otherwise the integer constant $no\_event$ is 
            returned. }*/


friend int read_mouse(window*& w, double& x, double& y)
{ return BASE_WINDOW::read_mouse((BASE_WINDOW*&)w,x,y); }
/*{\Mfuncl  waits for mouse input, assigns a pointer to the 
            corresponding window to $w$ and the position in 
            $*w$ to $(x,y)$ and returns the pressed button. }*/

friend int get_mouse(window*& w, double& x, double& y)
{ return BASE_WINDOW::get_mouse((BASE_WINDOW*&)w,x,y); }
/*{\Mfuncl  non-blocking variant of |read_mouse|, returns
            NO\_BUTTON if no button was pressed. }*/
 
        
friend void put_back_event();
/*{\Mfunc   puts last handled event back to the event queue. }*/




//------------------------------------------------------------------------------
// panel operations
//------------------------------------------------------------------------------

/*{\Mtext
\bigskip
{\bf 3.9 Panel Operations \label{panel-operations} }

The panel section of a window is used for displaying text messages and 
for updating the values of variables. It consists of a list of panel items and 
a list of buttons.  
The operations in this section add panel items or buttons to the panel section 
of $W$. Note that they have to be called before the window is displayed the
first time.

In general, a panel item consists of a string label and an associated variable
of a certain type (int, bool, string, double, color). The value of this variable
can be manipulated through the item.
Each button has a label (displayed on the button) and an associated number.
The number of a button is either defined by the user or is the rank of the
button in the list of all buttons. If a button is pressed  (i.e. selected
by a mouse click) during a $read\_mouse$ operation its number is returned.

{\em Action functions} can be associated with buttons and some items 
(e.g. slider items) whenever a button with an associated action function
is pressed this function is called with the number of the button as
actual parameter. Action functions of items are called whenever the
value of the corresponding variable is changed with the new value
as actual parameter. All action functions must have the type
$void\ \ func(int)$. 

Another way to define a button is to associate another window with it.
In this case the button will have a menu sign and
as soon as it is pressed the attached window will open. 
This method can be used to implement pop-up menues.
The return value of the current $read\_mouse$ operation will be 
the number associated with the button in the menu.
}*/

/*{\Mtext
\bigskip
{\bf 3.9.1 General Settings}
}*/

void set_panel_bg_color(color c)  { BASE_WINDOW::set_panel_bg_color(c); }
/*{\Mop  sets the background color of the panel area to $c$. }*/

void buttons_per_line(int n) { BASE_WINDOW::buttons_per_line(n); }
/*{\Mop  defines the maximal number $n$ of buttons per line. }*/

void set_button_space(int s) { BASE_WINDOW::set_button_space(s); }
/*{\Mop  sets the space between to adjacent buttons to $s$ pixels. }*/

void set_item_height(int h) { BASE_WINDOW::set_item_height(h); }
/*{\Mop  sets the vertical size of all items to $h$ pixels. }*/

void set_item_width(int w) { BASE_WINDOW::set_item_width(w); }
/*{\Mop  sets the horizontal size of all {\em slider} and {\em string} items
         to $w$ pixels. }*/

void set_bitmap_colors(int c0, int c1) { BASE_WINDOW::set_bitmap_colors(c0,c1);}
/*{\Mop  sets the unpressed/pressed colors used for drawing the pixels in 
         bitmap buttons to $c_0$ and $c_1$. }*/





/*{\Mtext
\bigskip
{\bf 3.9.2 Simple Panel Items}
}*/

panel_item text_item(string s);
/*{\Mop      adds a text\_item $s$ to $W$.}*/

panel_item bool_item(string s, bool& x, const char* hlp=0);
/*{\Mopl     adds a boolean item with label $s$ and variable $x$ to $W$.}*/

panel_item bool_item(string s, bool& x, void (*F)(int), const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item bool_item(string s, bool& x, const window_handler& obj, const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/


panel_item int_item(string s, int& x, const char* hlp=0);
/*{\Mopl     adds an integer item with label $s$ and variable $x$ to $W$.}*/


panel_item string_item(string s, string& x, void (*F)(char*), 
                                            const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item string_item(string s, string& x, const window_handler& obj, 
                                            const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

panel_item string_item(string s, string& x, const char* hlp=0);
/*{\Mopl     adds a string item with label $s$ and variable $x$ to $W$.}*/


panel_item real_item(string s, double& x, const char* hlp=0);

panel_item double_item(string s, double& x, const char* hlp=0);
/*{\Mopl     adds a real item with label $s$ and variable $x$ to $W$.}*/


panel_item color_item(string s, color& x, const char* hlp=0);
/*{\Mopl     adds a color item with label $s$ and variable $x$ to $W$.}*/

panel_item color_item(string s, color& x, void (*F)(int), const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item color_item(string s, color& x, const window_handler& obj, const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

panel_item pstyle_item(string s, point_style& x, const char* hlp=0);
/*{\Mopl     adds a point style item with label $s$ and variable $x$ to $W$.}*/

panel_item pstyle_item(string s,point_style& x,void(*F)(int),const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item pstyle_item(string s,point_style& x,const window_handler& obj,const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/


panel_item lstyle_item(string s, line_style& x, const char* hlp=0);
/*{\Mopl     adds a line style item with label $s$ and variable $x$ to $W$.}*/

panel_item lstyle_item(string s, line_style& x,void(*F)(int),const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item lstyle_item(string s, line_style& x, const window_handler& obj,const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

panel_item lwidth_item(string s, int& x, const char* hlp=0);
/*{\Mopl     adds a line width item with label $s$ and variable $x$ to $W$.}*/

panel_item lwidth_item(string s, int& x, void(*F)(int), const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item lwidth_item(string s, int& x, const window_handler& obj, const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

/*{\Mtext
\bigskip
{\bf 3.9.3 Integer Choice Items}
}*/

panel_item int_item(string s, int& x, int l, int h, int step,const char* hlp=0);
/*{\Mopl     adds an integer choice item with label $s$, variable $x$, 
	     range $l$,\dots, $h$, and step size $step$ to $W$.}*/

panel_item int_item(string s, int& x, int l, int h, int step,void (*F)(int),
                                                             const char* hlp=0);
/*{\Mopl     adds an integer choice item with label $s$, variable $x$, 
	     range $l$,\dots, $h$, and step size $step$ to $W$. Function 
             $F(x)$ is executed whenever the value of $x$ is changed. }*/

panel_item int_item(string s, int& x, int l, int h, int step,const window_handler& obj,
                                                             const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

panel_item int_item(string s, int& x, int l, int h, const char* hlp=0);
/*{\Mopl     adds an integer slider item with label $s$, variable $x$, and 
	     range $l$,\dots,$h$ to $W$.}*/

panel_item int_item(string s, int& x, int l, int h, void (*F)(int),const char* hlp=0);
/*{\Mopl     adds an integer slider item with label $s$, variable $x$, and 
	     range $l$,\dots,$h$ to $W$. Function $F(x)$ is executed whenever
             the value of $x$ has changed by moving the slider. }*/

panel_item int_item(string s, int& x, int l, int h, const window_handler& obj,const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

/*{\Mtext
\bigskip
{\bf 3.9.4 String Menu Items}
}*/


panel_item string_item(string s, string& x, const std::list<string>& L, const char* hlp=0);
/*{\Mopl     adds a string item with label $s$, variable $x$, and menu $L$ 
	     to $W$.}*/


panel_item string_item(string s, string& x, const std::list<string>& L, 
                                            void (*F)(char*),const char* hlp=0);
					    
panel_item string_item(string s, string& x, const std::list<string>& L, 
                                            const window_handler& obj,const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

panel_item string_item(string s, string& x, const std::list<string>& L, int sz,const char* hlp=0);
/*{\Mopl  menu |L| is displayed in a scroll box of height |sz|. }*/


panel_item string_item(string s, string& x, const std::list<string>& L, int sz,
                                            void (*F)(char*),const char* hlp=0);
/*{\Mopl     as above with action function |F|.}*/

panel_item string_item(string s, string& x, const std::list<string>& L, int sz,
                                            const window_handler& obj,const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

void set_menu(panel_item it, const std::list<string>& L, int sz=0);
/*{\Mopl     replaces the menu of string menu item |it| by a menu for 
             list |L| (table style if $sz = 0$ and scroll box with
             $sz$ entries otherwise).}*/

// synonym
void add_menu(panel_item it, const std::list<string>& L, int sz=0) 
{ set_menu(it,L,sz); }




/*{\Mtext
\bigskip
{\bf 3.9.5 Choice Items}
}*/

panel_item choice_item(string s, int& x, const std::list<string>& L, 
                       void (*F)(int)=0, const char* hlp=0);
/*{\Mopl     adds an integer item with label $s$, variable $x$, and choices 
	     from $L$ to $W$.}*/
	     
panel_item choice_item(string s, int& x, const std::list<string>& L, 
                       const window_handler& obj, const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

panel_item  choice_item(string label, int& x, const char* hlp, int n, ... );

panel_item choice_item(string, int& x, const char*, const char*);
panel_item choice_item(string, int& x, const char*, const char*, const char*);
panel_item choice_item(string, int& x, const char*, const char*, const char*, 
                                       const char*);
panel_item choice_item(string, int& x, const char*, const char*, const char*, 
                                       const char*, const char*);
panel_item choice_item(string, int& x, const char*, const char*, const char*, 
                                       const char*, const char*, const char*);
panel_item choice_item(string, int& x, const char*, const char*, const char*, 
                                       const char*, const char*, const char*, 
                                       const char*);
panel_item choice_item(string, int& x, const char*, const char*, const char*, 
                                       const char*, const char*, const char*, 
                                       const char*, const char*);

/*{\Moptions nextwarning=no}*/
/*
panel_item choice_item(string s, int& x, string s1,...,string sk);
*/
/*{\Mopl     adds an integer item with label $s$, variable $x$, and choices 
	     $s_1$, \dots, $s_k$ to $W$ ($k \le 8$).}*/


panel_item choice_item(string s, int& x, int n, int w, int h,unsigned char** bm,
                                                             const char* hlp=0);
/*{\Mopl     adds an integer item with label $s$, variable $x$, and $n$
             bitmaps |bm[0]|, \dots, |bm[n-1]| each of width $w$ 

             and height $h$.}*/ 


panel_item choice_item(string s, int& x, int n, int w, int h,unsigned char** bm,
                                                             void (*F)(int), 
                                                             const char* hlp=0);
/*{\Mopl  }*/


panel_item choice_item(string s, int& x, int n, int w, int h,unsigned char** bm,
                                                             const window_handler& obj, 
                                                             const char* hlp=0);
/*{\Mopl     as above with handler object |obj|.}*/

/*{\Mtext
\bigskip
{\bf 3.9.6 Multiple Choice Items}
}*/

panel_item choice_mult_item(string s,int& x,const std::list<string>& L, const char* hlp=0);
/*{\Mopl  }*/

panel_item choice_mult_item(string s,int& x,const std::list<string>& L,void (*F)(int), 
                                                            const char* hlp=0);
/*{\Mopl  }*/

panel_item choice_mult_item(string s,int& x,const std::list<string>& L,const window_handler& obj, 
                                                            const char* hlp=0);
/*{\Mopl  }*/

panel_item choice_mult_item(string s, int& x, int n, int w, int h,
                                                            unsigned char** bm,
                                                            const char* hlp=0);
/*{\Mopl  }*/

panel_item choice_mult_item(string s, int& x, int n, int w, int h, 
                                                            unsigned char** bm, 
                                                            void (*F)(int), 
                                                            const char* hlp=0);
/*{\Mopl  }*/

panel_item choice_mult_item(string s, int& x, int n, int w, int h, 
                                                            unsigned char** bm, 
                                                            const window_handler& obj, 
                                                            const char* hlp=0);
/*{\Mopl  }*/


/*{\Mtext
\bigskip
{\bf 3.9.7 Buttons}

The first occurence of character |'&'| in a button label makes the
following character |c| an {\em accelerator character}, i.e., the
button can be selected by typing ALT-|c| from the keyboard. 
}*/

int  button(string s, int n, const char* hlp=0);
/*{\Mop      adds a button with label $s$ and number $n$ to $W$.}*/

int  fbutton(string s, int n, const char* hlp=0);
/*{\Mop     as above but makes this button the focus button of $W$, i.e., 
            this button can be selected by pressing the return key. }*/

int  button(string s, const char* hlp=0);
/*{\Mop   adds a new button to $W$ with label $s$ and number equal to its
          position in the list of all buttons (starting with $0$).}*/

int  fbutton(string s, const char* hlp=0);
/*{\Mop     as above but makes this button the focus button. }*/ 


int  button(int w, int h, unsigned char* bm, string s,int n,const char* hlp=0);
/*{\Mop    adds a button with bitmap $bm$, label $s$, and number $n$ to $W$.}*/

int  button(char* pr1, char* pr2, string s, int n, const char* hlp=0);
/*{\Mop    adds a button with pixrects $pr1$ and $pr2$, label $s$, and number 
           $n$ to $W$. }*/

int  button(int w, int h, unsigned char* bm, string s,const char* hlp=0);
/*{\Mop   adds a new button to $W$ with bitmap $bm$, label $s$, and number equal
          to its position in the list of all buttons (starting with $0$).}*/

// leads to ambiguity in current release 
// int  button(char* pr1, char* pr2, string s, const char* hlp=0);
/*{\Xopl  as above, but with pixrects $pr1$ and $pr2$. }*/


int  button(string s, int n, void (*F)(int), const char* hlp=0);
/*{\Mopl  adds a button with label $s$, number $n$ and action 
          function $F$ to $W$. Function $F$ is called with actual
          parameter $n$  whenever the button is pressed. }*/
	  
int  button(string s, int n, const window_handler& obj, const char* hlp=0);
/*{\Mopl  as above with handler object |obj|.}*/

int  fbutton(string s, int n, void (*F)(int), const char* hlp=0);
/*{\Mop     as above but makes this button the focus button. }*/ 

int  fbutton(string s, int n, const window_handler& obj, const char* hlp=0);
/*{\Mopl  as above with handler object |obj|.}*/

int  button(int w, int h, unsigned char* bm, string s, int n, void (*F)(int), 
                                                       const char* hlp=0);
/*{\Mopl  adds a button with bitmap $bm$, label $s$, number $n$ and action 
          function $F$ to $W$. Function $F$ is called with actual
          parameter $n$  whenever the button is pressed. }*/
	  
int  button(int w, int h, unsigned char* bm, string s, int n, const window_handler& obj, 
                                                       const char* hlp=0);
/*{\Mopl }*/

int  button(char* pr1, char* pr2, string s, int n, void (*F)(int), const char* hlp=0);
/*{\Mopl  as above, but with pixrect $pr1$ and $pr2$. }*/

int  button(char* pr1, char* pr2, string s, int n, const window_handler& obj, const char* hlp=0);
/*{\Mopl }*/

int  button(string s, void (*F)(int), const char* hlp=0);
/*{\Mopl  adds a button with label $s$, number equal to its rank and action 
          function $F$ to $W$. Function $F$ is called with the value of the
          button as argument whenever the button is pressed. }*/
	  
int  button(string s, const window_handler& obj, const char* hlp=0);
/*{\Mopl }*/

int  button(int w, int h, unsigned char* bm, string s, void (*F)(int),
                                                       const char* hlp=0);
/*{\Mopl  adds a button with bitmap $bm$, label $s$, number equal to its rank 
          and action function $F$ to $W$. Function $F$ is called with the 
          value of the button as argument whenever the button is pressed. }*/
	  
int  button(int w, int h, unsigned char* bm, string s, const window_handler& obj,
                                                       const char* hlp=0);
/*{\Mopl }*/

int  button(char* pr1, char* pr2, string s, void (*F)(int),const char* hlp=0);
/*{\Mopl  as above, but with pixrect $pr1$ and $pr2$. }*/

int  button(char* pr1, char* pr2, string s, const window_handler& obj,const char* hlp=0);
/*{\Mopl }*/

int  button(string s, int n, window& M, const char* hlp=0);
/*{\Mopl  adds a button with label $s$, number $n$  and attached sub-window
          (menu) $M$ to $W$. Window $M$ is opened whenever the button is
          pressed. }*/

int  button(int w, int h, unsigned char* bm, string s, int n, window& M,
                                                     const char* hlp=0);
/*{\Mopl  adds a button with bitmap $bm$, label $s$, number $n$  and attached 
          sub-window (menu) $M$ to $W$. Window $M$ is opened whenever the 
          button is pressed. }*/

int  button(char* pr1, char* pr2, string s, int n, window& M,const char* hlp=0);
/*{\Mopl  as above, but with pixrect $pr1$ and $pr2$. }*/


int  button(string s, window& M, const char* hlp=0);
/*{\Mopl  adds a button with label $s$ and attached sub-window $M$ to $W$. 
          The number returned by $read\_mouse$ is the number of the
          button selected in sub-window $M$. }*/

int  button(int w, int h, unsigned char* bm, string s, window& M, 
                                                       const char* hlp=0);
/*{\Mopl  adds a button with bitmap $bm$, label $s$ and attached sub-window 
          $M$ to $W$. The number returned by $read\_mouse$ is the number of 
          the button selected in sub-window $M$. }*/

int  button(char* pr1, char* pr2, string s, window& M, const char* hlp=0);
/*{\Mopl  as above, but with pixrect $pr1$ and $pr2$. }*/



void make_menu_bar(int k) { BASE_WINDOW::make_menu_bar(k); }

void make_menu_bar() { BASE_WINDOW::make_menu_bar(); }
/*{\Mop  inserts a menu bar at the top of the panel section that 
         contains all previously added menu buttons (buttons with 
         an attached subwindow). }*/ 


int menu_button(string s, int n, window& M, const char* hlp=0);
/*{\Xopl  adds a {\em menu} button with label $s$, number $n$  and 
          attached sub-window (menu) $M$ to the panel bar of $W$. }*/


int menu_button(string s, window& M, const char* hlp=0);
/*{\Xopl  adds a {\em menu} button with label $s$ and attached sub-window $M$ 
          to $W$. The number returned by $read\_mouse$ is the number of the
          button selected in sub-window the panel bar of $M$. }*/



static window* get_call_window() 
{ return (window*)BASE_WINDOW::get_call_window(); } 
/*{\Mstatic  A static function that can be called in action functions attached 
             to panel items or buttons to retrieve a pointer to the window 
             containing the corresponding item or button. }*/

static panel_item get_call_item()   { return BASE_WINDOW::get_call_item(); } 
/*{\Mstatic  A static function that can be called in action functions attached
             to panel items to retrieve the corresponding item. }*/


static int   get_call_button() { return BASE_WINDOW::get_call_button(); } 
/*{\Mstatic  A static function that can be called in action functions attached
             to panel buttons to retrieve the number of the corresponding 
             button. }*/



/*{\Mtext
\bigskip
{\bf 3.9.8. Manipulating Panel Items and Buttons}
}*/


/*{\Mtext
{\bf Disabling and Enabling Items or buttons}
}*/

void disable_item(panel_item it) { BASE_WINDOW::disable_item(it); }
/*{\Mopl  disables panel item $it$. }*/

void enable_item(panel_item it)  { BASE_WINDOW::enable_item(it); }
/*{\Mopl  enables panel item $it$. }*/

bool is_enabled(panel_item it)   { return BASE_WINDOW::is_enabled(it); }
/*{\Mopl  tests whether item $it$ is enabled or not. }*/


void disable_button(int b) { disable_item(get_button_item(b)); }
/*{\Mopl  disables button $b$. }*/

void enable_button(int b)  { enable_item(get_button_item(b)); }
/*{\Mopl  enables button $b$. }*/

bool is_enabled(int b)   { return is_enabled(get_button_item(b)); }
/*{\Mopl  tests whether button $b$ is enabled or not. }*/


void disable_panel(bool disable_every_item=true) 
{ BASE_WINDOW::disable_panel(disable_every_item); }
/*{\Mopl  disables the entire panel section of $W$. }*/

void enable_panel()  { BASE_WINDOW::enable_panel(); }
/*{\Mopl  enables the entire panel section of $W$. }*/



/*{\Mtext
{\bf Accessing and Updating Item Data}
}*/

void set_text(panel_item it, string s);
/*{\Mopl  replaces the text of text item |it| by $s$. }*/


panel_item  get_item(string s) { return BASE_WINDOW::get_item(s.c_str()); }
/*{\Mopl  returns the item with label $s$ and |NULL| if no such
          item exists in $W$. }*/


int  get_button(string s) { return BASE_WINDOW::get_button(s.c_str()); }
/*{\Mopl  returns the button with label $s$ and $-1$ if no such
          button exists in $W$. }*/


string  get_button_label(int but) 
{ return BASE_WINDOW::get_button_label(but); }
/*{\Mopl  returns the label of button $but$. }*/

void    set_button_label(int but, string s) 
{ BASE_WINDOW::set_button_label(but,s.c_str()); }
/*{\Mopl  sets the label of button $but$ to $s$. }*/

void    set_button_pixrects(int but, char *pr1, char* pr2)
{ BASE_WINDOW::set_button_pixrects(but,pr1,pr2); }
/*{\Mopl  sets the pixrects of button $but$ to $pr1$ and $pr2$. }*/


window* get_window(panel_item it)
{ return (window*)BASE_WINDOW::get_window(it); }


window* get_window(int but);
/*{\Mopl     returns a pointer to the subwindow attached to button $but$ 
             (|NULL| if $but$ has no subwindow) }*/

window* set_window(int but, window* M);
/*{\Mopl     associates subwindow (menu) $*M$ with button $but$. Returns 
             a pointer to the window previously attached to $but$. }*/

void set_function(int but, void (*F)(int) )
{ BASE_WINDOW::set_item_func(get_button_item(but),F); }
/*{\Mopl     assign action function $F$ to button $but$. }*/

void set_object(int but, const window_handler& obj)
{ BASE_WINDOW::set_item_object(get_button_item(but),obj); }
/*{\Mopl     assign handler object $obj$ to button $but$. }*/


void redraw_panel() { BASE_WINDOW::draw_frame(); }
/*{\Mopl     redraw the panel area of $W$. }*/

void redraw_panel(panel_item it) { BASE_WINDOW::redraw_panel(it); }
/*{\Mopl     redraw item $i$ in the panel area of $W$. }*/



void display_help_text(string fname);
/*{\Mop  displays the help text contained in |name.hlp|. The file 
         |name.hlp| must exist either in the current working directory or in
         |\$LEDAROOT/incl/Help|. }*/



void string_edit(window_point p, string& s);

int string_edit(double x, double y, string& s, int& curs_pos);
int string_edit(double x, double y, string& s, int& curs_pos, int& timeout);


/*
// status_window
void    create_status_window()      { BASE_WINDOW::create_status_window(); }
void    open_status_window()        { BASE_WINDOW::open_status_window(); }
void    close_status_window()       { BASE_WINDOW::close_status_window(); }
void    destroy_status_window()     { BASE_WINDOW::destroy_status_window(); }
void    set_status_string(string s) { BASE_WINDOW::set_satus_string(s); }
window* get_status_window()         { return BASE_WINDOW::get_status_window(); }
*/



// zooming

bool read_zoom_rect(const window_point& center,window_point& p, window_point& q );
bool read_zoom_rect(window_point& p, window_point& q);

void adjust_zoom_rect(double& x0, double& y0, double& x1, double& y1);

void zoom_area(double x0, double y0, double x1, double y1, int steps = 16, 
                                     void (*redraw)(window*) = 0);

void zoom_area(const window_point& p, const window_point& q, int steps = 16, 
                                     void (*redraw)(window*) = 0)
{ zoom_area(p.x(),p.y(),q.x(),q.y(),steps,redraw); }

void zoom(double f, int steps = 16, void(*redraw)(window*) = 0);

void scroll_window(const window_point& p, 
                   void (*redraw_func)(window*,double,double,double,double)=0);


};

inline void put_back_event()      { BASE_WINDOW::put_back_event(); }


}




/*{\Moptions usesubscripts=no}*/

#include <CGAL/LEDA/panel.h>
#include <CGAL/LEDA/menu.h>


#if defined(WINMAIN)
#if !defined(_WINDOWS_)
extern "C" int   __stdcall WinMain(void*,void*,char*,int);
int main_tmp();
int __stdcall WinMain(void*,void*,char*,int) { return main_tmp(); }
#define main  main_tmp
#endif
#endif



#endif
