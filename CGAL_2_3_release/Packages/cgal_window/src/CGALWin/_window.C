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
// file          : src/CGALWin/_window.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

#include <CGAL/LEDA/basic.h>
#include <CGAL/LEDA/window.h>
#include <CGAL/LEDA/string_manip.h>
#include <CGAL/LEDA/file.h>

#include <cassert>
#include <cstdio>
#include <cstdarg>
#include <string>
#include <cctype>

#if defined(__BORLANDC__)
using std::isspace;
using std::isprint;
#endif

//------------------------------------------------------------------------------
// LEDA WINDOW BASICS
//------------------------------------------------------------------------------



namespace CGAL {


color window::fgcol = DEF_COLOR;
color window::bgcol = DEF_COLOR;


window::window() : BASE_WINDOW(LEDA::copyright_string)
{ 
  std_buttons();
  status_win = 0;
  point_buffer = 0;
  normalize_rat = false;
 }

window::window(int w, int h) : BASE_WINDOW(w,h,LEDA::copyright_string)
{  
  std_buttons();
  status_win = 0;
  point_buffer = 0;
  normalize_rat = false;
 }

window::window(int w, int h, const char* label) : BASE_WINDOW(w,h,label)
{ 
  std_buttons();
  status_win = 0;
  point_buffer = 0;
  normalize_rat = false;
 }

window::window(const char* label) : BASE_WINDOW(label)
{ 
  std_buttons();
  status_win = 0;
  point_buffer = 0;
  normalize_rat = false;
 }



//------------------------------------------------------------------------------
// WINDOW OUTPUT
//------------------------------------------------------------------------------


// pixels

void window::draw_pix(double x, double y, color c ) 
{ BASE_WINDOW::draw_pix(x,y,c); }



void window::draw_pixels(int n, double* xcoord, double* ycoord, color c ) 
{ BASE_WINDOW::draw_pixels(n,xcoord,ycoord,c); }





void window::draw_text(double x, double y, string s, color c)
{ BASE_WINDOW::draw_text(x,y,s.c_str(),c); }

void window::draw_text(const window_point& p, string s, color c)
{ draw_text(p.x(),p.y(),s.c_str(),c); }

void window::draw_ctext(double x, double y, string s, color c)
{ BASE_WINDOW::draw_ctext(x,y,s.c_str(),c); }

void window::draw_ctext(const window_point& p, string s, color c)
{ draw_ctext(p.x(),p.y(),s.c_str(),c); }

void window::draw_ctext(string s, color c)
{ BASE_WINDOW::draw_ctext(s.c_str(),c); }

double window::text_box(double x0, double x1, double y0, string s, bool draw)
{ return BASE_WINDOW::text_box(x0,x1,y0,s.c_str(),draw); }

void window::text_box(string s)
{ double d = pix_to_real(1);
  BASE_WINDOW::text_box(xmin()+d,xmax()-d,ymax()-d,s.c_str(),true); 
 }

// points

void window::draw_point(double x0,double y0,color c)
{ BASE_WINDOW::draw_point(x0,y0,c); }


// segments

void window::draw_segment(double x1, double y1, double x2, double y2, color c )
{ BASE_WINDOW::draw_segment(x1,y1,x2,y2,c); }



// lines

void window::draw_line(double x1, double y1, double x2, double y2, color c )
{ BASE_WINDOW::draw_line(x1,y1,x2,y2,c); }


void window::draw_hline(double y, color c )
{ BASE_WINDOW::draw_segment(xmin(),y,xmax(),y,c); }

void window::draw_vline(double x, color c )
{ BASE_WINDOW::draw_segment(x,ymin(),x,ymax(),c); }


// rays

void window::draw_ray(double x1, double y1, double x2, double y2, color c )
{ BASE_WINDOW::draw_ray(x1,y1,x2,y2,c); }




// nodes

void window::draw_node(double x0,double y0,color c) 
{ BASE_WINDOW::draw_node(x0,y0,c); }

void window::draw_node(const window_point& p, color c)
{ window::draw_node(p.x(),p.y(),c); }

void window::draw_filled_node(double x0,double y0,color c)
{ BASE_WINDOW::draw_filled_node(x0,y0,c); }

void window::draw_filled_node(const window_point& p, color c)
{ window::draw_filled_node(p.x(),p.y(),c); }

void window::draw_text_node(double x,double y,string s,color c)
{ BASE_WINDOW::draw_text_node(x,y,s.c_str(),c); }

void window::draw_text_node(const window_point& p ,string s,color c)
{ window::draw_text_node(p.x(),p.y(),s.c_str(),c); }

void window::draw_int_node(double x,double y,int i,color c)
{ BASE_WINDOW::draw_int_node(x,y,i,c); }

void window::draw_int_node(const window_point& p ,int i,color c)
{ window::draw_int_node(p.x(),p.y(),i,c); }


//circles

void window::draw_circle(double x,double y,double r,color c)
{ BASE_WINDOW::draw_circle(x,y,r,c); }

void window::draw_circle(const window_point& p,double r,color c)
{ BASE_WINDOW::draw_circle(p.x(),p.y(),r,c); }


// discs

void window::draw_disc(double x,double y,double r,color c)
{ BASE_WINDOW::draw_filled_circle(x,y,r,c); }

void window::draw_disc(const window_point& p,double r,color c)
{ window::draw_disc(p.x(),p.y(),r,c); }



//ellipses

void window::draw_ellipse(double x,double y,double r1,double r2,color c)
{ BASE_WINDOW::draw_ellipse(x,y,r1,r2,c); }

void window::draw_ellipse(const window_point& p, double r1, double r2, color c)
{ BASE_WINDOW::draw_ellipse(p.x(),p.y(),r1,r2,c); }

void window::draw_filled_ellipse(double x,double y,double r1,double r2,color c)
{ BASE_WINDOW::draw_filled_ellipse(x,y,r1,r2,c); }

void window::draw_filled_ellipse(const window_point& p, double r1, double r2, color c)
{ BASE_WINDOW::draw_filled_ellipse(p.x(),p.y(),r1,r2,c); }



// segment.angle() ; p1 - source, p2 target point ...

static double seg_angle(const window_point& p1, const window_point& p2)
{ 
 double dx = p2.x() - p1.x();
 double dy = p2.y() - p1.y();
 
 if (dx != 0 || dy != 0)  
      return CGAL_CLIB_STD::atan2(dy,dx); 
 else
      return 0;
}

// segment.y_proj and x_proj replaced ...

static double x_proj(const window_point& p1, const window_point& p2, double y)
{
  double dx  = p2.x() - p1.x(); 
  double dy  = p2.y() - p1.y(); 
  
  //compute slope value ...
  double slope;
  
  if (dx==0) slope = MAXDOUBLE;
  else slope = dy/dx;
  
  if (p1.x() == p2.x())  // vertical case ...
       return  p1.x();
  else
       return  p1.x() - (p1.y() - y)/slope;   
}

static double y_proj(const window_point& p1, const window_point& p2, double x)
{
  double dx  = p2.x() - p1.x(); 
  double dy  = p2.y() - p1.y(); 
  
  //compute slope value ...
  double slope;
  
  if (dx==0) slope = MAXDOUBLE;
  else slope = dy/dx;
  
  return  p1.y() - slope * (p1.x() - x);
}

// midpoint replaced ...

static window_point midpoint(const window_point& p1, const window_point& p2)
{
  double xm = (p1.x() + p2.x())/2;
  double ym = (p1.y() + p2.y())/2;
  
  return window_point(xm,ym);
}

// orientation replaced
static int compute_orientation(const window_point& a, const window_point& b, const window_point& c)
{
  double d1 = (a.xcoord() - b.xcoord()) * (a.ycoord() - c.ycoord());
  double d2 = (a.ycoord() - b.ycoord()) * (a.xcoord() - c.xcoord()); 
  if (d1 == d2) return 0; else return (d1 > d2) ? +1 : -1;
}


// polygons

void window::draw_polyline(const std::list<window_point>& lp, int arrow, double d, color c)
{ 
  int n = lp.size();

  if (n < 2) 
  { 
    return;
   }

  if ((arrow & 4) || (arrow & 8))
  { int m = n/2;
    std::list<window_point> L = lp;
    
    std::list<window_point>::iterator it = L.begin();
    
    while (--m) it++;
    it++; //one more now ...
    
    std::list<window_point> lp1,lp2;
    std::copy(L.begin(), it, std::back_inserter(lp1)); // first part
    std::copy(it, L.end(), std::back_inserter(lp2));   // second part
    
    // was L.split(it,lp1,lp2,LEDA::after);
    
    window_point s = midpoint(lp1.back(),lp2.front());
    lp1.reverse();
    
    lp1.push_front(s);
    lp2.push_front(s);
    
    int arr1 = (arrow & 4) ? 1 : 0;
    int arr2 = (arrow & 8) ? 1 : 0;
    draw_polyline(lp1,arr1,d,c);
    draw_polyline(lp2,arr2,d,c);
    return;
  }

  if (arrow == 2)
  { std::list<window_point> lpr = lp;
    lpr.reverse();
    draw_polyline(lpr,1,d,c);
    return;
   }

  double* X = new double[n];
  double* Y = new double[n];
  n = 0;
  window_point p;
  
  std::list<window_point>::const_iterator lit = lp.begin();
  
  for(;lit!=lp.end();lit++) 
  { 
    p = *lit;
    X[n] = p.x();
    Y[n] = p.y();
    n++;
   }

  BASE_WINDOW::adjust_polyline(n,X,Y);

  if (arrow & 1)
  { 
    window_point pa(X[1],Y[1]);
    window_point pb(X[0],Y[0]);
    window_point q = draw_arrow_head(pb,seg_angle(pa,pb),d,c);
    X[0] = q.x();
    Y[0] = q.y();
   }

  if (arrow & 2)
  { 
    window_point pa(X[n-2],Y[n-2]);
    window_point pb(X[n-1],Y[n-1]);
    window_point q = draw_arrow_head(pb,seg_angle(pa,pb),d,c);
    X[n-1] = q.x();
    Y[n-1] = q.y();
   }

  BASE_WINDOW::draw_polyline(n,X,Y,c);

  delete[] X;
  delete[] Y;    
}

void window::draw_polyline(const std::list<window_point>& lp, color c)
{ draw_polyline(lp,0,0,c); }

void window::draw_polyline_arrow(const std::list<window_point>& lp, color c)
{ double d = pix_to_real(get_line_width());
  draw_polyline(lp,2,d,c); 
 }


void window::draw_polygon(const std::list<window_point>& lp, color c)
{ int n = lp.size();
  double* X = new double[n];
  double* Y = new double[n];
  n = 0;
  window_point p;
  std::list<window_point>::const_iterator it = lp.begin();
   
  for(;it != lp.end();it++) 
  { p = *it;
    X[n] = p.x();
    Y[n] = p.y();
    n++;
   }
  BASE_WINDOW::draw_polygon(n,X,Y,c);
  delete[] X;
  delete[] Y;
}

void window::draw_filled_polygon(const std::list<window_point>& lp, color c)
{ int n = lp.size();
  double* X = new double[n];
  double* Y = new double[n];
  n = 0;
  window_point p;
  std::list<window_point>::const_iterator it = lp.begin();
  
  for(;it != lp.end();it++) 
  { p = *it;
    X[n] = p.x();
    Y[n] = p.y();
    n++;
   }
  BASE_WINDOW::draw_filled_polygon(n,X,Y,c);
  delete[] X;
  delete[] Y;
}


void window::draw_oriented_polygon(const std::list<window_point>& lp, color c)
{ draw_polygon(lp,c);

  std::list<window_point>::const_iterator it = lp.begin();
    
  window_point p = *it;
  it++;
  window_point q = *it;
  
  draw_arrow(p,midpoint(p,q));
}



void window::draw_rectangle(double x0,double y0,double x1,double y1, color col)
{ BASE_WINDOW::draw_rectangle(x0,y0,x1,y1,col); }

void window::draw_rectangle(window_point p, window_point q, color col)
{ BASE_WINDOW::draw_rectangle(p.x(),p.y(),q.x(),q.y(),col); }


void window::draw_filled_rectangle(double x0,double y0,double x1,double y1, color col)
{ BASE_WINDOW::draw_filled_rectangle(x0,y0,x1,y1,col); }

void window::draw_filled_rectangle(window_point p, window_point q, color col)
{ BASE_WINDOW::draw_filled_rectangle(p.x(),p.y(),q.x(),q.y(),col); }



void window::draw_triangle(window_point a, window_point b, window_point c, color col)
{ std::list<window_point> poly;
  poly.push_back(a);
  poly.push_back(b);
  poly.push_back(c);
  draw_polygon(poly,col);
}


void window::draw_filled_triangle(window_point a, window_point b, window_point c, color col)
{ std::list<window_point> poly;
  poly.push_back(a);
  poly.push_back(b);
  poly.push_back(c);
  draw_filled_polygon(poly,col);
}
   

// functions

void window::plot_xy(double x0, double x1, win_draw_func f, color c)
{ BASE_WINDOW::plot_xy(x0,x1,f,c); }

void window::plot_yx(double y0, double y1, win_draw_func f, color c)
{ BASE_WINDOW::plot_yx(y0,y1,f,c); }



// arrows
// translate_by_angle

static window_point translate_point_by_angle(const window_point& p, double phi, double d)
{
  double dx = CGAL_CLIB_STD::cos(phi) * d;
  double dy = CGAL_CLIB_STD::sin(phi) * d;
  if (CGAL_CLIB_STD::fabs(dx) < 1e-12) dx = 0; 
  if (CGAL_CLIB_STD::fabs(dy) < 1e-12) dy = 0; 
  return window_point(p.x()+dx,p.y()+dy);
}


window_point window::arrow_head(const window_point& q,double a,double d,double* X,double* Y)
{ 
  double alpha = a-LEDA_PI; 

  window_point l = translate_point_by_angle(q, alpha+LEDA_PI/6, 3.5*d);
  window_point m = translate_point_by_angle(q, alpha,           2.0*d);
  window_point r = translate_point_by_angle(q, alpha-LEDA_PI/6, 3.5*d);

  X[0] = q.x();
  Y[0] = q.y();
  X[1] = l.x();
  Y[1] = l.y();
  X[2] = m.x();
  Y[2] = m.y();
  X[3] = r.x();
  Y[3] = r.y();

  //return m;
  return translate_point_by_angle(q,alpha,d);
}

window_point window::draw_arrow_head(const window_point& q, double a, double d, color c)
{ double X[4];
  double Y[4];
  double pix = pix_to_real(1);

  if (d < pix)
     d = 1.5*d + pix;
  else
     d =  d + 1.5*pix;

  window_point m = arrow_head(q,a,d,X,Y);
  BASE_WINDOW::draw_filled_polygon(4,X,Y,c);
  BASE_WINDOW::draw_pix(q.x(),q.y(),c);

  return m;
}


window_point window::draw_arrow_head(const window_point& q, double a, color c)
{ double d = pix_to_real(get_line_width());
  return draw_arrow_head(q,a,d,c);
}



void window::draw_arrow(const window_point& p, const window_point& q, color c)
{ 
  window_point hp = draw_arrow_head(q,seg_angle(p,q),c);
  draw_segment(p.x(),p.y(),hp.x(),hp.y(),c);
}

void window::draw_arrow(double x0, double y0, double x1, double y1, color c)
{ draw_arrow(window_point(x0,y0),window_point(x1,y1),c); }




// edges

void window::draw_edge(double x0, double y0, double x1, double y1, color c)
{ BASE_WINDOW::draw_edge(x0,y0,x1,y1,c); }

void window::draw_edge(const window_point& p, const window_point& q, color c)
{ draw_edge(p.x(),p.y(),q.x(),q.y(),c); }


void window::draw_edge_arrow(double x0,double y0,double x1,double y1,color c)
{ draw_edge_arrow(window_point(x0,y0), window_point(x1,y1),c); }

void window::draw_edge_arrow(const window_point& p, const window_point& q, color c)
{ double  alpha = seg_angle(p,q);
  window_point hp = translate_point_by_angle(p, alpha,get_node_width()/scale());
  window_point hq = translate_point_by_angle(q, alpha,-get_node_width()/scale());
  draw_arrow(hp,hq,c);
}


//------------------------------------------------------------------------------
// WINDOW INPUT
//------------------------------------------------------------------------------


int window::get_mouse() 
{ return BASE_WINDOW::get_mouse(); }

int window::get_mouse(double& x, double& y) 
{ return BASE_WINDOW::get_mouse(x,y); }

int  window::get_mouse(window_point& q)
{ double x,y;
  int but = BASE_WINDOW::get_mouse(x,y);
  q = window_point(x,y);
  return but;
 }


int window::read_mouse() 
{ return BASE_WINDOW::read_mouse(); }

int  window::read_mouse(double& x, double& y)
{ return BASE_WINDOW::read_mouse(0,0.0,0.0,x,y); }

int  window::read_mouse(window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse(0,0.0,0.0,x,y);
  q = window_point(x,y);
  return but;
 }


int  window::read_mouse_seg(double x0, double y0, double& x, double& y)
{ return BASE_WINDOW::read_mouse(1,x0,y0,x,y); }

int  window::read_mouse_seg(const window_point& p, window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse(1,p.x(),p.y(),x,y);
  q = window_point(x,y);
  return but;
}


int  window::read_mouse_ray(double x0, double y0, double& x, double& y)
{ return BASE_WINDOW::read_mouse(2,x0,y0,x,y); }

int  window::read_mouse_ray(const window_point& p, window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse(2,p.x(),p.y(),x,y);
  q = window_point(x,y);
  return but;
}


int  window::read_mouse_line(double x0, double y0, double& x, double& y)
{ return BASE_WINDOW::read_mouse(3,x0,y0,x,y); }

int  window::read_mouse_line(const window_point& p, window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse(3,p.x(),p.y(),x,y);
  q = window_point(x,y);
  return but;
}


int  window::read_mouse_circle(double x0, double y0, double& x, double& y)
{ return BASE_WINDOW::read_mouse(4,x0,y0,x,y); }

int  window::read_mouse_circle(const window_point& p, window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse(4,p.x(),p.y(),x,y);
  q = window_point(x,y);
  return but;
}


int  window::read_mouse_rect(double x0, double y0, double& x, double& y)
{ return BASE_WINDOW::read_mouse(5,x0,y0,x,y); }

int  window::read_mouse_rect(const window_point& p, window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse(5,p.x(),p.y(),x,y);
  q = window_point(x,y);
  return but;
}



int window::read_mouse_action(mouse_action_func f, double& x, double& y)
{ return BASE_WINDOW::read_mouse_action(f,0.0,0.0,x,y); }

int window::read_mouse_action(mouse_action_func f, window_point& q)
{ double x,y;
  int but = BASE_WINDOW::read_mouse_action(f,0.0,0.0,x,y);
  q = window_point(x,y);
  return but;
}



window& window::draw(const window_point& p, color c)
{ draw_point(p.x(),p.y(),c); 
  return *this;
}

window& window::draw(const window_point& ps, const window_point& pt, color c)
{ 
  draw_segment(ps.x(),ps.y(),pt.x(),pt.y(),c);
  return *this;
}


window& window::draw(const window_point& m, double rad, color c)
{ draw_disc(m.x(),m.y(),rad,get_fill_color());
  draw_circle(m.x(),m.y(),rad,c); 
  return *this;
}




void window::set_point_buffer(const window_point p)
{ if (point_buffer) delete point_buffer;
  point_buffer = new window_point(p);
 }


window& window::read(window_point& p)
{ 
  if (point_buffer)
  { p = *point_buffer;
    delete point_buffer;
    point_buffer = 0;
    return *this;
   }

  double x,y;
  state = 1;
  int save_but[8];
  std_buttons(save_but);
  if (read_mouse(x,y) == MOUSE_BUTTON(1)) 
     p = window_point(x,y);
  else
     state = 0;
  set_buttons(save_but);
  return *this;
 }

// reading segments ...
window& window::read(window_point& ps, window_point& pt)
{ double x,y;
  window_point p;
  int key = 0;
  state = 1;

  int save_but[8];
  std_buttons(save_but);

  if (!read(p).state) return *this;

  for(;;)
  { key = read_mouse_seg(p.x(),p.y(),x,y);
    if (key == MOUSE_BUTTON(3))  
     { state = 0;
       break; 
      }
    if (key == MOUSE_BUTTON(1) && !(shift_key_down() && read(p).state)) break; 
   }

  if (state) {
    ps = window_point(p.x(), p.y());
    pt = window_point(x,y);    
  }  

  set_buttons(save_but);
  return *this;
}

// reading circles ...
window& window::read(window_point& m, double& rad)
{ double x,y;
  window_point p;
  int key = 0;
  state = 1;

  if (!read(p).state) return *this;

  int save_but[8];
  std_buttons(save_but);

  while ((key=read_mouse_circle(p.x(),p.y(),x,y)) != MOUSE_BUTTON(1))
  { if (key == MOUSE_BUTTON(3))  
     { state = 0;
       break; 
      }
    if (key == MOUSE_BUTTON(1) && shift_key_down() &&!read(p).state) break;
   }

  if (state) 
  { double dx = x-p.x();
    double dy = y-p.y();
    m   = p;
    rad = sqrt(dx*dx+dy*dy);
   }

  set_buttons(save_but);
  return *this;
}


bool window::confirm(string s)
{ panel p;
  p.text_item("");
  p.text_item("\\bf\\blue " + s);
  p.text_item("");
  p.button("YES",1);
  p.button("NO",0);
  return p.open(*this) == 1;
}


void window::acknowledge(string s)
{ panel p;
  p.text_item("");
  p.text_item("\\bf\\blue " + s);
  p.text_item("");
  p.fbutton("OK");
  p.open(*this);
}

void window::notice(string s) { acknowledge(s); }


int  window::read_panel(string header, int n, string* L)
{ panel P("LEDA PANEL");
  P.text_item(header);
  for(int i = 0; i < n; i++) P.button(L[i]);
  return P.open(*this);
 }


int  window::read_vpanel(string header, int n, string* L)
{ panel P("LEDA PANEL");
  P.buttons_per_line(1);
  P.text_item(header);
  for(int i = 0; i < n; i++)  P.button(L[i]);
  return P.open(*this);
 }

void  window::panel_read(string prompt, string& x)
{ panel P("STRING INPUT PANEL");
  P.string_item(prompt,x);
  P.fbutton("continue");
  P.open(*this);
 }

void  window::panel_read(string prompt, int& x)
{ panel P("INT INPUT PANEL");
  P.int_item(prompt,x);
  P.fbutton("continue");
  P.open(*this);
}

void  window::panel_read(string prompt, double& x)
{ panel P("DOUBLE INPUT PANEL");
  P.real_item(prompt,x);
  P.fbutton("continue");
  P.open(*this);
 }



string  window::read_string(string prompt)
{ string s;
  panel_read(prompt,s);
  return s;
 }

int  window::read_int(string prompt)
{ int i = 0;
  panel_read(prompt,i);
  return i;
}

double  window::read_real(string prompt)
{ double x = 0;
  panel_read(prompt,x);
  return x;
 }



//------------------------------------------------------------------------------
//   PANEL OPERATIONS
//------------------------------------------------------------------------------


int window::button(string s, const char* hlp)   
{ panel_action_func f= 0;
  return BASE_WINDOW::button(s.c_str(),-1,f,hlp); }

int window::button(int w, int h, unsigned char* bm, string s, const char* hlp)
{ panel_action_func f= 0;
  return BASE_WINDOW::button(w,h,bm,s.c_str(),-1,f,hlp); }


int window::button(string s,int v, const char* hlp)  
{ panel_action_func f= 0;
  return BASE_WINDOW::button(s.c_str(),v,f,hlp); }

int window::button(int w, int h, unsigned char* bm, string s, int v, 
                                                              const char* hlp)
{ panel_action_func f= 0;
  return BASE_WINDOW::button(w,h,bm,s.c_str(),v,f,hlp); }


int window::button(string s, panel_action_func F, const char* hlp)   
{ return BASE_WINDOW::button(s.c_str(),-1,F,hlp); }

int window::button(string s, const window_handler& obj, const char* hlp)   
{ panel_action_func f= 0;
  int bt = BASE_WINDOW::button(s.c_str(),-1,f,hlp);
  BASE_WINDOW::set_action(bt,obj);
  return bt; 
}



int window::button(int w, int h, unsigned char* bm, string s, 
                                                    panel_action_func F, 
                                                    const char* hlp)   
{ return BASE_WINDOW::button(w,h,bm,s.c_str(),-1,F,hlp); }

int window::button(int w, int h, unsigned char* bm, string s, 
                                                    const window_handler& obj, 
                                                    const char* hlp)   
{ 
  panel_action_func f= 0;
  int bt = BASE_WINDOW::button(w,h,bm,s.c_str(),-1,f,hlp); 
  BASE_WINDOW::set_action(bt,obj);
  return bt;
}




int window::button(string s, int v, panel_action_func F, const char* hlp)   
{ return BASE_WINDOW::button(s.c_str(),v,F,hlp); }

int window::button(string s, int v, const window_handler& obj, const char* hlp)   
{ panel_action_func f= 0;
  int bt = BASE_WINDOW::button(s.c_str(),v,f,hlp); 
  BASE_WINDOW::set_action(bt,obj);
  return bt;
}


int window::button(int w, int h, unsigned char* bm, string s, int v, 
                                                         panel_action_func F,
                                                         const char* hlp)
{ return BASE_WINDOW::button(w,h,bm,s.c_str(),v,F,hlp); }

int window::button(int w, int h, unsigned char* bm, string s, int v, 
                                                         const window_handler& obj,
                                                         const char* hlp)
{ panel_action_func f= 0;
  int bt = BASE_WINDOW::button(w,h,bm,s.c_str(),v,f,hlp); 
  BASE_WINDOW::set_action(bt,obj);
  return bt;
}


int window::button(string s, window& p, const char* hlp)   
{ return BASE_WINDOW::button(s.c_str(),-1,(BASE_WINDOW*)&p,hlp); }

int window::button(int w, int h, unsigned char* bm, string s, window& p, 
                                                              const char* hlp)
{ return BASE_WINDOW::button(w,h,bm,s.c_str(),-1,(BASE_WINDOW*)&p,hlp); }


int window::button(string s, int val, window& p, const char* hlp)   
{ return BASE_WINDOW::button(s.c_str(),val,(BASE_WINDOW*)&p,hlp); }

int window::button(int w, int h, unsigned char* bm, string s, int val, 
                                                              window& p,
                                                              const char* hlp)
{ return BASE_WINDOW::button(w,h,bm,s.c_str(),val,(BASE_WINDOW*)&p,hlp); }


int window::menu_button(string s, int val, window& p, const char* hlp)   
{ return BASE_WINDOW::menu_button(s.c_str(),val,(BASE_WINDOW*)&p,hlp); }

int window::menu_button(string s, window& p, const char* hlp)   
{ return BASE_WINDOW::menu_button(s.c_str(),-1,(BASE_WINDOW*)&p,hlp); }


// pixmaps

/*
int window::button(char* pr0, char* pr1, string s, const char* hlp)
{ panel_action_func F=0;
  return BASE_WINDOW::button(pr0,pr1,s,-1,F,hlp); }
*/

int window::button(char* pr0, char* pr1, string s, int v, const char* hlp)
{ panel_action_func F=0;
  return BASE_WINDOW::button(pr0,pr1,s.c_str(),v,F,hlp); }

int window::button(char* pr0, char* pr1, string s, panel_action_func F,
                                                          const char* hlp)
{ return BASE_WINDOW::button(pr0,pr1,s.c_str(),-1,F,hlp); }

int window::button(char* pr0, char* pr1, string s, const window_handler& obj,
                                                          const char* hlp)
{ panel_action_func f= 0;
  int bt = BASE_WINDOW::button(pr0,pr1,s.c_str(),-1,f,hlp); 
  BASE_WINDOW::set_action(bt, obj);
  return bt;
}



int window::button(char* pr0, char* pr1, string s, int v, panel_action_func F, 
                                                          const char* hlp)
{ return BASE_WINDOW::button(pr0,pr1,s.c_str(),v,F,hlp); }

int window::button(char* pr0, char* pr1, string s, int v, const window_handler& obj, 
                                                          const char* hlp)
{ panel_action_func f= 0;
  int bt = BASE_WINDOW::button(pr0,pr1,s.c_str(),v,f,hlp); 
  BASE_WINDOW::set_action(bt, obj);
  return bt;
}



int window::button(char* pr0, char* pr1, string s, window& p, const char* hlp)
{ return BASE_WINDOW::button(pr0,pr1,s.c_str(),-1,(BASE_WINDOW*)&p,hlp); }


int window::button(char* pr0, char* pr1, string s, int val, window& p, 
                                                            const char* hlp)
{ return BASE_WINDOW::button(pr0,pr1,s.c_str(),val,(BASE_WINDOW*)&p,hlp); }


// focus buttons

int window::fbutton(string s, int n, const char* hlp)
{ int b = button(s,n,hlp);
  set_focus_button(b);
  return b;
}

int window::fbutton(string s, const char* hlp)
{ int b = button(s,hlp);
  set_focus_button(b);
  return b;
}

int window::fbutton(string s, int n, void (*F)(int), const char* hlp)
{ int b = button(s,n,F,hlp);
  set_focus_button(b);
  return b;
}

int window::fbutton(string s, int n, const window_handler& obj, const char* hlp)
{ int b = button(s,n,obj,hlp);
  set_focus_button(b);
  return b;
}


window* window::get_window(int but)
{ return (window*)BASE_WINDOW::get_window(but); }

window* window::set_window(int but, window* M)
{ return (window*)BASE_WINDOW::set_window(but, M); }



panel_item window::text_item(string s)   
{ return BASE_WINDOW::text_item(s.c_str()); }

void window::set_text(panel_item it, string s)   
{ BASE_WINDOW::set_text(it,s.c_str()); }

panel_item window::int_item(string s,int& x,const char* hlp) 
{ return BASE_WINDOW::int_item(s.c_str(),&x,hlp);}

panel_item window::int_item(string s,int& x, int l, int h,const char* hlp) 
{ return BASE_WINDOW::slider_item(s.c_str(),&x,l,h,0,hlp);}

panel_item window::int_item(string s,int& x, int l, int h, void (*F)(int),const char* hlp) 
{ return BASE_WINDOW::slider_item(s.c_str(),&x,l,h,F,hlp);}

panel_item window::int_item(string s,int& x, int l, int h, const window_handler& obj,const char* hlp) 
{ panel_item it = BASE_WINDOW::slider_item(s.c_str(),&x,l,h,0,hlp);
  BASE_WINDOW::set_item_object(it,obj);
  return it;
}


panel_item window::double_item(string s, double& x,const char* hlp) 
{ return BASE_WINDOW::float_item(s.c_str(),&x,hlp);}

panel_item window::real_item(string s, double& x,const char* hlp)
{ return BASE_WINDOW::float_item(s.c_str(),&x,hlp); }

panel_item window::string_item(string label, string& x,const char* hlp)
{ //x = ~x; disconnect
  return BASE_WINDOW::string_item(label.c_str(),&x,0,hlp);
 }

panel_item window::string_item(string label, string& x, void (*F)(char*), 
                                                        const char* hlp)
{ //x = ~x; disconnect
  return BASE_WINDOW::string_item(label.c_str(),&x,F,hlp);
 }

panel_item window::string_item(string label, string& x, const window_handler& obj, 
                                                        const char* hlp)
{ //x = ~x; disconnect
  panel_item it = BASE_WINDOW::string_item(label.c_str(),&x,0,hlp);
  BASE_WINDOW::set_item_object_str(it, obj);
  return it;
 } 



panel_item  window::string_item(string label,string& x, const std::list<string>& L,
                                const char* hlp)
{ return string_item(label,x,L,0,0,hlp); }


panel_item  window::string_item(string label,string& x, const std::list<string>& L,
                                void (*F)(char*), const char* hlp)
{ return string_item(label,x,L,0,F,hlp); }


panel_item  window::string_item(string label,string& x, const std::list<string>& L,
                                const window_handler& obj, const char* hlp)
{ panel_item it = string_item(label,x,L,0,0,hlp); 
  BASE_WINDOW::set_item_object_str(it, obj);
  return it;
}



// with sz param

panel_item  window::string_item(string label,string& x, const std::list<string>& L,
                                int sz, const char* hlp)
{ return string_item(label,x,L,sz,0,hlp); }


panel_item  window::string_item(string label,string& x, const std::list<string>& L,
                                int sz, void (*F)(char*), const char* hlp)
{ //x = x.cstring(); disconnect
  const char** p = new const char*[L.size()];
  int    i = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();
  
  for(;iter!=L.end();iter++) {
     s = *iter;
     if (s.length() > 0) p[i++] = s.c_str();
  }   
     
  panel_item it = BASE_WINDOW::string_menu_item(label.c_str(),&x,"",i,p,sz,F,hlp); 
  delete[] p;
  return it;
}

panel_item  window::string_item(string label,string& x, const std::list<string>& L,
                                int sz, const window_handler& obj, const char* hlp)
{ //x = x.cstring(); disconnect
  const char** p = new const char*[L.size()];
  int    i = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();  
  
  for(;iter!=L.end();iter++){ 
     s = *iter;
     if (s.length() > 0) p[i++] = s.c_str();
  }   
     
  panel_item it = BASE_WINDOW::string_menu_item(label.c_str(),&x,"",i,p,sz,0,hlp); 
  delete[] p;
  BASE_WINDOW::set_item_object_str(it, obj);
  return it;
}


void window::set_menu(panel_item it, const std::list<string>& L, int sz)
{ const char** argv = new const char*[L.size()];
  int argc = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();    
  
  for(;iter!=L.end();iter++){ 
     s = *iter;
     if (s.length() > 0) argv[argc++] = s.c_str();
  }
  
  BASE_WINDOW::set_menu(it,argc,argv,sz); 
  delete[] argv;
}



// choice items

panel_item  window::choice_item(string label,int& x, const std::list<string>& L,
                                void (*F)(int), const char* hlp)
{ const char** p = new const char*[L.size()];
  int    i = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();       
  
  for(;iter!=L.end();iter++){
   s = *iter;
   p[i++] = s.c_str();
  }
  
  panel_item it = BASE_WINDOW::choice_item(label.c_str(),&x,i,p,1,0,F,hlp); 
  delete[] p;
  return it;
}

panel_item  window::choice_item(string label,int& x, const std::list<string>& L,
                                const window_handler& obj, const char* hlp)
{ const char** p = new const char*[L.size()];
  int    i = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();    
  
  for(;iter!=L.end();iter++){
   s=*iter;
   p[i++] = s.c_str();
  }
  
  panel_item it = BASE_WINDOW::choice_item(label.c_str(),&x,i,p,1,0,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  delete[] p;
  return it;
}


panel_item  window::choice_mult_item(string label,int& x, const std::list<string>& L,const char* hlp)
{ return choice_mult_item(label,x,L,0,hlp); }

panel_item  window::choice_mult_item(string label,int& x, const std::list<string>& L,
                                                    void (*F)(int),const char* hlp)
{ const char** p = new const char*[L.size()];
  int    i = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();    
  
  for(;iter!=L.end();iter++){
    s = *iter;
    p[i++] = s.c_str();
  }  
    
  panel_item it = BASE_WINDOW::choice_mult_item(label.c_str(),&x,i,p,F,hlp); 
  delete[] p;
  return it;
}

panel_item  window::choice_mult_item(string label,int& x, const std::list<string>& L,
                                                    const window_handler& obj,const char* hlp)
{ const char** p = new const char*[L.size()];
  int    i = 0;
  string s;
  std::list<string>::const_iterator iter = L.begin();    
  
  for(;iter!=L.end();iter++){
    s = *iter;
    p[i++] = s.c_str();
  }  
  
  panel_item it = BASE_WINDOW::choice_mult_item(label.c_str(),&x,i,p,0,hlp); 
  delete[] p;
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}


// bitmap choice item

panel_item window::choice_item(string label, int& x, int n, int w, int h, 
                               unsigned char** bm,const char* hlp) 
{ return BASE_WINDOW::bitmap_choice_item(label.c_str(),&x,n,w,h,bm,0,hlp); }

panel_item window::choice_item(string label, int& x, int n, int w, int h, 
                               unsigned char** bm, void (*F)(int),
                                                   const char* hlp) 
{ return BASE_WINDOW::bitmap_choice_item(label.c_str(),&x,n,w,h,bm,F,hlp); }

panel_item window::choice_item(string label, int& x, int n, int w, int h, 
                               unsigned char** bm, const window_handler& obj,
                                                   const char* hlp) 
{ panel_item it = BASE_WINDOW::bitmap_choice_item(label.c_str(),&x,n,w,h,bm,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}



panel_item window::choice_mult_item(string label, int& x, int n, int w, int h, 
                                    unsigned char** bm, const char* hlp) 
{ return BASE_WINDOW::bitmap_choice_mult_item(label.c_str(),&x,n,w,h,bm,0,hlp); }

panel_item window::choice_mult_item(string label, int& x, int n, int w, int h, 
                                    unsigned char** bm, void (*F)(int),
                                                        const char* hlp) 
{ return BASE_WINDOW::bitmap_choice_mult_item(label.c_str(),&x,n,w,h,bm,F,hlp); }

panel_item window::choice_mult_item(string label, int& x, int n, int w, int h, 
                                    unsigned char** bm, const window_handler& obj,
                                                        const char* hlp) 
{ panel_item it = BASE_WINDOW::bitmap_choice_mult_item(label.c_str(),&x,n,w,h,bm,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}


panel_item  window::choice_item(string label,int& x,const char* hlp,int n, ...)
{ const char* p[32];
  va_list arg_list;
  va_start(arg_list,n);
  for(int i=0; i<n; i++)
    p[i] = va_arg(arg_list,char*);
  va_end(arg_list);
  return BASE_WINDOW::choice_item(label.c_str(),&x,n,p,1,0,0,hlp);
 }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2)
{ return choice_item(label,x,0,2,s1,s2); }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2, const char* s3)
{ return choice_item(label,x,0,3,s1,s2,s3); }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2, const char* s3, const char* s4)
{ return choice_item(label,x,0,4,s1,s2,s3,s4); }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2, const char* s3, const char* s4, const char* s5)
{ return choice_item(label,x,0,5,s1,s2,s3,s4,s5); }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6)
{ return choice_item(label,x,0,6,s1,s2,s3,s4,s5,s6); }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6, const char* s7)
{ return choice_item(label,x,0,7,s1,s2,s3,s4,s5,s6,s7); }

panel_item  window::choice_item(string label,int& x,const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6, const char* s7, const char* s8)
{ return choice_item(label,x,0,8,s1,s2,s3,s4,s5,s6,s7,s8); }



panel_item window::int_item(string s,int& x,int low, int high, int step,const char* hlp)   
{ return int_item(s,x,low,high,step,0,hlp); }


panel_item window::int_item(string s,int& x,int low, int high, int step,
                                            void (*F)(int), const char* hlp)   
{ int n = (high-low)/step +1;
  char** p = new char*[n];
  for(int i = 0; i < n; i++) 
  { p[i] =  new char[16];
    CGAL_CLIB_STD::sprintf(p[i],"%d",low+i*step);
   }
  panel_item it = 
   BASE_WINDOW::choice_item(s.c_str(),&x,n,(const char**)p,step,low,F,hlp);
  for(int j = 0; j < n; j++)  delete[] p[j];
  delete[] p;
  return it;
 }
 
panel_item window::int_item(string s,int& x,int low, int high, int step,
                                            const window_handler& obj, const char* hlp)   
{ int n = (high-low)/step +1;
  char** p = new char*[n];
  for(int i = 0; i < n; i++) 
  { p[i] =  new char[16];
    CGAL_CLIB_STD::sprintf(p[i],"%d",low+i*step);
   }
  panel_item it = 
   BASE_WINDOW::choice_item(s.c_str(),&x,n,(const char**)p,step,low,0,hlp);
  for(int j = 0; j < n; j++)  delete[] p[j];
  delete[] p;
  BASE_WINDOW::set_item_object(it, obj);
  return it;
 } 
 


panel_item window::bool_item(string s, bool& x,void (*F)(int),const char* hlp)
{ return BASE_WINDOW::bool_item(s.c_str(),&x,F,hlp); }

panel_item window::bool_item(string s, bool& x,const char* hlp)
{ return BASE_WINDOW::bool_item(s.c_str(),&x,0,hlp); }

panel_item window::bool_item(string s, bool& x,const window_handler& obj,const char* hlp)
{ panel_item it = BASE_WINDOW::bool_item(s.c_str(),&x,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}


panel_item window::color_item(string s, color& x,const char* hlp)
{ return BASE_WINDOW::color_item(s.c_str(),(int*)&x,0,hlp); }

panel_item window::color_item(string s, color& x,void (*F)(int),const char* hlp)
{ return BASE_WINDOW::color_item(s.c_str(),(int*)&x,F,hlp); }

panel_item window::color_item(string s, color& x,const window_handler& obj,const char* hlp)
{ panel_item it = BASE_WINDOW::color_item(s.c_str(),(int*)&x,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}



panel_item window::pstyle_item(string s, point_style& x, const char* hlp)
{ return BASE_WINDOW::pstyle_item(s.c_str(),&x,0,hlp); }

panel_item window::pstyle_item(string s, point_style& x, void (*F)(int),
                                                        const char* hlp)
{ return BASE_WINDOW::pstyle_item(s.c_str(),&x,F,hlp); }

panel_item window::pstyle_item(string s, point_style& x, const window_handler& obj,
                                                        const char* hlp)
{ panel_item it = BASE_WINDOW::pstyle_item(s.c_str(),&x,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}




panel_item window::lstyle_item(string s, line_style& x, const char* hlp)
{ return BASE_WINDOW::lstyle_item(s.c_str(),&x,0,hlp); }

panel_item window::lstyle_item(string s, line_style& x, void (*F)(int),
                                                        const char* hlp)
{ return BASE_WINDOW::lstyle_item(s.c_str(),&x,F,hlp); }

panel_item window::lstyle_item(string s, line_style& x, const window_handler& obj,
                                                        const char* hlp)
{ panel_item it = BASE_WINDOW::lstyle_item(s.c_str(),&x,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}



panel_item window::lwidth_item(string s, int& x, const char* hlp)
{ return BASE_WINDOW::lwidth_item(s.c_str(),&x,0,hlp); }

panel_item window::lwidth_item(string s, int& x, void (*F)(int),
                                                        const char* hlp)
{ return BASE_WINDOW::lwidth_item(s.c_str(),&x,F,hlp); }

panel_item window::lwidth_item(string s, int& x, const window_handler& obj,
                                                        const char* hlp)
{ panel_item it = BASE_WINDOW::lwidth_item(s.c_str(),&x,0,hlp); 
  BASE_WINDOW::set_item_object(it, obj);
  return it;
}


static void casteljau(double t, int n, double* ax, double* ay,
                                       double& x, double& y)
{ double* bx = new double[n];
  double* by = new double[n];

  for(int j=0; j<n; j++)
  { bx[j] = ax[j];
    by[j] = ay[j];
   }

  for(int k=1; k<n; k++)
  { for(int i=n-1; i>=k; i--)
    { bx[i] = (1-t)*bx[i-1] + t*bx[i];
      by[i] = (1-t)*by[i-1] + t*by[i];
     }
  }

  x = bx[n-1];
  y = by[n-1];

  delete[] bx;
  delete[] by;
}


void window::compute_bezier(const std::list<window_point>& L, int m, double* X, double* Y)
{ 
  int n = L.size();

  double* bx = new double[n];
  double* by = new double[n];

  int i = 0;
  window_point p;
  std::list<window_point>::const_iterator it = L.begin();
  
  for(;it != L.end();it++)
  { 
    p = *it;
    bx[i] = p.x();
    by[i] = p.y();
    i++;
   }

  double d = 1.0/(m-1);

  X[0] = bx[0];
  Y[0] = by[0];

  for(int j=1; j<m-1; j++) casteljau(j*d,n,bx,by,X[j],Y[j]);

  X[m-1] = bx[n-1];
  Y[m-1] = by[n-1];

  delete[] bx;
  delete[] by;
}


void window::draw_bezier(const std::list<window_point>& L, int m, int arrow, double d, color c)
{ 
  double* X = new double[m];
  double* Y = new double[m];

  compute_bezier(L,m,X,Y);


  if (arrow & 1)
  { 
    window_point pa(X[1],Y[1]), pb(X[0],Y[0]);
    window_point q = draw_arrow_head(pb,seg_angle(pa,pb),d,c);
    X[0] = q.x();
    Y[0] = q.y();
   }

  if (arrow & 2)
  { 
    window_point pa(X[1],Y[1]), pb(X[0],Y[0]);  
    window_point q = draw_arrow_head(pb,seg_angle(pa,pb),d,c);
    X[m-1] = q.x();
    Y[m-1] = q.y();
   }

  BASE_WINDOW::draw_polyline(m,X,Y,c);

  delete[] X;
  delete[] Y;
}


void window::bezier_segments(const std::list<window_point>& L, int m, window_point& ps1, window_point& pt1, 
                                                          window_point& ps2, window_point& pt2)
{ double* X = new double[m];
  double* Y = new double[m];
  compute_bezier(L,m,X,Y);
  
  window_point  pa(X[1],Y[1]), pb(X[0],Y[0]);
  ps1 = pa;
  pt1 = pb;
  
  window_point  pc(X[m-2],Y[m-2]), pd(X[m-1],Y[m-1]);
  ps2 = pc;
  pt2 = pd;
  
  delete[] X;
  delete[] Y;
}



void window::draw_bezier(const std::list<window_point>& L, int m, color c)
{ draw_bezier(L,m,0,0,c); }

void window::draw_bezier_arrow(const std::list<window_point>& L, int m, color c)
{ double d = pix_to_real(get_line_width());
  draw_bezier(L,m,2,d,c); 
 }

  

// pixmaps

static void read_string_from_istream(string& sin, std::istream& s, char delim)
{ char buf[1024];
  char* q = buf+1023;
  bool go = true;

  sin = string(); // clear string

  while (s && go)
  { char* p;
    for(p = buf; p < q && s.get(*p); p++)
    { if (*p == delim || (delim == 0) && isspace(*p)) 
      { if (delim != '\n') s.putback(*p);
        go = false;
        break;
       }
     } 
    *p = '\0';
    //operator+=(buf);
    sin.append(buf);
   }

}


char* window::create_pixrect(string fname)
{
   std::ifstream in(fname.c_str());

   if (!in.good())
   { 
     string serr("window::create_pixrect: Cannot open ");
     serr.append(fname);
     
     std::cerr << serr.c_str() << "\n";
     return 0;
    }

   std::list<string> L;
   string line;

   char c = ' ';

   while (in)
   { while (in && c != '"') in >> c;
     if (c == '"')
     { read_string_from_istream(line,in,'"');
       in >> c;
       in >> c;
       L.push_back(line);
      }
    }

   int n = L.size();

   char** xpm = new char*[n];

   char** p = xpm;

   string s;
   std::list<string>::const_iterator it = L.begin();
   
   for(;it!=L.end();it++) {
     s = *it;
     *p++ = (char*) s.c_str();
   }

   return create_pixrect((const char**)xpm);
}


// edit string

void window::string_edit(window_point p, string& s)
{ BASE_WINDOW::string_edit(p.x(),p.y(),&s); }


void window::draw_text_with_cursor(double x, double y, string s, int cursor, 
                                                                 color col)
{  draw_text(x,y,s,col);

  if (s == "" && cursor >= 0)
  { double dx =  text_width("|");
    double dy =  text_height("|");
    draw_text(x,y,"| |",col);
    BASE_WINDOW::draw_text_cursor(x+dx,y-dy,col);
    return;
   }

  if (cursor < 0) return;

  std::list<string> L = break_into_lines(s);
  
  string st2 = s;

  string::size_type ST;
  ST = st2.find('\n');
    
  while (ST != string::npos){
     // replace ...
     st2.replace(ST,ST+1," ");
     ST = st2.find('\n');
  }  
  
  //was s.replace_all('\n',' ')
  double tht = text_height(st2.c_str());

  int len = 0;
  while (!L.empty())
  { s = L.front(); L.pop_front();
    y -= tht; 
    if (cursor >= len && 
        cursor <= len + ((int)s.length())) 
    { cursor -= len;
      x += text_width((s.substr(0,cursor-1)).c_str());
      BASE_WINDOW::draw_text_cursor(x,y,col);
      break;
     }
    len += s.length()+1;
   }

}


int window::string_edit(double x0, double y1, string& s, int& curs_pos)
{ int timeout = 0;
  return string_edit(x0,y1,s,curs_pos,timeout);
 }


int window::string_edit(double x0, double y1, string& s, int& curs_pos,
                                                          int& timeout)
{ 
  int val;
  unsigned long t0 = 0;
  unsigned long t  = 0;
  double x,y;
  int k;

  do { if (timeout > 0)
         { k = read_event(val,x,y,t,timeout);
           if (t0 == 0) 
             t0 = t;
           else
             if (int(t-t0) > timeout) k = no_event;
          }
       else
          k = read_event(val,x,y);
   } while (k != no_event && k != button_press_event && k != key_press_event);

  if (k == no_event)
  { timeout = 0;
    return 1;
   }

  if (k == button_press_event)
  { double w = text_width(s.c_str());
    double h = text_height(s.c_str());
    double x1 = x0 + w;
    double y0 = y1 - h;
    if (y > y0 && y < y1 && x > x0 && x < x1) 
    { double l = x - x0;
      int j = s.length()-1;
      while (text_width((s.substr(0,j)).c_str()) > l) j--;
      curs_pos = j+1;
      return 0;
     }
    return -1;
   }

  //if (val == KEY_RETURN) return -1;
  if (val == KEY_ESCAPE) return -1;

  char c    = (char)val;
  int  j    = curs_pos;
            
  if (isprint(c))
  { s = s.insert(j,1,c);
    j++;
   }
  else
  { 
    switch (c) {

     case KEY_RETURN:    s = s.insert(j++,"\n");
                         break;

     case KEY_BACKSPACE: if (j > 0) s = s.erase(--j,1);
                         break;
  
     case KEY_LEFT:      if (j > 0) j--;
                         break;
  
     case KEY_RIGHT:     if (j < ((int)s.length())) j++;
                         break;
  
     case KEY_HOME:      j = 0;
                         break;
  
     case KEY_END:       j = s.length();
                         break;
     }
   }

  curs_pos = j;
  return c;
}


void window::set_clip_polygon(const std::list<window_point>& lp)
{ int n = lp.size();
  double* X = new double[n];
  double* Y = new double[n];
  n = 0;
  window_point p;
  std::list<window_point>::const_iterator it = lp.begin();
  
  for(;it != lp.end(); it++) 
  { p= *it;
    X[n] = p.x();
    Y[n] = p.y();
    n++;
   }
  BASE_WINDOW::set_clip_polygon(n,X,Y);
  delete[] X;
  delete[] Y;
}


void window::screenshot(string fname, bool full_color)
{ BASE_WINDOW::screenshot(fname.c_str(),full_color); }


bool window::contains(double xc, double yc) const
{ return ( xc >= xmin() && xc <= xmax() &&
           yc >= ymin() && yc <= ymax() );
 }



void window::adjust_zoom_rect(double& x0, double& y0, double& x1, double& y1) 
{ 
  double wdx = xmax() - xmin();
  double wdy = ymax() - ymin();

  if (x0 > x1) wdy = -wdy;
  if (y0 > y1) wdx = -wdx;

  window_point lp1(x0,y0), lp2(x0+wdx,y0+wdy);

  if (compute_orientation(lp1,lp2,window_point(x1,y1)) == -1)  // was CGAL::RIGHTTURN
    y1 = y_proj(lp1,lp2,x1);
  else
    x1 = x_proj(lp1,lp2,y1);
 }


static void adjust_zoom_rect2(window& W, 
                              double x0, double y0, double x, double y, 
                              double& x1, double& y1, double& x2, double& y2)
{ 
  W.adjust_zoom_rect(x0,y0,x,y);
  
  //was p2 = p1.reflect(p0);  
  double xm = x-x0, ym = y-y0;
  window_point p2(x-xm, y-ym); // was + ....
  
  x1 = p2.x();
  y1 = p2.y();
  x2 = x;
  y2 = y;
  if (x1 > x2) { double tmp = x1; x1 = x2; x2 = tmp; } 
  if (y1 > y2) { double tmp = y1; y1 = y2; y2 = tmp; } 
}





static void draw_zoom_rect(window& W, char* buf, 
                           double x1, double y1, double x2, double y2, 
                           double x01,double y01,double x02,double y02) 
{ double d = W.pix_to_real(2);
  if (x01 > x1) x01 = x1;
  if (y01 > y1) y01 = y1;
  if (x02 < x2) x02 = x2;
  if (y02 < y2) y02 = y2;
  W.put_pixrect(buf);
  W.draw_point((x1+x2)/2,(y1+y2)/2);
  W.draw_rectangle(x1,y1,x2,y2);
  W.flush_buffer(x01-d,y01-d,x02+d,y02+d);
}



bool window::read_zoom_rect(const window_point& cent, window_point& p, window_point& q)
{

  point_style p_style = set_point_style(cross_point); 

  void (*c_handler)(window*,double,double) = get_show_coord_handler();
  set_show_coord_handler(nil);

  double x0 = cent.x(); 
  double y0 = cent.y(); 

  const char* msg2 = "~~~\\bf size:    \\rm move ~~~   \\bf position:\\rm click & drag~~~   \\bf zoom:    \\rm click & release~~~ \\bf cancel:  \\rm click right";
			 
  string st_str = get_status_string();
  set_status_string(msg2);


  char* buf = get_window_pixrect();

  start_buffering();

  bool canceled = false;

  double x1 = x0;
  double y1 = y0;
  double x2 = x0;
  double y2 = y0;

  double xpos = x0 + pix_to_real(8);
  double ypos = y0 - pix_to_real(8);

  adjust_zoom_rect2(*this,x0,y0,xpos,ypos,x1,y1,x2,y2);
  draw_zoom_rect(*this,buf,x1,y1,x2,y2,x1,y1,x2,y2);

  draw_point(x0,y0);
  flush_buffer();

  for(;;)
  { double x01 = x1;
    double y01 = y1;
    double x02 = x2;
    double y02 = y2;

    int event,but;
    while ( (event = read_event(but,xpos,ypos)) != button_press_event)
    { if (event != motion_event) continue;
      adjust_zoom_rect2(*this,x0,y0,xpos,ypos,x1,y1,x2,y2);
      draw_zoom_rect(*this,buf,x1,y1,x2,y2,x01,y01,x02,y02);
      x01 = x1;
      y01 = y1;
      x02 = x2;
      y02 = y2;
     }
    if (but == MOUSE_BUTTON(3)) 
    { canceled = true;
      break;
     }

    //if (but == MOUSE_BUTTON(1)) 
    { double xpos0 = xpos;
      double ypos0 = ypos;
      double x,y;
      while (read_event(but,x,y) != button_release_event)
      { double dx = x - xpos;
        double dy = y - ypos;
        draw_zoom_rect(*this,buf,x1+dx,y1+dy,x2+dx,y2+dy,x1,y1,x2,y2);
        x0 += dx;
        y0 += dy;
        x1 += dx;
        y1 += dy;
        x2 += dx;
        y2 += dy;
        xpos = x;
        ypos = y;
       }
       x -= xpos0;
       y -= ypos0;
       if (real_to_pix(x*x + y*y) < 16) break; // no dragging
     }
  }

  stop_buffering();

  set_point_style(p_style); 
  set_show_coord_handler(c_handler);
  set_status_string(st_str.c_str());

  if (canceled) return false;
  
  p = window_point(x1,y1);
  q = window_point(x2,y2);

  return true;
}



bool window::read_zoom_rect(window_point& p, window_point& q )
{
  set_cursor(XC_dotbox);

  const char* msg1 = "~ZOOM AREA~~~~~\\bf select center: \\rm click left~~~~~                                     \\bf cancel: \\rm click right";

  void (*c_handler)(window*,double,double) = get_show_coord_handler();
  set_show_coord_handler(nil);

  string st_str = get_status_string();
  set_status_string(msg1);

  window_point cent;
  int but = read_mouse(cent);

  set_cursor(-1);

  set_status_string(st_str.c_str());
  set_show_coord_handler(c_handler);

  if (but == MOUSE_BUTTON(3)) return false;

  return read_zoom_rect(cent,p,q);
}
   

  

void window::zoom_area(double x0, double y0, double x1, double y1, int steps,
                       void (*redraw_func)(window*))
{  
   if (!is_open()) {
     std::cerr << "zoom: window must be displayed first.\n";
   }
   
   double wx0 = xmin();
   double wy0 = ymin();
   double wx1 = xmax();
   double wy1 = ymax();

   adjust_zoom_rect(x0,y0,x1,y1);

   if (x0 > x1) { double tmp = x0; x0 = x1; x1 = tmp; }
   if (y0 > y1) { double tmp = y0; y0 = y1; y1 = tmp; }

   if (status_win)
   { double f = (x1-x0)/(wx1-wx0);
     y0 -= f*pix_to_real(status_win->height());
    }

   double dx0 = (x0 - wx0)/steps;
   double dy0 = (y0 - wy0)/steps;
   double dx1 = (x1 - wx1)/steps;
   double dy1 = (y1 - wy1)/steps;

   double g_dist = get_grid_dist();

   //grid_style g_style = get_grid_style();

   char* bg_pr = get_bg_pixrect();

   if (g_dist == 0 && bg_pr == 0 && real_to_pix(20) >= 8)
   { double gd = -20;
     //double gd = pix_to_real(25);
     set_grid_dist(-gd);
     set_grid_style(line_grid);
    }

   start_buffering();

   while (steps--)
   { wx0 += dx0;
     wy0 += dy0;
     wx1 += dx1;
     wy1 += dy1;
     init0(wx0,wx1,wy0,get_bg_pixrect() ? 1:0);
     if (redraw_func) 
        redraw_func(this); 
     else
        redraw();
     flush_buffer();
     clear();
    }

   stop_buffering();

   if (g_dist == 0)
   { leda_wait(0.75);
     set_grid_dist(0);
    }

    if (redraw_func) 
      redraw_func(this); 
    else
      redraw();
}



void window::zoom(double f, int steps, void (*redraw_func)(window*))
{
  if (f == 0) // default window
  { double x0 = 0;
    double x1 = 100;
    double y0 = 0;
    double ratio = (ymax() - ymin())/(xmax() - xmin());
    double y1 = y0 + ratio * (x1-x0);

    if (status_win) 
    { double f = (x1-x0)/(xmax()-xmin());
      y0 += f*pix_to_real(status_win->height());
     }

    zoom_area(x0,y0,x1,y1,steps,redraw_func);
    return;
   }


  double xmi=xmin();
  double xma=xmax();
  double ymi=ymin();
  double yma=ymax();
  double xd=(xma-xmi)*(1/f-1)/2;
  double yd=(yma-ymi)*(1/f-1)/2;
  zoom_area(xmi-xd,ymi-yd,xma+xd,yma+yd,steps,redraw_func);
}





void window::scroll_window(const window_point& p, 
                           void (*redraw_func)(window*,double,double,double,double)) 
{
  window& W = *this;

  W.grab_mouse();

  unsigned long t = W.button_press_time();

  double xmin = W.xmin();
  double ymin = W.ymin();
  double xmax = W.xmax();
  double ymax = W.ymax();

  double x = p.x();
  double y = p.y();
  

  W.start_buffering(2*W.width(), 2*W.height());

  double x_off = xmax - x;
  double y_off = y - ymin;

  double x0 = xmin - xmax + x;
  double x1 = xmax - xmin + x;
  double y0 = ymin - ymax + y;
  double y1 = ymax - ymin + y;


  // draw_background

  W.set_offset(x_off,y_off);

  double gd0 = W.get_grid_dist();
  char* bpr = W.get_bg_pixrect();
  grid_style gs = W.get_grid_style();
  double gd = -20;
  //double gd = -W.pix_to_real(25);
  y1 += W.pix_to_real(100);

  if (gd0 == 0 && bpr == 0 && W.real_to_pix(gd) <= -8) 
  { W.set_grid_dist(gd);
    W.set_grid_style(line_grid);
   }

  W.clear(x0,y0,x1,y1,0,0);

  if (redraw_func)
    redraw_func(&W,x0,y0,x1,y1);
  else
    W.redraw();

  W.set_grid_dist(gd0);
  W.set_grid_style(gs);
  W.set_offset(0,0);


  double d = W.pix_to_real(30);

  int dir = 0;

  if (fabs(ymax-y) < d || fabs(ymin-y) < d) dir = 1;
  if (fabs(xmax-x) < d || fabs(xmin-x) < d) dir = 2;

  switch (dir) {
  case 0:  W.set_cursor(XC_fleur);
           break;
  case 1:  W.set_cursor(XC_sb_h_double_arrow);
           break;
  case 2:  W.set_cursor(XC_sb_v_double_arrow);
           break;
  }


  double dx = 0;
  double dy = 0;

  void (*c_handler)(window*,double,double) = get_show_coord_handler();
  set_show_coord_handler(nil);

  string st_str = get_status_string();

  int event;
  do { //string s(" xmin = %6.2f  ymin = %6.2f", xmin-dx,ymin-dy);
       //W.set_status_string(s.c_str());

       unsigned long t1;
       int val;
       event = W.read_event(val,x1,y1,t1);

       if (event == motion_event && t1 - t < 50) continue;

       if (x1 < xmin || x1 > xmax || y1 < ymin || y1 > ymax) continue;

       t  = t1;
       dx = x1-x;
       dy = y1-y;

       if (dir == 1) dy = 0;
       if (dir == 2) dx = 0;

       W.flush_buffer(x_off-dx,y_off+dy);

  } while (event != button_release_event);


  W.stop_buffering();
  W.delete_buffer();
  W.set_cursor(-1);

  W.init0(xmin-dx,xmax-dx,ymin-dy,0);
  W.set_grid_dist(gd0);

  W.start_buffering();

  if (redraw_func)
    redraw_func(&W,W.xmin(),W.ymin(),W.xmax(),W.ymax());
  else
    W.redraw();

  W.flush_buffer();
  W.stop_buffering();

  W.set_status_string(st_str.c_str());
  W.set_show_coord_handler(c_handler);

  W.ungrab_mouse();
}


int window::read_mouse(double x0, double y0, int timeout1, int timeout2, 
                       bool& double_click, bool& drag) 
{
  int timeout3 = timeout1/2;

  double x,y;
  int key,event,val;

  unsigned long t0 = button_press_time();
  unsigned long t,t1;

  drag = false;
  double_click = false;

  do { event = read_event(val,x,y,t1);
       int dt = int(t1-t0);
       if (dt > timeout1) break;
       if (dt > timeout3)
       { int dx = real_to_pix(x-x0);
         int dy = real_to_pix(y-y0);
         //if (event == motion_event && (dx*dx + dy*dy) >= 8) break;
         if (event == motion_event && (dx*dx + dy*dy) >= 12) break;
        }
  } while (event != button_release_event);

  if (event != button_release_event)
     drag = true;
  else
    { double_click = ( int(t1-t0) < timeout1 && 
                      read_event(key,x,y,t,timeout2)==button_press_event );
      if (double_click)
      { do { event = read_event(val,x,y,t1);
             if (int(t1-t0) > timeout1) break;
            } while (event != button_release_event);
        drag = (event != button_release_event);
      }
    }

  return val;
}


// help texts

static void close_about_win(int but)
{ window* wp = window::get_call_window();
  wp->close();
  delete wp;
}


void window::display_help_text(string fname)
{
  fname += ".hlp";

  if (!is_file(fname))
  { char* lroot = getenv("LEDAROOT");
    if (lroot)
       fname = string(lroot) + "/incl/HELP/" + fname;
    if (!is_file(fname))
    { 
      string sn = string("Cannot open file ").append(fname);
      std::cerr << sn.c_str() << "\n";
      return;
     }
   }

  std::ifstream in(fname.c_str());
  string msg;
  while (in)
  { string str;
    read_string_from_istream(str,in,'\n');
    
    // was str = str.replace_all("\\\\","\\n");
    string::size_type ST;
    ST = str.find("\\\\");
    
    while (ST != string::npos){
     // replace ...
     str.replace(ST,ST+4,"\\n");
     // new find ...
     ST = str.find("\\\\");
    } 
  
    if (str == "") str.append("\\6");
    msg.append(str);
    msg.append(" ");
   }

  
  double l = real_to_pix(text_width(msg.c_str()));
  double h = real_to_pix(text_height(msg.c_str()));

  int w = int(1.1*sqrt(l*h));

  panel* pp = new panel(w,-1,fname);
  pp->text_item(msg);
  pp->button("close",close_about_win);
  pp->display(window::center,window::center);
}


} // end namespace 
