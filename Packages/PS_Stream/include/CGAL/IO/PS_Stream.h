// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/IO/PS_Stream.h
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IO_PS_STREAM_H
#define CGAL_IO_PS_STREAM_H

#include <CGAL/basic.h>

#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <fcntl.h>
#include <fstream>
#include <strstream>
#include <iterator>
#include <list>

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Point_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

class PS_Stream;

template <class T>
class PS_Manipulator {
  friend PS_Stream& operator<< CGAL_NULL_TMPL_ARGS (PS_Stream& ,
                                    const PS_Manipulator<T> &);
public:
  PS_Manipulator(PS_Stream& (PS_Stream::*f)(T),T v):
    _PS_func(f), param(v) {}
protected:
  PS_Stream& (PS_Stream::*_PS_func)(T);
  T param;
};

template <class T>
PS_Stream& operator << (PS_Stream& pss, const PS_Manipulator<T> &m)
{
  (pss.*m._PS_func)(m.param);
  return pss;
}

template<class T>
class PS_Manipulator_creator {
public:
  PS_Manipulator_creator(PS_Stream& (PS_Stream::*f)(T)):
    _PS_func(f) {}
  PS_Manipulator<T> operator() (T param)
  {
    return PS_Manipulator<T>(_PS_func,param);
  }
protected:
  PS_Stream& (PS_Stream::*_PS_func)(T);
};


typedef const char *DashStyle;
      
class PS_Stream {
public:
  
  typedef Bbox_2 PS_BBox;

  static const DashStyle SOLID;
  static const DashStyle DASH1;
  static const DashStyle DASH2;
  static const DashStyle DASH3;
  static const DashStyle DASH4;
  static const DashStyle DASH5;
  static const DashStyle DASH6;
  static const float CM;
  static const float INCH;
  static const float POINT;
  
  enum OutputMode {READABLE, QUIET, READABLE_EPS, QUIET_EPS, GS_VIEW};
  enum DotStyle {NONE, XCROSS, ICROSS, EDOT, FDOT, EBOX, FBOX};
  
  class Axis {
    
    friend class PS_Stream;
    friend PS_Stream & operator << (PS_Stream& ps, const Axis& g);
    
  public:
    Axis(double x, double y, unsigned int t=0)
      : _stepx(x), _stepy(y), _thick(t) {}
    Axis(double xy, unsigned int t=0)
      : _stepx(xy),_stepy(xy), _thick(t) {}
    Axis(): _stepx(1.0), _stepy(1.0), _thick(0) {}
    
    double stepx()     const { return _stepx;}
    double stepy()     const { return _stepy;}
    unsigned int thickness() const { return _thick;}
  protected:
    double _stepx;
    double _stepy;
    unsigned int _thick;
  };
  
class Grid {
  
  friend class PS_Stream;
  friend PS_Stream & operator << (PS_Stream& ps, const Grid& g);
  
public:
  Grid(double x, double y, DashStyle str="[1 5] 0")
    : _stepx(x), _stepy(y) {_style=strdup(str);}
  Grid(double xy, DashStyle str="[1 5] 0")
    : _stepx(xy), _stepy(xy) {_style=strdup(str);}
  Grid() : _stepx(1.0), _stepy(1.0), _style("[1 5] 0") {}
  
  double    stepx() const { return _stepx;}
  double    stepy() const { return _stepy;}
  DashStyle style() const { return _style;}
  
protected:
  double _stepx;
  double _stepy;
  DashStyle _style;
};

  class Label {
    
    friend class PS_Stream;
    friend PS_Stream & operator << (PS_Stream& ps, const Label& txt);
    
  public:
    Label(const char* txt) { _text=strdup(txt);}
    Label() { _text="";}
    const char* text() const { return _text;}
    
  protected:
    const char* _text;
    
  };
  
  class Latex_Label {
    
    friend class PS_Stream;
    friend PS_Stream& operator << (PS_Stream& ps, Latex_Label& txt);
  public:
    // Only text can be define by user
    Latex_Label(const char* txt) { _text=strdup(txt);posx=0;posy=0;}
    Latex_Label() { _text="";posx=0;posy=0;}
    
    float xpos() const {return posx;}
    float ypos() const {return posy;}
    const char* text() const { return _text;}
    
  protected:
    // These functions are private because the position of the string
    // must never appears to the user
    // Only stream modifiers will access these data
    void setposition(float x, float y) { posx=x;posy=y;}
    Latex_Label(const char* txt,float x, float y)
    {_text=strdup(txt);posx=x;posy=y;}
    
    // Slots
    const char* _text;
    float posx;
    float posy;
  };
  
  typedef std::list<Latex_Label> List_Label;
  
  class Border {
    
    friend class PS_Stream;
    friend PS_Stream& operator << (PS_Stream &ps, const Border &b);
    
  public:
    Border(int s=0) { _size=s;}
    int size() const { return _size;}
    
  protected:
    
    int _size;
  };

  class Context {
    
    friend class PS_Stream;
    
  public:
    
    Context() : _border_color(Color(0,0,0)),
      _fill_color(Color(0,0,0)),_dot_style(XCROSS),_dot_size(5),
      _thickness(0),_line_style(SOLID),_fill(false),
      _font("Helvetica"),_font_size(12),
      _anchor_point(Point_2<Cartesian <double> > (0,0))
      {}
    
    Context(const Context& c)
      :  _border_color(c.get_border_color()),_fill_color(c.get_fill_color()),
	 _dot_style(c.get_dot_style()), _dot_size(c.get_dot_size()),
	 _thickness(c.get_thickness()),
         _line_style(strdup(c.get_line_style())),
	 _fill(c.get_fill()), _font(strdup(c.get_font())),
         _font_size(c.get_font_size()),
	 _anchor_point(c.get_pos())
      {};

// Accessor

Color   get_border_color()  const {return _border_color;}
Color   get_fill_color() const {return _fill_color;}
DotStyle     get_dot_style() const {return _dot_style;}
unsigned int get_dot_size()  const {return _dot_size;}
unsigned int get_thickness()     const {return _thickness;}
unsigned int get_font_size()     const {return _font_size;}
DashStyle    get_line_style()    const {return _line_style;}
const char*  get_font()          const {return _font;}
bool         get_fill()          const {return _fill;}
Point_2<Cartesian <double> > get_pos() const {return _anchor_point;}

void set_border_color(const Color& c) {_border_color=c;}
void set_fill_color(const Color& c)   {_fill_color=c;}
void set_dot_style(DotStyle& s)    {_dot_style=s;}
void set_dot_size(unsigned int s)     {_dot_size=s;}
void set_thickness(unsigned int t)    {_thickness=t;}
void set_font_size(unsigned int s)    {_font_size=s;}
void set_fill(bool& b)     {_fill=b;}
void set_current_pos(const Point_2<Cartesian <double> >& p) {_anchor_point=p;}
void set_line_style(DashStyle style) {_line_style=strdup(style);}
void set_font(const char *font) {_font=strdup(font);}

protected:

// Store the current border color
Color _border_color;
// Store the current fill color
Color _fill_color;
// Store the current dot style
enum DotStyle _dot_style;
// Store the current dot size
unsigned int _dot_size;
// Store the current line thickness
unsigned int _thickness;
// Store the current line style
DashStyle _line_style;
// Define if a polygone must be fill or not
bool _fill;
// Define if direction must be shown.
bool _dir;
// Store the name of the font to use. It is only used  for standart Label,
// not for Latex Label.
const char* _font;
// Store the size of the font. It is only used  for standart Label,
// not for Latex Label.
unsigned int _font_size;
//Anchor point:
Point_2<Cartesian <double> > _anchor_point;

};

//Constructors used for PS_Stream 3D
PS_Stream(std::ostream& os,OutputMode = QUIET);
PS_Stream(const char* fname, OutputMode = QUIET);
PS_Stream(float H,std::ostream& os,OutputMode = QUIET);
PS_Stream(float H, const char* fname, OutputMode = QUIET);

  //Constructors
  PS_Stream(const PS_BBox& bb, std::ostream& os,
	    OutputMode = QUIET);
  PS_Stream(const PS_BBox& bb, const char* fname,
	    OutputMode = QUIET);
  PS_Stream(const PS_BBox& bb,float H, std::ostream& os,
	    OutputMode = QUIET);
  PS_Stream(const PS_BBox& bb,float H, const char* fname,
	    OutputMode = QUIET);
  PS_Stream(const PS_BBox& bb,float L, float H, std::ostream& os,
	    OutputMode = QUIET);
  PS_Stream(const PS_BBox& bb,float L, float H, const char* fname,
	    OutputMode = QUIET);
  PS_Stream(const PS_BBox& bb,float L, float H);
  
  ~PS_Stream();
  
  PS_Stream& set_border_color(const Color&);
  PS_Stream& set_fill_color(const Color&);
  PS_Stream& set_point_size(unsigned int);
  PS_Stream& set_line_width(unsigned int);
  PS_Stream& set_point_style(enum DotStyle);
  PS_Stream& set_line_style(DashStyle);
  PS_Stream& set_fill(bool);
  PS_Stream& set_default_context(void);
  PS_Stream& set_current_context(const Context&);
  PS_Stream& set_point(Point_2< Cartesian <double> >);
  PS_Stream& set_axis(Axis&);
  PS_Stream& set_grid(Grid&);
  PS_Stream& put_ps_label(const char*);
  PS_Stream& put_latex_label(const char*);
  PS_Stream& put_border(unsigned int);
  PS_Stream& set_font(const char*);
  PS_Stream& set_font_size(unsigned int);
  

  //Accessors
  std::ostream&        os()   {return _os;}
  List_Label&   list()        {return _ll;}
  const Context context() const {return ctxt;}
  bool         gs_output()     const {return (bool)(mode()==GS_VIEW);}
  PS_BBox      bbox()          const {return _bbox;}
  int          width()         const {return _width;}
  int          height()        const {return _height;}
  OutputMode   mode()          const {return _mode;}

  void set_scale(const PS_BBox& bb)
  {
      _xratio=_width/(bb.xmax()-bb.xmin());
      _yratio=_height/(bb.ymax()-bb.ymin());
  }
  void set_window(PS_BBox bb,float H)
  {
      _width=(int)((bb.xmax()-bb.xmin())*H/(bb.ymax()-bb.ymin()));
  }
      
  // Utils
  double xratio() { return _xratio;}
  double yratio() { return _yratio;}
  double x2ps(double x) { return (x-_bbox.xmin())*xratio();}
  double y2ps(double y) { return (y-_bbox.ymin())*yratio();}
  bool is_eps();
  bool is_readable();

protected:
  //   PS_Stream(const PS_BBox& bb);
  
  // PS_Stream(const PS_BBox& bb,float H);
  
  // PS_Stream(const PS_BBox& bb,float L, float H);
  
  
  // Manipulation du contexte.
  void setdefault();
  void setcontext();
  
// Pour inserer l'entete
  void insert_catalogue();
  
  // Define the scale.
  double _xratio;
  double _yratio;
  
// Define the boounding box
  PS_BBox _bbox;
  
  // OutputMode
  OutputMode _mode;

  // Size of output.
  int _width;
  int _height;
  
  // Graphical Context
  Context ctxt;

  // In case it's a file, we need to store the object to be able to reference
  // it, hence the of.
  std::ofstream of;
  std::ostream& _os;

  //List of Latex Labels. They will be inserted at the end of the file.
  List_Label _ll;
};

extern const PS_Stream::Context CTXT_DEFAULT;

PS_Stream & operator <<(PS_Stream& , const PS_Stream::Border& );

PS_Stream & operator <<(PS_Stream& , const PS_Stream::Axis& );

PS_Stream & operator <<(PS_Stream& , const PS_Stream::Grid& );

PS_Stream & operator <<(PS_Stream& , const PS_Stream::Label&);

PS_Stream& operator << (PS_Stream& , const PS_Stream::Latex_Label& txt);

extern PS_Manipulator_creator<const Color&> fill_color;

extern PS_Manipulator_creator<const Color&> border_color;

extern PS_Manipulator_creator<unsigned int> point_size;

extern PS_Manipulator_creator<PS_Stream::DotStyle> point_style;

extern PS_Manipulator_creator<DashStyle> line_style;

extern PS_Manipulator_creator<unsigned int> line_width;

extern PS_Manipulator_creator<bool> fill;

extern PS_Manipulator_creator<const PS_Stream::Context&> current_context;

extern PS_Manipulator_creator<Point_2< Cartesian <double> > > move_to;

extern PS_Manipulator_creator<PS_Stream::Axis&> show_axis;

extern PS_Manipulator_creator<PS_Stream::Grid&> show_grid;

extern PS_Manipulator_creator<const char*> ps_label;

extern PS_Manipulator_creator<const char*> latex_label;

extern PS_Manipulator_creator<unsigned int> border;

extern PS_Manipulator_creator<const char*> font;

extern PS_Manipulator_creator<unsigned int> font_size;



#ifdef CGAL_POINT_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps, const Point_2<R>& p)
{

  if (ps.is_readable())
    {
      ps.os() <<  "%CGAL% Point" << std::endl;
      ps.os() << "%CGAL% "<<p.x()<<" "<<p.y()<<std::endl;
    }
  if (ps.context().get_dot_style()!=PS_Stream::NONE)
    {
  ps.os() << ps.x2ps(p.x()) << " "
            << ps.y2ps(p.y()) << " "
            << ps.context().get_dot_size() << " ";
  
  switch (ps.context().get_dot_style())
    {
    case PS_Stream::EBOX:
      ps.os() << "eb" << std::endl;
      break;
    case PS_Stream::FBOX:
      ps.os() << "fb" << std::endl;
      break;
    case PS_Stream::EDOT:
      ps.os() << "ec" << std::endl;
      break;
    case PS_Stream::FDOT:
      ps.os() << "fc" << std::endl;
      break;
    case PS_Stream::ICROSS:
      ps.os() << "ic" << std::endl;
      break;
    default :
      ps.os() << "xc" << std::endl;
      break;
    }
    }  

else {
     ps.os() << ps.x2ps(p.x()) << " " 
	     << ps.y2ps(p.y())<< " " ;
     }
  return ps;
}

#endif // CGAL_POINT_2_H
#ifdef CGAL_SEGMENT_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps, const Segment_2<R>& s)
{
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Segment" << std::endl;
      ps.os() << "%CGAL% "<<s.source().x()<<" "<<s.source().y()<<" "
              <<s.target().x()<<" "<<s.target().y()<<std::endl;
    }
  ps.os() << ps.x2ps(s.source().x()) << " "
          << ps.y2ps(s.source().y()) << " mt ";
  ps.os() << ps.x2ps(s.target().x()) << " "
          << ps.y2ps(s.target().y())
          << " lt st" << std::endl;

  return ps;
}

#endif // CGAL_SEGMENT_2_H
#ifdef CGAL_LINE_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps, const Line_2<R>& l)
{
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Line" << std::endl;
      ps.os() << "%CGAL% "<<l.a()<<" "<<l.b()<<" "<<l.c()<<std::endl;
    }
  if (!l.is_degenerate())
    if (l.is_vertical())
      {
        double t=ps.x2ps(l.x_at_y(0));
        ps.os()<< t;
        ps.os()<< " 0 mt" << std::endl;
        ps.os()<< t << " " << ps.height() << " lt st" <<std::endl;
      }
    else
      {
        ps.os() << "0 "
                << ps.y2ps(l.y_at_x(ps.bbox().xmin()))
                << " mt" << std::endl;
        ps.os() << ps.width() << " "
                << ps.y2ps(l.y_at_x(ps.bbox().xmax()))
                << " lt st"<< std::endl;
      }
  return ps;
}

#endif // CGAL_LINE_2_H
#ifdef CGAL_RAY_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps, const Ray_2<R>& r)
{
  typedef Direction_2<Cartesian <double> > dir;
  Line_2<R> l=r.supporting_line();
  dir haut(0,1);
  dir bas(0,-1);
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Ray" << std::endl;
      ps.os() << "%CGAL% "<<r.source().x()<<" "<<r.source().y()<<" ";
      ps.os() << r.second_point().x() << " " << r.second_point().y();
      ps.os() << std::endl;
    }
  if (!r.is_degenerate())
    {
      ps.os()<< ps.x2ps(r.source().x()) << " "
             << ps.y2ps(r.source().y()) << " mt" << std::endl;
      if (r.is_vertical())
        {
          ps.os()<< ps.x2ps(r.source().x()) << " ";
          if (r.direction()==haut)
            ps.os() << ps.height();
          else
            ps.os() << "0 ";
          ps.os() << " lt st" << std::endl;
        }
      else
        if (r.direction()>bas || r.direction()<haut)
          ps.os() << ps.width()
                  << " "
                  <<ps.y2ps(l.y_at_x(ps.bbox().xmax()))
                  << " lt st" << std::endl;
        else
          ps.os() << "0 "
                  << ps.y2ps(l.y_at_x(ps.bbox().xmin()))
                  << " lt st" << std::endl;
    }
  return ps;
}

#endif // CGAL_RAY_2_H
#ifdef CGAL_PARABOLA_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps,const Parabola<R>& p)
{

  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Parabola" << std::endl;
      ps.os() << "%CGAL% Base "<<p.base().x()<<" "<<p.base().y()<<std::endl;
      ps.os() << "%CGAL% Vector " << p.vertor().x() << " " << p.vector().y()
	      << std::endl;
      ps.os() << "%CGAL% Curvature "<<p.curvature()<<std::endl;
    }
  return ps;
}

#endif // CGAL_PARABOLA_2_H
#ifdef CGAL_TRIANGLE_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps,const Triangle_2<R>& t)
{
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Triangle" << std::endl;
      for (int i=0;i<3;i++)
        ps.os() << "%CGAL " << t[i].x() << " " << t[i].y() << std::endl;
    }
  for (int i=0;i<4;i++)
    ps.os() << ps.x2ps(t[i].x())<< " " << ps.y2ps(t[i].y()) << " ";

  ps.os() << "tr ";
  if (ps.context().get_fill())
    {
      ps.os() << "gsave " << std::endl;
      ps.os() << ps.context().get_fill_color().r() << " "
              << ps.context().get_fill_color().g() << " "
              << ps.context().get_fill_color().b()
              << " setcolor fill grestore " <<std::endl;
    }
  ps.os() << "st" <<std::endl;
  return ps;
}

#endif // CGAL_TRIANGLE_2_H
#ifdef CGAL_ISO_RECTANGLE_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps,const Iso_rectangle_2<R>& r)
{
    if (ps.is_readable())
      {
        ps.os() << "%CGAL% Rectangle" << std::endl;
        for (int i=0;i<4;i++)
          ps.os() << "%CGAL " << r[i].x() << " " << r[i].y() << std::endl;
      }
    for (int i=0;i<5;i++)
      ps.os() << ps.x2ps(r[i].x()) << " " << ps.y2ps(r[i].y()) << " ";

    ps.os() << "re ";
    if (ps.context().get_fill())
      {
        ps.os() << "gsave " << std::endl;
        ps.os() << ps.context().get_fill_color().r() << " "
                << ps.context().get_fill_color().g() << " "
                << ps.context().get_fill_color().b()
                << " setcolor fill grestore " <<std::endl;
      }
    ps.os() << "st" <<std::endl;
    return ps;
  }

#endif // CGAL_ISO_RECTANGLE_2_H
#ifdef CGAL_CIRCLE_2_H

template < class R >
PS_Stream & operator <<(PS_Stream& ps, const Circle_2<R>& c)
{
    if (ps.is_readable())
      {
        ps.os() << "%CGAL% Circle" << std::endl;
        ps.os() << "%CGAL " << c.center().x()<<" "<<c.center().y()<< std::endl;
        ps.os() << "%CGAL " << c.squared_radius() << std::endl;
      }
    double ratio=ps.yratio()/ps.xratio();
    double radius=sqrt(to_double(c.squared_radius()));
    ps.os()<< "gsave 1 " << ratio << " scale" << std::endl;
    ps.os()<< ps.x2ps(c.center().x()) << " " << ps.y2ps(c.center().y())/ratio
           << " " << radius*ps.xratio()  << " 0 360 arc " << std::endl;
    if (ps.context().get_fill())
      {
        ps.os() << "gsave " << std::endl;
        ps.os() << ps.context().get_fill_color().r() << " "
                << ps.context().get_fill_color().g() << " "
                << ps.context().get_fill_color().b()
                << " setcolor fill grestore " <<std::endl;
      }
    ps.os() << "st grestore" <<std::endl;
    return ps;
  }

#endif // CGAL_CIRCLE_2_H

CGAL_END_NAMESPACE

#endif // CGAL_IO_PS_STREAM_H
