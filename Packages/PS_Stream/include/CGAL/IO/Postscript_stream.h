#ifndef CGAL_Postscript_STREAM
#define CGAL_Postscript_STREAM

// For g++ compiler... //

#include <CGAL/Cartesian.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <iomanip.h>
#include <fcntl.h>
#include <fstream.h>
#include <strstream.h>
#include <list.h>

#include <LEDA/basic.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Point_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Bbox_2.h>

#ifndef CGAL_PS_MANIP_H_
#define CGAL_PS_MANIP_H_
class CGAL_PS_Stream;

template <class T>
class CGAL_PS_Modifier {
  friend CGAL_PS_Stream& operator<<(CGAL_PS_Stream& pss,
                                    const CGAL_PS_Modifier<T> &);
public:
  CGAL_PS_Modifier(CGAL_PS_Stream& (CGAL_PS_Stream::*f)(T),T v):
    _PS_func(f), param(v) {}
private:
  CGAL_PS_Stream& (CGAL_PS_Stream::*_PS_func)(T);
  T param;
};

template <class T>
CGAL_PS_Stream& operator<<(CGAL_PS_Stream& pss, const CGAL_PS_Modifier<T> &m)
{
  (pss.*m._PS_func)(m.param);
  return pss;
};

template<class T>
class CGAL_PS_Modifier_creator {
public:
  CGAL_PS_Modifier_creator(CGAL_PS_Stream& (CGAL_PS_Stream::*f)(T)):
    _PS_func(f) {}
  CGAL_PS_Modifier<T> operator() (T param)
  {
    return CGAL_PS_Modifier<T>(_PS_func,param);
  }
private:
  CGAL_PS_Stream& (CGAL_PS_Stream::*_PS_func)(T);
};

#endif //CGAL_PS_MANIP_H_


typedef const char *DashStyle;
      
class CGAL_PS_Stream {
public:
  
  typedef CGAL_Bbox_2 PS_BBox;

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

friend CGAL_PS_Stream;
friend CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const Axis &g);

public:
  Axis(double x, double y, unsigned int t=0.0)
    : _stepx(x), _stepy(y), _thick(t) {}
  Axis(double xy, unsigned int t=0.0)
    : _stepx(xy),_stepy(xy), _thick(t) {}
  Axis(): _stepx(1.0), _stepy(1.0), _thick(0.0) {}

  double stepx()     const { return _stepx;}
  double stepy()     const { return _stepy;}
  double thickness() const { return _thick;}
private:
  double _stepx;
  double _stepy;
  bool _thick;
};

class Grid {

friend CGAL_PS_Stream;
friend CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const Grid &g);
      
public:
  Grid(double x, double y, DashStyle str="[1 5] 0")
    : _stepx(x), _stepy(y) {_style=strdup(str);}
  Grid(double xy, DashStyle str="[1 5] 0")
    : _stepx(xy), _stepy(xy) {_style=strdup(str);}
  Grid() : _stepx(1.0), _stepy(1.0), _style("[1 5] 0") {}

private:
  double    stepx() const { return _stepx;}
  double    stepy() const { return _stepy;}
  DashStyle style() const { return _style;}

  double _stepx;
  double _stepy;
  DashStyle _style;
};
class Label {
friend CGAL_PS_Stream;
  friend CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const Label &txt);
public:
  Label(const char* txt) { _text=strdup(txt);}
private:
  Label() { _text="";}
  const char* text() const { return _text;}
  const char* _text;
};
      
class Latex_Label {
  friend CGAL_PS_Stream;
  friend CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, Latex_Label &txt);
public:
  // Only text can be define by user
  Latex_Label(const char* txt) { _text=strdup(txt);posx=0;posy=0;}
  Latex_Label() { _text="";posx=0;posy=0;}
private:
  // These functions are private because the position of the string must never appears to the user
  // Only stream modifiers will access these data
  void setposition(float x, float y) { posx=x;posy=y;}
  Latex_Label(const char* txt,float x, float y) {_text=strdup(txt);posx=x;posy=y;}

  float xpos() const {return posx;};
  float ypos() const {return posy;};
  const char* text() const { return _text;};

  // Slots
  const char* _text;
  float posx;
  float posy;
};

typedef list<Latex_Label> List_Label;

class Border {
  friend CGAL_PS_Stream;
  friend CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const Border &b);
public:
  Border(int s=0) { _size=s;}
private:
  int size() const { return _size;}
  int _size;
};

class Context {

friend CGAL_PS_Stream;

public:

Context() : _border_color(CGAL_Color(0,0,0)),
  _fill_color(CGAL_Color(0,0,0)),_dot_style(XCROSS),_dot_size(5),
  _thickness(0),_line_style(SOLID),_fill(false),
  _font("Helvetica"),_font_size(12), _anchor_point(CGAL_Point_2<CGAL_Cartesian <double> > (0,0))
  {}

Context(const Context& c)
  :  _border_color(c.get_border_color()),_fill_color(c.get_fill_color()),
     _dot_style(c.get_dot_style()), _dot_size(c.get_dot_size()),
     _thickness(c.get_thickness()), _line_style(strdup(c.get_line_style())),
     _fill(c.get_fill()), _font(strdup(c.get_font())), _font_size(c.get_font_size()),
     _anchor_point(c.get_pos())
  {};

// Accessor

CGAL_Color   get_border_color()  const {return _border_color;}
CGAL_Color   get_fill_color() const {return _fill_color;}
DotStyle     get_dot_style() const {return _dot_style;}
unsigned int get_dot_size()  const {return _dot_size;}
unsigned int get_thickness()     const {return _thickness;}
unsigned int get_font_size()     const {return _font_size;}
DashStyle    get_line_style()    const {return _line_style;}
const char*  get_font()          const {return _font;}
bool         get_fill()          const {return _fill;}
CGAL_Point_2<CGAL_Cartesian <double> > get_pos() const {return _anchor_point;}

void set_border_color(CGAL_Color& c) {_border_color=c;}
void set_fill_color(CGAL_Color& c)   {_fill_color=c;}
void set_dot_style(DotStyle& s)    {_dot_style=s;}
void set_dot_size(unsigned int s)     {_dot_size=s;}
void set_thickness(unsigned int t)    {_thickness=t;}
void set_font_size(unsigned int s)    {_font_size=s;}
void set_fill(bool&b)     {_fill=b;}
void set_current_pos(const CGAL_Point_2<CGAL_Cartesian <double> >& p) {_anchor_point=p;}
void set_line_style(DashStyle style) {_line_style=strdup(style);}
void set_font(const char *font) {_font=strdup(font);}private:

// Store the current border color
CGAL_Color _border_color;
// Store the current fill color
CGAL_Color _fill_color;
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
// Store the name of the font to use. It is only used  for standart Label, not for Latex Label.
const char* _font;
// Store the size of the font. It is only used  for standart Label, not for Latex Label.
unsigned int _font_size;
//Anchor point:
CGAL_Point_2<CGAL_Cartesian <double> > _anchor_point;

};

  CGAL_PS_Stream(const PS_BBox &bb, ostream &os,
                 OutputMode = QUIET);
        CGAL_PS_Stream(const PS_BBox &bb, const char *fname,
                     OutputMode = QUIET);
    CGAL_PS_Stream(const PS_BBox &bb,float H, ostream &os,
                 OutputMode = QUIET);
    CGAL_PS_Stream(const PS_BBox &bb,float H, const char *fname,
                 OutputMode = QUIET);
    CGAL_PS_Stream(const PS_BBox &bb,float L, float H, ostream &os,
                 OutputMode = QUIET);
    CGAL_PS_Stream(const PS_BBox &bb,float L, float H, const char *fname,
                 OutputMode = QUIET);
    ~CGAL_PS_Stream();
    CGAL_PS_Stream& _SetBorderColor(CGAL_Color &);
  CGAL_PS_Stream& _SetFillColor(CGAL_Color &);
  CGAL_PS_Stream& _SetPointSize(unsigned int);
  CGAL_PS_Stream& _SetLineWidth(unsigned int);
  CGAL_PS_Stream& _SetPointStyle(enum DotStyle);
  CGAL_PS_Stream& _SetLineStyle(DashStyle);
  CGAL_PS_Stream& _SetFill(bool);
  CGAL_PS_Stream& _SetDefaultContext(void);
  CGAL_PS_Stream& _SetCurrentContext(const Context &);
  CGAL_PS_Stream& _ShowDirection(bool);
  CGAL_PS_Stream& _MoveTo(CGAL_Point_2< CGAL_Cartesian <double> >);
  CGAL_PS_Stream& _ShowAxis(Axis &);
  CGAL_PS_Stream& _ShowGrid(Grid &);
  CGAL_PS_Stream& _PutPsLabel(const char *);
  CGAL_PS_Stream& _PutLatexLabel(const char *);
  CGAL_PS_Stream& _PutBorder(unsigned int);
  CGAL_PS_Stream& _SetFont(const char *);
  CGAL_PS_Stream& _SetFontSize(unsigned int);


  //Accessors
  ostream&        os()   {return _os;}
  List_Label&   list()        {return _ll;}
  const Context context() const {return ctxt;}
  bool         gs_output()     const {return (bool)(mode()==GS_VIEW);}
  PS_BBox      bbox()          const {return _bbox;}
  int          width()         const {return _width;}
  int          height()        const {return _height;}
  OutputMode   mode()          const {return _mode;}
      
  // Utils
  double xratio() { return _xratio;}
  double yratio() { return _yratio;}
  double x2ps(double x) {
    return (x-_bbox.xmin())*xratio();
  }
  double y2ps(double y) {
    return (y-_bbox.ymin())*yratio();

  }

  bool is_eps();
  bool is_readable();

  private:
CGAL_PS_Stream(const PS_BBox &bb);

CGAL_PS_Stream(const PS_BBox &bb,float H);

CGAL_PS_Stream(const PS_BBox &bb,float L, float H);
CGAL_PS_Stream();

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
int _width,_height;
// Graphical Context
Context ctxt;
// The Output stream
#ifdef CGAL_WORKAROUND_016
_IO_ostream_withassign &_os;
#else
ostream_withassign &_os;
#endif

//List of Latex Labels. They will be inserted at the end of the file.
List_Label _ll;
};
extern const CGAL_PS_Stream::Context CGAL_CTXT_DEFAULT;

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Border &b);

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Axis &g);

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Grid &g);

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Label &txt);

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Latex_Label &txt);

extern CGAL_PS_Modifier_creator<CGAL_Color &> set_fill_color;

extern CGAL_PS_Modifier_creator<CGAL_Color &> set_border_color;

extern CGAL_PS_Modifier_creator<unsigned int> set_point_size;

extern CGAL_PS_Modifier_creator<CGAL_PS_Stream::DotStyle> set_point_style;

extern CGAL_PS_Modifier_creator<DashStyle> set_line_style;

extern CGAL_PS_Modifier_creator<unsigned int> set_line_width;

extern CGAL_PS_Modifier_creator<bool> set_fill;

extern CGAL_PS_Modifier_creator<const CGAL_PS_Stream::Context &> set_current_context;

extern CGAL_PS_Modifier_creator<bool> show_direction;

extern CGAL_PS_Modifier_creator<CGAL_Point_2< CGAL_Cartesian <double> > > move_to;

//extern CGAL_PS_Modifier_creator<CGAL_PS_Stream::Axis &> show_axis;

//extern CGAL_PS_Modifier_creator<CGAL_PS_Stream::Grid &> show_grid;

//extern CGAL_PS_Modifier_creator<const char *> put_ps_label;

//extern CGAL_PS_Modifier_creator<const char *> put_latex_label;

//extern CGAL_PS_Modifier_creator<unsigned int> put_border;

extern CGAL_PS_Modifier_creator<const char *> set_font;

extern CGAL_PS_Modifier_creator<unsigned int> set_font_size;


#ifdef CGAL_POINT_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_Point_2<R> &p)
{

  if (ps.is_readable())
    {
      ps.os() <<  "%CGAL% Point" << endl;
      ps.os() << "%CGAL% "<<p.x()<<" "<<p.y()<<endl;
    }
  if (ps.context().get_dot_style()!=CGAL_PS_Stream::NONE)
    ps.os() << ps.x2ps(p.x()) << " "
            << ps.y2ps(p.y()) << " "
            << ps.context().get_dot_size() << " ";
  switch (ps.context().get_dot_style())
    {
    case CGAL_PS_Stream::EBOX:
      ps.os() << "eb" << endl;
      break;
    case CGAL_PS_Stream::FBOX:
      ps.os() << "fb" << endl;
      break;
    case CGAL_PS_Stream::EDOT:
      ps.os() << "ec" << endl;
      break;
    case CGAL_PS_Stream::FDOT:
      ps.os() << "fc" << endl;
      break;
    case CGAL_PS_Stream::ICROSS:
      ps.os() << "ic" << endl;
      break;
    case CGAL_PS_Stream::XCROSS:
      ps.os() << "xc" << endl;
      break;
    }
  return ps;
}

#endif // CGAL_POINT_2_H
#ifdef CGAL_SEGMENT_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_Segment_2<R> &s)
{
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Segment" << endl;
      ps.os() << "%CGAL% "<<s.source().x()<<" "<<s.source().y()<<" "
              <<s.target().x()<<" "<<s.target().y()<<endl;
    }
  ps.os() << ps.x2ps(s.source().x()) << " "
          << ps.y2ps(s.source().y()) << " mt ";
  ps.os() << ps.x2ps(s.target().x()) << " "
          << ps.y2ps(s.target().y())
          << " lt st" << endl;

  return ps;
}

#endif // CGAL_SEGMENT_2_H
#ifdef CGAL_LINE_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_Line_2<R> &l)
{
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Line" << endl;
      ps.os() << "%CGAL% "<<l.a()<<" "<<l.b()<<" "<<l.c()<<endl;
    }
  if (!l.is_degenerate())
    if (l.is_vertical())
      {
        double t=ps.x2ps(l.x_at_y(0));
        ps.os()<< t;
        ps.os()<< " 0 mt" << endl;
        ps.os()<< t << " " << ps.height() << " lt st" <<endl;
      }
    else
      {
        ps.os() << "0 "
                << ps.y2ps(l.y_at_x(ps.bbox().xmin()))
                << " mt" << endl;
        ps.os() << ps.width() << " "
                << ps.y2ps(l.y_at_x(ps.bbox().xmax()))
                << " lt st"<< endl;
      }
  return ps;
}
#endif // CGAL_LINE_2_H
#ifdef CGAL_RAY_2_H
template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_Ray_2<R> &r)
{
  typedef CGAL_Direction_2<CGAL_Cartesian <double> > dir;
  CGAL_Line_2<R> l=r.supporting_line();
  dir haut(0,1);
  dir bas(0,-1);
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Ray" << endl;
      ps.os() << "%CGAL% "<<r.source().x()<<" "<<r.source().y()<<" ";
      ps.os() << r.second_point().x() << " " << r.second_point().y() <<endl;
    }
  if (!r.is_degenerate())
    {
      ps.os()<< ps.x2ps(r.source().x()) << " "
             << ps.y2ps(r.source().y()) << " mt" << endl;
      if (r.is_vertical())
        {
          ps.os()<< ps.x2ps(r.source().x()) << " ";
          if (r.direction()==haut)
            ps.os() << ps.height();
          else
            ps.os() << "0 ";
          ps.os() << " lt st" << endl;
        }
      else
        if (r.direction()>bas || r.direction()<haut)
          ps.os() << ps.width()
                  << " "
                  <<ps.y2ps(l.y_at_x(ps.bbox().xmax()))
                  << " lt st" << endl;
        else
          ps.os() << "0 "
                  << ps.y2ps(l.y_at_x(ps.bbox().xmin()))
                  << " lt st" << endl;
    }
  return ps;
}

#endif // CGAL_RAY_2_H
#ifdef CGAL_PARABOLA_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps,const CGAL_Parabola<R> &p)
{

  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Parabola" << endl;
      ps.os() << "%CGAL% Base "<<p.base().x()<<" "<<p.base().y()<<endl;
      ps.os() << "%CGAL% Vector "<<p.vertor().x()<<" "p.vector().y()<<endl;
      ps.os() << "%CGAL% Curvature "<<p.curvature()<<endl;
    }
  return ps;
}X


#endif // CGAL_PARABOLA_2_H
#ifdef CGAL_TRIANGLE_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps,const CGAL_Triangle_2<R> &t)
{
  if (ps.is_readable())
    {
      ps.os() << "%CGAL% Triangle" << endl;
      for (int i=0;i<3;i++)
        ps.os() << "%CGAL " << t[i].x() << " " << t[i].y() << endl;
    }
  for (int i=0;i<4;i++)
    ps.os() << ps.x2ps(t[i].x())<< " " << ps.y2ps(t[i].y()) << " ";

  ps.os() << "tr ";
  if (ps.context().get_fill())
    {
      ps.os() << "gsave " << endl;
      ps.os() << ps.context().get_fill_color().r() << " "
              << ps.context().get_fill_color().g() << " "
              << ps.context().get_fill_color().b()
              << " setcolor fill grestore " <<endl;
    }
  ps.os() << "st" <<endl;
  return ps;
}


#endif // CGAL_TRIANGLE_2_H
#ifdef CGAL_ISO_RECTANGLE_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps,const CGAL_Iso_rectangle_2<R> &r)
{
    if (ps.is_readable())
      {
        ps.os() << "%CGAL% Rectangle" << endl;
        for (int i=0;i<4;i++)
          ps.os() << "%CGAL " << r[i].x() << " " << r[i].y() << endl;
      }
    for (int i=0;i<5;i++)
      ps.os() << ps.x2ps(r[i].x()) << " " << ps.y2ps(r[i].y()) << " ";

    ps.os() << "re ";
    if (ps.context().get_fill())
      {
        ps.os() << "gsave " << endl;
        ps.os() << ps.context().get_fill_color().r() << " "
                << ps.context().get_fill_color().g() << " "
                << ps.context().get_fill_color().b()
                << " setcolor fill grestore " <<endl;
      }
    ps.os() << "st" <<endl;
    return ps;
  }
#endif // CGAL_ISO_RECTANGLE_2_H
#ifdef CGAL_CIRCLE_2_H

template < class R >
CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps,
                                   const CGAL_Circle_2<R> &c)
{
    if (ps.is_readable())
      {
        ps.os() << "%CGAL% Circle" << endl;
        ps.os() << "%CGAL " << c.center().x() << " " << c.center().y() << endl;
        ps.os() << "%CGAL " << c.squared_radius() << endl;
      }
    double ratio=ps.yratio()/ps.xratio();
    ps.os()<< "gsave 1 " << ratio << " scale" << endl;
    ps.os()<< ps.x2ps(c.center().x()) << " " << ps.y2ps(c.center().y())/ratio
           << " " << c.squared_radius()*ps.xratio()  << " 0 360 arc " << endl;
    if (ps.context().get_fill())
      {
        ps.os() << "gsave " << endl;
        ps.os() << ps.context().get_fill_color().r() << " "
                << ps.context().get_fill_color().g() << " "
                << ps.context().get_fill_color().b()
                << " setcolor fill grestore " <<endl;
      }
    ps.os() << "st grestore" <<endl;
    return ps;
  }
#endif // CGAL_CIRCLE_2_H


#endif  // CGAL_Postscript_STREAM
