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
// file          : src/PS_Stream.C
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Carine Bonetto
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#include <CGAL/basic.h>
#include <CGAL/IO/PS_Stream.h>

CGAL_BEGIN_NAMESPACE

const float PS_Stream::CM=28.4;
const float PS_Stream::INCH=72.0;
const float PS_Stream::POINT=1.0;

const DashStyle PS_Stream::SOLID="[] 0 ";
const DashStyle PS_Stream::DASH2="[5 5] 0 ";
const DashStyle PS_Stream::DASH3="[10 10] 0 ";
const DashStyle PS_Stream::DASH6="[10 6 4 6] 0 ";
const DashStyle PS_Stream::DASH4="[10 5] 0 ";
const DashStyle PS_Stream::DASH5="[5 10] 0 ";
const DashStyle PS_Stream::DASH1="[2 2] 0 ";
extern const PS_Stream::Context CTXT_DEFAULT=PS_Stream::Context();

PS_Stream::PS_Stream(std::ostream& os, OutputMode mode)
   :_bbox(PS_BBox(-2,-2,2,2)),_mode(mode),_width((int)(21*CM)),
    _height((int)(29.7*CM)),_os(os)
{
 insert_catalogue();
}

PS_Stream::PS_Stream(const char* fname, OutputMode mode)
  :_bbox(PS_BBox(-2,-2,2,2)), _mode(mode), _width((int)(21*CM)),
   _height((int)(29.7*CM)), of(fname,std::ios::out), _os(of)
{
  insert_catalogue();
}

PS_Stream::PS_Stream(float H, std::ostream& os, OutputMode mode)
   :_bbox(PS_BBox(-2,-2,2,2)),_mode(mode), _height((int)(H)),_os(os)
{}

PS_Stream::PS_Stream(float H, const char* fname, OutputMode mode)
  :_bbox(PS_BBox(-2,-2,2,2)),_mode(mode),_height((int)(H)),
  of(fname,std::ios::out), _os(of)
{}

PS_Stream::PS_Stream(const PS_BBox& bb, std::ostream& os, OutputMode mode)
   :_bbox(bb), _mode(mode),
   _width((int)(21*CM)), _height((int)(29.7*CM)), _os(os)
 {
   set_scale(bb);
   insert_catalogue();
 }

PS_Stream::PS_Stream(const PS_BBox& bb, const char* fname, OutputMode mode)
  : _bbox(bb),_mode(mode),_width((int)(21*CM)),
    _height((int)(29.7*CM)), of(fname,std::ios::out), _os(of)
{
  set_scale(bb);
  insert_catalogue();
}

PS_Stream::PS_Stream(const PS_BBox& bb,float H, std::ostream& os,
                               OutputMode mode)
  : _bbox(bb), _mode(mode), _height((int)H), _os(os)
{
  set_window(bb,H);
  set_scale(bb);
  insert_catalogue();
}

PS_Stream::PS_Stream(const PS_BBox& bb,float H, const char* fname,
                               OutputMode m)
  : _bbox(bb), _mode(m), _height((int)H), of(fname,std::ios::out), _os(of)
{
  set_window(bb,H);
  set_scale(bb);
  insert_catalogue();
}

PS_Stream::PS_Stream(const PS_BBox& bb,float L, float H,
                               std::ostream& os, OutputMode mode)
  : _bbox(bb), _mode(mode), _width((int)L), _height((int)H), _os(os)
{
  set_scale(bb);
  insert_catalogue();
}

PS_Stream::PS_Stream(const PS_BBox& bb,float L, float H,
                               const char* fname, OutputMode mode)
  : _bbox(bb), _mode(mode), _width((int)L), _height((int)H),
    of(fname,std::ios::out), _os(of)
{
  set_scale(bb);
  insert_catalogue();
}

PS_Stream::~PS_Stream()
{
  List_Label tmp=list();
  if (!list().empty())
    {
      os() << "%%}\\makeatletter\\let\\@notdefinable\\relax" <<std::endl;
      os() << "%%\\def\\IPEc#1[#2]#3{\\newcommand{#1}[#2]{#3}\\ignorespaces}\\@ifundefined" <<std::endl;
      os() << "%%{selectfont}{\\let\\selectfont\\relax\\def\\fontsize#1#2{}}{}\\makeatother" <<std::endl;
      os() << "%%\\IPEc\\IPEput[4]{\\put(0,0){\\special{psfile=\\IPEfile}}}" <<std::endl;
      os() << "%%\\IPEc\\IPEmp[2]{\\minipage[t]{#1bp}#2\\special{color pop}\\endminipage}" <<std::endl;
      os() << "%%\\IPEc\\IPEtext[1]{\\makebox(0,0)[lb]{#1\\special{color pop}}}" <<std::endl;
      os() << "%%\\IPEc\\IPEfs[1]{\\IPEcolfs{0 0 0}{#1}}" <<std::endl;
      os() << "%%\\IPEc\\IPEcolfs[2]{\\dimen0=#2pt\\fontsize{#2}{1.2\\dimen0}\\selectfont" <<std::endl;
      os() << "%%\\special{color push rgb #1}}" <<std::endl;
      os() << "%%\\IPEc\\IPEsize[2]{\\unitlength1bp\\ignorespaces}" <<std::endl;

      os() << "%%\\IPEsize{" << width() << "}{"<< height() << "}" <<std::endl;
      os() << "%%\\begin{picture}(" << width() << "," << height() << ")(0,0)" <<std::endl;
      os() << "%%\\IPEput{0}{0}{" << width() << "}{" << height() << "}" <<std::endl;

       while (!tmp.empty())
        {
          os() << "%%\\put(" << tmp.front().xpos() 
               << "," <<tmp.front().ypos() 
               <<"){\\IPEtext{\\IPEfs{10}\\rm "<< tmp.front().text() 
               <<" }}" << std::endl;
          tmp.pop_front();
        }
      os() << "%%\\end{picture}\\endinput}" <<std::endl;
    }
  os() << "showpage\nend" << std::endl;
  std::cout << std::flush;
}

#ifndef PS_MANIP_DEF
#define PS_MANIP_DEF

PS_Manipulator_creator<const Color&>
  border_color(&PS_Stream::set_border_color);

PS_Manipulator_creator<const Color&>
  fill_color(&PS_Stream::set_fill_color);

PS_Manipulator_creator<unsigned int>
  point_size(&PS_Stream::set_point_size);

PS_Manipulator_creator<PS_Stream::DotStyle>
 point_style(&PS_Stream::set_point_style);

PS_Manipulator_creator<DashStyle>
 line_style(&PS_Stream::set_line_style);

PS_Manipulator_creator<unsigned int>
 line_width(&PS_Stream::set_line_width);

PS_Manipulator_creator<bool>
 fill(&PS_Stream::set_fill);

PS_Manipulator_creator<const PS_Stream::Context&>
  current_context(&PS_Stream::set_current_context);


PS_Manipulator_creator<Point_2< Cartesian <double> > >
  move_to(&PS_Stream::set_point);

PS_Manipulator_creator<PS_Stream::Axis&>
  show_axis(&PS_Stream::set_axis);

PS_Manipulator_creator<PS_Stream::Grid&>
  show_grid(&PS_Stream::set_grid);

PS_Manipulator_creator<const char*>
 ps_label(&PS_Stream::put_ps_label);

PS_Manipulator_creator<const char*>
 latex_label(&PS_Stream::put_latex_label);

PS_Manipulator_creator<unsigned int>
  border(&PS_Stream::put_border);

PS_Manipulator_creator<const char*>
  font(&PS_Stream::set_font);

PS_Manipulator_creator<unsigned int>
  font_size(&PS_Stream::set_font_size);

#endif  //PS_MANIP_DEF

PS_Stream& PS_Stream::set_border_color(const Color& color)
{
  if (ctxt.get_border_color()!=color)
    {
      
      os() << color 
           << " setrgbcolor" <<std::endl;
      ctxt.set_border_color(color);
    }
  return *this;
}

PS_Stream& PS_Stream::set_fill_color (const Color& color)
{
  ctxt.set_fill_color(color);
  return *this;
}

PS_Stream& PS_Stream::set_point_size(unsigned int Size)
{
   ctxt.set_dot_size(Size);
  return *this;
}

PS_Stream& PS_Stream::set_line_width(unsigned int Width)
{
  if (ctxt.get_thickness()!=Width)
    {
      os()<< Width <<" setlinewidth" <<std::endl;
      ctxt.set_thickness(Width);
    }
  return *this;
}

PS_Stream& PS_Stream::set_point_style(DotStyle Style)
{
  ctxt.set_dot_style(Style);
  return *this;
}

PS_Stream& PS_Stream::set_line_style(DashStyle style)
{
  if (CGAL_CLIB_STD::strcmp(ctxt.get_line_style(),style))
    {
      ctxt.set_line_style(style);
      os() << style << " setdash" << std::endl;
    }
  return *this;
}

PS_Stream& PS_Stream::set_fill(bool test)
{
  ctxt.set_fill(test);
  return *this;
}

PS_Stream& PS_Stream::set_default_context(void)
{
  setdefault();
  return *this;
}

PS_Stream& PS_Stream::set_current_context(const PS_Stream::Context& c)
{
  if (ctxt.get_border_color()!=c.get_border_color())
    os()<< c.get_border_color().r() << " "
        << c.get_border_color().g() << " "
        << c.get_border_color().b() << " setrgbcolor"<<std::endl;
  if (CGAL_CLIB_STD::strcmp(ctxt.get_line_style(),c.get_line_style()))
    os() << c.get_line_style() << " setdash"<<std::endl;
  if (ctxt.get_thickness()!=c.get_thickness())
    os() << c.get_thickness() << " setlinewidth"<<std::endl;
  if (ctxt.get_font_size()!=c.get_font_size() ||
      CGAL_CLIB_STD::strcmp(ctxt.get_font(),c.get_font())!=0)
    {
      os() << "/" << c.get_font() << " findfont" <<std::endl;
      os() << c.get_font_size() << " scalefont setfont" << std::endl;
    }
  ctxt=c;
  return *this;
}

PS_Stream& PS_Stream::set_point(Point_2< Cartesian <double> > p)
{
  ctxt.set_current_pos(p);
  return *this;
}

PS_Stream& PS_Stream::set_axis(Axis& g)
{
  static bool test=false;
  double x0=x2ps(0);
  double y0=y2ps(0);
  double i;
  if (!test)
    {
      os() << "gsave 0 setgray " << g.thickness()
           << " setlinewidth" << std::endl;
      os() << "[] 0 setdash" << std::endl;
      os() << x0 << " " << 0 << " mt" <<std::endl;
      os() << x0 << " " << height() << " lt st" << std::endl;
      os() << 0 << " " << y0 << " mt" <<std::endl;
      os() << width() << " " << y0 << " lt st" <<std::endl;
      if (g.stepx())
        {
          for (i=((int) (bbox().xmin() / g.stepx())) *g.stepx();
               i<= bbox().xmax();
               i+=g.stepx()){
            double x=x2ps(i);
            os() << x << " " << y0 << " mt" <<std::endl;
            os() << x << " " << y0+2 << " lt st" << std::endl;
          }
        }
      if (g.stepy())
        {
          for (i=((int) (bbox().ymin() / g.stepy())) *g.stepy();
               i<=bbox().ymax();
               i+=g.stepy()){
            double y=y2ps(i);
            os() << x0 << " " << y << " mt" <<std::endl;
            os() << x0+2 << " " << y << " lt st" << std::endl;
          }
        }
      os() << "grestore" <<std::endl;
    }
  test=true;
  return *this;
}

PS_Stream& PS_Stream::set_grid(Grid& g)
{
  double i;
  os() << "gsave 0 setgray 0 setlinewidth" << std::endl;
  os() << g.style() << " setdash" << std::endl;
  if (g.stepx())
    {
      for (i=((int) (bbox().xmin() / g.stepx())) *g.stepx();
           i<=bbox().xmax();
           i+=g.stepx()){
        double x=x2ps(i);
        os() << x << " 0 mt" <<std::endl;
        os() << x << " " << height() << " lt st" << std::endl;
      }
    }
  if (g.stepy())
    {
      for (i=((int) (bbox().ymin() / g.stepy())) *g.stepy();
           i<=bbox().ymax();
           i+=g.stepy()){
        double y=y2ps(i);
        os() << " 0 " << y << " mt" <<std::endl;
        os() << width() << " " << y << " lt st" << std::endl;
      }
    }
  os() << "grestore" <<std::endl;
  return *this;
}

PS_Stream& PS_Stream::put_ps_label(const char* ch)
{
  os() << x2ps(context().get_pos().x()) << " "
       << y2ps(context().get_pos().y()) << " mt" <<std::endl;

  os() << "(" << ch << ") show" << std::endl;
  return *this;
}

PS_Stream& PS_Stream::put_latex_label(const char* ch)
{
  os() << "%% CGAL - LATEX : " << x2ps(context().get_pos().x()) << " "
    << y2ps(context().get_pos().y()) << " " << ch << std::endl;
  Latex_Label l= Latex_Label(ch,x2ps(context().get_pos().x()),y2ps(context().get_pos().y()));
  
 list().push_front(l);
  
  return *this;
}

PS_Stream& PS_Stream::put_border(unsigned int i)
{
  os() << "gsave" << std::endl;
  os() << "0 setgray [] 0 setdash" << std::endl;
  os() << i << " setlinewidth" << std::endl;
  os() << "0 0 " << width() << " " << height() << " rectstroke" << std::endl;
  os() << "grestore" << std::endl;
  return *this;
}

PS_Stream& PS_Stream::set_font(const char* ch)
{
  if (CGAL_CLIB_STD::strcmp(ch,context().get_font())!=0)
    {
      os() << "/" << ch << " findfont" << std::endl;
      os() << context().get_font_size() << " scalefont setfont" << std::endl;
      ctxt.set_font(ch);
    }
  return *this;
}

PS_Stream& PS_Stream::set_font_size(unsigned int i)
{
  if (context().get_font_size()!=i)
    {
      ctxt.set_font_size(i);
      os() << "/" << context().get_font() << " findfont" << std::endl;
      os() << i << " scalefont setfont" <<std::endl;
    }
  return *this;
}

void PS_Stream::setdefault()
{
  if (ctxt.get_border_color()!=CTXT_DEFAULT.get_border_color())
    os()<<"0 0 0 setrgbcolor"<<std::endl;
  if (ctxt.get_line_style()!=CTXT_DEFAULT.get_line_style())
    os() << PS_Stream::SOLID << " setdash"<<std::endl;
  if (ctxt.get_thickness()!=CTXT_DEFAULT.get_thickness())
    os() << 0 << " setlinewidth"<<std::endl;
  if (ctxt.get_font_size()!=CTXT_DEFAULT.get_font_size() ||
      CGAL_CLIB_STD::strcmp(ctxt.get_font(),CTXT_DEFAULT.get_font())!=0)
    {
      os() << "/Helvetica findfont" <<std::endl;
      os() << "12 scalefont setfont" << std::endl;
    }
  ctxt=CTXT_DEFAULT;
}

bool PS_Stream::is_eps()
{
  return (bool)(mode()==QUIET_EPS || mode()==READABLE_EPS);
}

bool PS_Stream::is_readable()
{
  return (bool)(mode()==READABLE || mode()==READABLE_EPS);
}

void PS_Stream::insert_catalogue()
{
  if (is_eps())
    {
      os() << "%!PS-Adobe-3.0 EPSF 3.0" << std::endl;
      os() << "%%BoundingBox: " << "0 0 "
           << width() << " " << height()<< std::endl;
      os() << "%%Creator: PS_Stream" << std::endl;
      os() << "%%Title: (CGAL Output)" << std::endl;
      os() << "%%CreationDate:" << std::endl;
    }
  else
    {
      os() << "%!PS-Adobe-3.0" << std::endl;
    }
  os() << "%%EndComments" <<std::endl<<std::endl;

  // The next line is used to include Latex commands in the file.
  // Thanks to this, it is possible to insert labels in the latex style.
  os() << "{\\catcode37=9\\def\\IPEdummy{({{)}} pop" <<std::endl;
  os() << "%% Ipe postscript prologue" << std::endl<<std::endl;

  os() << "/PS_Dict 14 dict def" << std::endl;
  os() << "PS_Dict begin" << std::endl;
  os() << "/lt {lineto} bind def" << std::endl;
  os() << "/mt {moveto} bind def" << std::endl;
  os() << "/st {stroke} bind def" << std::endl;
  os() << "/slw {setlinewidth} bind def" << std::endl;
  os() << "/box {/siz exch def /yy "
       << "exch def /xx exch def xx siz sub yy siz sub siz 2 mul dup} bind def" << std::endl;
  os() << "/fb {box rectfill} bind def" << std::endl;
  os() << "/eb {gsave box 4 copy gsave 1 setgray rectfill grestore "
       << "[] 0 setdash 0 setlinewidth rectstroke grestore} bind def" << std::endl;
  os() << "/xc {gsave [] 0 setdash 0 setlinewidth "
       << "/siz exch def /yy exch def /xx exch def "
       << "xx siz sub yy siz sub mt xx siz add yy siz add lineto stroke "
       << "xx siz sub yy siz add mt xx siz add yy siz sub lineto "
       << "stroke grestore} bind def" << std::endl;
  os() << "/ic {gsave [] 0 setdash 0 setlinewidth "
       << "/siz exch def /yy exch def /xx exch def "
       << "xx siz sub yy mt xx siz add yy lineto stroke "
       << "xx yy siz add mt xx yy siz sub "
       <<"lineto stroke grestore} bind def" << std::endl;
  os() << "/cir {0 360 arc} bind def" << std::endl;
  os() << "/ec {gsave 3 copy gsave 1 setgray cir fill grestore "
       << "[] 0 setdash 0 setlinewidth cir stroke grestore} bind def" << std::endl;
  os() << "/fc {cir fill} bind def" << std::endl;
  os() << "/sc {setrgbcolor} bind def" << std::endl;
  os() << "/tr {mt lt lt lt} bind def" << std::endl;
  os() << "/re {mt lt lt lt lt} bind def" << std::endl;

/*************************************************************/
//Rajout pour dessiner les aretes 

//Le stroke prend la couleur courante donc a definir avant de
//dessiner l'arete 

os() << "%Syntaxe xa ya xb yb arete" << std::endl; 
os() <<	"/arete {" << std::endl;
os() << "gsave" << std::endl;
os() << "/yb exch def " << std::endl;
os() << "/xb exch def " << std::endl;
os() << "/ya exch def " << std::endl;
os() << "/xa exch def " << std::endl;
os() << "xa ya moveto" << std::endl;
os() << "xb yb lineto" << std::endl;
os() << "closepath  " << std::endl;
os() << "stroke" << std::endl;
os() << "grestore" << std::endl;
os() << "} def" << std::endl;

//Rajout pour dessiner les faces 

os() << "%Syntaxe pt1x pt1y pt2x pt2y .. ptnx ptny nb_points face" << std::endl; 
os() <<	"/face {" << std::endl;
os() << std::endl;
os() << "/nbiter exch def" << std::endl;
os() << "newpath" << std::endl;
os() << "" << std::endl;
os() << "/ptfinaly exch def" << std::endl;
os() << "/ptfinalx exch def" << std::endl;
os() << "ptfinalx ptfinaly moveto " << std::endl;
os() << "/nbiter nbiter 1 sub def" << std::endl;
os() << std::endl;
os() << "nbiter {" << std::endl;
os() << "/ptay exch def " << std::endl;
os() << "/ptax exch def " << std::endl;
os() << "ptax ptay lineto" << std::endl;
os() << "} repeat" << std::endl;
os() << "closepath" << std::endl;
os() << std::endl;
os() << "} def" << std::endl; 

/*************************************************************/ 

  os() << "0 0 0 setrgbcolor"<<std::endl;
  os() << PS_Stream::SOLID << " setdash"<<std::endl;
  os() << 0 << " setlinewidth"<<std::endl;
  os() << "/Helvetica findfont" <<std::endl;
  os() << "12 scalefont setfont" << std::endl;
  os() << "0 0 " << width() << " " << height() << " rectclip" << std::endl;

  setdefault();
 }


#ifndef _PS_LABEL_
#define _PS_LABEL_

PS_Stream& operator << (PS_Stream& ps, const PS_Stream::Border& b)
{
  ps.os() << "gsave" << std::endl;
  ps.os() << "0 setgray [] 0 setdash" << std::endl;
  ps.os() << b.size() << " setlinewidth" << std::endl;
  ps.os() << "0 0 " << ps.width() << " " << ps.height() << " rectstroke" << std::endl;
  ps.os() << "grestore" << std::endl;
  return ps;
}

PS_Stream& operator << (PS_Stream& ps, const PS_Stream::Label& txt)
{
  ps.os() << ps.x2ps(ps.context().get_pos().x()) << " "
       << ps.y2ps(ps.context().get_pos().y()) << " mt" <<std::endl;

  ps.os() << "(" << txt.text() << ") show" << std::endl;
  return ps;
}

PS_Stream& operator << (PS_Stream& ps, PS_Stream::Latex_Label& txt)
{
  txt.setposition(ps.x2ps(ps.context().get_pos().x()),
                  ps.y2ps(ps.context().get_pos().y()));
  txt.text();
  ps.list().push_front(txt);
  return ps;
}

PS_Stream& operator << (PS_Stream& ps, const PS_Stream::Grid& g)
{
  double i;
  ps.os() << "gsave 0 setgray 0 setlinewidth" << std::endl;
  ps.os() << g.style() << " setdash" << std::endl;
  if (g.stepx())
    {
      for (i=((int) (ps.bbox().xmin() / g.stepx())) *g.stepx();
           i<=ps.bbox().xmax();
           i+=g.stepx()){
        double x=ps.x2ps(i);
        ps.os() << x << " 0 mt" <<std::endl;
        ps.os() << x << " " << ps.height() << " lt st" << std::endl;
      }
    }
  if (g.stepy())
    {
      for (i=((int) (ps.bbox().ymin() / g.stepy())) *g.stepy();
           i<=ps.bbox().ymax();
           i+=g.stepy()){
        double y=ps.y2ps(i);
        ps.os() << " 0 " << y << " mt" <<std::endl;
        ps.os() << ps.width() << " " << y << " lt st" << std::endl;
      }
    }
  ps.os() << "grestore" <<std::endl;
  return ps;
}

PS_Stream& operator << (PS_Stream& ps,const PS_Stream::Axis& g)
{
  static bool test=false;
  double x0=ps.x2ps(0);
  double y0=ps.y2ps(0);
  double i;
  if (!test)
    {
      ps.os() << "gsave 0 setgray " << g.thickness()
           << " setlinewidth" << std::endl;
      ps.os() << "[] 0 setdash" << std::endl;
      ps.os() << x0 << " " << 0 << " mt" <<std::endl;
      ps.os() << x0 << " " << ps.height() << " lt st" << std::endl;
      ps.os() << 0 << " " << y0 << " mt" <<std::endl;
      ps.os() << ps.width() << " " << y0 << " lt st" <<std::endl;
      if (g.stepx())
        {
          for (i=((int) (ps.bbox().xmin() / g.stepx())) *g.stepx();
               i<= ps.bbox().xmax();
               i+=g.stepx()){
            double x=ps.x2ps(i);
            ps.os() << x << " " << y0 << " mt" <<std::endl;
            ps.os() << x << " " << y0+2 << " lt st" << std::endl;
          }
        }
      if (g.stepy())
        {
          for (i=((int) (ps.bbox().ymin() / g.stepy())) *g.stepy();
               i<=ps.bbox().ymax();
               i+=g.stepy()){
            double y=ps.y2ps(i);
            ps.os() << x0 << " " << y << " mt" <<std::endl;
            ps.os() << x0+2 << " " << y << " lt st" << std::endl;
          }
        }
      ps.os() << "grestore" <<std::endl;
    }
  test=true;
  return ps;
}

#endif

CGAL_END_NAMESPACE
