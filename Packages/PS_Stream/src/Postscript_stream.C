
#include <CGAL/IO/Postscript_stream.h>


const float CGAL_PS_Stream::CM=28.4;
const float CGAL_PS_Stream::INCH=72.0;
const float CGAL_PS_Stream::POINT=1.0;

const DashStyle CGAL_PS_Stream::SOLID="[] 0 ";
const DashStyle CGAL_PS_Stream::DASH2="[5 5] 0 ";
const DashStyle CGAL_PS_Stream::DASH3="[10 10] 0 ";
const DashStyle CGAL_PS_Stream::DASH6="[10 6 4 6] 0 ";
const DashStyle CGAL_PS_Stream::DASH4="[10 5] 0 ";
const DashStyle CGAL_PS_Stream::DASH5="[5 10] 0 ";
const DashStyle CGAL_PS_Stream::DASH1="[2 2] 0 ";
extern const CGAL_PS_Stream::Context CGAL_CTXT_DEFAULT=CGAL_PS_Stream::Context();

/*
  CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb)
  : _bbox(bb), _os(cerr), _mode(GS_VIEW)
  {
  FILE *fp = popen("gs -","w");
  if(!fp){
  cerr << "Could not open pipe to gs" << endl ;
  exit(-1);
  }
  os().attach(fileno(fp));
  insert_catalogue();
  }
  CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb,float L, float H)
  : _bbox(bb), _os(cerr), _mode(GS_VIEW), _width(L), _height(H)
  {
  insert_catalogue();
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  };

*/  


CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb, ostream &os,
                               OutputMode mode)
  :_bbox(bb), _mode(mode),
  _width((int)21*CM), _height((int)29.7*CM), _os(cerr)
{
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  _os=os;
  insert_catalogue();
}

CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb, const char *fname,
                               OutputMode mode)
  : _bbox(bb),_mode(mode),_width((int)21*CM), _height((int)29.7*CM),_os(cerr) 
{
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  static ofstream os(fname);
  _os=os;
  insert_catalogue();
}

CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb,float H, ostream &os,
                               OutputMode mode)
  : _bbox(bb), _mode(mode), _height(H), _os(cerr)
{
  _width=(bb.xmax()-bb.xmin())*H/(bb.ymax()-bb.ymin());
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  _os=os;
  insert_catalogue();
}

CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb,float H, const char *fname,
                               OutputMode m)
  : _bbox(bb), _mode(m), _height(H), _os(cerr)
{
  static ofstream os(fname);
  _os=os;
  _width=(bb.xmax()-bb.xmin())*H/(bb.ymax()-bb.ymin());
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  insert_catalogue();
}

CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb,float L, float H,
                               ostream &os, OutputMode mode)
  : _bbox(bb), _mode(mode), _width(L), _height(H), _os(cerr)
{
  _os=os;
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  insert_catalogue();
}

CGAL_PS_Stream::CGAL_PS_Stream(const PS_BBox &bb,float L, float H,
                               const char *fname, OutputMode mode)
  : _bbox(bb), _mode(mode),_width(L),_height(H), _os(cerr)
{
  static ofstream os(fname);
  _os=os;
  _xratio=_width/(_bbox.xmax()-_bbox.xmin());
  _yratio=_height/(_bbox.ymax()-_bbox.ymin());
  insert_catalogue();
}

CGAL_PS_Stream::~CGAL_PS_Stream()
{
  List_Label tmp=list();
  if (!list().empty())
    {
      os() << "%%}\\makeatletter\\let\\@notdefinable\\relax" <<endl;
      os() << "%%\\def\\IPEc#1[#2]#3{\\newcommand{#1}[#2]{#3}\\ignorespaces}\\@ifundefined" <<endl;
      os() << "%%{selectfont}{\\let\\selectfont\\relax\\def\\fontsize#1#2{}}{}\\makeatother" <<endl;
      os() << "%%\\IPEc\\IPEput[4]{\\put(0,0){\\special{psfile=\\IPEfile}}}" <<endl;
      os() << "%%\\IPEc\\IPEmp[2]{\\minipage[t]{#1bp}#2\\special{color pop}\\endminipage}" <<endl;
      os() << "%%\\IPEc\\IPEtext[1]{\\makebox(0,0)[lb]{#1\\special{color pop}}}" <<endl;
      os() << "%%\\IPEc\\IPEfs[1]{\\IPEcolfs{0 0 0}{#1}}" <<endl;
      os() << "%%\\IPEc\\IPEcolfs[2]{\\dimen0=#2pt\\fontsize{#2}{1.2\\dimen0}\\selectfont" <<endl;
      os() << "%%\\special{color push rgb #1}}" <<endl;
      os() << "%%\\IPEc\\IPEsize[2]{\\unitlength1bp\\ignorespaces}" <<endl;

      os() << "%%\\IPEsize{" << width() << "}{"<< height() << "}" <<endl;
      os() << "%%\\begin{picture}(" << width() << "," << height() << ")(0,0)" <<endl;
      os() << "%%\\IPEput{0}{0}{" << width() << "}{" << height() << "}" <<endl;

      while (!tmp.empty())
        {
          os() << "%%\\put(" << tmp.front().xpos() 
               << "," <<tmp.front().ypos() 
               <<"){\\IPEtext{\\IPEfs{10}\\rm "<< tmp.front().text() 
               <<" }}" << endl;
          tmp.pop_front();
        }
      os() << "%%\\end{picture}\\endinput}" <<endl;
    }
  os() << "showpage\nend" << endl ;
  cout << flush;
}



#ifndef CGAL_PS_MANIP_DEF
#define CGAL_PS_MANIP_DEF

CGAL_PS_Modifier_creator<CGAL_Color &>
  set_border_color(&CGAL_PS_Stream::_SetBorderColor);

CGAL_PS_Modifier_creator<CGAL_Color &>
  set_fill_color(&CGAL_PS_Stream::_SetFillColor);

CGAL_PS_Modifier_creator<unsigned int>
  set_point_size(&CGAL_PS_Stream::_SetPointSize);

CGAL_PS_Modifier_creator<CGAL_PS_Stream::DotStyle>
  set_point_style(&CGAL_PS_Stream::_SetPointStyle);

CGAL_PS_Modifier_creator<DashStyle>
  set_line_style(&CGAL_PS_Stream::_SetLineStyle);

CGAL_PS_Modifier_creator<unsigned int>
  set_line_width(&CGAL_PS_Stream::_SetLineWidth);

CGAL_PS_Modifier_creator<bool>
  set_fill(&CGAL_PS_Stream::_SetFill);

CGAL_PS_Modifier_creator<const CGAL_PS_Stream::Context &>
  set_current_context(&CGAL_PS_Stream::_SetCurrentContext);

CGAL_PS_Modifier_creator<bool>
  show_direction(&CGAL_PS_Stream::_ShowDirection);

CGAL_PS_Modifier_creator<CGAL_Point_2< CGAL_Cartesian <double> > >
  move_to(&CGAL_PS_Stream::_MoveTo);

CGAL_PS_Modifier_creator<CGAL_PS_Stream::Axis &>
  show_axis(&CGAL_PS_Stream::_ShowAxis);

CGAL_PS_Modifier_creator<CGAL_PS_Stream::Grid &>
  show_grid(&CGAL_PS_Stream::_ShowGrid);

CGAL_PS_Modifier_creator<const char *>
  put_ps_label(&CGAL_PS_Stream::_PutPsLabel);

CGAL_PS_Modifier_creator<const char *>
  put_latex_label(&CGAL_PS_Stream::_PutLatexLabel);

CGAL_PS_Modifier_creator<unsigned int>
  put_border(&CGAL_PS_Stream::_PutBorder);

CGAL_PS_Modifier_creator<const char *>
  set_font(&CGAL_PS_Stream::_SetFont);

CGAL_PS_Modifier_creator<unsigned int>
  set_font_size(&CGAL_PS_Stream::_SetFontSize);
#endif  //CGAL_PS_MANIP_DEF


CGAL_PS_Stream& CGAL_PS_Stream::_SetBorderColor(CGAL_Color &color)
{
  if (ctxt.get_border_color()!=color)
    {
      os() << color.r() << " " << color.g() << " " << color.b()
           << " setrgbcolor" <<endl;
      ctxt.set_border_color(color);
    }
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetFillColor(CGAL_Color &color)
{
  ctxt.set_fill_color(color);
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetPointSize(unsigned int Size)
{
  ctxt.set_dot_size(Size);
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetLineWidth(unsigned int Width)
{
  if (ctxt.get_thickness()!=Width)
    {
      os()<< Width <<" setlinewidth" <<endl;
      ctxt.set_thickness(Width);
    }
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetPointStyle(DotStyle Style)
{
  ctxt.set_dot_style(Style);
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetLineStyle(DashStyle style)
{
  if (strcmp(ctxt.get_line_style(),style))
    {
      ctxt.set_line_style(style);
      os() << style << " setdash" << endl;
    }
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetFill(bool test)
{
  ctxt.set_fill(test);
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetDefaultContext(void)
{
  setdefault();
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetCurrentContext(const CGAL_PS_Stream::Context &c)
{
  if (ctxt.get_border_color()!=c.get_border_color())
    os()<< c.get_border_color().r() << " "
        << c.get_border_color().g() << " "
        << c.get_border_color().b() << " setrgbcolor"<<endl;
  if (strcmp(ctxt.get_line_style(),c.get_line_style()))
    os() << c.get_line_style() << " setdash"<<endl;
  if (ctxt.get_thickness()!=c.get_thickness())
    os() << c.get_thickness() << " setlinewidth"<<endl;
  if (ctxt.get_font_size()!=c.get_font_size() ||
      strcmp(ctxt.get_font(),c.get_font())!=0)
    {
      os() << "/" << c.get_font() << " findfont" <<endl;
      os() << c.get_font_size() << " scalefont setfont" << endl;
    }
  ctxt=c;
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_ShowDirection(bool choice)
{
  if (choice){};
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_MoveTo(CGAL_Point_2< CGAL_Cartesian <double> > p)
{
  ctxt.set_current_pos(p);
  return *this;
}
CGAL_PS_Stream& CGAL_PS_Stream::_ShowAxis(Axis &g)
{
  static bool test=false;
  double x0=x2ps(0);
  double y0=y2ps(0);
  double i;
  if (!test)
    {
      os() << "gsave 0 setgray " << g.thickness()
           << " setlinewidth" << endl;
      os() << "[] 0 setdash" << endl;
      os() << x0 << " " << 0 << " mt" <<endl;
      os() << x0 << " " << height() << " lt st" << endl;
      os() << 0 << " " << y0 << " mt" <<endl;
      os() << width() << " " << y0 << " lt st" <<endl;
      if (g.stepx())
        {
          for (i=((int) (bbox().xmin() / g.stepx())) *g.stepx();
               i<= bbox().xmax();
               i+=g.stepx()){
            double x=x2ps(i);
            os() << x << " " << y0 << " mt" <<endl;
            os() << x << " " << y0+2 << " lt st" << endl;
          }
        }
      if (g.stepy())
        {
          for (i=((int) (bbox().ymin() / g.stepy())) *g.stepy();
               i<=bbox().ymax();
               i+=g.stepy()){
            double y=y2ps(i);
            os() << x0 << " " << y << " mt" <<endl;
            os() << x0+2 << " " << y << " lt st" << endl;
          }
        }
      os() << "grestore" <<endl;
    }
  test=true;
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_ShowGrid(Grid &g)
{
  double i;
  os() << "gsave 0 setgray 0 setlinewidth" << endl;
  os() << g.style() << " setdash" << endl;
  if (g.stepx())
    {
      for (i=((int) (bbox().xmin() / g.stepx())) *g.stepx();
           i<=bbox().xmax();
           i+=g.stepx()){
        double x=x2ps(i);
        os() << x << " 0 mt" <<endl;
        os() << x << " " << height() << " lt st" << endl;
      }
    }
  if (g.stepy())
    {
      for (i=((int) (bbox().ymin() / g.stepy())) *g.stepy();
           i<=bbox().ymax();
           i+=g.stepy()){
        double y=y2ps(i);
        os() << " 0 " << y << " mt" <<endl;
        os() << width() << " " << y << " lt st" << endl;
      }
    }
  os() << "grestore" <<endl;
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_PutPsLabel(const char *ch)
{
  os() << x2ps(context().get_pos().x()) << " "
       << y2ps(context().get_pos().y()) << " mt" <<endl;

  os() << "(" << ch << ") show" << endl;
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_PutLatexLabel(const char *ch)
{
  //os() << "%% CGAL - LATEX : " << x2ps(context().get_pos().x()) << " "
  //   << y2ps(context().get_pos().y()) << " " << ch << endl;
  
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_PutBorder(unsigned int i)
{
  os() << "gsave" << endl;
  os() << "0 setgray [] 0 setdash" << endl;
  os() << i << " setlinewidth" << endl;
  os() << "0 0 " << width() << " " << height() << " rectstroke" << endl;
  os() << "grestore" << endl;
  return *this;
}

CGAL_PS_Stream& CGAL_PS_Stream::_SetFont(const char *ch)
{
  if (strcmp(ch,context().get_font())!=0)
    {
      os() << "/" << ch << " findfont" << endl;
      os() << context().get_font_size() << " scalefont setfont" << endl;
      ctxt.set_font(ch);
    }
  return *this;
}
CGAL_PS_Stream& CGAL_PS_Stream::_SetFontSize(unsigned int i)
{
  if (context().get_font_size()!=i)
    {
      ctxt.set_font_size(i);
      os() << "/" << context().get_font() << " findfont" << endl;
      os() << i << " scalefont setfont" <<endl;
    }
  return *this;
}




void CGAL_PS_Stream::setdefault()
{
  if (ctxt.get_border_color()!=CGAL_CTXT_DEFAULT.get_border_color())
    os()<<"0 0 0 setrgbcolor"<<endl;
  if (ctxt.get_line_style()!=CGAL_CTXT_DEFAULT.get_line_style())
    os() << CGAL_PS_Stream::SOLID << " setdash"<<endl;
  if (ctxt.get_thickness()!=CGAL_CTXT_DEFAULT.get_thickness())
    os() << 0 << " setlinewidth"<<endl;
  if (ctxt.get_font_size()!=CGAL_CTXT_DEFAULT.get_font_size() ||
      strcmp(ctxt.get_font(),CGAL_CTXT_DEFAULT.get_font())!=0)
    {
      os() << "/Helvetica findfont" <<endl;
      os() << "12 scalefont setfont" << endl;
    }
  ctxt=CGAL_CTXT_DEFAULT;
}

bool CGAL_PS_Stream::is_eps()
{
  return (bool)(mode()==QUIET_EPS || mode()==READABLE_EPS);
}

bool CGAL_PS_Stream::is_readable()
{
  return (bool)(mode()==READABLE || mode()==READABLE_EPS);
}

void CGAL_PS_Stream::insert_catalogue()
{
  if (is_eps())
    {
      os() << "%!PS-Adobe-3.0 EPSF 3.0" << endl;
      os() << "%%BoundingBox: " << "0 0 "
           << width() << " " << height()<< endl;
      os() << "%%Creator: CGAL_PS_Stream" << endl;
      os() << "%%Title: (CGAL Output)" << endl;
      os() << "%%CreationDate:" << endl;
    }
  else
    {
      os() << "%!PS-Adobe-3.0" << endl;
    }
  os() << "%%EndComments" <<endl<<endl;

  // The next line is used to include Latex commands in the file.
  // Thanks to this, it is possible to insert labels in the latex style.
  os() << "{\\catcode37=9\\def\\IPEdummy{({{)}} pop" <<endl;
  os() << "%% Ipe postscript prologue" << endl<<endl;

  os() << "/CGAL_PS_Dict 14 dict def" << endl;
  os() << "CGAL_PS_Dict begin" << endl;
  os() << "/lt {lineto} bind def" << endl;
  os() << "/mt {moveto} bind def" << endl;
  os() << "/st {stroke} bind def" << endl;
  os() << "/slw {setlinewidth} bind def" << endl;
  os() << "/box {/siz exch def /yy "
       << "exch def /xx exch def xx siz sub yy siz sub siz 2 mul dup} bind def" << endl;
  os() << "/fb {box rectfill} bind def" << endl;
  os() << "/eb {gsave box 4 copy gsave 1 setgray rectfill grestore "
       << "[] 0 setdash 0 setlinewidth rectstroke grestore} bind def" << endl;
  os() << "/xc {gsave [] 0 setdash 0 setlinewidth "
       << "/siz exch def /yy exch def /xx exch def "
       << "xx siz sub yy siz sub mt xx siz add yy siz add lineto stroke "
       << "xx siz sub yy siz add mt xx siz add yy siz sub lineto "
       << "stroke grestore} bind def" << endl;
  os() << "/ic {gsave [] 0 setdash 0 setlinewidth "
       << "/siz exch def /yy exch def /xx exch def "
       << "xx siz sub yy mt xx siz add yy lineto stroke "
       << "xx yy siz add mt xx yy siz sub "
       <<"lineto stroke grestore} bind def" << endl;
  os() << "/cir {0 360 arc} bind def" << endl;
  os() << "/ec {gsave 3 copy gsave 1 setgray cir fill grestore "
       << "[] 0 setdash 0 setlinewidth cir stroke grestore} bind def" << endl;
  os() << "/fc {cir fill} bind def" << endl;
  os() << "/sc {setrgbcolor} bind def" << endl;
  os() << "/tr {mt lt lt lt} bind def" << endl;
  os() << "/re {mt lt lt lt lt} bind def" << endl;
  os() << "0 0 0 setrgbcolor"<<endl;
  os() << CGAL_PS_Stream::SOLID << " setdash"<<endl;
  os() << 0 << " setlinewidth"<<endl;
  os() << "/Helvetica findfont" <<endl;
  os() << "12 scalefont setfont" << endl;
  os() << "0 0 " << width() << " " << height() << " rectclip" << endl;

  setdefault();
}


#ifndef _PS_LABEL_
#define _PS_LABEL_

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Border &b)
{
  ps.os() << "gsave" << endl;
  ps.os() << "0 setgray [] 0 setdash" << endl;
  ps.os() << b.size() << " setlinewidth" << endl;
  ps.os() << "0 0 " << ps.width() << " " << ps.height() << " rectstroke" << endl;
  ps.os() << "grestore" << endl;
  return ps;
}

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Label &txt)
{
  ps.os() << ps.x2ps(ps.context().get_pos().x()) << " "
       << ps.y2ps(ps.context().get_pos().y()) << " mt" <<endl;

  ps.os() << "(" << txt.text() << ") show" << endl;
  return ps;
}

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, CGAL_PS_Stream::Latex_Label &txt)
{
  txt.setposition(ps.x2ps(ps.context().get_pos().x()),
                  ps.y2ps(ps.context().get_pos().y()));
  txt.text();
  ps.list().push_front(txt);
  return ps;
}

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Grid &g)
{
  double i;
  ps.os() << "gsave 0 setgray 0 setlinewidth" << endl;
  ps.os() << g.style() << " setdash" << endl;
  if (g.stepx())
    {
      for (i=((int) (ps.bbox().xmin() / g.stepx())) *g.stepx();
           i<=ps.bbox().xmax();
           i+=g.stepx()){
        double x=ps.x2ps(i);
        ps.os() << x << " 0 mt" <<endl;
        ps.os() << x << " " << ps.height() << " lt st" << endl;
      }
    }
  if (g.stepy())
    {
      for (i=((int) (ps.bbox().ymin() / g.stepy())) *g.stepy();
           i<=ps.bbox().ymax();
           i+=g.stepy()){
        double y=ps.y2ps(i);
        ps.os() << " 0 " << y << " mt" <<endl;
        ps.os() << ps.width() << " " << y << " lt st" << endl;
      }
    }
  ps.os() << "grestore" <<endl;
  return ps;
}

CGAL_PS_Stream &operator<<(CGAL_PS_Stream &ps, const CGAL_PS_Stream::Axis &g)
{
  static bool test=false;
  double x0=ps.x2ps(0);
  double y0=ps.y2ps(0);
  double i;
  if (!test)
    {
      ps.os() << "gsave 0 setgray " << g.thickness()
           << " setlinewidth" << endl;
      ps.os() << "[] 0 setdash" << endl;
      ps.os() << x0 << " " << 0 << " mt" <<endl;
      ps.os() << x0 << " " << ps.height() << " lt st" << endl;
      ps.os() << 0 << " " << y0 << " mt" <<endl;
      ps.os() << ps.width() << " " << y0 << " lt st" <<endl;
      if (g.stepx())
        {
          for (i=((int) (ps.bbox().xmin() / g.stepx())) *g.stepx();
               i<= ps.bbox().xmax();
               i+=g.stepx()){
            double x=ps.x2ps(i);
            ps.os() << x << " " << y0 << " mt" <<endl;
            ps.os() << x << " " << y0+2 << " lt st" << endl;
          }
        }
      if (g.stepy())
        {
          for (i=((int) (ps.bbox().ymin() / g.stepy())) *g.stepy();
               i<=ps.bbox().ymax();
               i+=g.stepy()){
            double y=ps.y2ps(i);
            ps.os() << x0 << " " << y << " mt" <<endl;
            ps.os() << x0+2 << " " << y << " lt st" << endl;
          }
        }
      ps.os() << "grestore" <<endl;
    }
  test=true;
  return ps;
}

#endif

