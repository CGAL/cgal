// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_visualizor_OGL.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_visualizor_OGL.h    visualization of SNCs in a OpenGLUT viewer
// ============================================================================
#ifndef CGAL_SNC_VISUALIZOR_OGL_H
#define CGAL_SNC_VISUALIZOR_OGL_H

#include <CGAL/basic.h>
#include <CGAL/IO/Color.h>
#include <GL/glut.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>

//#undef TRACEN
//#define TRACEN(t) std::cerr << t << std::endl

#define CGAL_NEF3_MARKED_VERTEX_COLOR 221,221,255
#define CGAL_NEF3_UNMARKED_VERTEX_COLOR 15,15,31
#define CGAL_NEF3_MARKED_EDGE_COLOR 221,221,255
#define CGAL_NEF3_UNMARKED_EDGE_COLOR 15,15,31
#define CGAL_NEF3_MARKED_FACET_COLOR 221,255,255
#define CGAL_NEF3_UNMARKED_FACET_COLOR 31,31,63

CGAL_BEGIN_NAMESPACE
namespace OGL {

// ----------------------------------------------------------------------------
// Drawable double types:
// ----------------------------------------------------------------------------

  typedef CGAL::Simple_cartesian<double> DKernel;  
  typedef DKernel::Point_3               Double_point;
  typedef DKernel::Segment_3             Double_segment;

  // DPoint = a double point including a mark
  class DPoint : public Double_point {
    bool m_;
  public:
    DPoint() {}
    DPoint(const Double_point& p, bool m) : Double_point(p) { m_ = m; }
    DPoint(const DPoint& p) : Double_point(p) { m_ = p.m_; }
    DPoint& operator=(const DPoint& p)
    { Double_point::operator=(p); m_ = p.m_; return *this; }
    bool mark() const { return m_; }
  };

  // DSegment = a double segment including a mark
  class DSegment : public Double_segment {
    bool m_;
  public:
    DSegment() {}
    DSegment(const Double_segment& s, bool m) : Double_segment(s) { m_ = m; }
    DSegment(const DSegment& s) : Double_segment(s) { m_ = s.m_; }
    DSegment& operator=(const DSegment& s)
    { Double_segment::operator=(s); m_ = s.m_; return *this; }
    bool mark() const { return m_; }
  };

  // Double_triple = a class that stores a triple of double
  // coordinates; we need a pointer to the coordinates in a C array
  // for OpenGL
  class Double_triple {
    typedef double*       double_ptr;
    typedef const double* const_double_ptr;
    double coords_[3];
  public:
    Double_triple() 
    { coords_[0]=coords_[1]=coords_[2]=0.0; }
    Double_triple(double x, double y, double z)
    { coords_[0]=x; coords_[1]=y; coords_[2]=z; }
    Double_triple(const Double_triple& t)
    { coords_[0]=t.coords_[0];
      coords_[1]=t.coords_[1];
      coords_[2]=t.coords_[2];
    }
    Double_triple& operator=(const Double_triple& t)
    { coords_[0]=t.coords_[0];
      coords_[1]=t.coords_[1];
      coords_[2]=t.coords_[2];
      return *this; }
    operator double_ptr() const 
    { return const_cast<Double_triple&>(*this).coords_; }
    double operator[](unsigned i) 
    { CGAL_assertion(i<3); return coords_[i]; }
  }; // Double_triple

  static std::ostream& operator << (std::ostream& os,
				    const Double_triple& t)
  { os << "(" << t[0] << "," << t[1] << "," << t[2] << ")";
    return os; }


  // DFacet stores the facet cycle vertices in a continuus C array
  // of three double components, this is necessary due to the OpenGL
  // tesselator input format !
  class DFacet {
    typedef std::vector<Double_triple>   Coord_vector;
    typedef std::vector<unsigned>        Cycle_vector;
    Coord_vector    coords_;  // stores all vertex coordinates
    Cycle_vector    fc_ends_; // stores entry points of facet cycles
    Double_triple   normal_; // stores normal and mark of facet
    bool            mark_;

  public:
    typedef Coord_vector::iterator Coord_iterator;
    typedef Coord_vector::const_iterator Coord_const_iterator;

    DFacet() {}

    void push_back_vertex(double x, double y, double z)
    { coords_.push_back(Double_triple(x,y,z)); }

    DFacet(const DFacet& f) 
    { coords_  = f.coords_;
      fc_ends_ = f.fc_ends_;
      normal_  = f.normal_;
      mark_    = f.mark_;
    }

    DFacet& operator=(const DFacet& f) 
    { coords_ =  f.coords_;
      fc_ends_ = f.fc_ends_;
      normal_ =  f.normal_;
      mark_    = f.mark_;
      return *this; 
    }

    ~DFacet()
    { coords_.clear(); fc_ends_.clear(); }

    void push_back_vertex(const Double_point& p)
    { push_back_vertex(p.x(),p.y(),p.z()); }
   
    void set_normal(double x, double y, double z, bool m)
    { double l = sqrt(x*x + y*y + z*z);
      normal_ = Double_triple(x/l,y/l,z/l); mark_ = m; }

    double dx() const { return normal_[0]; }
    double dy() const { return normal_[1]; }
    double dz() const { return normal_[2]; }
    bool mark() const { return mark_; }
    double* normal() const 
    { return static_cast<double*>(normal_); }

    void new_facet_cycle()
    { fc_ends_.push_back(coords_.size()); }
    
    unsigned number_of_facet_cycles() const
    { return fc_ends_.size(); }

    Coord_iterator facet_cycle_begin(unsigned i) 
    { CGAL_assertion(i<number_of_facet_cycles());
      if (i==0) return coords_.begin();
      else return coords_.begin()+fc_ends_[i-1]; }

    Coord_iterator facet_cycle_end(unsigned i) 
    { CGAL_assertion(i<number_of_facet_cycles());
      if (i<fc_ends_.size()-1) return coords_.begin()+fc_ends_[i];
      else return coords_.end(); }

    Coord_const_iterator facet_cycle_begin(unsigned i) const
    { CGAL_assertion(i<number_of_facet_cycles());
      if (i==0) return coords_.begin();
      else return coords_.begin()+fc_ends_[i-1]; }

    Coord_const_iterator facet_cycle_end(unsigned i) const
    { CGAL_assertion(i<number_of_facet_cycles());
      if (i<fc_ends_.size()-1) return coords_.begin()+fc_ends_[i];
      else return coords_.end(); }

    void debug(std::ostream& os = std::cerr) const
    { os << "DFacet, normal=" << normal_ << ", mark=" << mark() << std::endl;
      for(unsigned i=0; i<number_of_facet_cycles(); ++i) {
	os << "  facet cycle ";
	// put all vertices in facet cycle into contour:
	Coord_const_iterator cit;
	for(cit = facet_cycle_begin(i); cit != facet_cycle_end(i); ++cit)
	  os << *cit;
	os << std::endl;
      }
    }
    
  }; // DFacet


// ----------------------------------------------------------------------------
// OGL Drawable Polyhedron:
// ----------------------------------------------------------------------------

  inline void beginCallback(GLenum which)
  { glBegin(which); }

  inline void endCallback(void)
  { glEnd(); }

  inline void errorCallback(GLenum errorCode)
  { const GLubyte *estring;
    estring = gluErrorString(errorCode);
    fprintf(stderr, "Tessellation Error: %s\n", estring);
    exit (0);
  }

  inline void vertexCallback(GLvoid* vertex,
			     GLvoid* user)
  { GLdouble* pc(static_cast<GLdouble*>(vertex));
    GLdouble* pu(static_cast<GLdouble*>(user));
    TRACEN("vertexCallback coord  "<<pc[0]<<","<<pc[1]<<","<<pc[2]);
    TRACEN("vertexCallback normal "<<pu[0]<<","<<pu[1]<<","<<pu[2]);
    glVertex3dv(pc); 
    glNormal3dv(pu);
  }

  class Polyhedron {
    std::list<DPoint>    vertices_;
    std::list<DSegment>  edges_;
    std::list<DFacet>    halffacets_;

    GLuint         object_list_;
    bool init_, axes_, surface_;

    typedef std::list<DPoint>::const_iterator   Vertex_iterator;
    typedef std::list<DSegment>::const_iterator Edge_iterator;
    typedef std::list<DFacet>::const_iterator   Halffacet_iterator;

  public:
    Polyhedron() 
    { object_list_ = 0;
      init_ = axes_ = false; surface_ = true; }

    Polyhedron(const Polyhedron& p)
    { object_list_ = 0;
      init_ = axes_ = false; surface_ = true; }

    Polyhedron& operator=(const Polyhedron& p)
    { return *this; }

    ~Polyhedron() 
    { if (object_list_) glDeleteLists(object_list_, 4); }

    void push_back(const Double_point& p, bool m) 
    { vertices_.push_back(DPoint(p,m)); }
    void push_back(const Double_segment& s, bool m) 
    { edges_.push_back(DSegment(s,m)); }
    void push_back(const DFacet& f) 
    { halffacets_.push_back(f); }
 
    void toggle_axes() { axes_ = !axes_; }
    void skeleton_on() { surface_ = false; }
    void boundary_on() { surface_ = true; }
    bool is_initialized() const { return init_; }

    void draw(Vertex_iterator v) const
    { TRACEN("drawing vertex "<<*v);
      CGAL::Color cf(CGAL_NEF3_MARKED_VERTEX_COLOR), 
	ct(CGAL_NEF3_UNMARKED_VERTEX_COLOR); // more blue-ish
      CGAL::Color c = v->mark() ? ct : cf;
      glPointSize(10);
      glColor3ub(c.red(), c.green(), c.blue());
      glBegin(GL_POINTS);
      glVertex3d(v->x(),v->y(),v->z());
      glEnd();           
    }

    void draw(Edge_iterator e) const
    { TRACEN("drawing edge "<<*e);
      Double_point p = e->source(), q = e->target();
      CGAL::Color cf(CGAL_NEF3_MARKED_EDGE_COLOR), 
	ct(CGAL_NEF3_UNMARKED_EDGE_COLOR); // more blue-ish
      CGAL::Color c = e->mark() ? ct : cf;
      glLineWidth(5);
      glColor3ub(c.red(),c.green(),c.blue());
      glBegin(GL_LINE_STRIP);
      glVertex3d(p.x(), p.y(), p.z());
      glVertex3d(q.x(), q.y(), q.z());
      glEnd();
    }

    void draw(Halffacet_iterator f) const
    { TRACEN("drawing facet "<<(f->debug(),""));
      GLUtesselator* tess_ = gluNewTess();
      gluTessCallback(tess_, GLenum(GLU_TESS_VERTEX_DATA),
		      (GLvoid (*)(...)) &vertexCallback);
      gluTessCallback(tess_, GLenum(GLU_TESS_BEGIN),
		      (GLvoid (*)(...)) &beginCallback);
      gluTessCallback(tess_, GLenum(GLU_TESS_END),
		      (GLvoid (*)(...)) &endCallback);
      gluTessCallback(tess_, GLenum(GLU_TESS_ERROR),
		      (GLvoid (*)(...)) &errorCallback);
      gluTessProperty(tess_, GLenum(GLU_TESS_WINDING_RULE),
		      GLU_TESS_WINDING_POSITIVE);

      DFacet::Coord_const_iterator cit;
      CGAL::Color cf(CGAL_NEF3_MARKED_FACET_COLOR),
	ct(CGAL_NEF3_UNMARKED_FACET_COLOR); // more blue-ish
      CGAL::Color c = (f->mark() ? ct : cf);
      glColor3ub(c.red(),c.green(),c.blue());
      gluTessBeginPolygon(tess_,f->normal());
      gluTessNormal(tess_,f->dx(),f->dy(),f->dz());
      // forall facet cycles of f:
      for(unsigned i = 0; i < f->number_of_facet_cycles(); ++i) {
        gluTessBeginContour(tess_);
	// put all vertices in facet cycle into contour:
	for(cit = f->facet_cycle_begin(i); 
	    cit != f->facet_cycle_end(i); ++cit) {
	  gluTessVertex(tess_, *cit, *cit);
	}
        gluTessEndContour(tess_);
      }
      gluTessEndPolygon(tess_);
      TRACEN("END drawing facet ");
      gluDeleteTess(tess_);
    }

    void construct_axes() const
    { 
      glLineWidth(2.0);
      // red x-axis
      glColor3f(1.0,0.0,0.0);
      glBegin(GL_LINES);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(5.0,0.0,0.0);
      glEnd();
       // green y-axis 
      glColor3f(0.0,1.0,0.0);
      glBegin(GL_LINES);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(0.0,5.0,0.0);
      glEnd();
      // blue z-axis and equator
      glColor3f(0.0,0.0,1.0);
      glBegin(GL_LINES);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(0.0,0.0,5.0);
      glEnd();
      // six coordinate points in pink:
      glPointSize(10);
      glBegin(GL_POINTS);
      glColor3f(1.0,0.0,0.0);
      glVertex3d(5,0,0);
      glColor3f(0.0,1.0,0.0);
      glVertex3d(0,5,0);
      glColor3f(0.0,0.0,1.0);
      glVertex3d(0,0,5);
      glEnd();
    }


    void fill_display_lists()
    {
      glNewList(object_list_, GL_COMPILE);
      Vertex_iterator v;
      CGAL_nef3_forall_iterators(v,vertices_) draw(v);
      glEndList();     

      glNewList(object_list_+1, GL_COMPILE);
      Edge_iterator e;
      CGAL_nef3_forall_iterators(e,edges_) draw(e);
      glEndList();     

      glNewList(object_list_+2, GL_COMPILE);
      Halffacet_iterator f;
      CGAL_nef3_forall_iterators(f,halffacets_) draw(f);
      glEndList();

      glNewList(object_list_+3, GL_COMPILE); // axes:
      construct_axes();
      glEndList();

    }

    void initialize() 
    { if (init_) return;
      init_ = true;
      axes_ = false; 
      object_list_ = glGenLists(4); CGAL_assertion(object_list_); 
      fill_display_lists();
    }


    void draw() const
    { 
      if (!is_initialized()) const_cast<Polyhedron&>(*this).initialize();
      glCallList(object_list_);   // vertices
      glCallList(object_list_+1); // edges
      if ( surface_ ) {
	glEnable(GL_LIGHTING); 
	glCallList(object_list_+2); // facets
	glDisable(GL_LIGHTING);
      }
      if (axes_) glCallList(object_list_+3); // axis
   }

    void debug(std::ostream& os = std::cerr) const
    {
      os << "OGL::Polyhedron" << std::endl;
      os << "Vertices:" << std::endl;
      Vertex_iterator v;
      CGAL_nef3_forall_iterators(v,vertices_) 
	os << "  "<<*v<<", mark="<<v->mark()<<std::endl;
      os << "Edges:" << std::endl;
      Edge_iterator e;
      CGAL_nef3_forall_iterators(e,edges_) 
	os << "  "<<*e<<", mark="<<e->mark()<<std::endl;
      os << "Facets:" << std::endl;
      Halffacet_iterator f;
      CGAL_nef3_forall_iterators(f,halffacets_) f->debug(); os << std::endl;
      os << std::endl;
    }

  }; // Polyhedron

// ----------------------------------------------------------------------------
// Viewer configuration:
// ----------------------------------------------------------------------------

  enum MenuEntries { ROTATE, SCALE, TRANSLATE, RESET_CONTROL, 
		     AXES, BOUNDARY, SKELETON, QUIT };

  int window_width  = 400;           // Breite und
  int window_height = 400;           // Hoehe des Fensters

  int mouse_x, mouse_y;              // Mauskoordinaten linker button
  int mouse_left_button = false;     // Mouse1 gedrueckt
  int motion_mode       = ROTATE;    // Bewegen der Maus bei Mouse1 gedrueckt
  int submenu1, submenu2;
  long double dx = 0;                // Translation
  long double dy = 0;                // Translation
  long double wx = 0;                // Rotation
  long double wy = 0;                // Rotation
  long double s  = 0.5;              // Skalierung
                       
  long double factor_d;              // Umrechnungsfaktor fuer Translation
  long double factor_w;              // Umrechnungsfaktor fuer Rotation
  long double factor_s;              // Umrechnungsfaktor fuer Skalierung

  // our draw object:
  static std::vector<Polyhedron>  polyhedra_;
  static std::vector<std::string> titles_;

static Polyhedron& add_polyhedron() 
{ polyhedra_.push_back(Polyhedron()); 
  return polyhedra_.back(); }

static void show (int mode)
{
  std::vector<Polyhedron>::iterator it;
  switch(mode)
  {
    case ROTATE:
    case SCALE:
    case TRANSLATE:
      motion_mode = mode;
      break;
    case RESET_CONTROL:
      dx = dy = wx = wy = 0.0;
      s = 0.5;
      motion_mode = ROTATE;
      CGAL_nef3_forall_iterators(it,polyhedra_) it->initialize();
      glutPostRedisplay();
      break;
    case AXES:
      CGAL_nef3_forall_iterators(it,polyhedra_) it->toggle_axes();
      glutPostRedisplay();
      break;
    case BOUNDARY:
      CGAL_nef3_forall_iterators(it,polyhedra_) it->boundary_on();
      CGAL_nef3_forall_iterators(it,polyhedra_) it->draw();
      glutPostRedisplay();
      break;
    case SKELETON:
      CGAL_nef3_forall_iterators(it,polyhedra_) it->skeleton_on();
      CGAL_nef3_forall_iterators(it,polyhedra_) it->draw();
      glutPostRedisplay();
      break;
    case QUIT: 
      exit(0);
  }
}



// Mausknopf gedrueckt
static void mouse (int button, int state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    mouse_x = x;
    mouse_y = y;
    mouse_left_button = true;
  } else 
    mouse_left_button = false;
}



// Objekt rotieren, zoomen oder verschieben
static void motion (int x, int y)
{
  if (mouse_left_button) {
    if (motion_mode == ROTATE) {
      wx += (y - mouse_y) * factor_w;       
      // Mausbewegung in y-Richtung entspricht Rotation um die x-Achse
      wy += (x - mouse_x) * factor_w;       
      // Mausbewegung in x-Richtung entspricht Rotation um die y-Achse
    }
    else if (motion_mode == SCALE) {
      s *= exp( (y - mouse_y) * factor_s );
    }
    else if (motion_mode == TRANSLATE) {
      dx += (x - mouse_x) * factor_d / s;
      dy -= (y - mouse_y) * factor_d / s;
    }
    mouse_x = x;
    mouse_y = y;
    glutPostRedisplay();
  }
}

static void initialize_olg()
{
  GLfloat mat_diffuse[4] = { 0.7, 0.7, 0.7, 1.0 };
  //GLfloat mat_specular[4] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_specular[4] = { 0.3, 0.3, 0.3, 1.0 };
  GLfloat mat_shininess[] = { 100.0 };
  
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse );
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );

//#define SCREENSHOTS
#ifdef SCREENSHOTS
  GLfloat mat_emission[] = { 0.1, 0.1, 0.2, 0.0 };
  glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
  //for screenshots enable this section
#endif

  GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };     // red diffuse light 
  GLfloat light_position[] = {  0.0, 0.0, 10.0, 0.0 }; //  finite location
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);

  GLfloat ambient_light[] = {  0.8, 0.8, 0.8, 1.0 };
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient_light);
  glEnable(GL_LIGHT1);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);

  glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

}

static void enter_leave(int state)
{ glutPostRedisplay(); }

static void draw()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glPushMatrix();
  glRotated(wy,0,1,0);
  glRotated(wx,1,0,0);
  glScaled(s,s,s);
  glTranslated(dx,dy,0.0);
  int win = glutGetWindow();
  polyhedra_[win-1].draw();
  glPopMatrix();
  glutSwapBuffers();
}

static void reshape(int width, int height)
{
  window_width = width;
  window_height = height;

  glViewport(0, 0, (GLint)width, (GLint)height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (width>height) {
    long double w = (long double) width / (long double) height;
    glOrtho( -2*w, 2*w, -2.0, 2.0, -100.0, 100.0 );
    factor_d = 2.0 / (height/2.0);           
    // halbe Fensterhoehe soll 2 LE entsprechen
    factor_w = 90.0 / (height/2.0);           
    // halbe Fensterhoehe soll 90 Grad entsprechen
    factor_s = std::log(4.0) / (height/2.0);           
    // halbe Fensterhoehe soll Faktor 4 entsprechen
  } else {
    long double h = (long double) height / (long double) width;
    glOrtho( -2.0, 2.0, -2*h, 2*h, -100.0, 100.0 );
    factor_d =  2.0   / (width/2.0);            
    // halbe Fensterbreite soll 2 LE entsprechen
    factor_w = 90.0   / (width/2.0);            
    // halbe Fensterbreite soll 90 Grad entsprechen
    factor_s = std::log(4.0) / (height/2.0);           
    // halbe Fensterhoehe soll Faktor 4 entsprechen
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}


static void start_viewer()
{ 
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(0,0);

  int submenu1 = glutCreateMenu(show);
  glutAddMenuEntry("Reset",RESET_CONTROL);
  glutAddMenuEntry("Rotate",ROTATE);
  glutAddMenuEntry("Scale",SCALE);
  glutAddMenuEntry("Translate",TRANSLATE);

  int submenu2 = glutCreateMenu(show);
  glutAddMenuEntry("Boundary",BOUNDARY);
  glutAddMenuEntry("Skeleton",SKELETON);
  glutAddMenuEntry("Toggle Axes",AXES);

  for (unsigned i = 0; i < polyhedra_.size(); ++i) {
    if (i > 0 ) glutInitWindowPosition(i*(window_width+12),0);
    if ( i < titles_.size() ) glutCreateWindow(titles_[i].c_str());
    else                      glutCreateWindow(" Polyhedron ");
    glutEntryFunc(enter_leave);
    initialize_olg();
    glutDisplayFunc(draw);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutCreateMenu(show);
    glutAddSubMenu("Control",submenu1);
    glutAddSubMenu("Render",submenu2);
    glutAddMenuEntry("Quit",QUIT);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }
  glutMainLoop();
}

} // namespace OGL

// ----------------------------------------------------------------------------
// SNC_visualizor_OGL<SNC_>
// when you define a SNC_visualizor_OGL<SNC_> object then you create a
// viewer window that shows the corresponding SNC_ object. All such viewers
// follow the same visualization motion. Thus visualizing binary operations
// is easy: just display the two input SNCs and the output SNC
//
// for similar viewer code see Nef_S2/demo/Nef_S2/Nef_polyhedron_S2-demo.C
// ----------------------------------------------------------------------------



template <typename SNC_>
class SNC_visualizor_OGL : public SNC_decorator<SNC_>
{ typedef SNC_ SNC_structure;
  typedef CGAL::SNC_visualizor_OGL<SNC_> Self;
  typedef CGAL::SNC_decorator<SNC_>      Base;
  typedef CGAL::SNC_FM_decorator<SNC_>   FM_decorator;

  CGAL::OGL::Polyhedron* ppoly_;

public:
  #define USING(t) typedef typename SNC_::t t
  USING(Vertex_iterator); 
  USING(Halfedge_iterator); 
  USING(Halffacet_iterator); 
  USING(Halffacet_cycle_iterator);

  USING(Object_handle);
  USING(SObject_handle); 
  USING(SHalfedge_handle); 
  USING(SHalfloop_handle); 

  USING(Vertex_handle); 
  USING(Halfedge_handle); 
  USING(Halffacet_handle);

  USING(Point_3);
  USING(Vector_3);
  USING(Segment_3);
  USING(Plane_3);
  USING(Mark);
  USING(SHalfedge_around_facet_circulator);
  USING(SHalfedge_around_facet_const_circulator);
  #undef USING

  SNC_visualizor_OGL(const SNC_structure& N) : Base(const_cast<SNC_&>(N))
  { ppoly_ = & CGAL::OGL::add_polyhedron(); }

  OGL::Double_point double_point(const Point_3& p) const
  { return OGL::Double_point(CGAL::to_double(p.x()),
			     CGAL::to_double(p.y()),
			     CGAL::to_double(p.z())); }

  OGL::Double_segment double_segment(const Segment_3& s) const
  { return OGL::Double_segment(double_point(s.source()),
			       double_point(s.target())); }

  void draw(Vertex_handle v) const
  { ppoly_->push_back(double_point(point(v)), mark(v)); }

  void draw(Halfedge_handle e) const
  { ppoly_->push_back(double_segment(segment(e)), mark(e)); }

  void draw(Halffacet_handle f) const
  { OGL::DFacet g;
    FM_decorator D(*sncp(),f);
    Halffacet_cycle_iterator fc; // all facet cycles:
    for (fc = D.facet_cycles_begin(); fc != D.facet_cycles_end(); ++ fc)
      if ( fc.is_shalfedge() ) { // non-trivial facet cycle
	g.new_facet_cycle();
	SHalfedge_handle h = fc;
	SHalfedge_around_facet_circulator hc(h), he(hc);
	CGAL_For_all(hc,he) // all vertex coordinates in facet cycle
	  g.push_back_vertex(double_point(point(source(hc))));
      }
    Vector_3 v = orthogonal_vector(f);
    g.set_normal(CGAL::to_double(v.x()), 
		 CGAL::to_double(v.y()), 
		 CGAL::to_double(v.z()), 
		 mark(f));
    ppoly_->push_back(g);
  }

  void draw() const
  { 
    Vertex_iterator v;
    CGAL_nef3_forall_vertices(v,*sncp()) draw(v);
    Halfedge_iterator e;
    CGAL_nef3_forall_edges(e,*sncp()) draw(e);
    Halffacet_iterator f;
    CGAL_nef3_forall_facets(f,*sncp()) draw(f);
  }

}; // SNC_visualizor_OGL<SNC_>





CGAL_END_NAMESPACE
#endif //CGAL_SNC_VISUALIZOR_OGL_H

