#ifndef CGAL_SPHERE_GEOMETRY_OGL_H
#define CGAL_SPHERE_GEOMETRY_OGL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/IO/Color.h>
#include <GL/glut.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>

#undef _DEBUG
#define _DEBUG 151
#include <CGAL/Nef_S2/debug.h>

CGAL_BEGIN_NAMESPACE
namespace OGL {

struct Gen_object {
  Gen_object() {}
  virtual ~Gen_object() {}
  virtual void draw() const {}
  virtual Gen_object* clone() const { return 0; }
  virtual void print() const {}
};

typedef CGAL::Simple_cartesian<double>      VKernel;
typedef VKernel::Vector_3                   VVector;
typedef VKernel::Point_3                    VPoint;
typedef VKernel::Aff_transformation_3       VTrafo;
typedef std::vector<VPoint>                 VSegment;
typedef VKernel::Triangle_3                 DTriangle;
typedef std::vector<DTriangle>              VTriangle;

const double refinement_angle = 0.1;
const double shrink_fac = 0.995;

template <typename R>
VVector convert(const CGAL::Vector_3<R>& v)
{ return VVector(CGAL::to_double(v.x()),
                 CGAL::to_double(v.y()),
                 CGAL::to_double(v.z())); }

template <typename R>
VPoint approximate(const CGAL::Sphere_point<R>& p)
{
  VVector v = convert(p-CGAL::ORIGIN);
  v = v / CGAL_NTS sqrt(v*v) ; // normalize
  return CGAL::ORIGIN+v;
}

template <class R_>
class Sphere_point : public VPoint, public Gen_object {
  typedef R_ R;
  CGAL::Sphere_point<R> p_;
  CGAL::Color           c_;
  unsigned              w_;
public:
  Sphere_point() {}
  Sphere_point(const CGAL::Sphere_point<R>& p,
    CGAL::Color c = CGAL::BLACK, unsigned w = 10) : 
    VPoint(approximate(p)), p_(p), c_(c), w_(w) {}
  Sphere_point(const Sphere_point<R>& p) : VPoint(p)
  { p_ = p.p_; c_ = p.c_; w_ = p.w_; }
  Sphere_point<R>& operator=(const Sphere_point<R>& p)
  { VPoint::operator=(p); p_ = p.p_;  c_ = p.c_; w_ = p.w_;
    return *this; }

  virtual ~Sphere_point() {}

  const CGAL::Sphere_point<R>& original() const
  { return p_; }

  virtual Gen_object* clone() const 
  { return new Sphere_point<R>(*this); }
 
  virtual void draw() const
  { glPointSize(w_);
    glColor3ub(c_.red(),c_.green(),c_.blue());
    glBegin(GL_POINTS);
    glNormal3d(x(),y(),z());
    glVertex3d(x(),y(),z());
    glEnd();
  }

  virtual void print() const 
  { std::cerr << "point " << p_; }

};



template <typename R>
VSegment approximate(const CGAL::Sphere_segment<R>& s)
{
  /* we construct the rotation matrix that transfers the x-axis
     into |ps_|, the z-axis into h_.orthogonal_vector() and the
     y-axix into the corresponding crossproduct of the two.*/
  if ( s.is_degenerate() ) {
    VSegment S(1);
    S[0] = approximate(s.source());
    return S;
  }
  
  VVector v0 = convert(s.source()-CGAL::ORIGIN);
  VVector v2 = convert(s.sphere_circle().orthogonal_vector());
  VVector v1(-cross_product(v0,v2));
  VVector v3 = convert(s.target()-CGAL::ORIGIN);
  double v0l = CGAL_NTS sqrt(v0*v0);
  double v1l = CGAL_NTS sqrt(v1*v1);
  double v2l = CGAL_NTS sqrt(v2*v2);
  double v3l = CGAL_NTS sqrt(v3*v3);
  double cosalpha = v0*v3 / v0l / v3l; 
  double alpha = acos(cosalpha);
  const int units_per_halfcircle = 50;
  int units = int(units_per_halfcircle/M_PI * alpha);
  if (units == 0) ++units;
  bool seg_is_short = s.is_short();
  bool seg_is_halfcircle = s.is_halfcircle();
  if ( seg_is_halfcircle ) units = units_per_halfcircle;
  else if ( !seg_is_short ) {
    units = 2*units_per_halfcircle - (units+1);
  } TRACEV(units); TRACEV(cosalpha); TRACEV(alpha);

  v0 = v0 / v0l;
  v1 = v1 / v1l;
  v2 = v2 / v2l;
  v3 = v3 / v3l;
  VTrafo T(v0.x(),v1.x(),v2.x(),
           v0.y(),v1.y(),v2.y(),
           v0.z(),v1.z(),v2.z());
  VSegment S(units+1);
  for (int i=0; i<units; ++i) 
    S[i] = VPoint(cos(M_PI*i/double(units_per_halfcircle)),
                  sin(M_PI*i/double(units_per_halfcircle)),
                  0.0);
  double sinalpha = 1 - cosalpha*cosalpha;
  if (sinalpha <0) sinalpha = 0; 
  else             sinalpha = CGAL_NTS sqrt(sinalpha);
  if ( seg_is_short ) 
    S[units] = VPoint(cosalpha, sinalpha, 0);
  else
    S[units] = VPoint(cosalpha, -sinalpha, 0);
  VSegment::iterator it;
  for(it = S.begin(); it != S.end(); ++it) { TRACEN(*it<<" "<<T(*it));
    *it = T(*it); 
  } TRACEN("");
  return S;
}


template <class R_>
class Sphere_segment : public VSegment, public Gen_object { 
  typedef R_ R;
  CGAL::Sphere_segment<R> s_;
  CGAL::Color             c_;
  unsigned                w_;
public:
  Sphere_segment() {}
  Sphere_segment(const CGAL::Sphere_segment<R>& s,
    CGAL::Color c = CGAL::BLACK, unsigned w = 2) 
    : VSegment(approximate(s)), s_(s), c_(c), w_(w) {}
  Sphere_segment(const Sphere_segment<R>& s) : VSegment(s)
  { s_ = s.s_; c_ = s.c_; w_ = s.w_; }
  Sphere_segment<R>& operator=(const Sphere_segment<R>& s)
  { VSegment::operator=(s); s_ = s.s_; c_ = s.c_; w_ = s.w_;
    return *this; }
  virtual ~Sphere_segment() {}

  const CGAL::Sphere_segment<R>& original() const
  { return s_; }

  virtual Gen_object* clone() const 
  { return new Sphere_segment<R>(*this); }

  virtual void draw() const
  { TRACEN("draw "<<s_);
    if ( size() == 1 ) {
      glPointSize(5*w_);
      glColor3ub(c_.red(),c_.green(),c_.blue());
      glBegin(GL_POINTS);
      glNormal3d(begin()->x(),begin()->y(),begin()->z());
      glVertex3d(begin()->x(),begin()->y(),begin()->z());
      glEnd();
    } else {
      glLineWidth(w_);
      glColor3ub(c_.red(),c_.green(),c_.blue()); 
      glBegin(GL_LINE_STRIP);
      VSegment::const_iterator it;
      for(it = begin(); it != end(); ++it) {
        glNormal3d(it->x(),it->y(),it->z());
        glVertex3d(it->x(),it->y(),it->z());
      }
      glEnd();
    }
  }

  virtual void print() const 
  { std::cerr << "segment " << s_; }

};

template <typename R>
VSegment approximate(const CGAL::Sphere_circle<R>& s)
{
  /* we construct the rotation matrix that transfers the x-axis
     into |ps_|, the z-axis into h_.orthogonal_vector() and the
     y-axix into the corresponding crossproduct of the two.*/
  
  VVector v0 = convert(s.base1());
  VVector v1 = convert(s.base2());
  VVector v2 = convert(s.orthogonal_vector());
  double v0l = CGAL_NTS sqrt(v0*v0);
  double v1l = CGAL_NTS sqrt(v1*v1);
  double v2l = CGAL_NTS sqrt(v2*v2);
  const int units = 100;
  v0 = v0 / v0l;
  v1 = v1 / v1l;
  v2 = v2 / v2l;
  VTrafo T(v0.x(),v1.x(),v2.x(),
           v0.y(),v1.y(),v2.y(),
           v0.z(),v1.z(),v2.z());
  VSegment S(units);
  for (int i=0; i<units; ++i) {
    S[i] = VPoint(cos(2.0*M_PI*i/double(units)),
                  sin(2.0*M_PI*i/double(units)),
                  0.0);
  }
  VSegment::iterator it;
  for(it = S.begin(); it != S.end(); ++it) *it = T(*it); 
  return S;
}



template <class R_>
class Sphere_circle : public VSegment, public Gen_object { 
  typedef R_ R;
  CGAL::Sphere_circle<R> s_;
  CGAL::Color            c_;
  unsigned               w_;
public:
  Sphere_circle() {}
  Sphere_circle(const CGAL::Sphere_circle<R>& s,
    CGAL::Color c = CGAL::BLACK, unsigned w = 2) 
    : VSegment(approximate(s)), s_(s), c_(c), w_(w) {}
  Sphere_circle(const Sphere_circle<R>& s) : VSegment(s)
  { s_ = s.s_; c_ = s.c_; w_ = s.w_; }
  Sphere_circle<R>& operator=(const Sphere_circle<R>& s)
  { VSegment::operator=(s); s_ = s.s_; c_ = s.c_; w_ = s.w_;
    return *this; }
  virtual ~Sphere_circle() {}

  const CGAL::Sphere_circle<R>& original() const
  { return s_; }

  virtual Gen_object* clone() const 
  { return new Sphere_circle<R>(*this); }

  virtual void draw() const
  { TRACEN("draw "<<s_);
    glLineWidth(w_);
    glColor3ub(c_.red(),c_.green(),c_.blue()); 
    glBegin(GL_LINE_LOOP);
    VSegment::const_iterator it;
    for(it = begin(); it != end(); ++it) {
      glNormal3d(it->x(),it->y(),it->z());
      glVertex3d(it->x(),it->y(),it->z());
    }
    glEnd();
  }

  virtual void print() const 
  { std::cerr << "circle " << s_; }

};


/* The following class approximates a spherical triangle by a list
   of flat triangles */

template <class R_>
class Sphere_triangle : public VTriangle, public Gen_object { 
  typedef R_ R;
  CGAL::Sphere_triangle<R> t_;
  CGAL::Color              c_;
public:
  Sphere_triangle() {}

  Sphere_triangle(const CGAL::Sphere_triangle<R>& t,
    CGAL::Color c = CGAL::GREY) 
    : VTriangle(approximate(t)), t_(t), c_(c) {}

  Sphere_triangle(const Sphere_triangle<R>& t) : VTriangle(t)
  { t_ = t.t_; c_ = t.c_; }

  Sphere_triangle<R>& operator=(const Sphere_triangle<R>& t)
  { VTriangle::operator=(t); t_ = t.t_; c_ = s.c_; return *this; }

  virtual ~Sphere_triangle() {}

  const CGAL::Sphere_triangle<R>& original() const
  { return t_; }

  virtual Gen_object* clone() const 
  { return new Sphere_triangle<R>(*this); }

  virtual void draw() const
  { TRACEN("draw "<<t_);
    VTriangle::const_iterator it;
    VPoint p;
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glColor3ub(c_.red(),c_.green(),c_.blue());
    glBegin(GL_TRIANGLES);
    for(it = begin(); it != end(); ++it) {
        p = it->vertex(0); 
	glNormal3d(p.x(),p.y(),p.z()); 	//glVertex3d(p.x(),p.y(),p.z());
	glVertex3d(shrink_fac*p.x(),shrink_fac*p.y(),shrink_fac*p.z());
        p = it->vertex(1); 
	glNormal3d(p.x(),p.y(),p.z()); 
	glVertex3d(shrink_fac*p.x(),shrink_fac*p.y(),shrink_fac*p.z());
        p = it->vertex(2); 
	glNormal3d(p.x(),p.y(),p.z()); 
	glVertex3d(shrink_fac*p.x(),shrink_fac*p.y(),shrink_fac*p.z());
    }
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
  }

  virtual void print() const 
  { std::cerr << "triangle " << t_; }

};

/* the following operation refines a sphere triangle as a list of flat
   triangles in 3d. The refinement only works for triangles that are
   contained in a perfect hemisphere (no long sphere segments are
   allowed as triangle segments). We split triangles along at the
   midpoint of their longest side into two. */

void refine(const DTriangle& t, VTriangle& T)
{
  double angle[3]; int i(0);
  angle[0] = acos((t[0]-CGAL::ORIGIN)*(t[1]-CGAL::ORIGIN));
  angle[1] = acos((t[1]-CGAL::ORIGIN)*(t[2]-CGAL::ORIGIN));
  angle[2] = acos((t[2]-CGAL::ORIGIN)*(t[0]-CGAL::ORIGIN));
  TRACEN("refine "<<angle[0]<<" "<<angle[1]<<" "<<angle[2]);
  if ( angle[1] > angle[0] ) {
    if ( angle[2] > angle[1] ) i=2;
    else                       i=1;
  } else { // angle[0] >= angle[1]
    if ( angle[2] > angle[0] ) i=2;
    else                       i=0;
  }
  // now i references the side of maximal angle
  if ( angle[i] < refinement_angle ) // refinement threshhold
  { T.push_back(t); return; }
  VVector v;
  switch (i) {
    case 0: v = (t[0]-CGAL::ORIGIN)+(t[1]-CGAL::ORIGIN); break;
    case 1: v = (t[1]-CGAL::ORIGIN)+(t[2]-CGAL::ORIGIN); break;
    case 2: v = (t[2]-CGAL::ORIGIN)+(t[0]-CGAL::ORIGIN); break;
  }
  v = v / CGAL_NTS sqrt(v*v) ; // normalize
  VPoint p = CGAL::ORIGIN+v;
  DTriangle t1,t2;
  switch (i) {
    case 0: t1=DTriangle(t[0],p,t[2]); t2=DTriangle(p,t[1],t[2]); break;
    case 1: t1=DTriangle(t[1],p,t[0]); t2=DTriangle(p,t[2],t[0]); break;
    case 2: t1=DTriangle(t[2],p,t[1]); t2=DTriangle(p,t[0],t[1]); break;
  }
  refine(t1,T);
  refine(t2,T);
}

template <typename R>
VTriangle approximate(const CGAL::Sphere_triangle<R>& t)
{
  /* we subdivide the triangle into a list of triangles until
     we reach a fine resolution on the surface.*/
  
  VTriangle T;
  DTriangle td(approximate(t.point(0)), 
	       approximate(t.point(1)),
	       approximate(t.point(2)));
  TRACEN("approximate " << td);
  refine(td,T);
  return T;
}


//----------------------------------------------------------------------------
// the sphere:
//----------------------------------------------------------------------------

enum Sphere_map_mode {
  SM_FACES, SM_SKELETON, SM_TRIANGULATION
};

class Unit_sphere  {
typedef std::list<Gen_object*> Object_list;
typedef Object_list::const_iterator  Object_const_iterator;
typedef Object_list::iterator  Object_iterator;

GLUquadricObj*        sphere_;
Sphere_map_mode       style_;
bool                  axes_, cube_, initialized_;
Object_list           objects_,triangles_,triangle_edges_;
GLuint                sphere_list_;

public:
void init(Sphere_map_mode style = SM_FACES)
{ style_ = style; axes_ = true; cube_ = false;
  gluQuadricNormals(sphere_,GLenum(GLU_SMOOTH));
}

Unit_sphere(Sphere_map_mode style = SM_FACES) 
{ sphere_ = gluNewQuadric(); initialized_ = false; init(style); }

void clear_list() 
{ while ( objects_.begin() != objects_.end() ) {
    delete (*objects_.begin());
    objects_.pop_front();
  }
  while ( triangles_.begin() != triangles_.end() ) {
    delete (*triangles_.begin());
    triangles_.pop_front();
  }
  while ( triangle_edges_.begin() != triangle_edges_.end() ) {
    delete (*triangle_edges_.begin());
    triangle_edges_.pop_front();
  }
}

void copy_list(const Unit_sphere& S)
{ Object_const_iterator it;
  CGAL_forall_iterators (it,S.objects_)
    objects_.push_back( (*it)->clone() );
  CGAL_forall_iterators (it,S.triangles_)
    triangles_.push_back( (*it)->clone() );
  CGAL_forall_iterators (it,S.triangle_edges_)
    triangle_edges_.push_back( (*it)->clone() );
}

void print() const
{ std::cerr << "Dumping Unit_sphere:\n";
  for (Object_const_iterator it = objects_.begin();
       it != objects_.end(); ++it)
     (*it)->print();
  std::cerr << std::endl;
  for (Object_const_iterator it = triangles_.begin();
       it != triangles_.end(); ++it)
     (*it)->print();
  std::cerr << std::endl;
}

Unit_sphere(const Unit_sphere& S) 
{ TRACEN("copyconstruction");
  sphere_ = gluNewQuadric();
  initialized_ = S.initialized_;
  style_ = S.style_;
  axes_ = S.axes_;
  cube_ = S.cube_;
  copy_list(S);
}


Unit_sphere& operator=(const Unit_sphere& S)
{ TRACEN("assignment");
  initialized_ = S.initialized_;
  style_ = S.style_;
  axes_ = S.axes_;
  cube_ = S.cube_;
  clear_list(); copy_list(S);
  return *this;
}

~Unit_sphere() { clear_list(); gluDeleteQuadric(sphere_); }

template <typename R>
void push_back(const CGAL::Sphere_point<R>& p,
  CGAL::Color c = CGAL::YELLOW, unsigned w = 5)
{ objects_.push_back(new Sphere_point<R>(p,c,w)); }

template <typename R>
void push_back(const CGAL::Sphere_segment<R>& s,
  CGAL::Color c = CGAL::BLACK, unsigned w = 1)
{ objects_.push_back(new Sphere_segment<R>(s,c,w)); }

template <typename R>
void push_back(const CGAL::Sphere_circle<R>& s,
  CGAL::Color c = CGAL::BLACK, unsigned w = 1)
{ objects_.push_back(new Sphere_circle<R>(s,c,w)); }

template <typename R>
void push_back(const CGAL::Sphere_triangle<R>& t,
  CGAL::Color c = CGAL::WHITE)
{ triangles_.push_back(new Sphere_triangle<R>(t,c)); }

template <typename R>
void push_back_triangle_edge(const CGAL::Sphere_segment<R>& s,
  CGAL::Color c = CGAL::BLUE, unsigned w = 1)
{ triangle_edges_.push_back(new Sphere_segment<R>(s,c,w)); }

void set_style(Sphere_map_mode style) { style_ = style; }
void toggle_axes() { axes_ = !axes_; }
void toggle_cube() { cube_ = !cube_; }

private:

void construct_axes() const
{ int i;
  register double f(1.02);
  glNewList(sphere_list_+3, GL_COMPILE);
  glLineWidth(2.0);
  // red x-axis and equator
  glColor3f(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(1.2,0.0,0.0);
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(i=0;i<100;i++)
    glVertex3d(f*cos(2.0*M_PI*i/100.0),f*sin(2.0*M_PI*i/100.0),0.0);
  glEnd();
  // green y-axis and equator
  glColor3f(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,1.2,0.0);
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(i=0;i<100;i++)
    glVertex3d(0.0,f*cos(2.0*M_PI*i/100.0),f*sin(2.0*M_PI*i/100.0));
  glEnd();
  // blue z-axis and equator
  glColor3f(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,0.0,1.2);
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(i=0;i<100;i++)
    glVertex3d(f*cos(2.0*M_PI*i/100.0),0.0,f*sin(2.0*M_PI*i/100.0));
  glEnd();
  // six coordinate points in pink:
  glPointSize(10);
  glColor3f(1,0,1);
  glBegin(GL_POINTS);
  glVertex3d(1,0,0);
  glVertex3d(0,1,0);
  glVertex3d(0,0,1);
  glVertex3d(-1,0,0);
  glVertex3d(0,-1,0);
  glVertex3d(0,0,-1);
  glEnd(); glEndList();
}

void construct_cube() const
{ 
  glNewList(sphere_list_+4, GL_COMPILE);
  glColor3f(1,1,0); // yellow
  glLineWidth(2.0);
  glBegin(GL_LINE_LOOP);
  glVertex3f(-1.0,-1.0,-1.0);
  glVertex3f( 1.0,-1.0,-1.0);
  glVertex3f( 1.0, 1.0,-1.0);
  glVertex3f(-1.0, 1.0,-1.0);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f(-1.0,-1.0, 1.0);
  glVertex3f( 1.0,-1.0, 1.0);
  glVertex3f( 1.0, 1.0, 1.0);
  glVertex3f(-1.0, 1.0, 1.0);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(-1.0,-1.0,-1.0);
  glVertex3f(-1.0,-1.0, 1.0);
  glVertex3f( 1.0,-1.0,-1.0);
  glVertex3f( 1.0,-1.0, 1.0);
  glVertex3f( 1.0, 1.0,-1.0);
  glVertex3f( 1.0, 1.0, 1.0);
  glVertex3f(-1.0, 1.0,-1.0);
  glVertex3f(-1.0, 1.0, 1.0);
  glEnd(); glEndList();
}

void construct()
{ initialized_=true;
  sphere_list_ = glGenLists(5);
  CGAL_assertion_msg(sphere_list_!=0,"no display list.");
  // skeleton:
  glNewList(sphere_list_, GL_COMPILE);
  for (Object_const_iterator it = objects_.begin();
       it != objects_.end(); ++it) 
     (*it)->draw();
  glEndList();
  // triangles:
  glNewList(sphere_list_+1, GL_COMPILE);
  for (Object_const_iterator it = triangles_.begin();
       it != triangles_.end(); ++it) 
     (*it)->draw();
  glEndList();
  glNewList(sphere_list_+2, GL_COMPILE);
  for (Object_const_iterator it = triangle_edges_.begin();
       it != triangle_edges_.end(); ++it) 
     (*it)->draw();
  glEndList();
  // orientation features:
  construct_axes();
  construct_cube();
}

public:

void draw() const
{ //gluQuadricDrawStyle(sphere_,style_);
  gluQuadricDrawStyle(sphere_,GLenum(GLU_FILL));
  glEnable(GL_LIGHTING);
  if ( style_ == SM_SKELETON ) {
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(0.7,0.7,0.7);
    gluSphere(sphere_,shrink_fac,50,50);
    glDisable(GL_COLOR_MATERIAL);
  }
  glDisable(GL_LIGHTING);
  if ( !initialized_ ) const_cast<Unit_sphere*>(this)->construct();
  if ( style_ == SM_FACES || style_ == SM_TRIANGULATION ) {
    glEnable(GL_LIGHTING); 
    glCallList(sphere_list_+1); // triangle list
    glDisable(GL_LIGHTING);
  }
  if ( style_ == SM_TRIANGULATION ) {
    glCallList(sphere_list_+2); // triangle edge list
  }
  
  if ( style_ != SM_TRIANGULATION ) {
    glCallList(sphere_list_);
  }
  if ( axes_ ) glCallList(sphere_list_+3);
  if ( cube_ ) glCallList(sphere_list_+4);
}


};

  enum MenuEntries { ROTATE, SCALE, TRANSLATE, RESET_CONTROL, 
		     UNITY_CUBE, AXES, 
  		     FACES, SKELETON, TRIANGULATION, 
		     QUIT };

  int window_width  = 300;           // Breite und
  int window_height = 300;           // Hoehe des Fensters

  int mouse_x, mouse_y;              // Mauskoordinaten linker button
  int mouse_left_button = false;     // Mouse1 gedrueckt
  int motion_mode       = ROTATE;    // Bewegen der Maus bei Mouse1 gedrueckt
  int submenu1, submenu2;
  long double dx = 0;                // Translation
  long double dy = 0;                // Translation
  long double wx = 0;                // Rotation
  long double wy = 0;                // Rotation
  long double s  = 1.5;              // Skalierung
                       
  long double factor_d;              // Umrechnungsfaktor fuer Translation
  long double factor_w;              // Umrechnungsfaktor fuer Rotation
  long double factor_s;              // Umrechnungsfaktor fuer Skalierung

  // our draw object:
  static std::vector<Unit_sphere> spheres_;
  static std::vector<std::string> titles_;

void show (int mode)
{
  std::vector<Unit_sphere>::iterator it;
  switch(mode)
  {
    case ROTATE:
    case SCALE:
    case TRANSLATE:
      motion_mode = mode;
      break;
    case RESET_CONTROL:
      dx = dy = wx = wy = 0.0;
      s = 1.5;
      motion_mode = ROTATE;
      CGAL_forall_iterators(it,spheres_) it->init();
      glutPostRedisplay();
      break;
    case UNITY_CUBE:
      CGAL_forall_iterators(it,spheres_) it->toggle_cube();
      glutPostRedisplay();
      break;
    case AXES:
      CGAL_forall_iterators(it,spheres_) it->toggle_axes();
      glutPostRedisplay();
      break;
    case FACES:
      CGAL_forall_iterators(it,spheres_) it->set_style(SM_FACES);
      glutPostRedisplay();
      break;
    case SKELETON:
      CGAL_forall_iterators(it,spheres_) it->set_style(SM_SKELETON);
      glutPostRedisplay();
      break;
    case TRIANGULATION:
      CGAL_forall_iterators(it,spheres_) it->set_style(SM_TRIANGULATION);
      glutPostRedisplay();
      break;
    case QUIT: 
      exit(0);
  }
}



// Mausknopf gedrueckt
void mouse (int button, int state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
      mouse_x = x;
      mouse_y = y;
      mouse_left_button = true;
    }
  else
    mouse_left_button = false;
}



// Objekt rotieren, zoomen oder verschieben
void motion (int x, int y)
{
  if (mouse_left_button)
    {
      if (motion_mode == ROTATE)
        {
          wx += (y - mouse_y) * factor_w;       
          // Mausbewegung in y-Richtung entspricht Rotation um die x-Achse
          wy += (x - mouse_x) * factor_w;       
          // Mausbewegung in x-Richtung entspricht Rotation um die y-Achse
        }
      else if (motion_mode == SCALE)
        {
          s *= exp( (y - mouse_y) * factor_s );
        }
      else if (motion_mode == TRANSLATE)
        {
          dx += (x - mouse_x) * factor_d / s;
          dy -= (y - mouse_y) * factor_d / s;
        }
      mouse_x = x;
      mouse_y = y;
      glutPostRedisplay();
    }
}

void init()
{
  GLfloat mat_diffuse[4] = { 0.7, 0.7, 0.7, 1.0 };
  GLfloat mat_specular[4] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[] = { 100.0 };
  GLfloat ambient_light[] = { 0.2, 0.2, 0.2, 1.0 };
  
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse );
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
#define SCREENSHOTS
#ifdef SCREENSHOTS
  GLfloat mat_emission[] = { 0.1, 0.1, 0.2, 0.0 };
  glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
  //for screenshots enable this section
#endif

  GLfloat light0[4] = { 4.0, 4.0, 10.0, 1.0 };
  glLightfv (GL_LIGHT0, GL_POSITION, light0);
  glEnable (GL_LIGHT0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_POINT_SMOOTH);
}

void enter_leave(int state)
{ glutPostRedisplay(); }

void draw()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glPushMatrix();
  glRotated(wy,0,1,0);
  glRotated(wx,1,0,0);
  glScaled(s,s,s);
  glTranslated(dx,dy,0.0);
  int win = glutGetWindow();
  //std::cerr << "WINDOW" << win << std::endl;
  spheres_[win-1].draw();
  glPopMatrix();
  glutSwapBuffers();
}

void reshape(int width, int height)
{
  window_width = width;
  window_height = height;

  glViewport(0, 0, (GLint)width, (GLint)height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (width>height)
    {
      long double w = (long double) width / (long double) height;
      glOrtho( -2*w, 2*w, -2.0, 2.0, -100.0, 100.0 );
      factor_d =  2.0   / (height/2.0);           
      // halbe Fensterhoehe soll 2 LE entsprechen
      factor_w = 90.0   / (height/2.0);           
      // halbe Fensterhoehe soll 90 Grad entsprechen
      factor_s = std::log(4.0) / (height/2.0);           
      // halbe Fensterhoehe soll Faktor 4 entsprechen
    }  
  else
    {
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

void add_sphere()
{ spheres_.push_back(Unit_sphere()); }

void start_viewer()
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
  glutAddMenuEntry("Toggle Axes",AXES);
  glutAddMenuEntry("Toggle Unity Cube",UNITY_CUBE);

  int submenu3 = glutCreateMenu(show);
  glutAddMenuEntry("Faces",FACES);
  glutAddMenuEntry("Skeleton",SKELETON);
  glutAddMenuEntry("Triangulation",TRIANGULATION);

  //for (unsigned i = 0; i < spheres_.size(); ++i) spheres_[i].print();
  for (unsigned i = 0; i < spheres_.size(); ++i) {
    if (i > 0 ) glutInitWindowPosition(i*(window_width+12),0);
    if ( i < titles_.size() ) glutCreateWindow(titles_[i].c_str());
    else                      glutCreateWindow("Sphere Window");
    glutEntryFunc(enter_leave);
    init();
    glutDisplayFunc(draw);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutCreateMenu(show);
    glutAddSubMenu("Control",submenu1);
    glutAddSubMenu("Render",submenu3);
    glutAddSubMenu("Options",submenu2);
    glutAddMenuEntry("Quit",QUIT);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }

  glutMainLoop();
}


} // OGL
CGAL_END_NAMESPACE



template <class R>
std::ostream& operator<<(std::ostream& os, 
                         const CGAL::OGL::Sphere_segment<R>& s)
{ CGAL::OGL::VSegment::const_iterator it;
   os << s.original() << " ";
  for (it = s.begin(); it != s.end(); ++it) 
    os << *it;
  return os;
}

template <class R>
std::ostream& operator<<(std::ostream& os, 
                         const CGAL::OGL::Sphere_circle<R>& s)
{ CGAL::OGL::VSegment::const_iterator it;
   os << s.original() << " ";
  for (it = s.begin(); it != s.end(); ++it) 
    os << *it;
  return os;
}



template <class R>
std::ostream& operator<<(std::ostream& os, 
                         const CGAL::OGL::Sphere_point<R>& p)
{ os << p.original() << CGAL::OGL::VPoint(p); return os; }


#endif //CGAL_SPHERE_GEOMETRY_OGL_H

