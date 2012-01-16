// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_GEOMETRY_OGL_H
#define CGAL_SPHERE_GEOMETRY_OGL_H

#include <CGAL/Nef_S2/OGL_base_object.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_triangulator.h>
#include <CGAL/IO/Color.h>
#include <qgl.h>
#include <CGAL/glu.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 151
#include <CGAL/Nef_2/debug.h>

#define CGAL_NEF_LGREY CGAL::Color(170,170,200)
#define CGAL_NEF_DGREY CGAL::Color(30,30,50)

namespace CGAL {
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
typedef VKernel::Aff_transformation_3       Affine_3;
typedef std::vector<VPoint>                 VSegment;
typedef VKernel::Triangle_3                 DTriangle;
typedef std::vector<DTriangle>              VTriangle;

template <typename R>
VVector convert(const CGAL::Vector_3<R>& v)
{ return VVector(CGAL::to_double(v.x()),
                 CGAL::to_double(v.y()),
                 CGAL::to_double(v.z())); }


const double refinement_angle = 0.1;
const double shrink_fac = 0.995;

template <typename R>
class Approximator {

 public:
  static VPoint approximate(const CGAL::Sphere_point<R>& p) {
    VVector v = convert(p-CGAL::ORIGIN);
    v = v / CGAL_NTS sqrt(v*v) ; // normalize
    return CGAL::ORIGIN+v;
  }
   
  static VSegment approximate(const CGAL::Sphere_segment<R>& s) {
    
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
    double alpha = std::acos(cosalpha);
    const int units_per_halfcircle = 50;
    int units = int(units_per_halfcircle/CGAL_PI * alpha);
    if (units == 0) ++units;
    bool seg_is_short = s.is_short();
    bool seg_is_halfcircle = s.is_halfcircle();
    if ( seg_is_halfcircle ) units = units_per_halfcircle;
    else if ( !seg_is_short ) {
      units = 2*units_per_halfcircle - (units+1);
    } CGAL_NEF_TRACEV(units); CGAL_NEF_TRACEV(cosalpha); CGAL_NEF_TRACEV(alpha);
    
    v0 = v0 / v0l;
    v1 = v1 / v1l;
    v2 = v2 / v2l;
    v3 = v3 / v3l;
    VTrafo T(v0.x(),v1.x(),v2.x(),
	     v0.y(),v1.y(),v2.y(),
	     v0.z(),v1.z(),v2.z());
    VSegment S(units+1);
    for (int i=0; i<units; ++i) 
      S[i] = VPoint(std::cos(CGAL_PI*i/double(units_per_halfcircle)),
		    std::sin(CGAL_PI*i/double(units_per_halfcircle)),
		    0.0);
    double sinalpha = 1 - cosalpha*cosalpha;
    if (sinalpha <0) sinalpha = 0; 
    else             sinalpha = CGAL_NTS sqrt(sinalpha);
    if ( seg_is_short ) 
      S[units] = VPoint(cosalpha, sinalpha, 0);
    else
      S[units] = VPoint(cosalpha, -sinalpha, 0);
    VSegment::iterator it;
    for(it = S.begin(); it != S.end(); ++it) { CGAL_NEF_TRACEN(*it<<" "<<T(*it));
    *it = T(*it); 
    } CGAL_NEF_TRACEN("");
    return S;
  }

  static VSegment approximate(const CGAL::Sphere_circle<R>& s) {
    
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
      S[i] = VPoint(std::cos(2.0*CGAL_PI*i/double(units)),
		    std::sin(2.0*CGAL_PI*i/double(units)),
		    0.0);
    }
    VSegment::iterator it;
    for(it = S.begin(); it != S.end(); ++it) *it = T(*it); 
    return S;
  }


/* the following operation refines a sphere triangle as a list of flat
   triangles in 3d. The refinement only works for triangles that are
   contained in a perfect hemisphere (no long sphere segments are
   allowed as triangle segments). We split triangles along at the
   midpoint of their longest side into two. */

  static void refine(const DTriangle& t, VTriangle& T) {
    double angle[3]; int i(0);
    angle[0] = std::acos((t[0]-CGAL::ORIGIN)*(t[1]-CGAL::ORIGIN));
    angle[1] = std::acos((t[1]-CGAL::ORIGIN)*(t[2]-CGAL::ORIGIN));
    angle[2] = std::acos((t[2]-CGAL::ORIGIN)*(t[0]-CGAL::ORIGIN));
    CGAL_NEF_TRACEN("refine "<<angle[0]<<" "<<angle[1]<<" "<<angle[2]);
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

  static VTriangle approximate(const CGAL::Sphere_triangle<R>& t) {
    // we subdivide the triangle into a list of triangles until
    // we reach a fine resolution on the surface.
    
    VTriangle T;
    DTriangle td(approximate(t.point(0)), 
		 approximate(t.point(1)),
		 approximate(t.point(2)));
    CGAL_NEF_TRACEN("approximate " << td);
    refine(td,T);
    return T;
  }

};


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
    VPoint(Approximator<R>::approximate(p)), p_(p), c_(c), w_(w) {}
  Sphere_point(const Sphere_point<R>& p) : VPoint(p), Gen_object()
  { p_ = p.p_; c_ = p.c_; w_ = p.w_; }
  Sphere_point<R>& operator=(const Sphere_point<R>& p)
  { VPoint::operator=(p); p_ = p.p_;  c_ = p.c_; w_ = p.w_;
    return *this; }

  virtual ~Sphere_point() {}

  const CGAL::Sphere_point<R>& original() const
  { return p_; }

  virtual Gen_object* clone() const 
  { return new Sphere_point<R>(*this); }
 
  virtual void draw() const { 
    glPointSize(w_);
    glColor3ub(c_.red(),c_.green(),c_.blue());
    glBegin(GL_POINTS);
    glNormal3d(x(),y(),z());
    glVertex3d(x(),y(),z());
    glEnd();
  }

  virtual void print() const 
  { std::cerr << "point " << p_; }

};





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
    : VSegment(Approximator<R>::approximate(s)), s_(s), c_(c), w_(w) {}
  Sphere_segment(const Sphere_segment<R>& s) : VSegment(s), Gen_object()
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
  { CGAL_NEF_TRACEN("draw "<<s_);
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
    : VSegment(Approximator<R>::approximate(s)), s_(s), c_(c), w_(w) {}
  Sphere_circle(const Sphere_circle<R>& s) : VSegment(s), Gen_object()
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
  { CGAL_NEF_TRACEN("draw "<<s_);
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
    CGAL::Color c = CGAL::Color(100,100,120)) 
    : VTriangle(Approximator<R>::approximate(t)), t_(t), c_(c) {}

  Sphere_triangle(const Sphere_triangle<R>& t) : VTriangle(t), Gen_object()
  { t_ = t.t_; c_ = t.c_; }

  Sphere_triangle<R>& operator=(const Sphere_triangle<R>& t)
  { VTriangle::operator=(t); t_ = t.t_; c_ = t.c_; return *this; }

  virtual ~Sphere_triangle() {}

  const CGAL::Sphere_triangle<R>& original() const
  { return t_; }

  virtual Gen_object* clone() const 
  { return new Sphere_triangle<R>(*this); }

  virtual void draw() const { 
    CGAL_NEF_TRACEN("draw "<<t_ << " in " << c_);
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

//----------------------------------------------------------------------------
// the sphere:
//----------------------------------------------------------------------------

 enum { SM_FACES, SM_SKELETON, SM_TRIANGULATION };
 enum { SM_AXES, SM_CUBE };

class Unit_sphere : public OGL_base_object {
  typedef std::list<Gen_object*> Object_list;
  typedef Object_list::const_iterator  Object_const_iterator;
  typedef Object_list::iterator  Object_iterator;
  
  GLUquadricObj*        sphere_;
  int                   style_;
  bool                  initialized_;
  Object_list           objects_,triangles_,triangle_edges_;
  GLuint                sphere_list_;
  std::vector<bool>     switches;

public:
void init() { 
  style_ = SM_FACES; 
  switches[0] = true;
  switches[1] = false;
  gluQuadricNormals(sphere_,GLenum(GLU_SMOOTH));
}

Unit_sphere() : switches(2) { 
  sphere_ = gluNewQuadric(); 
  initialized_ = false; 
  init(); 
}

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

Unit_sphere(const Unit_sphere& S) : OGL_base_object(), switches(2)
{ CGAL_NEF_TRACEN("copyconstruction");
  sphere_ = gluNewQuadric();
  initialized_ = S.initialized_;
  style_ = S.style_;
  switches[0] = S.switches[0];
  switches[1] = S.switches[1];
  copy_list(S);
}


Unit_sphere& operator=(const Unit_sphere& S)
{ CGAL_NEF_TRACEN("assignment");
  initialized_ = S.initialized_;
  style_ = S.style_;
  switches[0] = S.switches[0];
  switches[1] = S.switches[1];
  clear_list(); copy_list(S);
  return *this;
}

~Unit_sphere() {
  clear_list(); 
  gluDeleteQuadric(sphere_); 
}

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

void set_style(int style) { 
  style_ = style; 
}

void toggle(int index) { 
  switches[index] = !switches[index]; 
}
private:

void construct_axes() const
{ int i;
  register double f(1.02);
  glNewList(sphere_list_+3, GL_COMPILE);
  glLineWidth(2.0);
  // red x-axis and equator
  glColor3f(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3d(0.0,0.0,0.0);
  glVertex3d(1.2,0.0,0.0);
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(i=0;i<100;i++)
    glVertex3d(f*std::cos(2.0*CGAL_PI*i/100.0),f*std::sin(2.0*CGAL_PI*i/100.0),0.0);
  glEnd();
  // green y-axis and equator
  glColor3f(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex3d(0.0,0.0,0.0);
  glVertex3d(0.0,1.2,0.0);
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(i=0;i<100;i++)
    glVertex3d(0.0,f*std::cos(2.0*CGAL_PI*i/100.0),f*std::sin(2.0*CGAL_PI*i/100.0));
  glEnd();
  // blue z-axis and equator
  glColor3f(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3d(0.0,0.0,0.0);
  glVertex3d(0.0,0.0,1.2);
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(i=0;i<100;i++)
    glVertex3d(f*std::cos(2.0*CGAL_PI*i/100.0),0.0,f*std::sin(2.0*CGAL_PI*i/100.0));
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
  glVertex3d(-1.0,-1.0,-1.0);
  glVertex3d( 1.0,-1.0,-1.0);
  glVertex3d( 1.0, 1.0,-1.0);
  glVertex3d(-1.0, 1.0,-1.0);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3d(-1.0,-1.0, 1.0);
  glVertex3d( 1.0,-1.0, 1.0);
  glVertex3d( 1.0, 1.0, 1.0);
  glVertex3d(-1.0, 1.0, 1.0);
  glEnd();
  glBegin(GL_LINES);
  glVertex3d(-1.0,-1.0,-1.0);
  glVertex3d(-1.0,-1.0, 1.0);
  glVertex3d( 1.0,-1.0,-1.0);
  glVertex3d( 1.0,-1.0, 1.0);
  glVertex3d( 1.0, 1.0,-1.0);
  glVertex3d( 1.0, 1.0, 1.0);
  glVertex3d(-1.0, 1.0,-1.0);
  glVertex3d(-1.0, 1.0, 1.0);
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

//void draw(GLdouble z_vec[3]) const
void draw() const
{ 
  gluQuadricDrawStyle(sphere_,GLenum(GLU_FILL));
  glEnable(GL_LIGHTING);
  if ( style_ == SM_SKELETON ) {
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(0.7f,0.7f,0.7f);
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
  if ( switches[0] ) glCallList(sphere_list_+3);
  if ( switches[1] ) glCallList(sphere_list_+4);
}


};

template <typename Map_>
class SM_BooleColor 
{
  typedef typename Map_::SVertex_const_handle   SVertex_const_handle;   
  typedef typename Map_::SHalfedge_const_handle SHalfedge_const_handle;   
  typedef typename Map_::SHalfloop_const_handle SHalfloop_const_handle;   
  typedef typename Map_::SFace_const_handle     SFace_const_handle;   
  typedef typename Map_::Mark Mark;   
public:
  Color color(SVertex_const_handle, Mark m) const
  { return ( m ? CGAL::BLACK : CGAL::WHITE ); }
  Color color(SHalfedge_const_handle, Mark m) const
  { return ( m ? CGAL::BLACK : CGAL::WHITE ); }
  Color color(SHalfloop_const_handle, Mark m) const
  { return ( m ? CGAL::BLACK : CGAL::WHITE ); }
  Color color(SFace_const_handle, Mark m) const
  { return ( m ? CGAL_NEF_DGREY : CGAL_NEF_LGREY ); }
};

template <typename Nef_polyhedron>
class NefS2_to_UnitSphere
{
  typedef typename Nef_polyhedron::Sphere_map          Sphere_map;
  typedef typename Nef_polyhedron::Const_decorator     SM_const_decorator;
  typedef CGAL::OGL::SM_BooleColor<Sphere_map>           Color_;
  typedef typename Sphere_map::Sphere_kernel        Sphere_kernel;
  typedef SM_decorator<Sphere_map>                  Decorator;
  typedef CGAL::SM_triangulator<Decorator>                Base;

  typedef typename Sphere_map::SVertex_const_handle SVertex_const_handle;   
  typedef typename Sphere_map::SHalfedge_const_handle SHalfedge_const_handle; 
  typedef typename Sphere_map::SFace_const_handle SFace_const_handle;     
  typedef typename Sphere_map::SVertex_const_iterator SVertex_const_iterator;
  typedef typename Sphere_map::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Sphere_map::SFace_const_iterator SFace_const_iterator;
  typedef typename Sphere_map::Mark Mark;

  typedef typename Sphere_kernel::Sphere_point    Sphere_point;
  typedef typename Sphere_kernel::Sphere_segment  Sphere_segment;
  typedef typename Sphere_kernel::Sphere_circle   Sphere_circle;
  typedef typename Sphere_kernel::Sphere_triangle Sphere_triangle;

  typedef Color_                                  Color_objects;  

 public:
  static CGAL::OGL::Unit_sphere convert(const SM_const_decorator& E_,
					const Color_objects& CO_ = Color_objects()) {
    Sphere_map MT_(true);
    Base T_(&MT_, &E_);
    CGAL::OGL::Unit_sphere S_;

    T_.triangulate();

    SHalfedge_const_iterator e;
    CGAL_forall_sedges(e,E_) {
      if ( e->source() == e->twin()->source() ) {
	S_.push_back(e->circle(), CO_.color(e,e->mark()));} else {
	S_.push_back(Sphere_segment(e->source()->point(),
				    e->twin()->source()->point(),
				    e->circle()),CO_.color(e,e->mark()));
      }
    }
    // draw sphere circles underlying loops of E_:
    
    if ( E_.has_shalfloop() )
      S_.push_back(
		   E_.shalfloop()->circle(),
		   CO_.color(E_.shalfloop(),E_.shalfloop()->mark()));
    
    // draw points underlying vertices of E_:
    SVertex_const_iterator v;
    CGAL_forall_svertices(v,E_)
      S_.push_back(v->point(),CO_.color(v,v->mark()));
    
    Unique_hash_map<SHalfedge_const_iterator,bool> Done(false);
    CGAL_forall_shalfedges(e,T_) {
      if ( Done[e] ) continue;
      SHalfedge_const_handle en(e->snext()),enn(en->snext());
      CGAL_NEF_TRACEV(T_.incident_triangle(e));
      CGAL_NEF_TRACEN(T_.incident_mark(e)<<T_.incident_mark(en)<<T_.incident_mark(enn));
      CGAL_assertion(T_.incident_mark(e)==T_.incident_mark(en) &&
		     T_.incident_mark(en)==T_.incident_mark(enn));
      Mark m = T_.incident_mark(e);
      Sphere_triangle t = T_.incident_triangle(e);
      S_.push_back(t, (m ? CGAL_NEF_DGREY : CGAL_NEF_LGREY) );
      Done[e]=Done[en]=Done[enn]=true;
    }
    
    Done.clear(false);
    CGAL_forall_shalfedges(e,T_) {
      if ( Done[e] ) continue;
      S_.push_back_triangle_edge(Sphere_segment(e->source()->point(),
						e->twin()->source()->point(),
						e->circle()));
      Done[e]=Done[e->twin()]=true;
    }    
    return S_;
  }
};

} // OGL
} //namespace CGAL



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
