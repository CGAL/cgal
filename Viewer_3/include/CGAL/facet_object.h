// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
#ifndef V_UTILS
#include <CGAL/v_utils.h>
#endif 

// enum Style {FILL=1, WIRE, RAW, FOS1, FOS2, FOS3, FOS4, FOS5};
CGAL_BEGIN_NAMESPACE


typedef  std::vector<double> vertex ;
typedef  std::vector<vertex> edge ;
typedef  std::list<vertex> facet;

typedef  std::list<vertex> list_vertices;
typedef  std::list<edge>   list_edges;
typedef  std::list<facet>   list_facets;


void draw_vertex(vertex v, Size s, Style sty, Precision prec)
{
  double x = v[0] ;
  double y = v[1] ;
  double z = v[2] ;
  	if (sty==FILL) {
	  GLUquadricObj *q= gluNewQuadric();
	  glPushMatrix();
	  glTranslatef(x, y, z);
	  gluQuadricNormals(q, (GLenum) GL_SMOOTH);
	  gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
	  gluSphere(q,s,prec,prec);
	  glPopMatrix();
	  gluDeleteQuadric(q);
        }
	else if (sty==WIRE) {
	   glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	  glBegin(GL_QUADS);
	  glNormal3f(0,0,1);
	  glVertex3f(x+s,y+s,z+s);
	  glVertex3f(x+s,y-s,z+s);
	  glVertex3f(x-s,y-s,z+s);
	  glVertex3f(x-s,y+s,z+s);

	  glNormal3f(0,0,-1);
	  glVertex3f(x+s,y+s,z-s);
	  glVertex3f(x+s,y-s,z-s);
	  glVertex3f(x-s,y-s,z-s);
	  glVertex3f(x-s,y+s,z-s);

	  glNormal3f(0,1,0);
	  glVertex3f(x+s,y+s,z+s);
	  glVertex3f(x+s,y+s,z-s);
	  glVertex3f(x-s,y+s,z-s);
	  glVertex3f(x-s,y+s,z+s);

	  glNormal3f(0,-1,0);
	  glVertex3f(x+s,y-s,z+s);
	  glVertex3f(x+s,y-s,z-s);
	  glVertex3f(x-s,y-s,z-s);
	  glVertex3f(x-s,y-s,z+s);

	  glNormal3f(1,0,0);
	  glVertex3f(x+s,y+s,z+s);
	  glVertex3f(x+s,y+s,z-s);
	  glVertex3f(x+s,y-s,z-s);
	  glVertex3f(x+s,y-s,z+s);

	  glNormal3f(-1,0,0);
	  glVertex3f(x-s,y+s,z+s);
	  glVertex3f(x-s,y+s,z-s);
	  glVertex3f(x-s,y-s,z-s);
	  glVertex3f(x-s,y-s,z+s);
	  glEnd();
	}
	  else if (sty==RAW) {
	    glPointSize(s);
	    glBegin(GL_POINTS);
	    glVertex3f(x,y,z);
	    glEnd();
	  }
}

void draw_edge(edge e, Size s, Style sty, Precision prec)
{
  double x1=e[0][0];
  double y1=e[0][1];
  double z1=e[0][2];
  double x2=e[1][0];
  double y2=e[1][1];
  double z2=e[1][2];
  if (sty==FILL) {
     GLUquadricObj *q= gluNewQuadric();
     std::pair<double, double> ag=get_angles(x1,y1,z1,x2,y2,z2);
     double l=sqrt(pow(x1-x2,2) +pow(y1-y2,2) + pow(z1-z2,2));
     glPushMatrix();
     glTranslatef(x1, y1, z1);
     glRotatef(ag.first,1,0,0);
     glRotatef(ag.second,0,1,0);
     gluQuadricNormals(q, (GLenum) GL_SMOOTH);
     gluCylinder(q, s, s, l, prec, 1);
     glPopMatrix();
     gluDeleteQuadric(q);
   }
   else {
     glLineWidth(s);
     glBegin(GL_LINES);
       glNormal3f(0,0,1);
       glVertex3f(x1,y1,z1);
       glVertex3f(x2,y2,z2);
     glEnd();
   }
}

void draw_facet(facet f, Size s, Style sty, Precision prec)
{
  facet::iterator fit=f.begin();
  vertex v1=*fit;
  vertex v2=*(++fit);
  vertex v3=*(++fit);
  std::vector<double> v(3);
  v=normal(v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],v3[0],v3[1],v3[2]);
  glNormal3d(v[0],v[1],v[2]);
  int l=f.size()/2;
  int i=0;
  glBegin(GL_POLYGON);
  for (fit=f.begin(); fit!=f.end();fit++) {
      glVertex3f((*fit)[0],(*fit)[1],(*fit)[2]);
      i++;
      if (i==l)
	glNormal3d(-v[0],-v[1],-v[2]);
  }
  glEnd();
}


void draw_list_vertices(list_vertices v, Size s, Style sty, Precision prec)
{
  list_vertices::iterator it;
  for (it=v.begin(); it!=v.end(); it++)
    draw_vertex(*it,s,sty,prec);
}

void draw_list_edges(list_edges v, Size s, Style sty, Precision prec)
{
  list_edges::iterator it;
  for (it=v.begin(); it!=v.end(); it++)
    draw_edge(*it,s,sty,prec);
}

void draw_list_facets(list_facets v, Size s, Style sty, Precision prec)
{
  list_facets::iterator it;
  for (it=v.begin(); it!=v.end(); it++)
    draw_facet(*it,s,sty,prec);
}

void draw_wire_hidden(list_facets v, Size s, Color fcol, Color bcol)
{
  glLineWidth(s);
  list_facets::iterator it;
  glEnable(GL_STENCIL_TEST);
  glClear(GL_STENCIL_BUFFER_BIT);
  glStencilFunc(GL_ALWAYS,0,1);
  glStencilOp(GL_INVERT,GL_INVERT,GL_INVERT);
  set_color(fcol);
  for (it=v.begin(); it!=v.end(); it++) {
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    draw_facet((*it),2,WIRE,0);

    set_color(bcol);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glStencilFunc(GL_EQUAL,0,1);
    glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
    draw_facet((*it),0,RAW,0);

    set_color(fcol);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glStencilFunc(GL_ALWAYS,0,1);
    glStencilOp(GL_INVERT,GL_INVERT,GL_INVERT);
    draw_facet((*it),1,WIRE,0);
  }
  glDisable(GL_STENCIL_TEST);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
}



template<class facets_object_3>
class Drawable_facets_object_3: public  Drawable_object
{
private:

list_vertices list_v;
list_edges    list_e;
list_facets   list_f;

public:
Drawable_facets_object_3(const facets_object_3 &obj ,Color c, Color
			 c2=BLACK, char* name="Facet Object", Style
			 sty=FOS1, Size s=2, Precision prec=10)
    {
      list_v=get_vertices(obj);
      list_e=get_edges(obj);
      list_f=get_facets(obj);
      set_center();
      lind=0;
      color=c; col2=c2;type=name;style=sty;size=s;precision=prec;
    }

void set_center()
    {
      list_vertices::iterator it;
      o_center[0]=0; o_center[1]=0;o_center[2]=0;
      for (it=list_v.begin(); it!=list_v.end();it++) {
	o_center[0]+=(*it)[0];o_center[1]+=(*it)[1];o_center[2]+=(*it)[2];
      }
      o_center[0] = o_center[0]/list_v.size();
      o_center[1] = o_center[1]/list_v.size();
      o_center[2] = o_center[2]/list_v.size();
    }

void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	switch(style) {
	case FOS1:
	  set_color(color);
	  draw_list_edges(list_e,size,WIRE,0);
	  break;
	case FOS2:
	  set_color(color);
	  draw_list_edges(list_e,size,WIRE,0);
	  set_color(col2);
	  draw_list_vertices(list_v,size*2,WIRE,0);
	  break;
	case FOS3:
	  draw_wire_hidden(list_f,size,color,col2);
	  break;
	case FOS4:
	  set_color(color);
	  draw_list_facets(list_f,0,RAW,0);
	  break;
	case FOS5:
	  set_color(color);
	  draw_list_edges(list_e,2*size,FILL,precision);
	  set_color(col2);
	  draw_list_vertices(list_v,3*size,FILL,precision);
	  break;
	default:
	  draw_wire_hidden(list_f,size,color,col2);
	  break;
	}
	glEndList();
      }
    }
};


CGAL_END_NAMESPACE
