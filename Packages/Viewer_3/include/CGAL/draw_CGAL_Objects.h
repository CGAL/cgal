// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/draw_CGAL_Objects.h
// revision      : $Revision$
//
// author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
#include <CGAL/Cartesian.h>
#include <CGAL/Aff_transformation_3.h>
#ifndef V_UTILS
#include <CGAL/v_utils.h>
#endif 
#ifndef DRAWABLE 
#include <CGAL/Drawable_object.h>
#endif

#ifndef DRAW_CGAL_OBJECTS_H
#define DRAW_CGAL_OBJECTS_H


typedef CGAL::Cartesian <double > CT;
typedef CGAL::Aff_transformation_3 < CT> Transformation;

// style for drawing objects
// enum Style {FILL=1, WIRE, RAW, FOS1, FOS2, FOS3, FOS4, FOS5};

CGAL_BEGIN_NAMESPACE


//#### DRAWABLE POINT #######
template<class Point3>
class Drawable_point_3: public  Drawable_object
{
private:
  Point3 pt;
  
public:

  ~Drawable_point_3(){}
  Drawable_point_3(){type="Point";}
  Drawable_point_3(const Point3 &p, Color c, Style sty=FILL, Size s =
		   5 , Precision prec=15)
    {
      pt=p;
      color = c; size = s;
      set_center(); precision=prec; lind=0; style=sty;type="Point";
    }

  Drawable_point_3(const double x,const double y,const double z,
		   Color c, Style sty=FILL, Size s = 5 , Precision
		   prec=15)
    {
      pt=Point3(x,y,z);
      color = c; size = s;
      set_center(); precision=prec; lind=0; style=sty;type="Point";
    }
  

  void set_center()
    {
      o_center[0]=to_double(pt.x()); o_center[1]=to_double(pt.y()); o_center[2]=to_double(pt.z());
    }




  void draw() 
    {
      
      if (lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	set_color(color);
	if (style==WIRE) {
	  GLUquadricObj *q= gluNewQuadric();
	  glPushMatrix();
	  glTranslatef(to_double(pt.x()), to_double(pt.y()), to_double(pt.z()));
	  gluQuadricNormals(q, (GLenum) GL_SMOOTH);
	  gluQuadricDrawStyle(q,(GLenum) GLU_LINE);
	  glLineWidth(1);
	  gluSphere(q,size,precision,precision);
	  glPopMatrix();
	  gluDeleteQuadric(q);
	}
	else if (style==FILL) {
	  GLUquadricObj *q= gluNewQuadric();
	  glPushMatrix();
	  glTranslatef(to_double(pt.x()), to_double(pt.y()), to_double(pt.z()));
	  gluQuadricNormals(q, (GLenum) GL_SMOOTH);
	  gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
	  gluSphere(q,size,precision,precision);
	  glPopMatrix();
	  gluDeleteQuadric(q);
        }
	else {
	  glPointSize(size);
	  glBegin(GL_POINTS);
	  glVertex3f(to_double(pt.x()),to_double(pt.y()),to_double(pt.z()));
	  glEnd();
	}
	glEndList();
      }
    }

  void to_ps(PS_Stream_3 &ps)
    {
       CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

     ps << pt;
    }


};


// #### DRAWABLE SEGMENT #####
// draw a segment as a tube.

template<class segment3>
class Drawable_segment_3: public  Drawable_object
{
private:
  segment3 seg;
public:

  Drawable_segment_3(){type="Segment";}
  Drawable_segment_3(const segment3 &sg,Color c, Style sty= RAW, Size
		     s=5, Precision prec=15)
    {
      seg=sg;
      color = c; size=s; style = sty;
      set_center(); precision=prec;lind=0;type="Segment";
    }

  void set_center()
    {
      o_center[0]=(to_double(seg.source().x())+to_double(seg.target().x()))/2;
      o_center[1]=(to_double(seg.source().y())+to_double(seg.target().y()))/2;
      o_center[2]=(to_double(seg.source().z())+to_double(seg.target().z()))/2;
    }
  
  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	double x1=to_double(seg.source().x());
	double y1=to_double(seg.source().y());
	double z1=to_double(seg.source().z());
	double x2=to_double(seg.target().x());
	double y2=to_double(seg.target().y());
	double z2=to_double(seg.target().z());
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	set_color(color);
	if (style==FILL) {
	  GLUquadricObj *q= gluNewQuadric();
	  std::pair<double, double> ag=get_angles(x1,y1,z1,x2,y2,z2);
	  double l=sqrt(pow(x1-x2,2) +pow(y1-y2,2) + pow(z1-z2,2));
	  glPushMatrix();
	  glTranslatef(x1, y1, z1);
	  glRotatef(ag.first,1,0,0);
	  glRotatef(ag.second,0,1,0);
	  gluQuadricNormals(q, (GLenum) GL_SMOOTH);
	  gluCylinder(q, size, size, l, precision, 1);
	  glPopMatrix();
	  gluDeleteQuadric(q);
	}
	else {
	  glLineWidth(size);
	  glBegin(GL_LINES);
	  glVertex3f(x1,y1,z1);
	  glVertex3f(x2,y2,z2);
	  glEnd();
	}
	
	glEndList();
      }
    }
  
  void to_ps(PS_Stream_3 &ps)
    {
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      ps << seg;
    }

};




//### DRAWABLE SET OF POINTS ######
// The type pointed by InputIterator must be Point
template < class InputIterator, class Point >
class Drawable_points_set_3: public  Drawable_object
{
private:
  std::list<Point> LP;

public:
  Drawable_points_set_3(){type="Point Set";}
  Drawable_points_set_3(InputIterator first, InputIterator
 			last,Color c, Style sty=FILL, Size s=5,
			Precision prec=15)
    {
      while(first!=last) {
        LP.push_back(*first);
        ++first;
      }
      color = c; size = s;
      set_center(); precision=prec; lind=0;
      style=sty;type="Point Set";
    }
  
  
  
  void draw() 
    {
      if (lind)
	glCallList(lind);
      else 
	{
  	  typename std::list<Point>::iterator it;
	  lind = glGenLists(1);
	  glNewList(lind,GL_COMPILE_AND_EXECUTE);
	  set_color(color);
	  for (it=LP.begin();it!=LP.end();it++) {
	    GLUquadricObj *q= gluNewQuadric();
	    glPushMatrix();
	    glTranslatef(to_double(it->x()), to_double(it->y()), to_double(it->z()));
	    gluQuadricNormals(q, (GLenum) GL_SMOOTH);
	    if (style != FILL) {
	      gluQuadricDrawStyle(q,(GLenum) GLU_LINE);
	      glLineWidth(1);
	    }
	    else
	      gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
	    gluSphere(q,size,precision,precision);
	    glPopMatrix();
	    gluDeleteQuadric(q);
	  }
	  glEndList();
	}
    }


  void set_center()
    {
      typename std::list<Point>::iterator it;
      int s=LP.size();
      double o1=0; double o2=0; double o3=0;
      for (it=LP.begin();it!=LP.end();it++) {
	o1+=to_double(it->x());
	o2+=to_double(it->y());
	o3+=to_double(it->z());
      }
      o_center[0]=o1/s; o_center[1]=o2/s; o_center[2]=o3/s;
    }

  void add_point(double x, double y, double z) 
    {
      glDeleteLists(lind,1);
      lind=0;
      LP.push_back(Point(x,y,z));
      set_center();
    }

  void to_ps(PS_Stream_3 &ps)
    {
      typename std::list<Point>::iterator it;
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      for (it=LP.begin();it!=LP.end();it++) {
	ps << (*it);
      }
    }
};






template<class triangle>
class Drawable_triangle_3: public  Drawable_object
{
private:
  triangle trg;

public:

  Drawable_triangle_3(){type="Triangle";}
  Drawable_triangle_3(const triangle &tri,Color c, Style sty=WIRE,
		      Size s=2, Precision prec = 0)
    {
      trg=tri;
      color = c; size=s;
      set_center();lind=0;type="Triangle";style=sty;
    }

  void set_center()
    {
      o_center[0]=to_double((trg[0].x()+trg[1].x()+trg[2].x()))/3; 
      o_center[1]=to_double((trg[0].y()+trg[1].y()+trg[2].y()))/3;
      o_center[2]=to_double((trg[0].z()+trg[1].z()+trg[2].z()))/3;
    }

  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	glLineWidth(size);
	set_color(color);
	if (style==WIRE) 
	  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	else
	  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    if ((style==WIRE) || (style==RAW))
      draw_triangle(to_double(trg[0].x()),to_double(trg[0].y()),to_double(trg[0].z()),to_double(trg[1].x()),to_double(trg[1].y()),to_double(trg[1].z()),to_double(trg[2].x()),to_double(trg[2].y()),to_double(trg[2].z()));
    
    else
      draw_triangle_nice(to_double(trg[0].x()),to_double(trg[0].y()),to_double(trg[0].z()),to_double(trg[1].x()),to_double(trg[1].y()),to_double(trg[1].z()),to_double(trg[2].x()),to_double(trg[2].y()),to_double(trg[2].z()));
    glEndList();
      }
    }

  void to_ps(PS_Stream_3 &ps)
    {
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      ps << trg;
    }
  
};

template<class tetrahedron>
class Drawable_tetrahedron_3: public  Drawable_object
{
private:
  tetrahedron tr;
  
public:
  
  Drawable_tetrahedron_3(){type="Tetrahedron";}
  Drawable_tetrahedron_3(const tetrahedron &tet,Color c, Style
			 sty=FILL, Size s=5, Precision prec=15)
    
    {
      tr=tet;
      color = c; size=s;
      set_center();lind=0;type="Tetrahedron";style=sty;precision=prec;
    }
  
  void set_center()
    {
      o_center[0]=to_double((tr[0].x()+tr[1].x()+tr[2].x()+tr[3].x())/3); 
      o_center[1]=to_double((tr[0].y()+tr[1].y()+tr[2].y()+tr[3].y())/3); 
      o_center[2]=to_double((tr[0].z()+tr[1].z()+tr[2].z()+tr[3].z())/3); 
    }
  
  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	glLineWidth(size);
	set_color(color);
	if (style==RAW) {
	  glBegin(GL_LINES);
	  glVertex3d(to_double(tr[0].x()),to_double(tr[0].y()),to_double(tr[0].z()));
	  glVertex3d(to_double(tr[1].x()),to_double(tr[1].y()),to_double(tr[1].z()));
	  glVertex3d(to_double(tr[0].x()),to_double(tr[0].y()),to_double(tr[0].z()));
	  glVertex3d(to_double(tr[3].x()),to_double(tr[3].y()),to_double(tr[3].z()));
	  glVertex3d(to_double(tr[0].x()),to_double(tr[0].y()),to_double(tr[0].z()));
	  glVertex3d(to_double(tr[2].x()),to_double(tr[2].y()),to_double(tr[2].z()));
	  glVertex3d(to_double(tr[1].x()),to_double(tr[1].y()),to_double(tr[1].z()));
	  glVertex3d(to_double(tr[2].x()),to_double(tr[2].y()),to_double(tr[2].z()));
	  glVertex3d(to_double(tr[1].x()),to_double(tr[1].y()),to_double(tr[1].z()));
	  glVertex3d(to_double(tr[3].x()),to_double(tr[3].y()),to_double(tr[3].z()));
	  glVertex3d(to_double(tr[2].x()),to_double(tr[2].y()),to_double(tr[2].z()));
	  glVertex3d(to_double(tr[3].x()),to_double(tr[3].y()),to_double(tr[3].z()));
	  
	  glEnd();
	}
	else if (style==WIRE){
	  draw_tube(to_double(tr[0].x()),to_double(tr[0].y()),
		    to_double(tr[0].z()),to_double(tr[1].x()),
		    to_double(tr[1].y()),to_double(tr[1].z()),size,
		    precision);
	  
	  draw_tube(to_double(tr[0].x()),to_double(tr[0].y()),
		    to_double(tr[0].z()),to_double(tr[2].x()),
		    to_double(tr[2].y()),to_double(tr[2].z()),size,
		    precision);
	  draw_tube(to_double(tr[0].x()),to_double(tr[0].y()),
		    to_double(tr[0].z()),to_double(tr[3].x()),
		    to_double(tr[3].y()),to_double(tr[3].z()),size,
		    precision);
	  draw_tube(to_double(tr[1].x()),to_double(tr[1].y()),
		    to_double(tr[1].z()),to_double(tr[2].x()),
		    to_double(tr[2].y()),to_double(tr[2].z()),size,
		    precision);
	  draw_tube(to_double(tr[1].x()),to_double(tr[1].y()),
		    to_double(tr[1].z()),to_double(tr[3].x()),
		    to_double(tr[3].y()),to_double(tr[3].z()),size,
		    precision);
	  draw_tube(to_double(tr[2].x()),to_double(tr[2].y()),
		    to_double(tr[2].z()),to_double(tr[3].x()),
		    to_double(tr[3].y()),to_double(tr[3].z()),size,
		    precision);
	  
	  draw_sphere(to_double(tr[0].x()),to_double(tr[0].y()),
		      to_double(tr[0].z()),size+1, precision);
	  draw_sphere(to_double(tr[1].x()),to_double(tr[1].y()),
		      to_double(tr[1].z()),size+1, precision);
	  draw_sphere(to_double(tr[2].x()),to_double(tr[2].y()),
		      to_double(tr[2].z()),size+1, precision);
	  draw_sphere(to_double(tr[3].x()),to_double(tr[3].y()),
		      to_double(tr[3].z()),size+1, precision);
	  
	  
	}
	else {
	  draw_triangle(to_double(tr[0].x()), to_double(tr[0].y()),
			to_double(tr[0].z()), to_double(tr[1].x()),
			to_double(tr[1].y()), to_double(tr[1].z()),
			to_double(tr[2].x()), to_double(tr[2].y()),
			to_double(tr[2].z()));
	  draw_triangle(to_double(tr[0].x()), to_double(tr[0].y()),
			to_double(tr[0].z()), to_double(tr[1].x()),
			to_double(tr[1].y()), to_double(tr[1].z()),
			to_double(tr[3].x()), to_double(tr[3].y()),
			to_double(tr[3].z()));
	  draw_triangle(to_double(tr[0].x()), to_double(tr[0].y()),
			to_double(tr[0].z()), to_double(tr[3].x()),
			to_double(tr[3].y()), to_double(tr[3].z()),
			to_double(tr[2].x()), to_double(tr[2].y()),
			to_double(tr[2].z()));
	  draw_triangle(to_double(tr[1].x()), to_double(tr[1].y()),
			to_double(tr[1].z()), to_double(tr[2].x()),
			to_double(tr[2].y()), to_double(tr[2].z()),
			to_double(tr[3].x()), to_double(tr[3].y()),
			to_double(tr[3].z()));
	  
	  
	}

	glEndList();
	
      }
    }
  
  void to_ps(PS_Stream_3 &ps)
    {
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      ps << tr;
    }
};





template<class line3>
class Drawable_line_3: public  Drawable_object
{
private:
  double x1,y1,z1,x2,y2,z2;
  line3 l;
public:

  Drawable_line_3(){type="Line";}
  Drawable_line_3(const line3 &ln, Color c, Style sty=RAW, Size s=2 , Precision
  		   prec=15)
    {
      l= ln;
      x1=to_double(l.point(5000).x());
      y1=to_double(l.point(5000).y());
      z1=to_double(l.point(5000).z());
      x2=to_double(l.point(-5000).x());
      y2=to_double(l.point(-5000).y());
      z2=to_double(l.point(-5000).z());
      color = c; size=s;style=sty;
      set_center(); precision=prec;lind=0;type="Line";
    }

  void set_center()
    {
      o_center[0]=(x1+x2)/2; o_center[1]=(y1+y2)/2;
      o_center[2]=(z1+z2)/2;
    }

void draw()
  {
  if(lind)
    glCallList(lind);
  else {
   lind = glGenLists(1);
   glNewList(lind,GL_COMPILE_AND_EXECUTE);
   set_color(color);
   if (style==FILL) {
     GLUquadricObj *q = gluNewQuadric();
     std::pair<double, double> ag=get_angles(x1,y1,z1,x2,y2,z2);
     double l=sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
     glPushMatrix();
     glTranslatef(x1, y1, z1);
     glRotatef(ag.first,1,0,0);
     glRotatef(ag.second,0,1,0);
     gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
     gluQuadricNormals(q, (GLenum) GL_SMOOTH);
     gluCylinder(q, size, size, l, precision, precision);
     glPopMatrix();
     gluDeleteQuadric(q);
   }
   else {
     glLineWidth(size);
     glBegin(GL_LINES);
       glVertex3f(x1,y1,z1);
       glVertex3f(x2,y2,z2);
     glEnd();
   }

   glEndList();
  }
  }

  
  void to_ps(PS_Stream_3 &ps)
    {
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      ps << l;
    }

};


template<class ray3>
class Drawable_ray_3: public  Drawable_object {
private:
  double x1 , y1, z1,x2,y2,z2 ;
  ray3 r;
public:

  Drawable_ray_3(){type="Ray";}
  Drawable_ray_3(const ray3 &ry, Color c, Style sty=RAW, Size s=2 , Precision
  		   prec=15)
    {
      r=ry;
      x1=to_double(r.source().x());
      y1=to_double(r.source().y());
      z1=to_double(r.source().z());
      x2=to_double(r.point(5000).x());
      y2=to_double(r.point(5000).y());
      z2=to_double(r.point(5000).z());
      color = c; size=s;style=sty;
      set_center(); precision=prec;lind=0;type="Ray";
    }

  void set_center()
    {
      o_center[0]=x1; o_center[1]=y1;
      o_center[2]=z1;
    }

void draw()
  {
  if(lind)
    glCallList(lind);
  else {
   lind = glGenLists(1);
   glNewList(lind,GL_COMPILE_AND_EXECUTE);
   set_color(color);
   if (style==FILL) {
     GLUquadricObj *q= gluNewQuadric();
     std::pair<double, double> ag=get_angles(x1,y1,z1,x2,y2,z2);
     double l=sqrt(pow(x1-x2,2) +pow(y1-y2,2) + pow(z1-z2,2));
     glPushMatrix();
     glTranslatef(x1, y1, z1);
     glRotatef(ag.first,1,0,0);
     glRotatef(ag.second,0,1,0);
     gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
     gluQuadricNormals(q, (GLenum) GL_SMOOTH);
     gluCylinder(q, size, size, l, precision, precision);
     glPopMatrix();
     gluDeleteQuadric(q);
   }
   else {
     glLineWidth(size);
     glBegin(GL_LINES);
       glVertex3f(x1,y1,z1);
       glVertex3f(x2,y2,z2);
     glEnd();
   }

   glEndList();
  }
  }
  
  void to_ps(PS_Stream_3 &ps)
    {
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      ps << r;
    }

};



// triangulation 2D avec des points 3D
template<class triangulation_3>
class Drawable_triangulation_3: public  CGAL::Drawable_object
{
private:
  triangulation_3 tr;
  
public:  
  
  Drawable_triangulation_3(){type="Triangulation_3";}
  
  
  Drawable_triangulation_3(const triangulation_3 &tet,Color c1, Color
			   c2=BLACK, Style sty=WIRE, Size s=3, Precision prec=10)
    
    {
      tr=tet;
      color = c1; col2 = c2; size=s;
      set_center();lind=0;type="Triangulation_3";style=sty;precision=prec;
    }
  
  void set_center()
    {
      
      typename triangulation_3::Vertex_iterator vit;
      o_center[0]=0;o_center[1]=0;o_center[2]=0;
      for (vit=tr.finite_vertices_begin() ; vit != tr.vertices_end(); vit++) {
	o_center[0]= o_center[0] + to_double(vit->point().x());
	o_center[1]= o_center[1] +  to_double(vit->point().y());
	o_center[2]= o_center[2] +  to_double(vit->point().z());
      }
      o_center[0] = o_center[0]/tr.number_of_vertices();
      o_center[1] = o_center[1]/tr.number_of_vertices();
      o_center[2] = o_center[2]/tr.number_of_vertices();
    }
  
  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	int no;
	typedef typename triangulation_3::Edge_iterator Edge_iterator;
	typedef typename triangulation_3::Vertex_handle Vertex_handle;
	typedef typename triangulation_3::Face_handle Face_handle;
	Vertex_handle v1, v2;
	Face_handle f;
	Edge_iterator it = tr.finite_edges_begin();
	Edge_iterator   beyond = tr.edges_end();
	
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
        if (style==WIRE) {
	  set_color(color);
	  for ( ;it != beyond; ++it) {
	    f = (*it).first;
	    no = (*it).second;
	    v1 = f->vertex(f->ccw(no));
	    v2 = f->vertex(f->cw(no));
	    draw_tube(to_double(v1->point().x()),to_double(v1->point().y()),to_double(v1->point().z()),to_double(v2->point().x()),to_double(v2->point().y()),to_double(v2->point().z()),size, precision);
	  }
	}
	else {
	  typename triangulation_3::Face_iterator fit;
	  set_color(color);

	  switch(style) {
	  case 1:
	    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	    break;
	  case 2:
	    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	    break;
	  case 3:
	    glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
	    break;
	  default:
	    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	  }

	  std::vector<float> vp1(3), vp2(3), vp3(3);
	  for (fit=tr.finite_faces_begin(); fit!=tr.finite_faces_end();fit++) {
	    vp1[0] =to_double(fit->vertex(0)->point().x());
	    vp1[1] =to_double(fit->vertex(0)->point().y());
	    vp1[2] =to_double(fit->vertex(0)->point().z());
	    vp2[0] =to_double(fit->vertex(1)->point().x());
	    vp2[1] =to_double(fit->vertex(1)->point().y());
	    vp2[2] =to_double(fit->vertex(1)->point().z());
	    vp3[0] =to_double(fit->vertex(2)->point().x());
	    vp3[1] =to_double(fit->vertex(2)->point().y());
	    vp3[2] =to_double(fit->vertex(2)->point().z());
	    draw_triangle(vp1[0],vp1[1],vp1[2]-1,vp2[0],vp2[1],vp2[2]-1,vp3[0],vp3[1],vp3[2]-1);
	  }
	  
	  it = tr.finite_edges_begin();
	  set_color(col2);
	  glLineWidth(2);
	  glBegin(GL_LINES);
	  for ( ;it != beyond; ++it) {
	    f = (*it).first;
	    no = (*it).second;
	    v1 = f->vertex(f->ccw(no));
	    v2 = f->vertex(f->cw(no));
	    glVertex3f(to_double(v1->point().x()),to_double(v1->point().y()),to_double(v1->point().z()));
	    glVertex3f(to_double(v2->point().x()),to_double(v2->point().y()),to_double(v2->point().z()));
	  }
	  glEnd();
	}
	glEndList();
      }
    }

 void add_point(double x, double y ,double z)
    {
      glDeleteLists(lind,1);
      lind=0;
      typedef typename triangulation_3::Point Point;
      tr.insert(Point(x,y,z));
      set_center();
    }

  void to_ps(PS_Stream_3 &ps)
    {
      CGAL::FILLING fill;

      switch(style) {
      default:
	fill = CGAL::NORMAL_FILL;
	break;
      case 2:
	fill = CGAL::WIRED_CULLBACK_FACING;
	break;
      }

      ps.set_current_filling(fill);
      ps.set_border_color(col2);
      ps.set_fill_color(color);

      ps.add_triangulation(tr);

    }
};





// ###################### 2-D drawable Objetct #######################


//#### DRAWABLE POINT #######
template<class Point2>
class Drawable_point_2: public  Drawable_object
{
private:
  Point2 pt;
public:

  Drawable_point_2(){type="Point";}
  Drawable_point_2(const Point2 &p, Color c, Style
		   sty=RAW, Size s = 5 , Precision
		   prec=5)
    {
     pt=p;
     x=to_double(p.x());
     y=to_double(p.y());
     color = c; size = s;
     set_center(); precision=prec; lind=0; style=sty;type="Point_2";
    }

   Drawable_point_2(const double x,const double y,
		    Color c, Style sty=RAW, Size s = 5 , Precision
		    prec=5)

    {
     pt=Point2(x,y);
     color = c; size = s;
     set_center(); precision=prec; lind=0; style=sty;type="Point_2";
    }


  void set_center()
    {
      o_center[0]=to_double(pt.x()); 
      o_center[1]=to_double(pt.y()); 
      o_center[2]=0;
    }

 void draw() 
    {
      if (lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	set_color(color);

	if (style==RAW) {
	  glPointSize(size);
	  glBegin(GL_POINTS);
	  glVertex3f(x,y,0);
	  glEnd();
	}
	else {
	  GLUquadricObj *q= gluNewQuadric();
	  glPushMatrix();
	  glTranslatef(to_double(pt.x()),to_double(pt.y()) , 0);
	  gluQuadricNormals(q, (GLenum) GL_NONE);
	  glLineWidth(1);
	  if (style==FILL) {
	    gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
	    gluDisk(q,0,size,precision,1);
	  }
	  else {
	    gluQuadricDrawStyle(q,(GLenum) GLU_LINE);
	    gluDisk(q,size, size,precision,1);
	  }
	  glPopMatrix();
	  gluDeleteQuadric(q);
	}
	glEndList();
      }
    }

};


template<class segment2>
class Drawable_segment_2: public  Drawable_object
{
private:

  segment2 seg;
public:

  Drawable_segment_2(){type="Segment";}
  Drawable_segment_2(const segment2 &sg,Color c, Style sty= RAW, Size
		     s=2, Precision prec=0)
    {
      seg=sg;
      color = c; size=s; style = sty;
      set_center(); precision=prec;lind=0;type="Segment_2";
    }

  void set_center()
    {
      o_center[0]=(to_double(seg.source().x())+to_double(seg.target().x()))/2;
      o_center[1]=(to_double(seg.source().y())+to_double(seg.target().y()))/2;
      o_center[2]=0;
    }

void draw()
  {
    if(lind)
      glCallList(lind);
    else {
      lind = glGenLists(1);
      glNewList(lind,GL_COMPILE_AND_EXECUTE);
      set_color(color);
      glLineWidth(size);
      glBegin(GL_LINES);
      glVertex2f(to_double(seg.source().x()),to_double(seg.source().y()));
      glVertex2f(to_double(seg.target().x()),to_double(seg.target().y()));
      glEnd();
    }
    glEndList();
  }
  
};


template<class line2>
class Drawable_line_2: public  Drawable_object
{
private:
  double x1 , y1, x2,y2;
  line2 l;
public:

  Drawable_line_2(){type="Line";}
  Drawable_line_2(const line2 &ln, Color c, Style sty=RAW, Size s=2 , Precision
  		   prec=15)
    {
      l=ln;
      x1=to_double(l.point(5000).x());
      y1=to_double(l.point(5000).y());
      x2=to_double(l.point(-5000).x());
      y2=to_double(l.point(-5000).y());
      color = c; size=s;style=sty;
      set_center(); precision=prec;lind=0;type="Line2";
    }

  void set_center()
    {
      o_center[0]=(x1+x2)/2; o_center[1]=(y1+y2)/2;
      o_center[2]=0;
    }
void draw()
  {
    if(lind)
      glCallList(lind);
    else {
      lind = glGenLists(1);
      glNewList(lind,GL_COMPILE_AND_EXECUTE);
      set_color(color);
      glLineWidth(size);
      glBegin(GL_LINES);
      glVertex2f(x1,y1);
      glVertex2f(x2,y2);
      glEnd();
    }
    glEndList();
  }

};



template<class ray2>
class Drawable_ray_2: public  Drawable_object
{
private:
  double x1 , y1,x2,y2;
  ray2 r;
public:

  Drawable_ray_2(){type="Ray2";}
  Drawable_ray_2(const ray2 &ry, Color c, Style sty=RAW, Size s=2 , Precision
  		   prec=15)
    {
      r=ry;
      x1=to_double(r.source().x());
      y1==to_double(r.source().y());
      x2==to_double(r.point(5000).x());
      y2==to_double(r.point(5000).y());
      color = c; size=s;style=sty;
      set_center(); precision=prec;lind=0;type="Ray2";
    }

  void set_center()
    {
      o_center[0]=x1; o_center[1]=y1;
      o_center[2]=0;
    }
void draw()
  {
    if(lind)
      glCallList(lind);
    else {
      lind = glGenLists(1);
      glNewList(lind,GL_COMPILE_AND_EXECUTE);
      set_color(color);
      glLineWidth(size);
      glBegin(GL_LINES);
      glVertex2f(x1,y1);
      glVertex2f(x2,y2);
      glEnd();
    }
    glEndList();
  }


};



template<class triangle2>
class Drawable_triangle_2: public  Drawable_object
{
private:
  triangle2 tr;

public:

  Drawable_triangle_2(){type="Triangle2";}
  Drawable_triangle_2(const triangle2 &tri,Color c, Color c2=BLACK, Style sty=FILL,
		      Size s=2, Precision prec = 0)
    {
      tr=tri;
      color = c; col2=c2; size=s;
      set_center();lind=0;type="Triangle2";style=sty;
    }

  void set_center()
    {
      o_center[0]=to_double((tr[0].x()+tr[1].x()+tr[2].x()))/3; 
      o_center[1]=to_double((tr[0].y()+tr[1].y()+tr[2].y()))/3;
      o_center[2]=0;
    }

void draw()
{
  if(lind)
    glCallList(lind);
  else {
    lind = glGenLists(1);
    glNewList(lind,GL_COMPILE_AND_EXECUTE);
    glLineWidth(size);
    set_color(color);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    draw_triangle_2(to_double(tr[0].x()),to_double(tr[0].y()),
		    to_double(tr[1].x()),to_double(tr[1].y()), 
		    to_double(tr[2].x()),to_double(tr[2].y()));

    if (style==FILL) {
      set_color(col2);
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      draw_shrink_triangle_2(to_double(tr[0].x()),to_double(tr[0].y()),
			     to_double(tr[1].x()),to_double(tr[1].y()),
			     to_double(tr[2].x()),to_double(tr[2].y()));
    }
    glEndList();
}
}
};

template<class circle2>
class Drawable_circle_2: public  Drawable_object
{
private:
  double x,y,r;
  circle2 crc;
public:

  Drawable_circle_2(){type="Circle2";}
  Drawable_circle_2(const circle2 &crl, Color c, Style
		   sty=RAW, Size s = 5 , Precision
		   prec=5)
    {
     crc=crl;
     x=to_double(crc.center().x());
     y=to_double(crc.center().y());
     r=sqrt(to_double(crc.squared_radius()));
     color = c; size = s;
     set_center(); precision=prec; lind=0; style=sty;type="Circle2";
    }



  void set_center()
    {
      o_center[0]=x; o_center[1]=y; o_center[2]=0;
    }

 void draw() 
    {
      if (lind)
	glCallList(lind);
      else {
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	set_color(color);
	GLUquadricObj *q= gluNewQuadric();
	glPushMatrix();
	glTranslatef(x, y, 0);
	gluQuadricNormals(q, (GLenum) GL_NONE);
	glLineWidth(size);
	if (style==FILL) {
	  gluQuadricDrawStyle(q,(GLenum) GLU_FILL);
	  gluDisk(q,0,r,precision,10);
	}
	else {
	  gluQuadricDrawStyle(q,(GLenum) GLU_LINE);
	  gluDisk(q,r, r,precision,1);
	  }
	glPopMatrix();
	gluDeleteQuadric(q);
      }
	glEndList();
      }

};

template<class triangulation_2>
class Drawable_triangulation_2: public  Drawable_object
{
private:
  triangulation_2 tr;
 
public:  
  

Drawable_triangulation_2(){type="Triangulation_2";}

Drawable_triangulation_2(const triangulation_2 &tet,Color c, Color c2, Style
			 sty=WIRE, Size s=5, Precision prec=15)
  {
    tr=tet;
    color = c; col2= c2; size=s;
      set_center();lind=0;type="Triangulation_2";style=sty;precision=prec;
    }

  void set_center()
    {

      typename triangulation_2::Vertex_iterator vit;
      o_center[0]=0;o_center[1]=0;o_center[2]=0;
      for (vit=tr.finite_vertices_begin() ; vit != tr.vertices_end(); vit++) {
	o_center[0]= o_center[0] + to_double(vit->point().x());
	o_center[1]= o_center[1] + to_double(vit->point().y());
      }
      o_center[0] = o_center[0]/tr.number_of_vertices();
      o_center[1] = o_center[1]/tr.number_of_vertices();
    }

  void draw()
    {
      if(lind)
	glCallList(lind);
      else {
	typename triangulation_2::Vertex_handle v1, v2;
	typename triangulation_2::Face_handle f;
	int no;
	lind = glGenLists(1);
	glNewList(lind,GL_COMPILE_AND_EXECUTE);
	set_color(color);
        glLineWidth(size);
	typename triangulation_2::Edge_iterator it;
	for (it = tr.edges_begin(); it != tr.edges_end(); it++) 
	  {
	    glBegin( GL_LINES );  
	    f = (*it).first;
	    no = (*it).second;
	    v1 = f->vertex(f->ccw(no));
	    v2 = f->vertex(f->cw(no));
	    glVertex2d(to_double(v1->point().x()),to_double(v1->point().y()));
	    glVertex2d(to_double(v2->point().x()),to_double(v2->point().y()));
	  }
	glEnd();
      }
      glEndList();

    }


 void add_point(double x, double y ,double z)
    {
      glDeleteLists(lind,1);
      lind=0;
      typedef typename triangulation_2::Point Point;
      tr.insert(Point(x,y));
      set_center();
    }
  
  void to_ps(PS_Stream_3 &ps)
    {
      ps.add_triangulation(tr);
    }


};



CGAL_END_NAMESPACE

#endif
