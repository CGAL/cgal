//Author: Sebastien Loriot sebastien.loriot@sophia.inria.fr

#ifndef _TOOL_H_
#define _TOOL_H_

#define NB_PTS 30
#define NB_PTS_SM 10
#define PTS_WIDTH  0.01f

#include <CGAL/glu.h>

template <class Kernel>
struct drawing_tools{
  typedef typename Kernel::Point_3 Point_3;
  
  static void drawSphere(GLUquadricObj *quadric,const Point_3 center,const double rad,const Point_3 color= Point_3(1,0,0) ){
    glPushMatrix();
    glColor3f(color.x(),color.y(),color.z());
    set_material();
    glTranslatef(center.x(),center.y(),center.z());
    gluSphere(quadric,rad,NB_PTS,NB_PTS);
    glPopMatrix();
  }
  
  static void draw_a_point_as_sphere(GLUquadricObj *quadric,const Point_3 center,const Point_3 color= Point_3(1,0,0) ){
    glPushMatrix();
    glColor3f(color.x(),color.y(),color.z());
    glTranslatef(center.x(),center.y(),center.z());
    gluSphere(quadric,PTS_WIDTH,NB_PTS_SM,NB_PTS_SM);
    glPopMatrix();
  }  

  static void drawCylinder(GLUquadricObj *quadric,Point_3& P1,Point_3& P2,double radius=PTS_WIDTH*2/3.,int subdivisions=NB_PTS_SM)
  {
    typename Kernel::Vector_3 V(P2.x()-P1.x(),P2.y()-P1.y(),P2.z()-P1.z());

    //handle the degenerate case of z1 == z2 with an approximation
    if(V.z() == 0)
      V=V+typename Kernel::Vector_3(0,0,.00000001);

    double v = sqrt( V.squared_length () );
    double ax = 57.2957795*acos( V.z()/v );
    if ( V.z() < 0.0 )
       ax = -ax;
    double rx = -V.y()*V.z();
    double ry = V.x()*V.z();
    glPushMatrix();

    //draw the cylinder body
    glTranslatef( P1.x(),P1.y(),P1.z() );
    glRotatef(ax, rx, ry, 0.0);
    gluQuadricOrientation(quadric,GLU_OUTSIDE);
    gluCylinder(quadric, radius, radius, v, subdivisions, 1);
    glPopMatrix();
  }


  
  template <class Edge>
  static void draw_edge(const Edge& e){
    glLineWidth(2);
    glBegin(GL_LINES);
    glColor3f(0,0,0);
    glVertex3f(e.first->vertex(e.second)->point().x(),e.first->vertex(e.second)->point().y(),e.first->vertex(e.second)->point().z());
    glVertex3f(e.first->vertex(e.third)->point().x(),e.first->vertex(e.third)->point().y(),e.first->vertex(e.third)->point().z());
    glEnd();
  }

  
  
  template<class Facet>
  static void draw_facet(const Facet& f){
#if 0    
  glBegin(GL_TRIANGLES);
  glColor3f(0,0,1); 
  glVertex3f(f[2].x(),f[2].y(),f[2].z());    
  glVertex3f(f[1].x(),f[1].y(),f[1].z());    
  glVertex3f(f[3].x(),f[3].y(),f[3].z());    
  glColor3f(0.5859,0.61719,0);
  glVertex3f(f[1].x(),f[1].y(),f[1].z());    
  glVertex3f(f[2].x(),f[2].y(),f[2].z());    
  glVertex3f(f[3].x(),f[3].y(),f[3].z());    
    
  glEnd();  
#else
    glDisable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    if ((f.second%2)==0)
      glColor3f(0,0,1);
    else
      glColor3f(0.5859,0.61719,0);
    glVertex3f(f.first->vertex((f.second+1)%4)->point().x(),f.first->vertex((f.second+1)%4)->point().y(),f.first->vertex((f.second+1)%4)->point().z());
    glVertex3f(f.first->vertex((f.second+2)%4)->point().x(),f.first->vertex((f.second+2)%4)->point().y(),f.first->vertex((f.second+2)%4)->point().z());
    glVertex3f(f.first->vertex((f.second+3)%4)->point().x(),f.first->vertex((f.second+3)%4)->point().y(),f.first->vertex((f.second+3)%4)->point().z());
    //Back in yellow
    if ((f.second%2)!=0)
      glColor3f(0,0,1);
    else
      glColor3f(0.5859,0.61719,0);
    glVertex3f(f.first->vertex((f.second+1)%4)->point().x(),f.first->vertex((f.second+1)%4)->point().y(),f.first->vertex((f.second+1)%4)->point().z());
    glVertex3f(f.first->vertex((f.second+3)%4)->point().x(),f.first->vertex((f.second+3)%4)->point().y(),f.first->vertex((f.second+3)%4)->point().z());
    glVertex3f(f.first->vertex((f.second+2)%4)->point().x(),f.first->vertex((f.second+2)%4)->point().y(),f.first->vertex((f.second+2)%4)->point().z());
    glEnd();
    glEnable(GL_LIGHTING);
#endif
  }
  
  template <class Cell_handle>
  static void draw_cell(Cell_handle c,GLUquadricObj* cylinder,GLUquadricObj* sphere){
    for (int i=0; i<4;++i){
      draw_a_point_as_sphere(sphere,c->vertex(i)->point(),Point_3(0,1,0));
      for (int j=i+1;j<4;++j)
        drawCylinder(cylinder,c->vertex(i)->point(),c->vertex(j)->point());
    }
  }
  
  // change material
  static void set_material() 
  {
    float ambient[]={0.0f,0.0f,0.0f,1.0f};
    float diffuse[]={0.0f,0.0f,0.0f,1.0f};
    float specular[]={0.0f,0.0f,0.0f,1.0f};
    float emission[]={0.3f,0.3f,0.3f,1.0f};
    float shininess[]={0.0f};

    // Ambient
    ambient[0] = 0.19225f;
    ambient[1] = 0.19225f;
    ambient[2] = 0.19225f;
    ambient[3] = 1.0f;
    // Diffuse
    diffuse[0] = 0.9f;
    diffuse[1] = 0.9f;
    diffuse[2] = 0.9f;
    diffuse[3] = 1.0f;
    // Specular
    specular[0] = 0.508273f;
    specular[1] = 0.508273f;
    specular[2] = 0.508273f;
    specular[3] = 1.0f;
    // Shininess
    shininess[0] = 50.0f;
    
    // apply
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,ambient);
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,diffuse);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,emission);
  }
  
};


#endif


