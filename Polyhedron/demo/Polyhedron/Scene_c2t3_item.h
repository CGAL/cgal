#ifndef SCENE_C2T3_ITEM_H
#define SCENE_C2T3_ITEM_H

#include "Scene_c2t3_item_config.h"
#include "Scene_item_with_display_list.h"
#include "C2t3_type.h"
#include <iostream>

#include <qgl.h>
#include <Qt/qglobal.h>
#include <CGAL/gl.h>

class Q_DECL_EXPORT Scene_c2t3_item : public Scene_item
{
  Q_OBJECT
public:
  Scene_c2t3_item(const C2t3& c2t3)
    : sphere_display_list(0), quadric(0), c2t3_(c2t3)
  {
  }

  ~Scene_c2t3_item()
  {
    if(quadric != 0)
      gluDeleteQuadric(quadric);
    if(sphere_display_list  != 0)
      glDeleteLists(sphere_display_list, 1);
  }

  C2t3& c2t3() {
    return c2t3_;
  }

  const C2t3& c2t3() const {
    return c2t3_;
  }
  bool isFinite() const { return true; }
  bool isEmpty() const {
    return c2t3().triangulation().number_of_vertices() == 0;
  }

  Bbox bbox() const {
    if(isEmpty())
      return Bbox();
    else {
      bool first = true;
      CGAL::Bbox_3 result;
      for(Tr::Finite_vertices_iterator
            vit = ++c2t3().triangulation().finite_vertices_begin(),
            end = c2t3().triangulation().finite_vertices_end();
          vit != end; ++vit)
      {
        if(vit->point().weight() > 0) {
          if(first) {
            result = vit->point().bbox();
            first = false;
          } else { 
            result = result + vit->point().bbox();
          }
        }
      }
      return Bbox(result.xmin(), result.ymin(), result.zmin(),
                  result.xmax(), result.ymax(), result.zmax());
    }
  }


  Scene_c2t3_item* clone() const {
    return 0;
  }

  QString toolTip() const {
    return tr("<p><b>2D complex in a 3D triangulation</b></p>"
              "<p>Number of vertices: %1<br />"
              "Number of surface facets: %2<br />")
      .arg(c2t3().triangulation().number_of_vertices())
      .arg(c2t3().number_of_facets());
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud); // CHECK THIS!
  }

  void draw() const {
    ::glBegin(GL_TRIANGLES);
    for(C2t3::Facet_iterator
          fit = c2t3().facets_begin(),
          end = c2t3().facets_end();
        fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Tr::Geom_traits::Point_3& pa = cell->vertex((index+1)&3)->point();
      const Tr::Geom_traits::Point_3& pb = cell->vertex((index+2)&3)->point();
      const Tr::Geom_traits::Point_3& pc = cell->vertex((index+3)&3)->point();
      draw_triangle(pa, pb, pc);
    }
    ::glEnd();
    
    GLenum gl_error = ::glGetError();
    if(gl_error != GL_NO_ERROR)
      std::cerr << "GL error: " << gluErrorString(gl_error) << std::endl;

    if(!draw_spheres)
      return;

    // force wireframe for protecting spheres
    GLint polygon_mode[2];
    ::glGetIntegerv(GL_POLYGON_MODE, &polygon_mode[0]);
    ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    for(Tr::Finite_vertices_iterator 
          vit = c2t3().triangulation().finite_vertices_begin(),
          end =  c2t3().triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      draw_sphere(vit->point());
    }
    ::glPolygonMode(GL_FRONT_AND_BACK, polygon_mode[0]);
  }

  void draw_sphere(const Tr::Point p) const 
  {
    if(p.weight() > 0) {
      if(sphere_display_list == 0) {
        sphere_display_list = glGenLists(1);
        if(sphere_display_list == 0)
          std::cerr << "ERROR: Cannot create display list!\n";
        if(quadric == 0)
          quadric = gluNewQuadric();
        if(quadric == 0)
          std::cerr << "ERROR: Cannot create GLU quadric!\n";
        glNewList(sphere_display_list, GL_COMPILE);
        gluSphere(quadric, 1., 10, 10);
        glEndList();
        if(glGetError() != GL_NO_ERROR)
          std::cerr << gluErrorString(glGetError());
      }
      glPushMatrix();
      glTranslated(CGAL::to_double(p.point().x()),
                   CGAL::to_double(p.point().y()),
                   CGAL::to_double(p.point().z()));
      const GLdouble r = CGAL::to_double(CGAL_NTS sqrt(p.weight()));
      glScaled(r, r, r);
      glCallList(sphere_display_list);
      glPopMatrix();
    }
  }

public slots:
  void show_spheres(bool b) {
    draw_spheres = b;
  }

private:
  static void draw_triangle(const Tr::Point& pa,
                            const Tr::Point& pb,
                            const Tr::Point& pc) {
    Tr::Geom_traits::Vector_3 n = cross_product(pb - pa, pc -pa);
    n = n / CGAL::sqrt(n*n);

    ::glNormal3d(n.x(),n.y(),n.z());

    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glVertex3d(pc.x(),pc.y(),pc.z());
  }

private:
  mutable GLuint sphere_display_list;
  mutable GLUquadric* quadric;
  C2t3 c2t3_;
  bool draw_spheres;
};

#endif // SCENE_C2T3_ITEM
