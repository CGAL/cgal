#ifndef CGAL_VOLUME_PLANE_INTERSECTION_H_
#define CGAL_VOLUME_PLANE_INTERSECTION_H_

#include <CGAL_demo/Scene_item.h>

#include <QColor>
#include <QString>
#include<QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

class Volume_plane_interface;

class Volume_plane_intersection
  : public Scene_item {
  typedef std::pair<Volume_plane_interface*, Volume_plane_interface*> Interface_pair;
Q_OBJECT
public:
  Volume_plane_intersection(float x, float y, float z)
    : a(NULL), b(NULL), c(NULL), x(x), y(y), z(z) {
    setColor(QColor(255, 0, 0));
    setName("Volume plane intersection");
    compile_shaders();
    compute_elements();
    init_buffers();
  }

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  bool manipulatable() const { return false; }
  Volume_plane_intersection* clone() const { return 0; }
  bool supportsRenderingMode(RenderingMode) const { return true; }
  QString toolTip() const { return "Tooling"; }

  void draw(Viewer*)const;

  void setX(Volume_plane_interface* x) { a = x; }
  void setY(Volume_plane_interface* x) { b = x; }
  void setZ(Volume_plane_interface* x) { c = x; }

public Q_SLOTS:
  void planeRemoved(Volume_plane_interface* i) {
    if(a == i) {
      a = NULL;
    } else if(b == i) {
      b = NULL;
    } else if(c == i) {
      c = NULL;
    }
  }

private:
  Volume_plane_interface *a, *b, *c;
  float x, y, z;

  static const int vaoSize = 3;
  static const int vboSize = 3;

  mutable int vertexLocation[1];
  mutable int mvpLocation[1];

  std::vector<float> a_vertex;
  std::vector<float> b_vertex;
  std::vector<float> c_vertex;

  mutable QOpenGLBuffer buffers[vboSize];
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  void compute_elements();
  void init_buffers();
  void attrib_buffers(Viewer*) const;
  void compile_shaders();
};

#endif /* CGAL_VOLUME_PLANE_INTERSECTION_H_ */

