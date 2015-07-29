#ifndef SCENE_POLYGON_SOUP_H
#define SCENE_POLYGON_SOUP_H

#include "Scene_polygon_soup_config.h"
#include <CGAL_demo/Scene_item.h>
#include <iostream>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

struct Polygon_soup;

class SCENE_POLYGON_SOUP_EXPORT Scene_polygon_soup 
  : public Scene_item
{
  Q_OBJECT
public:  
  Scene_polygon_soup();
  ~Scene_polygon_soup();

  Scene_polygon_soup* clone() const;
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return m != Gouraud; } // CHECK THIS!
  // OpenGL drawing in a display list
  void draw(Viewer*) const;
  void draw_edges(Viewer*) const;
  void draw_points(Viewer*) const;

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  void shuffle_orientations();
  bool orient();
  void inside_out();

  void setDisplayNonManifoldEdges(const bool);
  bool displayNonManifoldEdges() const;
private:
  Polygon_soup* soup;
  bool oriented;


  static const int vaoSize = 3;
  static const int vboSize = 6;

  mutable int poly_vertexLocation[1];
  mutable int normalsLocation[1];
  mutable int mvpLocation[1];
  mutable int mvLocation[1];
  mutable int colorLocation[1];
  mutable int lightLocation[5];
  mutable int twosideLocation;

  std::vector<float> v_poly;
  std::vector<float> v_edge;
  std::vector<float> v_points;
  std::vector<float> normal;
  std::vector<float> vertex_nm;

  mutable bool are_buffers_initialized;
  mutable QOpenGLBuffer buffers[vboSize];
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  void initialize_buffers() const;
  void compute_elements();
  void attrib_buffers(Viewer*) const;
  void compile_shaders();

public Q_SLOTS:
  void changed();

}; // end class Scene_polygon_soup

#endif // SCENE_POLYGON_SOUP_H
