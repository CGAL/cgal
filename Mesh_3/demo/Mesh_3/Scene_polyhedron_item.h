#ifndef SCENE_POLYHEDRON_ITEM_H
#define SCENE_POLYHEDRON_ITEM_H

#include "Scene_polyhedron_item_config.h"
#include <CGAL_demo/Scene_item.h>
#include "Polyhedron_type_fwd.h"
#include <iostream>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_ITEM_EXPORT Scene_polyhedron_item 
  : public Scene_item{
  Q_OBJECT
public:  
  Scene_polyhedron_item();
//   Scene_polyhedron_item(const Scene_polyhedron_item&);
  Scene_polyhedron_item(const Polyhedron& p);
  Scene_polyhedron_item(Polyhedron* const p);
  ~Scene_polyhedron_item();

  Scene_polyhedron_item* clone() const;
  
  // IO
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode) const { return true; }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void draw(Viewer*) const ;
  virtual void draw_edges(Viewer*) const ;
  virtual void draw_points(Viewer*) const ;

  // Get wrapped polyhedron
  Polyhedron*       polyhedron();
  const Polyhedron* polyhedron() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

private:
  Polyhedron* poly;
  mutable
  bool are_buffers_initialized;

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
   std::vector<float> normal_flat;
   std::vector<float> normal_smooth;
   std::vector<float> vertex_nm;


  mutable QOpenGLBuffer buffers[vboSize];
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  void initialize_buffers() const;
  void compute_elements();
  void attrib_buffers(Viewer*) const;
  void compile_shaders();

public Q_SLOTS:
    void changed();

}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_ITEM_H
