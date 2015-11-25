#ifndef SCENE_SEGMENTED_IMAGE_ITEM_H
#define SCENE_SEGMENTED_IMAGE_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include "Image_type_fwd.h"
#include "Scene_segmented_image_item_config.h"
#include <CGAL/gl.h>

#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>

typedef CGAL::Image_3 Image;
using namespace CGAL::Three;
class SCENE_SEGMENTED_IMAGE_ITEM_EXPORT Scene_segmented_image_item 
  : public Scene_item
{
  Q_OBJECT
public:

  Scene_segmented_image_item(Image* im, int drawing_scale);
  ~Scene_segmented_image_item();

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  void compute_bbox() const;

  Scene_segmented_image_item* clone() const { return NULL; }

  // rendering mode
  virtual bool supportsRenderingMode(RenderingMode m) const;

  // draw
  virtual void direct_draw(CGAL::Three::Viewer_interface* viewer) const
  { draw(viewer); }
  virtual void direct_draw_edges(CGAL::Three::Viewer_interface* viewer) const
  { draw_edges(viewer); }
  virtual void draw(CGAL::Three::Viewer_interface*) const;
  virtual void draw_edges(CGAL::Three::Viewer_interface* viewer) const
  { draw_gl(viewer); }
  
  virtual QString toolTip() const;
  
  const Image* image() const { return m_image; }

private:
  void draw_gl(CGAL::Three::Viewer_interface* viewer) const;
  
  void initialize_buffers();
  GLint ibo_size() const;
  
public:
  Image* m_image;

private:
  bool m_initialized;
#ifdef SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
  int m_voxel_scale;
  static const int vaoSize = 2;
  static const int vboSize = 6;
  mutable int poly_vertexLocation[1];
  mutable int normalsLocation[1];
  mutable int mvpLocation[1];
  mutable int mvLocation[1];
  mutable int colorLocation[1];
  mutable int lightLocation[5];
  mutable int twosideLocation;

   std::vector<float> *v_box;
   std::vector<float> color;


  mutable QOpenGLBuffer m_vbo[vboSize];
  mutable QOpenGLBuffer *m_ibo;
  mutable QOpenGLVertexArrayObject vao[vaoSize];
  mutable QOpenGLShaderProgram rendering_program;
  void draw_bbox();
  void attrib_buffers(CGAL::Three::Viewer_interface*) const;
  void compile_shaders();
  void draw_Bbox(Bbox bbox, std::vector<float> *vertices);
public Q_SLOTS:
    void changed();
#endif // SCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE
};

#endif // SCENE_SEGMENTED_IMAGE_ITEM_H
