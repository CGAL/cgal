#ifndef SCENE_SEGMENTED_IMAGE_ITEM_H
#define SCENE_SEGMENTED_IMAGE_ITEM_H

#include <CGAL_demo/Scene_item.h>
#include <CGAL_demo/Scene_interface.h>
#include "Image_type_fwd.h"
#include "Scene_segmented_image_item_config.h"
#include <CGAL/gl.h>

typedef CGAL::Image_3 Image;

class SCENE_SEGMENTED_IMAGE_ITEM_EXPORT Scene_segmented_image_item 
  : public Scene_item
{
  Q_OBJECT
public:

  Scene_segmented_image_item(Image* im, int drawing_scale);
  ~Scene_segmented_image_item();

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  Bbox bbox() const;

  Scene_segmented_image_item* clone() const { return NULL; }

  // rendering mode
  virtual bool supportsRenderingMode(RenderingMode m) const;

  // draw
  virtual void direct_draw() const { draw(); }
  virtual void direct_draw_edges() const { draw_edges(); }
  virtual void draw() const;
  virtual void draw_edges() const { draw_gl_edges(); }
  
  virtual QString toolTip() const;
  
  const Image* image() const { return m_image; }

private:
  void draw_gl() const;
  void draw_gl_edges() const;
  
  void initialize_buffers();
  GLint ibo_size() const;
  
public:
  Image* m_image;

private:
  bool m_initialized;
  bool m_draw_edges;
  int m_voxel_scale;
  GLuint m_vbo[3];
  GLuint m_ibo;
};

#endif // SCENE_SEGMENTED_IMAGE_ITEM_H
