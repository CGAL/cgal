#ifndef SCENE_C3T3_ITEM_H
#define SCENE_C3T3_ITEM_H

#include "Scene_c3t3_item_config.h"
#include "C3t3_type.h"
#include "Scene_item_with_display_list.h"

#include <QVector>
#include <QColor>
#include <map>

#include <Qt/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

struct Scene_c3t3_item_priv;

class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item
  : public Scene_item_with_display_list
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item(const C3t3& c3t3);
  ~Scene_c3t3_item();

  const C3t3& c3t3() const;

  bool manipulatable() const {
    return true;
  }

  ManipulatedFrame* manipulatedFrame() {
    return frame;
  }

  void setPosition(float x, float y, float z) {
    frame->setPosition(x, y, z);
  }

  void setNormal(float x, float y, float z) {
    frame->setOrientation(x, y, z, 0.f);
  }

  Kernel::Plane_3 plane() const;

  bool isFinite() const { return true; }
  bool isEmpty() const {
    return c3t3().triangulation().number_of_vertices() == 0;
  }

  Bbox bbox() const;

  Scene_c3t3_item* clone() const {
    return 0;
  }

  QString toolTip() const;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud); // CHECK THIS!
  }

  void direct_draw() const;
  void direct_draw_edges() const;

protected:
  void direct_draw(int) const;
  Scene_c3t3_item_priv* d;

  qglviewer::ManipulatedFrame* frame;
};

#endif // SCENE_C3T3_ITEM_H
