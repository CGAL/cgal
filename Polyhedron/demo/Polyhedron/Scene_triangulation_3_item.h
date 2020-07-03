#ifndef SCENE_TRIANGULATION_3_ITEM_H
#define SCENE_TRIANGULATION_3_ITEM_H

#include "T3_type.h"

#include <QVector>
#include <QColor>
#include <QMenu>

#include <QtCore/qglobal.h>
#include <CGAL/Qt/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>

#ifdef scene_triangulation_3_item_EXPORTS
#  define mesh_3_demo_scene_triangulation_3_item_EXPORTS 1
#endif

#ifdef mesh_3_demo_scene_triangulation_3_item_EXPORTS
#  define SCENE_TRIANGULATION_3_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_TRIANGULATION_3_ITEM_EXPORT Q_DECL_IMPORT
#endif

struct Scene_triangulation_3_item_priv;
using namespace CGAL::Three;
class SCENE_TRIANGULATION_3_ITEM_EXPORT Scene_triangulation_3_item
    : public Scene_item_rendering_helper
{
  Q_OBJECT
public:
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_triangulation_3_item();
  Scene_triangulation_3_item(T3* t3);
  ~Scene_triangulation_3_item();

  void common_constructor();
  void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;
  void changed();
  const T3& triangulation() const;
  T3& triangulation();
  bool isFinite() const Q_DECL_OVERRIDE { return true; }
  bool isEmpty() const Q_DECL_OVERRIDE {
    //todo
    return false;
  }

  void compute_bbox() const Q_DECL_OVERRIDE;
  Scene_item::Bbox bbox() const Q_DECL_OVERRIDE
  {
    //todo
    return Scene_item::bbox();
  }
  void frame_changed();
  Scene_triangulation_3_item* clone() const  Q_DECL_OVERRIDE;

  QString toolTip() const Q_DECL_OVERRIDE;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const  Q_DECL_OVERRIDE{
    return (m != Gouraud && m != PointsPlusNormals && m != Points && m != ShadedPoints);
  }

  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawPoints(CGAL::Three::Viewer_interface * viewer) const Q_DECL_OVERRIDE;

public Q_SLOTS:
  void initializeBuffers(Viewer_interface *) const Q_DECL_OVERRIDE;
  void computeElements() const Q_DECL_OVERRIDE;
  void reset_intersection_item();
  void show_intersection(bool b);

protected:
  friend struct Scene_triangulation_3_item_priv;
  Scene_triangulation_3_item_priv* d;
};

#endif // SCENE_TRIANGULATION_3_ITEM_H
