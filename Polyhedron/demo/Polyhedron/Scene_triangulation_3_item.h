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
#include <CGAL/Three/Scene_group_item.h>

#ifdef scene_triangulation_3_item_EXPORTS
#  define mesh_3_demo_scene_triangulation_3_item_EXPORTS 1
#endif

#ifdef mesh_3_demo_scene_triangulation_3_item_EXPORTS
#  define SCENE_TRIANGULATION_3_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_TRIANGULATION_3_ITEM_EXPORT Q_DECL_IMPORT
#endif

struct Scene_triangulation_3_item_priv;
class QSlider;
using namespace CGAL::Three;
class SCENE_TRIANGULATION_3_ITEM_EXPORT Scene_triangulation_3_item
    : public Scene_group_item
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
  void setAlpha(int alpha) Q_DECL_OVERRIDE;
  QSlider* alphaSlider();
  Geom_traits::Plane_3 plane(CGAL::qglviewer::Vec offset = CGAL::qglviewer::Vec(0,0,0)) const;
  void compute_bbox() const Q_DECL_OVERRIDE;
  Scene_item::Bbox bbox() const Q_DECL_OVERRIDE
  {
    //todo
    return Scene_item::bbox();
  }
  //Scene_triangulation_3_item* clone() const  Q_DECL_OVERRIDE;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const  Q_DECL_OVERRIDE{
    return (m != Gouraud && m != PointsPlusNormals && m != Points && m != ShadedPoints);
  }

  QMenu* contextMenu() Q_DECL_OVERRIDE;
  bool manipulatable() const  Q_DECL_OVERRIDE{
    return true;
  }
  void setColor(QColor c) Q_DECL_OVERRIDE;
  bool has_tets() const;
  bool has_grid() const;
  void resetCutPlane();
  void setNormal(float x, float y, float z) ;
  void setPosition(float x, float y, float z) ;
  bool eventFilter(QObject *, QEvent *) Q_DECL_OVERRIDE;
  ManipulatedFrame* manipulatedFrame() Q_DECL_OVERRIDE;

  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawPoints(CGAL::Three::Viewer_interface * viewer) const Q_DECL_OVERRIDE;
  QList<Scene_interface::Item_id> getChildrenForSelection() const Q_DECL_OVERRIDE { return QList<Scene_interface::Item_id>(); }
public Q_SLOTS:
  void itemAboutToBeDestroyed(Scene_item *) Q_DECL_OVERRIDE;
  void initializeBuffers(Viewer_interface *) const Q_DECL_OVERRIDE;
  void computeElements() const Q_DECL_OVERRIDE;
  void show_intersection(bool b);
  void show_grid(bool b);
  void newViewer(Viewer_interface *) Q_DECL_OVERRIDE;
  void updateCutPlane();
  void frame_changed();

protected:
  friend struct Scene_triangulation_3_item_priv;
  Scene_triangulation_3_item_priv* d;
};

#endif // SCENE_TRIANGULATION_3_ITEM_H
