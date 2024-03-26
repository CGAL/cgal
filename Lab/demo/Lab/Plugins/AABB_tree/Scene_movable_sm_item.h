#ifndef Scene_movable_sm_item_H
#define Scene_movable_sm_item_H

#include "SMesh_type.h"

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

#include <QKeyEvent>


#if defined( scene_movable_sm_item_EXPORTS)
#  define SCENE_MOVABLE_SM_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_MOVABLE_SM_ITEM_EXPORT Q_DECL_IMPORT
#endif

using namespace CGAL::Three;
struct Scene_movable_sm_item_priv;
// This class represents a polyhedron in the OpenGL scene
class SCENE_MOVABLE_SM_ITEM_EXPORT Scene_movable_sm_item
    : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT

public:
  Scene_movable_sm_item(const CGAL::qglviewer::Vec& pos, SMesh *sm,
                        const QString name);
  Scene_item* clone() const{return nullptr;}
  QString toolTip() const;
  void draw(CGAL::Three::Viewer_interface *) const;
  void drawEdges(CGAL::Three::Viewer_interface *) const;
  void compute_bbox() const;
  Scene_item::Bbox bbox() const ;
  ~Scene_movable_sm_item();
  bool manipulatable() const { return true;}
  ManipulatedFrame* manipulatedFrame();
  const CGAL::qglviewer::Vec& center() const;
  virtual bool supportsRenderingMode(RenderingMode m) const { return m==Flat
        || m==Wireframe; }
  virtual void invalidateOpenGLBuffers();
  void setFMatrix(const GLdouble matrix[16]);
  void computeElements() const;
  bool isEmpty() const {return false;}
  SMesh *getFaceGraph();
protected:
  friend struct Scene_movable_sm_item_priv;
  Scene_movable_sm_item_priv* d;

Q_SIGNALS:
  void stop();
  void killed();
}; // end class Scene_movable_sm_item

#endif // Scene_movable_sm_item_H
