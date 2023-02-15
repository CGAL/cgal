#ifndef SCENE_TETRAHEDRA_ITEM_H
#define SCENE_TETRAHEDRA_ITEM_H

#ifndef SCENE_TETRAHEDRA_ITEM_CONFIG_H
#define SCENE_TETRAHEDRA_ITEM_CONFIG_H

#ifdef scene_tetrahedra_item_EXPORTS
#  define mesh_3_demo_scene_tetrahedra_item_EXPORTS 1
#endif

#ifdef mesh_3_demo_scene_tetrahedra_item_EXPORTS
#  define SCENE_TETRAHEDRA_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_TETRAHEDRA_ITEM_EXPORT Q_DECL_IMPORT
#endif

#endif // SCENE_TETRAHEDRA_ITEM_CONFIG_H


#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include "Scene_c3t3_item.h"

struct tet_item_priv;
class QLabel;
class DoubleEdit;

class SCENE_TETRAHEDRA_ITEM_EXPORT Scene_tetrahedra_item : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public :
  Scene_tetrahedra_item(Scene_c3t3_item*);
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE {
    return (m == Flat);
  }
  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void invalidateOpenGLBuffers()Q_DECL_OVERRIDE;
  void computeElements() const Q_DECL_OVERRIDE;
  Scene_item* clone() const Q_DECL_OVERRIDE;
  QString toolTip() const Q_DECL_OVERRIDE{return QString();}
  bool isEmpty() const Q_DECL_OVERRIDE{ return false;}
  bool isFinite() const Q_DECL_OVERRIDE{ return true;}
  void initializeBuffers(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;
  void compute_bbox() const Q_DECL_OVERRIDE;
  void setMinMinLabelPointer(QLabel*);
  void setMinMaxLabelPointer(QLabel*);
  void setMaxMinLabelPointer(QLabel*);
  void setMaxMaxLabelPointer(QLabel*);
  void setMinEditPointer(DoubleEdit* );
  void setMaxEditPointer(DoubleEdit* );
  Scene_c3t3_item* c3t3_item();
  public Q_SLOTS:
  void setMinThreshold(int);
  void setMinThreshold(void);
  void setMaxThreshold(int);
  void setMaxThreshold(void);
  void setFilter(int);

private:
  friend struct tet_item_priv;
  tet_item_priv* d;
  void updateFilter() const;
  void updateThresholds(bool update = true);
}; //end of class Scene_tetrahedra_item

#endif // SCENE_TETRAHEDRA_ITEM_H
