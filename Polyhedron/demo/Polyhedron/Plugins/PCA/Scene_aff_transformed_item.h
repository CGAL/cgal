#ifndef SCENE_AFF_TRANSFORMED_ITEM_H
#define SCENE_AFF_TRANSFORMED_ITEM_H

#if defined( scene_aff_transformed_item_EXPORTS)
#  define SCENE_AFF_TRANSFORMED_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_AFF_TRANSFORMED_ITEM_EXPORT Q_DECL_IMPORT
#endif

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Three.h>

#include <QKeyEvent>

using namespace CGAL::Three;

struct Scene_aff_transformed_item_priv
{
  bool manipulable;
  CGAL::Three::Scene_item::ManipulatedFrame* frame;
  QMatrix4x4 f_matrix;

  Scene_aff_transformed_item_priv(const CGAL::qglviewer::Vec& pos)
    : frame(new CGAL::Three::Scene_item::ManipulatedFrame())
  {
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    frame->setPosition(pos + offset);
  }

  ~Scene_aff_transformed_item_priv()
  {
    if(frame)
    {
      delete frame;
      frame = nullptr;
    }
  }
};

// This is an abstract representation of a scene item (points, polygon mesh, ...)
// that is transformed by an affine transformation.
// In Affine_transform_plugin.cpp, type erasure is used to factorize most of the code,
// using this base class.
class SCENE_AFF_TRANSFORMED_ITEM_EXPORT Scene_aff_transformed_item
  : public Scene_item_rendering_helper
{
  Q_OBJECT

protected:
  friend Scene_aff_transformed_item_priv;
  Scene_aff_transformed_item_priv* d;

public:
  Scene_aff_transformed_item(const CGAL::qglviewer::Vec& pos);

  ~Scene_aff_transformed_item();

  void itemAboutToBeDestroyed(Scene_item *item) Q_DECL_OVERRIDE;

  void setManipulatable(bool b = true) { d->manipulable = b;}
  bool manipulatable() const Q_DECL_OVERRIDE { return d->manipulable; }
  CGAL::Three::Scene_item::ManipulatedFrame* manipulatedFrame() Q_DECL_OVERRIDE { return d->frame; }
  void setFMatrix(double matrix[16]);
  const QMatrix4x4& getFMatrix() const { return d->f_matrix; }

  // below is defined in the specific aff_transformed items
  virtual void compute_bbox() const Q_DECL_OVERRIDE = 0;
  virtual const CGAL::qglviewer::Vec& center() const = 0;

  virtual bool keyPressEvent(QKeyEvent* e) Q_DECL_OVERRIDE;

Q_SIGNALS:
  void applyTransformation();
};

#endif // SCENE_AFF_TRANSFORMED_ITEM_H
