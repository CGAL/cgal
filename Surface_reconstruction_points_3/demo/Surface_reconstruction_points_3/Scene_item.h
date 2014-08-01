#ifndef SCENE_ITEM_H
#define SCENE_ITEM_H

#include "Scene_item_config.h"
#include "Scene_interface.h"
#include <QString>
#include <QFont>

namespace qglviewer {
  class ManipulatedFrame;
}

// This class represents an object in the OpenGL scene
class SCENE_ITEM_EXPORT Scene_item : public QObject {
  Q_OBJECT
  Q_PROPERTY(QColor color READ color WRITE setColor)
  Q_PROPERTY(QString name READ name WRITE setName)
  Q_PROPERTY(bool visible READ visible WRITE setVisible)
  Q_ENUMS(RenderingMode)
  Q_PROPERTY(RenderingMode renderingMode READ renderingMode WRITE setRenderingMode)
public:
  typedef Scene_interface::Bbox Bbox;
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  static const QColor defaultColor; // defined in Scene_item.cpp

  static QString renderingModeNameList(); // Rendering mode list as a human readable string

  Scene_item()
    : name_("unamed"),
      color_(defaultColor),
      visible_(true),
      rendering_mode(FlatPlusEdges)
  {};
  virtual ~Scene_item();
  virtual Scene_item* clone() const = 0;

  // Is selection empty?
  virtual bool isSelectionEmpty() const { return true; }
  // Delete selection
  virtual void deleteSelection() {}
  // Reset selection mark
  virtual void resetSelection() {}

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const = 0;
  // Flat/Gouraud OpenGL drawing
  virtual void draw() const = 0;
  // Wireframe OpenGL drawing
  virtual void draw_edges() const { draw(); }
  // Points OpenGL drawing
  virtual void draw_points() const { draw(); }
  // Normals OpenGL drawing
  virtual void draw_normals() const {}
  // Draws oriented points with radius
  virtual void draw_splats() const {}

  // Functions for displaying meta-data of the item
  virtual QString toolTip() const = 0;
  virtual QFont font() const { return QFont(); }

  // Functions that help the Scene to compute its bbox
  virtual bool isFinite() const { return true; }
  virtual bool isEmpty() const { return true; }
  virtual Bbox bbox() const { return Bbox(); }

  // Function about manipulation
  virtual bool manipulatable() const { return false; }
  virtual ManipulatedFrame* manipulatedFrame() { return 0; }

  // Getters for the four basic properties
  virtual QColor color() const { return color_; }
  virtual QString name() const { return name_; }
  virtual bool visible() const { return visible_; }
  virtual RenderingMode renderingMode() const { return rendering_mode; }
  virtual QString renderingModeName() const; // Rendering mode as a human readable string

public slots:
  // Call that once you have finished changing something in the item
  // (either the properties or internal data)
  virtual void changed() {}

  // Setters for the four basic properties
  virtual void setColor(QColor c) { color_ = c; }
  virtual void setName(QString n) { name_ = n; }
  virtual void setVisible(bool b) { visible_ = b; }
  virtual void setRenderingMode(RenderingMode m) {
    if (supportsRenderingMode(m))
      rendering_mode = m;
  }

  virtual void itemAboutToBeDestroyed(Scene_item*);

signals:
  void aboutToBeDestroyed();

protected:
  // The four basic properties
  QString name_;
  QColor color_;
  bool visible_;
  RenderingMode rendering_mode;

}; // end class Scene_item

#endif // SCENE_ITEM_H
