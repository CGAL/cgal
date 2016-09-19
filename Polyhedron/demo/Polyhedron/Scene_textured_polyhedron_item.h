#ifndef SCENE_TEXTURED_POLYHEDRON_ITEM_H
#define SCENE_TEXTURED_POLYHEDRON_ITEM_H
#include "Scene_textured_polyhedron_item_config.h"
#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include "Textured_polyhedron_type_fwd.h"
#include <iostream>
#include "texture.h"

struct Scene_textured_polyhedron_item_priv;
// This class represents a textured polyhedron in the OpenGL scene
class SCENE_TEXTURED_POLYHEDRON_ITEM_EXPORT Scene_textured_polyhedron_item 
  : public CGAL::Three::Scene_item {
  Q_OBJECT
public:  
  Scene_textured_polyhedron_item();
//   Scene_textured_polyhedron_item(const Scene_textured_polyhedron_item&);
  Scene_textured_polyhedron_item(const Textured_polyhedron& p);
  Scene_textured_polyhedron_item(Textured_polyhedron* const p);
  ~Scene_textured_polyhedron_item();

  Scene_textured_polyhedron_item* clone() const;
  
  // IO
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return (m != Splatting && m != PointsPlusNormals && m != Points && m != Gouraud && m != ShadedPoints); }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
   void draw() const {}
  virtual void draw(CGAL::Three::Viewer_interface*) const;
   virtual void drawEdges() const {}
   virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const;

  // Get wrapped textured_polyhedron
  Textured_polyhedron*       textured_polyhedron();
  const Textured_polyhedron* textured_polyhedron() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  void compute_bbox() const;
  virtual void invalidateOpenGLBuffers();
  virtual void selection_changed(bool);
  void add_border_edges(std::vector<float> border_edges);

Q_SIGNALS:
  void selectionChanged();

protected:
  friend struct Scene_textured_polyhedron_item_priv;
  Scene_textured_polyhedron_item_priv* d;

}; // end class Scene_textured_polyhedron_item

#endif // SCENE_TEXTURED_POLYHEDRON_ITEM_H
