#ifndef SCENE_POLYGON_SOUP_H
#define SCENE_POLYGON_SOUP_H

#include "Scene_polygon_soup_config.h"
#include <CGAL_demo/Scene_item_with_display_list.h>
#include <iostream>

struct Polygon_soup;

class SCENE_POLYGON_SOUP_EXPORT Scene_polygon_soup 
  : public Scene_item_with_display_list 
{
  Q_OBJECT
public:  
  Scene_polygon_soup();
  ~Scene_polygon_soup();

  Scene_polygon_soup* clone() const;
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return m != Gouraud; } // CHECK THIS!
  // OpenGL drawing in a display list
  void direct_draw() const;

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  void shuffle_orientations();
  bool orient();
  void inside_out();

  void setDisplayNonManifoldEdges(const bool);
  bool displayNonManifoldEdges() const;
private:
  Polygon_soup* soup;
  bool oriented;
}; // end class Scene_polygon_soup

#endif // SCENE_POLYGON_SOUP_H
