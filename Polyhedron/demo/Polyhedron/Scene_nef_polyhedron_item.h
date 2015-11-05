#ifndef SCENE_NEF_POLYHEDRON_ITEM_H
#define SCENE_NEF_POLYHEDRON_ITEM_H
#include "Scene_item.h"
#include "Scene_nef_polyhedron_item_config.h"
#include "Nef_type_fwd.h"
#include <iostream>
#include <queue>
class Scene_polyhedron_item;

class SCENE_NEF_POLYHEDRON_ITEM_EXPORT Scene_nef_polyhedron_item
 : public Scene_item
{
  Q_OBJECT
public:
  Scene_nef_polyhedron_item();
//   Scene_nef_polyhedron_item(const Scene_nef_polyhedron_item&);
  Scene_nef_polyhedron_item(const Nef_polyhedron& p);
  Scene_nef_polyhedron_item(Nef_polyhedron* const p);
  ~Scene_nef_polyhedron_item();

  Scene_nef_polyhedron_item* clone() const;
  bool load_from_off(std::istream& in);
  bool load(std::istream& in);
  bool save(std::ostream& in) const;

  QFont font() const;
  QString toolTip() const;

  virtual void invalidate_buffers();
  virtual void selection_changed(bool);
  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return m != Gouraud && m!=Splatting; } // CHECK THIS!
  // OpenGL drawing in a display list
  void direct_draw() const;

  virtual void draw(Viewer_interface*) const;
  virtual void draw_edges() const {}
  virtual void draw_edges(Viewer_interface* viewer) const;
  virtual void draw_points(Viewer_interface*) const;
  // Wireframe OpenGL drawing

  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  Nef_polyhedron* nef_polyhedron();
  const Nef_polyhedron* nef_polyhedron() const;

  bool is_simple() const;
  bool is_Triangle;
  // conversion operations
  static Scene_nef_polyhedron_item* from_polyhedron(Scene_polyhedron_item*);
  Scene_polyhedron_item* convert_to_polyhedron() const;

  // Nef boolean operations
  Scene_nef_polyhedron_item&
  operator+=(const Scene_nef_polyhedron_item&); // union

  Scene_nef_polyhedron_item&
  operator*=(const Scene_nef_polyhedron_item&); // intersection

  Scene_nef_polyhedron_item&
  operator-=(const Scene_nef_polyhedron_item&); // difference

  static Scene_nef_polyhedron_item*
  sum(const Scene_nef_polyhedron_item&, 
      const Scene_nef_polyhedron_item&);

  void convex_decomposition(std::list< Scene_polyhedron_item*>&);
  
private:
  typedef Scene_item Base;
  typedef std::vector<QColor> Color_vector;

  Nef_polyhedron* nef_poly;


  mutable std::vector<double> positions_lines;
  mutable std::vector<double> positions_facets;
  mutable std::vector<double> positions_points;
  mutable std::vector<double> normals;
  mutable std::vector<double> color_lines;
  mutable std::vector<double> color_facets;
  mutable std::vector<double> color_points;
  mutable std::size_t nb_points;
  mutable std::size_t nb_lines;
  mutable std::size_t nb_facets;

  mutable QOpenGLShaderProgram *program;

  using Scene_item::initialize_buffers;
  void initialize_buffers(Viewer_interface *viewer) const;
  void compute_normals_and_vertices(void);

  void triangulate_facet();
  void triangulate_facet_color();
}; // end class Scene_nef_polyhedron_item

#endif // SCENE_NEF_POLYHEDRON_ITEM_H
