#ifndef CGAL_SCENE_SURFACE_MESH_ITEM_H
#define CGAL_SCENE_SURFACE_MESH_ITEM_H
//Defines the precision of the positions (for performance/precision sake)
#define CGAL_GL_DATA GL_FLOAT
#define cgal_gl_data float
#define CGAL_IS_FLOAT 1

#include "Scene_surface_mesh_item_config.h"
#include "SMesh_type.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>
#include <QColor>

#include "properties.h"


struct Scene_surface_mesh_item_priv;

class SCENE_SURFACE_MESH_ITEM_EXPORT Scene_surface_mesh_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  typedef SMesh Face_graph;
  typedef SMesh::Property_map<vertex_descriptor,int> Vertex_selection_map;
  typedef SMesh::Property_map<face_descriptor,int> Face_selection_map;
  Scene_surface_mesh_item();
  // Takes ownership of the argument.
  Scene_surface_mesh_item(SMesh*);
  Scene_surface_mesh_item(SMesh);
  Scene_surface_mesh_item(const Scene_surface_mesh_item& other);

  ~Scene_surface_mesh_item();


  Scene_surface_mesh_item* clone() const Q_DECL_OVERRIDE;
  void draw(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;
  void drawEdges(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;
  void drawPoints(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;

  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE;
  bool isFinite() const Q_DECL_OVERRIDE { return true; }
  bool isEmpty() const Q_DECL_OVERRIDE;
  Bbox bbox() const Q_DECL_OVERRIDE;
  QString toolTip() const Q_DECL_OVERRIDE;

  // Only needed for Scene_polyhedron_item
  void setItemIsMulticolor(bool);
  void update_vertex_indices(){}
  void update_halfedge_indices(){}
  void update_facet_indices(){}
  Vertex_selection_map vertex_selection_map();
  Face_selection_map face_selection_map();

  std::vector<QColor>& color_vector();
  void show_feature_edges(bool);
  SMesh* polyhedron();
  const SMesh* polyhedron() const;

  Face_graph*       face_graph() { return polyhedron(); }
  const Face_graph* face_graph() const { return polyhedron(); }

  void invalidate_aabb_tree();
  void invalidateOpenGLBuffers()Q_DECL_OVERRIDE;


  void compute_bbox()const Q_DECL_OVERRIDE;
  void standard_constructor(SMesh *sm);

Q_SIGNALS:
  void item_is_about_to_be_changed();
  void selection_done();
  void selected_vertex(void*);
  void selected_facet(void*);
  void selected_edge(void*);
  void selected_halfedge(void*);


public Q_SLOTS:
  void itemAboutToBeDestroyed(Scene_item *) Q_DECL_OVERRIDE;
  virtual void selection_changed(bool) Q_DECL_OVERRIDE;
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z) Q_DECL_OVERRIDE;
  bool intersect_face(double orig_x,
                      double orig_y,
                      double orig_z,
                      double dir_x,
                      double dir_y,
                      double dir_z,
                      const face_descriptor &f);
protected:
  friend struct Scene_surface_mesh_item_priv;
  Scene_surface_mesh_item_priv* d;
};

#endif /* CGAL_SCENE_SURFACE_MESH_ITEM_H */
