#ifndef CGAL_SCENE_SURFACE_MESH_ITEM_H
#define CGAL_SCENE_SURFACE_MESH_ITEM_H
//Defines the precision of the positions (for performance/precision sake)
#define CGAL_GL_DATA GL_FLOAT
#define cgal_gl_data float
#define CGAL_IS_FLOAT 1

#include "Scene_surface_mesh_item_config.h"
#include <CGAL/Three/Scene_zoomable_item_interface.h>
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
  : public CGAL::Three::Scene_item,
    public CGAL::Three::Scene_zoomable_item_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Scene_zoomable_item_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.ZoomInterface/1.0")
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
  bool isItemMulticolor();
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
  bool save(std::ostream& out) const;
  bool save_obj(std::ostream& out) const;
  bool load_obj(std::istream& in);

  enum STATS {
    NB_VERTICES = 0,
    NB_CONNECTED_COMPOS,
    NB_BORDER_EDGES,
    IS_PURE_TRIANGLE,
    NB_DEGENERATED_FACES,
    HOLES,
    AREA,
    VOLUME,
    SELFINTER,
    NB_FACETS,
    MIN_AREA,
    MAX_AREA,
    MED_AREA,
    MEAN_AREA,
    MIN_ALTITUDE,
    MIN_ASPECT_RATIO,
    MAX_ASPECT_RATIO,
    MEAN_ASPECT_RATIO,
    GENUS,
    NB_EDGES,
    MIN_LENGTH,
    MAX_LENGTH,
    MID_LENGTH,
    MEAN_LENGTH,
    NB_NULL_LENGTH,
    MIN_ANGLE,
    MAX_ANGLE,
    MEAN_ANGLE
  };

  bool has_stats()const Q_DECL_OVERRIDE{return true;}
  QString computeStats(int type)Q_DECL_OVERRIDE;
  CGAL::Three::Scene_item::Header_data header() const Q_DECL_OVERRIDE;
  void zoomToPosition(const QPoint &point, CGAL::Three::Viewer_interface *)const Q_DECL_OVERRIDE;

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
