#ifndef CGAL_SCENE_SURFACE_MESH_ITEM_H
#define CGAL_SCENE_SURFACE_MESH_ITEM_H
//Defines the precision of the positions (for performance/precision sake)
#define CGAL_GL_DATA GL_FLOAT
#define cgal_gl_data float
#define CGAL_IS_FLOAT 1

#include "Scene_surface_mesh_item_config.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <QColor>

#include <CGAL/Polygon_mesh_processing/properties.h>

struct Scene_surface_mesh_item_priv;

class SCENE_SURFACE_MESH_ITEM_EXPORT Scene_surface_mesh_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef SMesh FaceGraph;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  typedef SMesh::Property_map<vertex_descriptor,int> Vertex_selection_map;
  typedef SMesh::Property_map<face_descriptor,int> Face_selection_map;
  typedef SMesh::Property_map<halfedge_descriptor, bool> Halfedge_is_feature_map;
  typedef SMesh::Property_map<face_descriptor, int> Face_patch_id_map;
  typedef SMesh::Property_map<vertex_descriptor,int> Vertex_num_feature_edges_map;


  Scene_surface_mesh_item();
  // Takes ownership of the argument.
  Scene_surface_mesh_item(SMesh*);
  Scene_surface_mesh_item(const Scene_surface_mesh_item& other);

  ~Scene_surface_mesh_item();

  // Only needed for Scene_polyhedron_item
  void setItemIsMulticolor(bool){} 
  void update_vertex_indices(){}
  void update_halfedge_indices(){}
  void update_facet_indices(){}

  Scene_surface_mesh_item* clone() const;

  Vertex_selection_map vertex_selection_map();
  Face_selection_map face_selection_map();

  std::vector<QColor>& color_vector();
  void set_patch_id(SMesh::Face_index f,int i)const;
  int patch_id(SMesh::Face_index f)const;
  void draw() const {}
  void draw(CGAL::Three::Viewer_interface *) const;

  virtual void drawEdges() const {}
  void drawEdges(CGAL::Three::Viewer_interface *) const;
  void drawPoints(CGAL::Three::Viewer_interface *) const;

  bool supportsRenderingMode(RenderingMode m) const;
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;
  QString toolTip() const;

  SMesh* polyhedron();
  const SMesh* polyhedron() const;
  void compute_bbox()const;
  void invalidate_aabb_tree();
  void invalidateOpenGLBuffers();

Q_SIGNALS:
  void item_is_about_to_be_changed();
  void selection_done();
  void selected_vertex(void*);
  void selected_facet(void*);
  void selected_edge(void*);
  void selected_halfedge(void*);

public Q_SLOTS:
  virtual void selection_changed(bool);
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);
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

typedef Scene_surface_mesh_item::SMesh SMesh;
namespace CGAL{
SMesh::Property_map<boost::graph_traits<SMesh>::halfedge_descriptor,bool>
inline get(CGAL::halfedge_is_feature_t, SMesh& smesh)
{
  typedef boost::graph_traits<SMesh >::halfedge_descriptor halfedge_descriptor;
  return smesh.add_property_map<halfedge_descriptor,bool>("h:is_feature").first;
}


 SMesh::Property_map< boost::graph_traits<SMesh >::face_descriptor,int>
inline get(CGAL::face_patch_id_t, SMesh& smesh)
{
  typedef boost::graph_traits<SMesh >::face_descriptor face_descriptor;
  return smesh.add_property_map<face_descriptor,int>("f:patch_id").first;
}



 SMesh::Property_map< boost::graph_traits<SMesh >::face_descriptor,int>
 inline get(CGAL::face_selection_t, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::face_descriptor face_descriptor;
  return smesh.add_property_map<face_descriptor,int>("f:selection").first;
}



 SMesh::Property_map< boost::graph_traits<SMesh >::vertex_descriptor,int>
inline get(CGAL::vertex_selection_t, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,int>("v:selection").first;
}


 SMesh::Property_map< boost::graph_traits<SMesh >::vertex_descriptor,int>
 inline get(CGAL::vertex_num_feature_edges_t, SMesh& smesh)
{
  typedef  boost::graph_traits<SMesh >::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,int>("v:nfe").first;
}
}
namespace boost
{
template<>
struct property_map<SMesh, CGAL::halfedge_is_feature_t>
{
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;

  typedef SMesh::Property_map<halfedge_descriptor, bool> type;
};

template<>
struct property_map<SMesh, CGAL::face_patch_id_t>
{
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

  typedef SMesh::Property_map<face_descriptor, int> type;
};

template<>
struct property_map<SMesh, CGAL::vertex_selection_t>
{

  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  typedef SMesh::Property_map<vertex_descriptor, int> type;
};

template<>
struct property_map<SMesh, CGAL::face_selection_t>
{

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

  typedef SMesh::Property_map<face_descriptor, int> type;
};

template<>
struct property_map<SMesh, CGAL::vertex_num_feature_edges_t>
{

  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  typedef SMesh::Property_map<vertex_descriptor, int> type;
};
} //boost
#endif /* CGAL_SCENE_SURFACE_MESH_ITEM_H */
