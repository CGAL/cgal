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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <QColor>

#include "properties.h"

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
  void selected_vertex(std::size_t);
  void selected_facet(std::size_t);
  void selected_edge(std::size_t);
  void selected_halfedge(std::size_t);

public Q_SLOTS:
  virtual void selection_changed(bool);
  void select(double orig_x,
              double orig_y,
              double orig_z,
              double dir_x,
              double dir_y,
              double dir_z);
protected:
  friend struct Scene_surface_mesh_item_priv;
  Scene_surface_mesh_item_priv* d;
};
#if 0
Scene_surface_mesh_item::Halfedge_is_feature_map
get(halfedge_is_feature_t, Scene_surface_mesh_item::SMesh& smesh);

Scene_surface_mesh_item::Face_patch_id_map
get(face_patch_id_t, Scene_surface_mesh_item::SMesh& smesh);

Scene_surface_mesh_item::Vertex_selection_map
get(vertex_selection_t, Scene_surface_mesh_item::SMesh& smesh);

Scene_surface_mesh_item::Face_selection_map
get(face_selection_t, Scene_surface_mesh_item::SMesh& smesh);

Scene_surface_mesh_item::Vertex_num_feature_edges_map
get(vertex_num_feature_edges_t, Scene_surface_mesh_item::SMesh& smesh);
#else
Scene_surface_mesh_item::Halfedge_is_feature_map
inline get(halfedge_is_feature_t, Scene_surface_mesh_item::SMesh& smesh)
{
  typedef boost::graph_traits<Scene_surface_mesh_item::SMesh>::halfedge_descriptor halfedge_descriptor;
  return smesh.add_property_map<halfedge_descriptor,bool>("h:is_feature").first;
}


Scene_surface_mesh_item::Face_patch_id_map
inline get(face_patch_id_t, Scene_surface_mesh_item::SMesh& smesh)
{
  typedef  boost::graph_traits<Scene_surface_mesh_item::SMesh>::face_descriptor face_descriptor;
  return smesh.add_property_map<face_descriptor,int>("f:patch_id").first;
}


Scene_surface_mesh_item::Face_selection_map
inline get(face_selection_t, Scene_surface_mesh_item::SMesh& smesh)
{
  typedef  boost::graph_traits<Scene_surface_mesh_item::SMesh>::face_descriptor face_descriptor;
  return smesh.add_property_map<face_descriptor,int>("f:selection").first;
}


Scene_surface_mesh_item::Vertex_selection_map
inline get(vertex_selection_t, Scene_surface_mesh_item::SMesh& smesh)
{
  typedef  boost::graph_traits<Scene_surface_mesh_item::SMesh>::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,int>("v:selection").first;
}

Scene_surface_mesh_item::Vertex_num_feature_edges_map
inline get(vertex_num_feature_edges_t, Scene_surface_mesh_item::SMesh& smesh)
{
  typedef  boost::graph_traits<Scene_surface_mesh_item::SMesh>::vertex_descriptor vertex_descriptor;
  return smesh.add_property_map<vertex_descriptor,int>("v:nfe").first;
}
#endif

namespace boost {
  
  template <>
  struct property_map<Scene_surface_mesh_item::SMesh, halfedge_is_feature_t>
  {
    typedef Scene_surface_mesh_item::SMesh SMesh;
    typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  
    typedef SMesh::Property_map<halfedge_descriptor, bool> type;
  };
  
  template <>
  struct property_map<Scene_surface_mesh_item::SMesh, face_patch_id_t>
  {
    typedef Scene_surface_mesh_item::SMesh SMesh;
    typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  
    typedef SMesh::Property_map<face_descriptor, int> type;
  }; 

  template <>
  struct property_map<Scene_surface_mesh_item::SMesh, vertex_selection_t>
  {
    typedef Scene_surface_mesh_item::SMesh SMesh;
    typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  
    typedef SMesh::Property_map<vertex_descriptor, int> type;
  };

  template <>
  struct property_map<Scene_surface_mesh_item::SMesh, face_selection_t>
  {
    typedef Scene_surface_mesh_item::SMesh SMesh;
    typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  
    typedef SMesh::Property_map<face_descriptor, int> type;
  };

  template <>
  struct property_map<Scene_surface_mesh_item::SMesh, vertex_num_feature_edges_t>
  {
    typedef Scene_surface_mesh_item::SMesh SMesh;
    typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
  
    typedef SMesh::Property_map<vertex_descriptor, int> type;
  };

}



#endif /* CGAL_SCENE_SURFACE_MESH_ITEM_H */
