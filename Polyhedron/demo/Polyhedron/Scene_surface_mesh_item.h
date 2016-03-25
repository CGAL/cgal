#ifndef CGAL_SCENE_SURFACE_MESH_ITEM_H
#define CGAL_SCENE_SURFACE_MESH_ITEM_H
//Defines the precision of the positions (for performance/precision sake)
#define GL_DATA GL_FLOAT
#define gl_data float

#include "Scene_surface_mesh_item_config.h"
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/array.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>

class SCENE_SURFACE_MESH_ITEM_EXPORT Scene_surface_mesh_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT
  typedef CGAL::Simple_cartesian<gl_data> Kernel;
  typedef Kernel::Point_3 Point;
public:
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;

  // Takes ownership of the argument.
  Scene_surface_mesh_item(SMesh*);

  Scene_surface_mesh_item(const Scene_surface_mesh_item& other);

  ~Scene_surface_mesh_item(){delete smesh_;}

  Scene_surface_mesh_item* clone() const;
  void draw(CGAL::Three::Viewer_interface *) const;
  void draw_edges(CGAL::Three::Viewer_interface *) const;
  void draw_points(CGAL::Three::Viewer_interface *) const;

  bool supportsRenderingMode(RenderingMode m) const;
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

  enum VAOs {
   Flat_facets = 0,
   Smooth_facets,
   Edges,
   NbOfVaos
  };
  enum VBOs {
    Flat_vertices = 0,
    Smooth_vertices,
    Flat_normals,
    Smooth_normals,
    Colors,
    NbOfVbos
  };

  QString toolTip() const;

  SMesh* polyhedron() { return smesh_; }
  const SMesh* polyhedron() const { return smesh_; }
private:
  SMesh* smesh_;
  void initializeBuffers(CGAL::Three::Viewer_interface *) const;
  std::vector<unsigned int> idx_data_;
  std::vector<unsigned int> idx_edge_data_;
  mutable std::vector<gl_data> flat_vertices;
  mutable std::vector<gl_data> flat_normals;
  mutable QOpenGLShaderProgram *program;

};


#endif /* CGAL_SCENE_SURFACE_MESH_ITEM_H */
