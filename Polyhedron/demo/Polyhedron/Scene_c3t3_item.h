#ifndef SCENE_C3T3_ITEM_H
#define SCENE_C3T3_ITEM_H

#include "Scene_c3t3_item_config.h"
#include "C3t3_type.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QMenu>
#include <set>

#include <QtCore/qglobal.h>
#include <CGAL/gl.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Three/Viewer_interface.h>

#include <CGAL/Three/Scene_item.h>
#include <Scene_polyhedron_item.h>
#include <Scene_polygon_soup_item.h>
#include <CGAL/IO/File_binary_mesh_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_C3T3_triangle_primitive.h>

struct Scene_c3t3_item_priv;

using namespace CGAL::Three;
  class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item
  : public Scene_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item();
  Scene_c3t3_item(const C3t3& c3t3);
  ~Scene_c3t3_item();

  bool save_binary(std::ostream& os) const
  {
    return CGAL::Mesh_3::save_binary_file(os, c3t3());
  }
  bool save_ascii(std::ostream& os) const
  {
      os << "ascii CGAL c3t3 " << CGAL::Get_io_signature<C3t3>()() << "\n";
      CGAL::set_ascii_mode(os);
      return !!(os << c3t3());
  }

  void invalidateOpenGLBuffers()
  {
    are_buffers_filled = false;
    compute_bbox();
  }

  void c3t3_changed();

  const C3t3& c3t3() const;
  C3t3& c3t3();

  bool manipulatable() const {
    return true;
  }
  ManipulatedFrame* manipulatedFrame() {
    return frame;
  }

  void setPosition(float x, float y, float z) {
    frame->setPosition(x, y, z);
  }

  void setNormal(float x, float y, float z) {
    frame->setOrientation(x, y, z, 0.f);
  }

  Kernel::Plane_3 plane() const;

  bool isFinite() const { return true; }
  bool isEmpty() const {
    return c3t3().triangulation().number_of_vertices() == 0
      || c3t3().number_of_cells() == 0;
  }

  void compute_bbox() const;

  Scene_c3t3_item* clone() const {
    return 0;
  }

  bool load_binary(std::istream& is);

  // data item
  const Scene_item* data_item() const;
  void set_data_item(const Scene_item* data_item);

  QString toolTip() const;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m != Gouraud && m != PointsPlusNormals && m != Splatting);
  }

  void draw(CGAL::Three::Viewer_interface* viewer) const;
  void draw_edges(CGAL::Three::Viewer_interface* viewer) const;
  void draw_points(CGAL::Three::Viewer_interface * viewer) const;
private:

  typedef CGAL::AABB_C3T3_triangle_primitive<Kernel,C3t3> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;
  typedef Tree::Point_and_primitive_id Point_and_primitive_id;

  Tree tree;

  bool need_changed;
  void reset_cut_plane();
  void draw_triangle(const Kernel::Point_3& pa,
    const Kernel::Point_3& pb,
    const Kernel::Point_3& pc, bool /* is_cut */) const;

  void draw_triangle_edges(const Kernel::Point_3& pa,
    const Kernel::Point_3& pb,
    const Kernel::Point_3& pc)const;

  double complex_diag() const;

  void compute_color_map(const QColor& c);

  public Q_SLOTS:
  void export_facets_in_complex();

  void data_item_destroyed();

  void show_spheres(bool b)
  {
    spheres_are_shown = b;
    Q_EMIT redraw();

  }
  virtual QPixmap graphicalToolTip() const;

  void update_histogram();

  void changed();

  void timerEvent(QTimerEvent*);

private:
  void build_histogram();

  QColor get_histogram_color(const double v) const;

public:
  QMenu* contextMenu();

  void set_scene(CGAL::Three::Scene_interface* scene){ last_known_scene = scene; }

protected:
  Scene_c3t3_item_priv* d;

private:

  
  
  mutable bool are_intersection_buffers_filled;
  enum Buffer
  {
      Facet_vertices =0,
      Facet_normals,
      Facet_colors,
      Edges_vertices,
      Grid_vertices,
      Sphere_vertices,
      Sphere_normals,
      Sphere_colors,
      Sphere_radius,
      Sphere_center,
      Wired_spheres_vertices,
      iEdges_vertices,
      iFacet_vertices,
      iFacet_normals,
      iFacet_colors,
      NumberOfBuffers
  };
  enum Vao
  {
      Facets=0,
      Edges,
      Grid,
      Spheres,
      Wired_spheres,
      iEdges,
      iFacets,
      NumberOfVaos
  };
  qglviewer::ManipulatedFrame* frame;
  CGAL::Three::Scene_interface* last_known_scene;

  bool spheres_are_shown;
  const Scene_item* data_item_;
  QPixmap histogram_;

  typedef std::set<int> Indices;
  Indices indices_;

  //!Allows OpenGL 2.1 context to get access to glDrawArraysInstanced.
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  //!Allows OpenGL 2.1 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  //!Allows OpenGL 2.1 context to get access to gkFrameBufferTexture2D.
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  //!Allows OpenGL 2.1 context to get access to glVertexAttribDivisor.
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;

  mutable std::size_t positions_poly_size;
  mutable std::size_t positions_lines_size;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_grid;
  mutable std::vector<float> positions_poly;
  mutable std::vector<float> normals;
  mutable std::vector<float> f_colors;
  mutable std::vector<float> s_normals;
  mutable std::vector<float> s_colors;
  mutable std::vector<float> s_vertex;
  mutable std::vector<float> ws_vertex;
  mutable std::vector<float> s_radius;
  mutable std::vector<float> s_center;
  mutable QOpenGLShaderProgram *program;
  mutable QOpenGLShaderProgram *program_sphere;

  using Scene_item::initialize_buffers;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer);
  void initialize_intersection_buffers(CGAL::Three::Viewer_interface *viewer);
  void compute_elements();
  void compute_intersections();
  void compute_intersection(const Primitive& facet);
  void compile_shaders();

struct Compute_intersection {
  Scene_c3t3_item& item;

  Compute_intersection(Scene_c3t3_item& item)
    : item(item)
  {}

  void operator()(const Primitive& facet) const
  {
    item.compute_intersection(facet);
  }
};


};

#endif // SCENE_C3T3_ITEM_H
