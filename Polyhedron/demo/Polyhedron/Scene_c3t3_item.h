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
#include <CGAL/Three/Scene_group_item.h>
#include <Scene_polyhedron_item.h>
#include <Scene_polygon_soup_item.h>
#include <CGAL/IO/File_binary_mesh_3.h>

struct Scene_c3t3_item_priv;
class Scene_spheres_item;
class Scene_intersection_item;
using namespace CGAL::Three;
  class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item
  : public Scene_group_item
{
  Q_OBJECT
public:
  typedef qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item();
  Scene_c3t3_item(const C3t3& c3t3);
  ~Scene_c3t3_item();
  void setColor(QColor c);
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
      || (    c3t3().number_of_vertices_in_complex() == 0
           && c3t3().number_of_facets_in_complex()   == 0
           && c3t3().number_of_cells_in_complex()    == 0  );
  }


  void compute_bbox() const;
  Scene_item::Bbox bbox() const
  {
      return Scene_item::bbox();
  }
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
    return (m != Gouraud && m != PointsPlusNormals && m != Splatting && m != Points);
  }

  void draw(CGAL::Three::Viewer_interface* viewer) const;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
  void drawPoints(CGAL::Three::Viewer_interface * viewer) const;
  public:
    QMenu* contextMenu();

  protected:
    friend struct Scene_c3t3_item_priv;

    Scene_c3t3_item_priv* d;

};

#endif // SCENE_C3T3_ITEM_H
