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
#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_group_item.h>
#include <Scene_polygon_soup_item.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/Three/Scene_item_with_properties.h>

struct Scene_c3t3_item_priv;
class Scene_spheres_item;
class Scene_intersection_item;
using namespace CGAL::Three;
  class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item
  : public Scene_group_item, public Scene_item_with_properties
{
  Q_OBJECT
public:
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

  Scene_c3t3_item(bool is_surface = false);
  Scene_c3t3_item(const C3t3& c3t3, bool is_surface = false);
  ~Scene_c3t3_item();

  void common_constructor(bool is_surface);
  bool has_stats()const  Q_DECL_OVERRIDE {return true;}
  QString computeStats(int type)  Q_DECL_OVERRIDE;
  CGAL::Three::Scene_item::Header_data header() const Q_DECL_OVERRIDE;


  void setColor(QColor c) Q_DECL_OVERRIDE;
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

  void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;

  void c3t3_changed();

  void resetCutPlane();

  void set_valid(bool);

  const C3t3& c3t3() const;
  C3t3& c3t3();

  bool manipulatable() const  Q_DECL_OVERRIDE{
    return true;
  }

  bool has_spheres() const;
  bool has_grid() const;
  bool has_cnc() const;
  bool has_tets() const;
  bool is_valid() const;//true if the c3t3 is correct, false if it was made from a .mesh, for example
  float alpha() const Q_DECL_OVERRIDE;
  void setAlpha(int alpha) Q_DECL_OVERRIDE;
  QSlider* alphaSlider();
  ManipulatedFrame* manipulatedFrame() Q_DECL_OVERRIDE;

  void setPosition(float x, float y, float z) ;

  void setNormal(float x, float y, float z) ;

  Geom_traits::Plane_3 plane(CGAL::qglviewer::Vec offset = CGAL::qglviewer::Vec(0,0,0)) const;

  bool isFinite() const Q_DECL_OVERRIDE { return true; }
  bool isEmpty() const Q_DECL_OVERRIDE {
    return c3t3().triangulation().number_of_vertices() == 0
      || (    c3t3().number_of_vertices_in_complex() == 0
           && c3t3().number_of_facets_in_complex()   == 0
           && c3t3().number_of_cells_in_complex()    == 0  );
  }


  void compute_bbox() const Q_DECL_OVERRIDE;
  Scene_item::Bbox bbox() const Q_DECL_OVERRIDE
  {
      return Scene_item::bbox();
  }

  Scene_c3t3_item* clone() const  Q_DECL_OVERRIDE;

  bool load_binary(std::istream& is);

  // data item
  const Scene_item* data_item() const;
  void set_data_item(const Scene_item* data_item);

  QString toolTip() const Q_DECL_OVERRIDE;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const  Q_DECL_OVERRIDE{
    return (m != Gouraud && m != PointsPlusNormals && m != Points && m != ShadedPoints);
  }

  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
  void drawPoints(CGAL::Three::Viewer_interface * viewer) const Q_DECL_OVERRIDE;
   //When selecting a c3t3 item, we don't want to select its children, so we can still apply Operations to it
  QList<Scene_interface::Item_id> getChildrenForSelection() const Q_DECL_OVERRIDE { return QList<Scene_interface::Item_id>(); }
  public:
    QMenu* contextMenu() Q_DECL_OVERRIDE;
    void copyProperties(Scene_item *) Q_DECL_OVERRIDE;
    float getShrinkFactor() const;
    bool keyPressEvent(QKeyEvent *) Q_DECL_OVERRIDE;
    bool eventFilter(QObject *, QEvent *) Q_DECL_OVERRIDE;
  public Q_SLOTS:

  void on_spheres_color_changed();
  void export_facets_in_complex();

  void data_item_destroyed();

  void reset_spheres();

  void reset_intersection_item();
  void show_spheres(bool b);
  void show_intersection(bool b);
  void show_grid(bool b);
  void show_cnc(bool b);

  virtual QPixmap graphicalToolTip() const Q_DECL_OVERRIDE;

  void update_histogram();

  void changed();

  void updateCutPlane();

  void build_histogram();

  QColor get_histogram_color(const double v) const;

  void set_sharp_edges_angle(double d);
  double get_sharp_edges_angle();

  void set_detect_borders(bool b);
  bool get_detect_borders();

  void itemAboutToBeDestroyed(Scene_item *) Q_DECL_OVERRIDE;

  void initializeBuffers(Viewer_interface *) const Q_DECL_OVERRIDE;
  void computeElements() const Q_DECL_OVERRIDE;
  void newViewer(Viewer_interface *viewer) Q_DECL_OVERRIDE;

  protected:
    friend struct Scene_c3t3_item_priv;
    Scene_c3t3_item_priv* d;
    enum Face_Containers{
      C3t3_faces = 0
    };
    enum Edge_Containers{
      C3t3_edges = 0,
      Grid_edges,
      CNC
    };
    enum Point_Container{
      C3t3_points = 0
    };

};

#endif // SCENE_C3T3_ITEM_H
