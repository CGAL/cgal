#ifndef SCENE_TRIANGULATION_3_ITEM_H
#define SCENE_TRIANGULATION_3_ITEM_H

#include "Scene_triangulation_3_item_config.h"
#include "T3_type.h"

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

struct Scene_triangulation_3_item_priv;
using namespace CGAL::Three;
  class SCENE_TRIANGULATION_3_ITEM_EXPORT Scene_triangulation_3_item
  : public Scene_group_item, public Scene_item_with_properties
{
  Q_OBJECT
public:
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;
  static const int number_of_bitset = 4; // also defined in "shader_c3t3.frag" and "shader_c3t3_edges.frag"

  Scene_triangulation_3_item(bool display_elements = true);
  Scene_triangulation_3_item(const T3 t3, bool display_elements = true);
  ~Scene_triangulation_3_item();

  void common_constructor(bool);
  bool has_stats()const  override {return true;}
  QString computeStats(int type)  override;
  CGAL::Three::Scene_item::Header_data header() const override;


  void setColor(QColor c) override;
  void setUseSubdomainColors(bool use_subdomain_colors);

  void invalidateOpenGLBuffers() override;

  void triangulation_changed();

  void resetCutPlane();

  void set_valid(bool);

  virtual const T3& triangulation() const;
  virtual T3& triangulation();

  bool manipulatable() const  override{
    return true;
  }

  Scene_item::Bbox bbox() const override;

  bool has_spheres() const;
  bool has_grid() const;
  bool has_tets() const;

  float alpha() const override;
  void setAlpha(int alpha) override;
  QSlider* alphaSlider();
  ManipulatedFrame* manipulatedFrame() override;

  void setPosition(float x, float y, float z) ;

  void setNormal(float x, float y, float z) ;

  Geom_traits::Plane_3 plane(CGAL::qglviewer::Vec offset = CGAL::qglviewer::Vec(0,0,0)) const;

  bool isFinite() const override { return true; }
  bool isEmpty() const override {
    return triangulation().number_of_vertices() == 0
      || ( triangulation().number_of_finite_facets()   == 0
           && triangulation().number_of_finite_cells()    == 0  );
    return false;
  }


  void compute_bbox() const override;

  Scene_triangulation_3_item* clone() const  override;

  virtual bool load_binary(std::istream& is);

  // data item
  const Scene_item* data_item() const;
  void set_data_item(const Scene_item* data_item);

  QString toolTip() const override;

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const  override{
    return (m != Gouraud && m != PointsPlusNormals && m != Points && m != ShadedPoints);
  }

  void draw(CGAL::Three::Viewer_interface* viewer) const override;
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const override;
  void drawPoints(CGAL::Three::Viewer_interface * viewer) const override;
   //When selecting a t3 item, we don't want to select its children, so we can still apply Operations to it
  QList<Scene_interface::Item_id> getChildrenForSelection() const override { return QList<Scene_interface::Item_id>(); }
  public:
    QMenu* contextMenu() override;
    void copyProperties(Scene_item *) override;
    float getShrinkFactor() const;
    bool keyPressEvent(QKeyEvent *) override;
    bool eventFilter(QObject *, QEvent *) override;
    const std::set<int> &subdomain_indices() const;
    QColor getSubdomainIndexColor(int i) const;
  public Q_SLOTS:

  void on_spheres_color_changed();

  void data_item_destroyed();

  void check_new_item(int id);
  void check_deleted_item(Scene_interface::Item_id id);

  void reset_spheres();

  void reset_intersection_item();
  void show_intersection(bool b);
  void show_grid(bool b);
  void show_spheres(bool b);
  void computeIntersection();
  void set_cut_edge(bool b);

  virtual QPixmap graphicalToolTip() const override;

  void update_histogram();

  void changed();

  void updateCutPlane();

  void build_histogram();

  QColor get_histogram_color(const double v) const;

  void resetVisibleSubdomain();
  void switchVisibleSubdomain(int);
  bool isVisibleSubdomain(int) const;

  void itemAboutToBeDestroyed(Scene_item *) override;

  void initializeBuffers(Viewer_interface *) const override;
  void computeElements() const override;
  void newViewer(Viewer_interface *viewer) override;

  protected:
    friend struct Scene_triangulation_3_item_priv;
    mutable Scene_triangulation_3_item_priv* d;
    enum Face_Containers{
      T3_faces = 0
    };
    enum Edge_Containers{
      T3_edges = 0,
      Grid_edges
    };

    virtual bool do_take_cell(const T3::Cell_handle&) const { return true; }
    virtual bool do_take_facet(const T3::Facet&)const { return true; }
    virtual bool do_take_vertex(const T3::Vertex_handle&)const { return true; }
    virtual bool is_facet_oriented(const T3::Facet&)const { return true; }
    virtual bool is_surface() const { return false; }
};

#endif // SCENE_TRIANGULATION_3_ITEM_H
