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
#include "Scene_triangulation_3_item.h"

struct Scene_c3t3_item_priv;
using namespace CGAL::Three;
  class SCENE_C3T3_ITEM_EXPORT Scene_c3t3_item
  : public Scene_triangulation_3_item
  {
    Q_OBJECT
  public:
    typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;

    Scene_c3t3_item(bool is_surface = false);
    Scene_c3t3_item(const C3t3& c3t3, bool is_surface = false);
    ~Scene_c3t3_item();

    Scene_c3t3_item* clone() const  override;

    const C3t3& c3t3() const;
    C3t3& c3t3();

    void c3t3_changed();

    const T3& triangulation() const override;
    T3& triangulation() override;

  bool save_binary(std::ostream& os) const
  {
    return CGAL::IO::save_binary_file(os, c3t3());
  }

  bool save_ascii(std::ostream& os) const
  {
      os << "ascii CGAL c3t3 " << CGAL::Get_io_signature<C3t3>()() << "\n";
      CGAL::IO::set_ascii_mode(os);
      return !!(os << c3t3());
    }
    bool load_binary(std::istream& is) override;

    bool is_valid() const;//true if the c3t3 is correct, false if it was made from a .mesh, for example
    void set_valid(bool);
    QMenu* contextMenu() override;

    void drawEdges(Viewer_interface *viewer) const override;

    //stats
    QString computeStats(int type)  override;

    void copyProperties(Scene_item *) override;

    CGAL::Three::Scene_item::Header_data header() const override;
    bool has_cnc() const;
  public Q_SLOTS:
    void show_cnc(bool);
    void export_facets_in_complex();
    void initializeBuffers(Viewer_interface *) const override;
    void computeElements() const override;
    void set_sharp_edges_angle(double d);
    double get_sharp_edges_angle();

    void set_detect_borders(bool b);
    bool get_detect_borders();

  protected:
    friend struct Scene_c3t3_item_priv;
    mutable Scene_c3t3_item_priv* d;

    enum Edge_Containers{
      CNC = 2
    };
    bool do_take_cell(const T3::Cell_handle&) const override;
    bool do_take_facet(const T3::Facet&)const override;
    bool do_take_vertex(const T3::Vertex_handle &)const override;
    bool is_facet_oriented(const T3::Facet&)const override;
    bool is_surface()const override;
    void common_constructor(bool is_surface);
  };

#endif // SCENE_C3T3_ITEM_H
