#include "config.h"
#include "Scene_spheres_item.h"
#include "Scene_c3t3_item.h"
#include "Scene_surface_mesh_item.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QApplication>
#include <QPainter>
#include <QtCore/qglobal.h>
#include <QGuiApplication>
#include <QSlider>
#include <QWidgetAction>
#include <QKeyEvent>
#include <QMouseEvent>

#include <map>
#include <vector>

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

#include <CGAL/Real_timer.h>

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangulation_3_cell_primitive.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include "Scene_polygon_soup_item.h"


typedef CGAL::AABB_triangulation_3_cell_primitive<EPICK,
                                                  C3t3::Triangulation> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
using namespace CGAL::Three;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Point_container Pc;
typedef Viewer_interface Vi;


struct Scene_c3t3_item_priv {

  void init_default_values() {
    sharp_edges_angle = -1;
    detect_borders = false;
  }

  Scene_c3t3_item_priv(Scene_c3t3_item* item)
    : item(item), c3t3(), is_valid(true){
    init_default_values();
  }

  Scene_c3t3_item_priv(const C3t3& c3t3_, Scene_c3t3_item* item)
    : item(item), c3t3(c3t3_), is_valid(true){
    init_default_values();
  }

  Scene_c3t3_item* item;
  C3t3 c3t3;
  bool is_surface;
  //only for optimizers
  double sharp_edges_angle;
  bool detect_borders;
  bool is_valid;
};


Scene_c3t3_item::Scene_c3t3_item(bool is_surface)
  : Scene_triangulation_3_item()
{
  d = new Scene_c3t3_item_priv(this);
}

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3, bool is_surface)
  : Scene_triangulation_3_item(c3t3.triangulation())
{
  d = new Scene_c3t3_item_priv(c3t3, this);
}

Scene_c3t3_item::~Scene_c3t3_item()
{
}

const C3t3&
Scene_c3t3_item::c3t3() const {
  return d->c3t3;
}

C3t3&
Scene_c3t3_item::c3t3()
{
  return d->c3t3;
}


const T3&
Scene_c3t3_item::triangulation() const {
  return d->c3t3.triangulation();
}

T3&
Scene_c3t3_item::triangulation()
{
  return d->c3t3.triangulation();
}

bool Scene_c3t3_item::do_take_cell(const T3::Cell_handle& c) const
{
  return d->c3t3.is_in_complex(c);
}

bool Scene_c3t3_item::do_take_facet(const T3::Facet& f)const
{
  return (d->c3t3.is_in_complex(f));
}

bool Scene_c3t3_item::do_take_vertex(const T3::Vertex& )const
{
  return true;
}

bool Scene_c3t3_item::is_facet_oriented(const T3::Facet& f)const
{
  const Tr::Cell_handle& cell = f.first;
  const int& index = f.second;
  return (index % 2 == 1) == d->c3t3.is_in_complex(cell);
}

QMenu* Scene_c3t3_item::contextMenu()
{

  const char* prop_name = "Menu modified by Scene_c3t3_item.";

  QMenu* menu = Scene_triangulation_3_item::contextMenu();

  // Use dynamic properties:
  // https://doc.qt.io/qt-5/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if (!menuChanged) {

    QAction* actionExportFacetsInComplex =
      menu->addAction(tr("Export facets in complex"));
    actionExportFacetsInComplex->setObjectName("actionExportFacetsInComplex");
    connect(actionExportFacetsInComplex,
      SIGNAL(triggered()), this,
      SLOT(export_facets_in_complex()));

    if(is_valid())
    {
      QAction* actionShowSpheres =
          menu->addAction(tr("Show protecting &spheres"));
      actionShowSpheres->setCheckable(true);
      actionShowSpheres->setObjectName("actionShowSpheres");
      connect(actionShowSpheres, SIGNAL(toggled(bool)),
              this, SLOT(show_spheres(bool)));

      QAction* actionShowCNC =
          menu->addAction(tr("Show cells not in complex"));
      actionShowCNC->setCheckable(true);
      actionShowCNC->setObjectName("actionShowCNC");
      connect(actionShowCNC, SIGNAL(toggled(bool)),
              this, SLOT(show_cnc(bool)));
    }
    menu->setProperty(prop_name, true);
  }
  return menu;
}


bool Scene_c3t3_item::is_valid() const
{
  return d->is_valid;
}
void Scene_c3t3_item::set_valid(bool b)
{
  d->is_valid = b;
}

void Scene_c3t3_item::show_spheres(bool b)
{
  if(is_valid())
  {
    d->spheres_are_shown = b;
    contextMenu()->findChild<QAction*>("actionShowSpheres")->setChecked(b);
    if(b && !d->spheres)
    {
      d->spheres = new Scene_spheres_item(this, triangulation().number_of_vertices(), true);
      connect(d->spheres, &Scene_spheres_item::picked,
              this, [this](std::size_t id)
      {
        if(id == (std::size_t)(-1))
          return;
        QString msg = QString("Vertex's index : %1; Vertex's in dimension: %2.").arg(d->tr_vertices[id].index()).arg(d->tr_vertices[id].in_dimension());
        CGAL::Three::Three::information(msg);
        CGAL::Three::Three::mainViewer()->displayMessage(msg, 5000);

      });
      d->spheres->setName("Protecting spheres");
      d->spheres->setRenderingMode(Gouraud);
      connect(d->spheres, SIGNAL(destroyed()), this, SLOT(reset_spheres()));
      connect(d->spheres, SIGNAL(on_color_changed()), this, SLOT(on_spheres_color_changed()));
      d->computeSpheres();
      lockChild(d->spheres);
      scene->addItem(d->spheres);
      scene->changeGroup(d->spheres, this);
    }
    else if (!b && d->spheres!=NULL)
    {
      unlockChild(d->spheres);
      scene->erase(scene->item_id(d->spheres));
    }
    Q_EMIT redraw();
  }
}

void Scene_c3t3_item::set_sharp_edges_angle(double a) { d->sharp_edges_angle = a; }
double Scene_c3t3_item::get_sharp_edges_angle() { return d->sharp_edges_angle; }

void Scene_c3t3_item::set_detect_borders(bool b) { d->detect_borders = b;}
bool Scene_c3t3_item::get_detect_borders() { return d->detect_borders; }

std::size_t Scene_c3t3_item::number_of_patches() const
{
  return d->surface_patch_indices_.size();
}
#include "Scene_c3t3_item.moc"

