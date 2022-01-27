#include "config.h"
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
    : item(item), c3t3(), is_valid(true), computed_stats(false){
  }

  Scene_c3t3_item_priv(const C3t3& c3t3_, Scene_c3t3_item* item)
    : item(item), c3t3(c3t3_), is_valid(true), computed_stats(false){
  }

  ~Scene_c3t3_item_priv()
  {
    c3t3.clear();
  }

  void draw_triangle_edges_cnc(const Tr::Bare_point& pa,
                               const Tr::Bare_point& pb,
                               const Tr::Bare_point& pc) const;
  Scene_c3t3_item* item;
  C3t3 c3t3;
  bool is_surface;
  //only for optimizers
  double sharp_edges_angle;
  bool detect_borders;
  bool is_valid;
  bool cnc_are_shown;

  enum STATS {
    MIN_EDGES_LENGTH = 0,
    MAX_EDGES_LENGTH,
    MEAN_EDGES_LENGTH,
    MIN_DIHEDRAL_ANGLE,
    MAX_DIHEDRAL_ANGLE,
    MEAN_DIHEDRAL_ANGLE,
    NB_SPHERES,
    NB_VERTICES,
    NB_TETS,
    SMALLEST_RAD_RAD,
    SMALLEST_EDGE_RAD,
    BIGGEST_VL3_CUBE,
    NB_SUBDOMAINS,
    NB_CNC
  };

  mutable std::vector<float> positions_lines_not_in_complex;
  mutable std::size_t positions_lines_not_in_complex_size;
  mutable std::size_t nb_cnc;
  mutable bool computed_stats;

  void push_point(std::vector<float>& points, const EPICK::Point_3& p,
                  const CGAL::qglviewer::Vec& offset) const
  {
    points.push_back(static_cast<float>(p.x()+offset.x));
    points.push_back(static_cast<float>(p.y()+offset.y));
    points.push_back(static_cast<float>(p.z()+offset.z));
  }

  void push_edge(std::vector<float>& edges,
                 const EPICK::Point_3& pa,
                 const EPICK::Point_3& pb,
                 const CGAL::qglviewer::Vec& offset) const
  {
    push_point(edges, pa, offset);
    push_point(edges, pb, offset);
  }
};

void Scene_c3t3_item_priv::draw_triangle_edges_cnc(const Tr::Bare_point& pa,
                                                   const Tr::Bare_point& pb,
                                                   const Tr::Bare_point& pc) const
{
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  push_edge(positions_lines_not_in_complex, pa, pb, offset);
  push_edge(positions_lines_not_in_complex, pb, pc, offset);
  push_edge(positions_lines_not_in_complex, pc, pa, offset);
}


void Scene_c3t3_item::common_constructor(bool is_surface)
{
  d->is_surface = is_surface;
  d->cnc_are_shown = false;
  setEdgeContainer(CNC, new Ec(Vi::PROGRAM_NO_SELECTION, false));
}

Scene_c3t3_item::Scene_c3t3_item(bool is_surface)
  : Scene_triangulation_3_item(!is_surface)
{
  d = new Scene_c3t3_item_priv(this);
  common_constructor(is_surface);

}

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3, bool is_surface)
  : Scene_triangulation_3_item(c3t3.triangulation(), !is_surface)
{
  d = new Scene_c3t3_item_priv(c3t3, this);
  common_constructor(is_surface);
}

Scene_c3t3_item::~Scene_c3t3_item()
{
  delete d;
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

void
Scene_c3t3_item::c3t3_changed()
{
  triangulation_changed();
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

bool Scene_c3t3_item::do_take_vertex(const T3::Vertex_handle& v)const
{
  return d->c3t3.is_in_complex(v);
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

    if(is_valid())
    {
      QAction* actionShowCNC =
          menu->addAction(tr("Show cells not in complex"));
      actionShowCNC->setCheckable(true);
      actionShowCNC->setObjectName("actionShowCNC");
      connect(actionShowCNC, SIGNAL(toggled(bool)),
              this, SLOT(show_cnc(bool)));
    }

    QAction* actionExportFacetsInComplex =
      menu->addAction(tr("Export facets in complex"));
    actionExportFacetsInComplex->setObjectName("actionExportFacetsInComplex");
    connect(actionExportFacetsInComplex,
      SIGNAL(triggered()), this,
      SLOT(export_facets_in_complex()));

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

void Scene_c3t3_item::show_cnc(bool b)
{
  if(is_valid())
  {
    d->cnc_are_shown = b;
    contextMenu()->findChild<QAction*>("actionShowCNC")->setChecked(b);
    Q_EMIT redraw();
  }
}

bool Scene_c3t3_item::load_binary(std::istream& is)
{
  if(!CGAL::Mesh_3::load_binary_file(is, c3t3())) return false;
  resetCutPlane();
  if(is.good()) {
    c3t3_changed();
    changed();
    return true;
  }
  else
    return false;
}

void Scene_c3t3_item::export_facets_in_complex()
{
  SMesh outmesh;
  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3(), outmesh);
  Scene_surface_mesh_item* item = new Scene_surface_mesh_item(std::move(outmesh));
  item->setName(QString("%1_%2").arg(this->name()).arg("facets"));
  scene->addItem(item);
  this->setVisible(false);
}

void Scene_c3t3_item::drawEdges(Viewer_interface *viewer) const
{
  Scene_triangulation_3_item::drawEdges(viewer);
//add cnc
  if(d->cnc_are_shown)
  {
    getEdgeContainer(CNC)->setColor(QColor(Qt::black));
    getEdgeContainer(CNC)->draw(viewer, true);
  }
}

void Scene_c3t3_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
{
  Scene_triangulation_3_item::initializeBuffers(viewer);
  // add cnc
  {
    getEdgeContainer(Scene_c3t3_item::CNC)->initializeBuffers(viewer);
    getEdgeContainer(Scene_c3t3_item::CNC)->setFlatDataSize(
          d->positions_lines_not_in_complex_size);
    d->positions_lines_not_in_complex.clear();
    d->positions_lines_not_in_complex.shrink_to_fit();
  }
}
void Scene_c3t3_item::computeElements()const
{
  Scene_triangulation_3_item::computeElements();
  QApplication::setOverrideCursor(Qt::WaitCursor);
  //add cnc
  //the cells not in the complex
  Geom_traits::Construct_point_3 wp2p
    = c3t3().triangulation().geom_traits().construct_point_3_object();
  for(C3t3::Triangulation::Cell_iterator
      cit = c3t3().triangulation().finite_cells_begin(),
      end = c3t3().triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    if(!c3t3().is_in_complex(cit))
    {

      bool has_far_point = false;
      for(int i=0; i<4; i++)
        if(c3t3().in_dimension(cit->vertex(i)) == -1)
        {
          has_far_point = true;
          break;
        }
      if(!has_far_point)
      {
        const Tr::Bare_point& p1 = wp2p(cit->vertex(0)->point());
        const Tr::Bare_point& p2 = wp2p(cit->vertex(1)->point());
        const Tr::Bare_point& p3 = wp2p(cit->vertex(2)->point());
        const Tr::Bare_point& p4 = wp2p(cit->vertex(3)->point());
        d->draw_triangle_edges_cnc(p1, p2, p4);
        d->draw_triangle_edges_cnc(p1, p3, p4);
        d->draw_triangle_edges_cnc(p2, p3, p4);
        d->draw_triangle_edges_cnc(p1, p2, p3);
      }
    }
  }
  getEdgeContainer(CNC)->allocate(
        Ec::Vertices,
        d->positions_lines_not_in_complex.data(),
        static_cast<int>(d->positions_lines_not_in_complex.size()*sizeof(float)));
  d->positions_lines_not_in_complex_size = d->positions_lines_not_in_complex.size();
  QApplication::restoreOverrideCursor();
}

QString Scene_c3t3_item::computeStats(int type)
{
  if(type != Scene_c3t3_item_priv::NB_CNC)
    return Scene_triangulation_3_item::computeStats(type);
  if(!d->computed_stats)
  {
    d->nb_cnc = 0;
    for(C3t3::Triangulation::Cell_iterator
        cit = d->c3t3.triangulation().finite_cells_begin(),
        end = d->c3t3.triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      if(!d->c3t3.is_in_complex(cit))
      {

        bool has_far_point = false;
        for(int i=0; i<4; i++)
          if(d->c3t3.in_dimension(cit->vertex(i)) == -1)
          {
            has_far_point = true;
            break;
          }
        if(!has_far_point)
          ++d->nb_cnc;
      }
    }
    d->computed_stats = true;
  }
  return QString::number(d->nb_cnc);
}

CGAL::Three::Scene_item::Header_data Scene_c3t3_item::header() const
{
  CGAL::Three::Scene_item::Header_data data = Scene_triangulation_3_item::header();
  data.categories[0] = std::pair<QString,int>(QString("Properties"),14);
  data.titles.append(QString("#Cells not in Complex"));
  return data;
}

bool Scene_c3t3_item::has_cnc()const { return d->cnc_are_shown;}

bool Scene_c3t3_item::is_surface() const
{
  return d->is_surface;
}

void Scene_c3t3_item::set_sharp_edges_angle(double a) { d->sharp_edges_angle = a; }
double Scene_c3t3_item::get_sharp_edges_angle() { return d->sharp_edges_angle; }

void Scene_c3t3_item::set_detect_borders(bool b) { d->detect_borders = b;}
bool Scene_c3t3_item::get_detect_borders() { return d->detect_borders; }

Scene_c3t3_item* Scene_c3t3_item::clone() const
{
  return new Scene_c3t3_item(d->c3t3, d->is_surface);
}

void Scene_c3t3_item::copyProperties(Scene_item *item)
{
  Scene_triangulation_3_item::copyProperties(item);
  Scene_c3t3_item* c3t3_item = qobject_cast<Scene_c3t3_item*>(item);
  if(!c3t3_item)
    return;
  show_cnc(c3t3_item->has_cnc());

}
