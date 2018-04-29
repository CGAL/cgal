#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "ui_Remeshing_dialog.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMenu>
#include <QApplication>
#include <QtPlugin>
#include "Scene_polyhedron_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include <QInputDialog>
#include <QStringList>

#include "C2t3_type.h"
#include "Scene_c2t3_item.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_surface_3.h>

namespace CGAL{

template <class TriangleMesh, class GeomTraits>
class Offset_function
{
  typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef AABB_traits<GeomTraits, Primitive> Traits;
  typedef AABB_tree<Traits> Tree;
  typedef Side_of_triangle_mesh<TriangleMesh, GeomTraits> Side_of;

public:

  Offset_function(TriangleMesh& tm, double offset_distance)
    : m_tree_ptr(new Tree(boost::begin(faces(tm)),
                          boost::end(faces(tm)),
                          tm) )
    , m_side_of_ptr( new Side_of(*m_tree_ptr) )
    , m_offset_distance(offset_distance)
    , m_is_closed( is_closed(tm) )
  {
    m_tree_ptr->accelerate_distance_queries();
  }

  double operator()(const typename GeomTraits::Point_3& p) const
  {
    using CGAL::sqrt;

    Bounded_side side = m_is_closed?m_side_of_ptr->operator()(p):ON_UNBOUNDED_SIDE;
    if (side==ON_BOUNDARY) return m_offset_distance;

    typename GeomTraits::Point_3 closest_point = m_tree_ptr->closest_point(p);
    double distance = sqrt(squared_distance(p, closest_point));

    return (side == ON_UNBOUNDED_SIDE ? -distance : distance) + m_offset_distance;
  }

private:
  boost::shared_ptr<Tree> m_tree_ptr;
  boost::shared_ptr<Side_of> m_side_of_ptr;
  double m_offset_distance;
  bool m_is_closed;

};

} //end of CGAL namespace


Scene_polyhedron_item* make_item(Polyhedron* poly)
{
  return new Scene_polyhedron_item(poly);
}

Scene_surface_mesh_item* make_item(SMesh* sm)
{
  return new Scene_surface_mesh_item(sm);
}

// declare the CGAL function
template<class Mesh>
CGAL::Three::Scene_item* cgal_off_meshing(QWidget*,
                                          Mesh* tm_ptr,
                                          const double offset_value,
                                          const double angle,
                                          const double sizing,
                                          const double approx,
                                          int tag)
{
  typedef Tr::Geom_traits GT;
  typedef CGAL::Offset_function<Mesh, GT> Offset_function;
  typedef CGAL::Implicit_surface_3<GT, Offset_function> Surface_3;
  typedef GT::Sphere_3 Sphere_3;

  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(*tm_ptr);

  GT::Point_3 center((bbox.xmax()+bbox.xmin())/2,
                     (bbox.ymax()+bbox.ymin())/2,
                     (bbox.zmax()+bbox.zmin())/2);
  double sqrad = 0.6 * std::sqrt( CGAL::square(bbox.xmax()-bbox.xmin())+
                                  CGAL::square(bbox.ymax()-bbox.ymin())+
                                  CGAL::square(bbox.zmax()-bbox.zmin()) )
                + offset_value;
  sqrad=CGAL::square(sqrad);

  CGAL::Timer timer;
  timer.start();

  Offset_function offset_function(*tm_ptr, offset_value);

  Tr tr;
  C2t3 c2t3 (tr);

  // defining the surface
  Surface_3 surface(offset_function,
                    Sphere_3(center, sqrad)); // bounding sphere

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle, sizing, approx);

  // meshing surface
  switch(tag) {
  case 0:
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    break;
  case 1:
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
    break;
  default:
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());
  }

  timer.stop();
  std::cerr << "done (" << timer.time() << " ms, " << c2t3.triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(c2t3.triangulation().number_of_vertices() > 0)
  {
    // add remesh as new polyhedron
    Mesh *pRemesh = new Mesh;
    CGAL::facets_in_complex_2_to_triangle_mesh<C2t3, Mesh>(c2t3, *pRemesh);
    if(c2t3.number_of_facets() != num_faces(*pRemesh))
    {
      delete pRemesh;
      std::stringstream temp_file;
      if(!CGAL::output_surface_facets_to_off(temp_file, c2t3))
      {
        std::cerr << "Cannot write the mesh to an off file!\n";
        return 0;
      }
      Scene_polygon_soup_item* soup = new Scene_polygon_soup_item();
      if(!soup->load(temp_file))
      {
        std::cerr << "Cannot reload the mesh from an off file!\n";
        return 0;
      }
      else
        return soup;
    } else {
      return make_item(pRemesh);
    }
  }
  else
    return 0;
}

using namespace CGAL::Three;
class Polyhedron_demo_offset_meshing_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionOffsetMeshing = new QAction(tr("Offset meshing"), mw);
    actionOffsetMeshing->setProperty("subMenuName", "3D Surface Mesh Generation");
    if(actionOffsetMeshing) {
      connect(actionOffsetMeshing, SIGNAL(triggered()),
              this, SLOT(offset_meshing()));
    }
  }

  bool applicable(QAction*) const {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    return qobject_cast<Scene_polyhedron_item*>(item) ||
        qobject_cast<Scene_surface_mesh_item*>(item);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionOffsetMeshing;
  }
public Q_SLOTS:
  void offset_meshing();

private:
  QAction* actionOffsetMeshing;
  Scene_interface *scene;
  QMainWindow *mw;
}; // end class Polyhedron_demo_offset_meshing_plugin

void Polyhedron_demo_offset_meshing_plugin::offset_meshing()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_item* item = scene->item(index);

  Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(item);
  Scene_surface_mesh_item* sm_item =
      qobject_cast<Scene_surface_mesh_item*>(item);

  Polyhedron* pMesh = NULL;
  SMesh* sMesh = NULL;
  if(poly_item)
  {
    pMesh = poly_item->polyhedron();

    if(!pMesh)
      return;
  }
  else if(sm_item)
  {
    sMesh = sm_item->face_graph();
    if(!sMesh)
      return;
  }
  else
    return;

  double diag = scene->len_diagonal();
  double offset_value = QInputDialog::getDouble(mw,
                                                QString("Choose Offset Value"),
                                                QString("Offset Value (use negative number for inset)"),
                                                0.1*diag,
                                                -(std::numeric_limits<double>::max)(),
                                                (std::numeric_limits<double>::max)(), 10);

  QDialog dialog(mw);
  Ui::Remeshing_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  ui.sizing->setDecimals(4);
  ui.sizing->setRange(diag * 10e-6, // min
                      diag); // max
  ui.sizing->setValue(diag * 0.05); // default value

  ui.approx->setDecimals(6);
  ui.approx->setRange(diag * 10e-7, // min
                      diag); // max
  ui.approx->setValue(diag * 0.005);


  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;

  const double angle = ui.angle->value();
  const double approx = ui.approx->value();
  const double sizing = ui.sizing->value();
  const int tag_index = ui.tags->currentIndex();

  if(tag_index < 0) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cerr << "mesh with:"
            << "\n  angle=" << angle
            << "\n  sizing=" << sizing
            << "\n  approx=" << approx
            << "\n  tag=" << tag_index
            << std::boolalpha
            << std::endl;
  CGAL::Three::Scene_item* new_item;
  if(pMesh)
    new_item = cgal_off_meshing(mw,
                                pMesh,
                                offset_value,
                                angle,
                                sizing,
                                approx,
                                tag_index);
  else
    new_item = cgal_off_meshing(mw,
                                sMesh,
                                offset_value,
                                angle,
                                sizing,
                                approx,
                                tag_index);

  if(new_item)
  {
    new_item->setName(tr("%1 offset %5 (%2 %3 %4)")
                      .arg(item->name())
                      .arg(angle)
                      .arg(sizing)
                      .arg(approx)
                      .arg(offset_value));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(item->renderingMode());
    scene->addItem(new_item);
    item->setVisible(false);
    scene->itemChanged(index);
  }

  QApplication::restoreOverrideCursor();
}

#include "Offset_meshing_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
