#include <QtCore/qglobal.h>

#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include "Variational_medial_axis_skeleton_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_polylines_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/boost/graph/helpers.h>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QInputDialog>

#include <CGAL/make_mesh_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/variational_medial_axis_sampling.h>

#include "C3t3_type.h"

using namespace CGAL::Three;
class Variational_medial_axis_skeleton_plugin :
    public QObject,
    public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTraceSkeleton << actionTraceSurfaceFromSkeleton;
  }

  bool applicable(QAction*) const {
    if(scene->numberOfEntries() > 0)
    {
      int item_id = scene->mainSelectionIndex();
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(item_id));
    }
    return false;
  }
public Q_SLOTS:

  void trace_skeleton();
  void trace_surface();
  void enableAction();
  void connectNewViewer(QObject* o)
  {
    for(int i=0; i<scene->numberOfEntries(); ++i)
    {
      Variational_medial_axis_skeleton_item* item = qobject_cast<Variational_medial_axis_skeleton_item*>(
            scene->item(i));
      if(item)
        o->installEventFilter(item);
    }
  }

private:
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  QAction* actionTraceSkeleton;
  QAction* actionTraceSurfaceFromSkeleton;
}; // end Variational_medial_axis_skeleton_plugin

void Variational_medial_axis_skeleton_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionTraceSkeleton = new QAction(tr("Variational Medial Axis Skeleton"), mainWindow);
  connect(actionTraceSkeleton, SIGNAL(triggered()),
          this, SLOT(trace_skeleton()));
  actionTraceSurfaceFromSkeleton = new QAction(tr("Surface from Variational Medial Axis Skeleton"), mainWindow);
  connect(actionTraceSurfaceFromSkeleton, SIGNAL(triggered()),
          this, SLOT(trace_surface()));
  connect(mw, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));
}

void Variational_medial_axis_skeleton_plugin::trace_skeleton()
{
  bool ok=true;
  const int nbs = QInputDialog::getInt(mw,
                                       QString("Expected number of spheres"),
                                       QString("Expected number of spheres"),
                                       200,
                                       0,
                                       (std::numeric_limits<int>::max)(),
                                       1,
                                       &ok);

  if (!ok) return;


  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

  //todo: add dialog

  Variational_medial_axis_skeleton_item* item = new Variational_medial_axis_skeleton_item(scene, sm_item, nbs);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  item->setName(sm_item->name() + " medial axis");
  item->setRenderingMode(FlatPlusEdges);
  for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(item);

  scene->addItem(item);
  item->fill_subitems();

  sm_item->setVisible(false);

  QApplication::restoreOverrideCursor();
}

void Variational_medial_axis_skeleton_plugin::trace_surface()
{
  Scene_surface_mesh_item* source_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  const double angle=25;

  bool ok=true;
  const int nbs = QInputDialog::getInt(mw,
                                       QString("Expected number of spheres used for the skeleton"),
                                       QString("Expected number of spheres used for the skeleton"),
                                       200,
                                       0,
                                       (std::numeric_limits<int>::max)(),
                                       1,
                                       &ok);

  if (!ok) return;

  double diag = source_item->bboxDiagonal();
  const double sizing = QInputDialog::getDouble(mw,
                                                QString("Facet sizing"),
                                                QString("Facet sizing"),
                                                diag * 0.05,
                                                diag * 10e-6,
                                                diag,
                                                10,
                                                &ok);

  if (!ok) return;

  const double approx = QInputDialog::getDouble(mw,
                                                QString("Approximation distance"),
                                                QString("Approximation distance"),
                                                diag * 0.005,
                                                diag * 10e-7,
                                                diag,
                                                10,
                                                &ok);

  if (!ok) return;


  using GT = EPICK;
  using Sphere_3 = GT::Sphere_3;

  using Mesh_domain_base = CGAL::Labeled_mesh_domain_3<GT, int, int>;
  using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain_base>;
  using Tr = C3t3::Triangulation;
  using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

  using Mesh = Scene_surface_mesh_item::Face_graph;

  const CGAL::Bbox_3 bbox = source_item->bbox();



  QApplication::setOverrideCursor(Qt::WaitCursor);


  const GT::Point_3 center((bbox.xmax() + bbox.xmin()) / 2,
                           (bbox.ymax() + bbox.ymin()) / 2,
                           (bbox.zmax() + bbox.zmin()) / 2);
  const double rad = 0.6 * std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                     CGAL::square(bbox.ymax() - bbox.ymin()) +
                                     CGAL::square(bbox.zmax() - bbox.zmin()));
  const double sqrad = CGAL::square(rad);

  CGAL::Timer timer;
  timer.start();

  CGAL::Variational_medial_axis<Mesh, CGAL::Parallel_if_available_tag> vmas(*source_item->face_graph());
  vmas.compute_variational_medial_axis_sampling(CGAL::parameters::number_of_spheres(nbs));
  auto skeleton = vmas.export_skeleton();

  CGAL::Medial_skeleton_offset_function<Mesh> implicit_function(skeleton);

  namespace p = CGAL::parameters;

  Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain
    (p::function = implicit_function,
     p::bounding_object = Sphere_3(center, sqrad),
     p::relative_error_bound = 1e-7,
     p::construct_surface_patch_index = [](int i, int j) { return (i * 1000 + j); });

  const CGAL::Mesh_facet_topology topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH;
  auto manifold_option = p::non_manifold();
//  if(tag == 1)
    manifold_option = p::manifold_with_boundary();
//  if(tag == 2)
//    manifold_option = p::manifold();

  Mesh_criteria criteria(p::facet_angle = angle,
                         p::facet_size = sizing,
                         p::facet_distance = approx,
                         p::facet_topology = topology);

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      p::no_perturb(),
                                      p::no_exude(),
                                      manifold_option);

  timer.stop();
  std::cerr << "done (" << timer.time() << " ms, " << c3t3.triangulation().number_of_vertices() << " vertices)" << std::endl;

  if(c3t3.number_of_facets_in_complex() > 0)
  {
    Mesh* smesh = new Mesh();
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, *smesh);

    Scene_surface_mesh_item* skel_mesh_item = new Scene_surface_mesh_item(smesh);
    skel_mesh_item->setName(tr("%1 skeleton mesh (%2 %3 %4)").arg(source_item->name())
                                                             .arg(angle)
                                                             .arg(sizing)
                                                             .arg(approx));
    skel_mesh_item->setColor(Qt::magenta);
    skel_mesh_item->setRenderingMode(Wireframe);
    Three::scene()->addItem(skel_mesh_item);
  }

  QApplication::restoreOverrideCursor();
}

void Variational_medial_axis_skeleton_plugin::enableAction() {
  actionTraceSkeleton->setEnabled(true);
  actionTraceSurfaceFromSkeleton->setEnabled(true);
}

#include "Variational_medial_axis_skeleton_plugin.moc"
