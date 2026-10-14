#include <QElapsedTimer>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QStringList>
#include <QInputDialog>
#include "Scene_surface_mesh_item.h"
#include "Scene_lcc_item.h"
#include "SMesh_type.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>

#include <CGAL/generate_hexahedral_mesh_using_two_refinement.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Mesh_smoothing_3/mesh_representations.h>
#include <CGAL/Mesh_smoothing_3/Mesh_smoothing_3.h>
#include <CGAL/Mesh_smoothing_3/Projectors.h>

using namespace CGAL::Three;
class CGAL_Lab_hexmeshing_lcc_plugin :
  public QObject,
  public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")
public:
  void init(QMainWindow*mw,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    scene = scene_interface;
    this->mw = mw;
    actionHexmeshing = new QAction("Hexmeshing using 2-refinement", mw);
    connect(actionHexmeshing, SIGNAL(triggered()), this, SLOT(on_actionHexmeshing_triggered()));
    actionHexmeshingSmoothing = new QAction("Hexmeshing Smoothing", mw);
    connect(actionHexmeshing, SIGNAL(triggered()), this, SLOT(on_actionHexmeshingSmoothing_triggered()));
    _actions << actionHexmeshingSmoothing;
  }

  QList<QAction*> actions()const {return _actions;}

  bool applicable(QAction* action) const {
    if (action == actionHexmeshing)
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    else
    {
      // smoothing
      if (scene->selectionIndices().size()==2)
      {
        return ( qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->selectionIndices()[0])) && qobject_cast<Scene_lcc_item*>(scene->item(scene->selectionIndices()[1])) ) ||
               ( qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->selectionIndices()[1])) && qobject_cast<Scene_lcc_item*>(scene->item(scene->selectionIndices()[0])) );
      }
      else
        return false;
    }

  }

public Q_SLOTS:
  void on_actionHexmeshing_triggered();
  void on_actionHexmeshingSmoothing_triggered();
private:
  QList<QAction*> _actions;
  QAction* actionHexmeshing;
  QAction* actionHexmeshingSmoothing;
  Scene_interface* scene;
  QMainWindow* mw;
}; // end CGAL_Lab_hexmeshing_lcc_plugin


void CGAL_Lab_hexmeshing_lcc_plugin::on_actionHexmeshing_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(sm_item)
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    QElapsedTimer time;
    time.start();
    std::cout << "Hexmeshing using 2-refinement...";

    // add convex hull as new polyhedron
    SMesh* pMesh = sm_item->polyhedron();

    bool ok = true;
    int cube_cells_per_dim = QInputDialog::getInt(mw,
                                                  QString("Initial grid resolution per dimension"),
                                                  QString("Choose initial grid resolution"),
                                                  10 /* value */, 2 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                                  &ok);
    if (!ok)
      cube_cells_per_dim=10;

    int nb_levels = QInputDialog::getInt(mw,
                                                  QString("Initial grid resolution per dimension"),
                                                  QString("Choose initial grid resolution"),
                                                  2 /* value */, 1 /* min */, (std::numeric_limits<int>::max)() /* max */, 1 /*step*/,
                                                  &ok);
    if (!ok)
      nb_levels=2;

    bool trim=true;
    bool smooth=true;

    using LCC=Scene_lcc_item::LCC;
    LCC lcc=CGAL::generate_hexahedral_mesh_using_two_refinement<LCC>
      (*pMesh, cube_cells_per_dim, nb_levels, CGAL::parameters::use_trimming(trim).use_smoothing(smooth));

    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    sm_item->setVisible(false);
    Scene_lcc_item* new_item = new Scene_lcc_item(lcc);
    new_item->invalidateOpenGLBuffers();
    new_item->setName(tr("%1 (hexmeshing-%2-%3)").arg(scene->item(index)->name())
                                                 .arg(cube_cells_per_dim)
                                                 .arg(nb_levels));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(FlatPlusEdges);
    scene->addItem(new_item);

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

void CGAL_Lab_hexmeshing_lcc_plugin::on_actionHexmeshingSmoothing_triggered()
{
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->selectionIndices()[0]));
  if (!sm_item) sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->selectionIndices()[1]));
  Scene_lcc_item* lcc_item = qobject_cast<Scene_lcc_item*>(scene->item(scene->selectionIndices()[0]));
  if (!lcc_item) lcc_item = qobject_cast<Scene_lcc_item*>(scene->item(scene->selectionIndices()[1]));

  using LCC=Scene_lcc_item::LCC;
  LCC* pLcc = lcc_item->lcc();
  SMesh* pMesh = sm_item->polyhedron();


  QApplication::setOverrideCursor(Qt::WaitCursor);

  SMesh poly = *pMesh;
  CGAL::Polygon_mesh_processing::triangulate_faces(poly);
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<SMesh, K>;
  Mesh_domain mesh_domain(poly);
  CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector<Mesh_domain> projector(mesh_domain);

  using ConcurrencyTag = CGAL::Mesh_smoothing_3::Parallel_if_available_tag;
  unsigned max_nb_metric_evaluations = 10000;
  double time_limit = 0.;
  bool verbose = true;

  LCC& lcc = *pLcc;

  using LCC_mixed_mesh = CGAL::Mesh_smoothing_3::cgal_types::LCC_mixed_mesh<LCC>;
  using Mixed_mesh = CGAL::Mesh_smoothing_3::helper_structures::Mixed_mesh_wrapper<LCC_mixed_mesh>;

  LCC_mixed_mesh lcc_mixed_mesh(lcc);
  Mixed_mesh mixed_mesh(lcc_mixed_mesh);

  using Lcc_surface = CGAL::Mesh_smoothing_3::cgal_types::LCC_surface_mesh<LCC>;
  Lcc_surface surface(lcc);

  using Curve_network = CGAL::Mesh_smoothing_3::default_structures::Empty_edge_network<Mixed_mesh::Vertex_descriptor>;

  CGAL::Mesh_smoothing_3::Mesh_smoother<Mixed_mesh, Lcc_surface, Curve_network, ConcurrencyTag>
  smoother(mixed_mesh, surface);

  smoother.set_predicates_mode(CGAL::Mesh_smoothing_3::Parameters::STRONG_ENFORCEMENT);
  smoother.set_verbose(verbose);
  smoother.set_maximum_running_time(time_limit);
  smoother.set_maximum_number_of_metric_evaluations(max_nb_metric_evaluations);

  smoother.set_boundary_query([&](std::vector<K::Point_3> const& pts, Lcc_surface::Surface_patch_index patch_face) {
      auto proj = projector.patch_face_projection_plane(patch_face, pts);
      return std::make_tuple(proj.origin(), proj.vector(), 1.);
  });

  auto result = smoother.run();

  std::cout << "Number of inverted elements: " << result.nb_invalid_elements << std::endl;
  std::cout << "Number of vertex updates: " << result.nb_vertex_updates << std::endl;
  std::cout << "Number of metric evaluations: " << result.nb_metric_evaluations << std::endl;
  std::cout << "Pre-processing time: " << result.pre_processing_time << std::endl;
  std::cout << "Smoothing time: " << result.optimization_time << std::endl;

  QApplication::restoreOverrideCursor();

  lcc_item->invalidateOpenGLBuffers();
  lcc_item->itemChanged();
}

#include "Hexmeshing_lcc_plugin.moc"
