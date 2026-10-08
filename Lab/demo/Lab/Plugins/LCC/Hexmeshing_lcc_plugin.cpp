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
    QAction *actionHexmeshing = new QAction("Hexmeshing using 2-refinement", mw);
    connect(actionHexmeshing, SIGNAL(triggered()), this, SLOT(on_actionHexmeshing_triggered()));
    _actions <<actionHexmeshing;
  }

  QList<QAction*> actions()const {return _actions;}

  bool applicable(QAction*) const {
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionHexmeshing_triggered();
private:
  QList<QAction*> _actions;
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

#include "Hexmeshing_lcc_plugin.moc"
