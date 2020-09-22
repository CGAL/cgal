#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_DUMP_REMESHING_STEPS
//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
//#define CGAL_TETRAHEDRAL_REMESHING_PROFILE

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_c3t3_item.h"
#include "C3t3_type.h"

#include <CGAL/tetrahedral_remeshing.h>

#include <unordered_map>
#include <memory>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Tetrahedral_remeshing_dialog.h"

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
std::size_t nb_topology_test = 0;
std::size_t nb_impossible = 0;
std::size_t nb_valid_collapse = 0;
std::size_t nb_invalid_collapse = 0;
std::size_t nb_invalid_lengths = 0;
std::size_t nb_invalid_collapse_short = 0;
std::size_t nb_orientation_v0 = 0;
std::size_t nb_orientation_v1 = 0;
std::size_t nb_orientation_midpoint = 0;
std::size_t nb_test_v0 = 0;
std::size_t nb_test_v1 = 0;
std::size_t nb_test_midpoint = 0;
#endif

using namespace CGAL::Three;
class Polyhedron_demo_tetrahedral_remeshing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "tetrahedral_remeshing_plugin.json")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionTetrahedralRemeshing_ = new QAction("Tetrahedral Remeshing", mw);
    if (actionTetrahedralRemeshing_) {
      connect(actionTetrahedralRemeshing_, SIGNAL(triggered()),
        this, SLOT(tetrahedral_remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTetrahedralRemeshing_;
  }

  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_c3t3_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void tetrahedral_remeshing()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_c3t3_item* c3t3_item =
      qobject_cast<Scene_c3t3_item*>(scene->item(index));

    if (c3t3_item)
    {
      // Create dialog box
      QDialog dialog(mw);
      Ui::Tetrahedral_remeshing_dialog ui
        = tet_remeshing_dialog(&dialog, c3t3_item);

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }
      double target_length = ui.edgeLength_dspinbox->value();
      unsigned int nb_iter = ui.nbIterations_spinbox->value();
      bool protect = ui.protect_checkbox->isChecked();
      bool smooth_edges = ui.smoothEdges_checkBox->isChecked();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QElapsedTimer time;
      time.start();

      CGAL::tetrahedral_isotropic_remeshing(
          c3t3_item->c3t3(),
          target_length,
          CGAL::parameters::remesh_boundaries(!protect)
          .number_of_iterations(nb_iter)
          .smooth_constrained_edges(smooth_edges));

      std::cout << "Remeshing done (" << time.elapsed() << " ms)" << std::endl;

      c3t3_item->invalidateOpenGLBuffers();
      c3t3_item->c3t3_changed();
      this->scene->itemChanged(index);

      // default cursor
      QApplication::restoreOverrideCursor();

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      std::cout << "nb_topology_test = " << nb_topology_test << std::endl;
      std::cout << "nb_impossible = " << nb_impossible << std::endl;

      std::cout << "nb_invalid_collapse_short = " << nb_invalid_collapse_short << std::endl;

      std::cout << "nb_orientation_fail_v0 = " << nb_orientation_v0 << std::endl;
      std::cout << "nb_orientation_fail_v1 = " << nb_orientation_v1 << std::endl;
      std::cout << "nb_orientation_fail_midpoint = " << nb_orientation_midpoint << std::endl;

      std::cout << "nb_test_v0 = " << nb_test_v0 << std::endl;
      std::cout << "nb_test_v1 = " << nb_test_v1 << std::endl;
      std::cout << "nb_test_midpoint = " << nb_test_midpoint << std::endl;
      std::cout << std::endl;

      if (nb_test_v0 > 0)
        std::cout << "nb_orientation_v0 / nb_test_v0  = "
        << ((float)nb_orientation_v0 / (float)nb_test_v0) << std::endl;
      if (nb_test_v1 > 1)
        std::cout << "nb_orientation_v1 / nb_test_v1  = "
        << ((float)nb_orientation_v1 / (float)nb_test_v1) << std::endl;
      if (nb_test_midpoint > 0)
        std::cout << "nb_orientation_midpoint / nb_test_midpoint  = "
        << ((float)nb_orientation_midpoint / (float)nb_test_midpoint) << std::endl;

      std::cout << std::endl;
      std::cout << "nb_valid_collapse = " << nb_valid_collapse << std::endl;
      std::cout << "nb_invalid_collapse = " << nb_invalid_collapse << std::endl;
      std::cout << "nb_invalid_lengths = " << nb_invalid_lengths << std::endl;
#endif
    }
    else
    {
      std::cout << "Can't remesh that type of thing" << std::endl;
    }
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;

  Ui::Tetrahedral_remeshing_dialog
  tet_remeshing_dialog(QDialog* dialog,
                       Scene_c3t3_item* c3t3_item)
  {
    Ui::Tetrahedral_remeshing_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    //Set default parameters
    Scene_interface::Bbox bbox = c3t3_item->bbox();
    ui.objectName->setText(c3t3_item->name());

    ui.objectNameSize->setText(
      tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
      .arg(bbox.xmax()-bbox.xmin(), 0, 'g', 3)
      .arg(bbox.ymax()-bbox.ymin(), 0, 'g', 3)
      .arg(bbox.zmax()-bbox.zmin(), 0, 'g', 3));

    double diago_length = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
                                   + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
                                   + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
    double log = std::log10(diago_length);
    unsigned int nb_decimals = (log > 0) ? 5 : (std::ceil(-log)+3);

    ui.edgeLength_dspinbox->setDecimals(nb_decimals);
    ui.edgeLength_dspinbox->setSingleStep(1e-3);
    ui.edgeLength_dspinbox->setRange(1e-6 * diago_length, //min
                                     2.   * diago_length);//max
    ui.edgeLength_dspinbox->setValue(0.05 * diago_length);

    std::ostringstream oss;
    oss << "Diagonal length of the Bbox of the triangulation to remesh is ";
    oss << diago_length << "." << std::endl;
    oss << "Default is 5% of it" << std::endl;
    ui.edgeLength_dspinbox->setToolTip(QString::fromStdString(oss.str()));

    ui.nbIterations_spinbox->setSingleStep(1);
    ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
    ui.nbIterations_spinbox->setValue(1);

    ui.protect_checkbox->setChecked(false);
    ui.smoothEdges_checkBox->setChecked(false);

    connect(ui.protect_checkbox, SIGNAL(toggled(bool)),
            ui.smoothEdges_checkBox, SLOT(setDisabled(bool)));

    return ui;
  }


private:
  QAction* actionTetrahedralRemeshing_;

}; // end Polyhedron_demo_tetrahedral_remeshing_plugin

#include "Tetrahedral_remeshing_plugin.moc"
