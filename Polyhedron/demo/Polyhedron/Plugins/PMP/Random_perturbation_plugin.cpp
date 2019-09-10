
#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>

#include <CGAL/Polygon_mesh_processing/random_perturbation.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Random_perturbation_dialog.h"

using namespace CGAL::Three;
class Polyhedron_demo_random_perturbation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionRandomPerturbation_ = new QAction("Random perturbation", mw);
    actionRandomPerturbation_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionRandomPerturbation_) {
      connect(actionRandomPerturbation_, SIGNAL(triggered()),
        this, SLOT(random_perturb()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRandomPerturbation_;
  }

  bool applicable(QAction*) const
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
      return true;
    else if (qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
      return true;
    else
      return false;
  }

public Q_SLOTS:
  void random_perturb()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    //Create dialog box
    QDialog dialog(mw);
    Ui::Random_perturbation_dialog ui
      = perturb_dialog(&dialog, poly_item, selection_item);

    //Get values
    int i = dialog.exec();
    if (i == QDialog::Rejected)
    {
      std::cout << "Perturbation aborted" << std::endl;
      return;
    }

    double max_move = ui.moveSize_dspinbox->value();
    bool project = ui.project_checkbox->isChecked();

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();

    std::cout << "Perturbation..." << std::endl;

    namespace PMP = CGAL::Polygon_mesh_processing;

    if (poly_item)
    {
      Polyhedron& pmesh = *poly_item->polyhedron();
      if(ui.deterministic_checkbox->isChecked())
      {
        unsigned int seed = static_cast<unsigned int>(ui.seed_spinbox->value());
        PMP::random_perturbation(pmesh, max_move,
            PMP::parameters::do_project(project)
            .random_seed(seed));
      }
      else
      {
        PMP::random_perturbation(pmesh, max_move,
            PMP::parameters::do_project(project));
      }

      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }
    else if (selection_item)
    {
      if (selection_item->selected_vertices.empty())
      {
        QMessageBox msg(QMessageBox::Warning,
          QString("Empty selection"),
          QString("Selection of vertices is empty.\nPerturbation aborted."),
          QMessageBox::Ok,
          this->mw);
        msg.exec();
        QApplication::restoreOverrideCursor();
        return;
      }
      Polyhedron& pmesh = *selection_item->polyhedron();
      if (ui.deterministic_checkbox->isChecked())
      {
        unsigned int seed = static_cast<unsigned int>(ui.seed_spinbox->value());
        PMP::random_perturbation(
            selection_item->selected_vertices,
            pmesh,
            max_move,
            PMP::parameters::do_project(project)
            .random_seed(seed));
      }
      else
      {
        std::cout << "selection_item->selected_vertices : "
          << selection_item->selected_vertices.size() << std::endl;
        PMP::random_perturbation(
          selection_item->selected_vertices,
          pmesh,
          max_move,
          PMP::parameters::do_project(project));
      }

      selection_item->invalidateOpenGLBuffers();
      Q_EMIT selection_item->itemChanged();
      selection_item->polyhedron_item()->invalidateOpenGLBuffers();
      Q_EMIT selection_item->polyhedron_item()->itemChanged();
    }
    else
      std::cerr << "Can't perturb that type of item" << std::endl;

    std::cout << " ok (" << time.elapsed() << " ms)" << std::endl;

    // default cursor
    QApplication::restoreOverrideCursor();
  }

  Ui::Random_perturbation_dialog
  perturb_dialog(QDialog* dialog,
                 Scene_polyhedron_item* poly_item,
                 Scene_polyhedron_selection_item* selection_item)
  {
    Ui::Random_perturbation_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    //Set default parameters
    Scene_interface::Bbox bbox = poly_item != NULL ? poly_item->bbox()
      : (selection_item != NULL ? selection_item->bbox()
      : scene->bbox());
    ui.objectName->setText(poly_item != NULL ? poly_item->name()
      : (selection_item != NULL ? selection_item->name()
      : QString("Remeshing parameters")));

    ui.objectNameSize->setText(
      tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
      .arg(bbox.xmax() - bbox.xmin(), 0, 'g', 3)
      .arg(bbox.ymax() - bbox.ymin(), 0, 'g', 3)
      .arg(bbox.zmax() - bbox.zmin(), 0, 'g', 3));

    double diago_length = CGAL::sqrt((bbox.xmax() - bbox.xmin())*(bbox.xmax() - bbox.xmin())
      + (bbox.ymax() - bbox.ymin())*(bbox.ymax() - bbox.ymin())
      + (bbox.zmax() - bbox.zmin())*(bbox.zmax() - bbox.zmin()));
    double log = std::log10(diago_length);
    unsigned int nb_decimals = (log > 0) ? 5 : (std::ceil(-log) + 3);

    //parameters
    ui.moveSize_dspinbox->setDecimals(nb_decimals);
    ui.moveSize_dspinbox->setSingleStep(1e-3);
    ui.moveSize_dspinbox->setRange(0., diago_length);
    ui.moveSize_dspinbox->setValue(0.05 * diago_length);

    ui.project_checkbox->setChecked(true);

    //advanced
    ui.deterministic_checkbox->setChecked(false);
    ui.seed_spinbox->setEnabled(false);
    ui.seed_label->setEnabled(false);
    ui.seed_spinbox->setValue(0);
    ui.seed_spinbox->setMinimum(0);

    connect(ui.deterministic_checkbox, SIGNAL(toggled(bool)),
            ui.seed_spinbox, SLOT(setEnabled(bool)));
    connect(ui.deterministic_checkbox, SIGNAL(toggled(bool)),
            ui.seed_label, SLOT(setEnabled(bool)));

    return ui;
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;
  QAction* actionRandomPerturbation_;

}; // end Polyhedron_demo_random_perturbation_plugin


#include "Random_perturbation_plugin.moc"
