#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_c3t3_item.h"
#include "C3t3_type.h"

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

//#include "ui_Tetrahedral_remeshing_dialog.h"

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

    actionTetrahedralRemeshing_ = new QAction("Tetrehedral Remeshing", mw);
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
//      QDialog dialog(mw);
//      Ui::Isotropic_remeshing_dialog ui
//        = remeshing_dialog(&dialog, poly_item, selection_item);
//
//      // Get values
//      int i = dialog.exec();
//      if (i == QDialog::Rejected)
//      {
//        std::cout << "Remeshing aborted" << std::endl;
//        return;
//      }
//      bool edges_only = ui.splitEdgesOnly_checkbox->isChecked();
//      bool preserve_duplicates = ui.preserveDuplicates_checkbox->isChecked();
//      double target_length = ui.edgeLength_dspinbox->value();
//      unsigned int nb_iter = ui.nbIterations_spinbox->value();
//      unsigned int nb_smooth = ui.nbSmoothing_spinbox->value();
//      bool protect = ui.protect_checkbox->isChecked();
//      bool smooth_features = ui.smooth1D_checkbox->isChecked();

      bool ok;
      double target_edge_length = QInputDialog::getDouble(mw,
        tr("Tetrahedral remeshing"),
        tr("target edge length = "),
        0.1,    //value
        1e-10,  //min
        2147483647,//max
        10,//decimals
        &ok);
      if (!ok)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      time.start();


      Tr& tr = c3t3_item->c3t3().triangulation();


      CGAL::tetrahedral_adaptive_remeshing(tr, target_edge_length);

      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      c3t3_item->invalidateOpenGLBuffers();

      Q_EMIT c3t3_item->itemChanged();

    }
    else
    {
      std::cout << "Can't remesh that type of thing" << std::endl;
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;

  //Ui::Isotropic_remeshing_dialog
  //remeshing_dialog(QDialog* dialog,
  //                 Scene_facegraph_item* poly_item,
  //                 Scene_polyhedron_selection_item* selection_item = NULL)
  //{
  //  Ui::Isotropic_remeshing_dialog ui;
  //  ui.setupUi(dialog);
  //  connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
  //  connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

  //  //connect checkbox to spinbox
  //  connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
  //          ui.nbIterations_spinbox, SLOT(setDisabled(bool)));
  //  connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
  //          ui.protect_checkbox, SLOT(setDisabled(bool)));
  //  connect(ui.protect_checkbox, SIGNAL(toggled(bool)),
  //          ui.smooth1D_checkbox, SLOT(setDisabled(bool)));
  //  connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
  //          ui.smooth1D_checkbox, SLOT(setDisabled(bool)));
  //  connect(ui.preserveDuplicates_checkbox, SIGNAL(toggled(bool)),
  //          ui.protect_checkbox, SLOT(setChecked(bool)));
  //  connect(ui.preserveDuplicates_checkbox, SIGNAL(toggled(bool)),
  //          ui.protect_checkbox, SLOT(setDisabled(bool)));

  //  //Set default parameters
  //  Scene_interface::Bbox bbox = poly_item != NULL ? poly_item->bbox()
  //    : (selection_item != NULL ? selection_item->bbox()
  //      : scene->bbox());
  //  ui.objectName->setText(poly_item != NULL ? poly_item->name()
  //    : (selection_item != NULL ? selection_item->name()
  //      : QString("Remeshing parameters")));

  //  ui.objectNameSize->setText(
  //    tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
  //    .arg(bbox.xmax()-bbox.xmin(), 0, 'g', 3)
  //    .arg(bbox.ymax()-bbox.ymin(), 0, 'g', 3)
  //    .arg(bbox.zmax()-bbox.zmin(), 0, 'g', 3));

  //  double diago_length = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
  //                                 + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
  //                                 + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
  //  double log = std::log10(diago_length);
  //  unsigned int nb_decimals = (log > 0) ? 5 : (std::ceil(-log)+3);

  //  ui.edgeLength_dspinbox->setDecimals(nb_decimals);
  //  ui.edgeLength_dspinbox->setSingleStep(1e-3);
  //  ui.edgeLength_dspinbox->setRange(1e-6 * diago_length, //min
  //                                   2.   * diago_length);//max
  //  ui.edgeLength_dspinbox->setValue(0.05 * diago_length);

  //  std::ostringstream oss;
  //  oss << "Diagonal length of the Bbox of the selection to remesh is ";
  //  oss << diago_length << "." << std::endl;
  //  oss << "Default is 5% of it" << std::endl;
  //  ui.edgeLength_dspinbox->setToolTip(QString::fromStdString(oss.str()));

  //  ui.nbIterations_spinbox->setSingleStep(1);
  //  ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
  //  ui.nbIterations_spinbox->setValue(1);

  //  ui.protect_checkbox->setChecked(false);
  //  ui.smooth1D_checkbox->setChecked(true);

  //  if (NULL != selection_item)
  //  {
  //    //do not preserve duplicates in selection mode
  //    ui.preserveDuplicates_checkbox->setDisabled(true);
  //    ui.preserveDuplicates_checkbox->setChecked(false);
  //  }

  //  return ui;
  //}


private:
  QAction* actionTetrahedralRemeshing_;

}; // end Polyhedron_demo_isotropic_remeshing_plugin

#include "Tetrahedral_remeshing_plugin.moc"
