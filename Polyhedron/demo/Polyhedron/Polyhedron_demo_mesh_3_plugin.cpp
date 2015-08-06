#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include "Scene_polyhedron_item.h"
#include <QInputDialog>
#include <QFileDialog>
#include <fstream>
#include "ui_Meshing_dialog.h"

// declare the CGAL function
Scene_item* cgal_code_mesh_3(const Polyhedron*,
                             QString filename,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing,
                             const double tet_shape,
                             const bool protect_features,
                             Scene_interface* scene);

class Polyhedron_demo_mesh_3_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionMesh_3 = new QAction("Create a tetrahedral mesh", mw);
    if(actionMesh_3) {
      connect(actionMesh_3, SIGNAL(triggered()),
              this, SLOT(mesh_3()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMesh_3;
  }


  bool applicable(QAction*) const {
      return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void mesh_3();

private:
  QAction* actionMesh_3;
}; // end class Polyhedron_demo_mesh_3_plugin

double
get_approximate(double d, int precision, int& decimals)
{
    if ( d<0 ) { return 0; }

    double i = std::pow(10.,precision-1);

    decimals = 0;
    while ( d > i*10 ) { d = d/10.; ++decimals; }
    while ( d < i ) { d = d*10.; --decimals; }

    return std::floor(d)*std::pow(10.,decimals);
}

void Polyhedron_demo_mesh_3_plugin::mesh_3()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
  qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(!item) return;

  Polyhedron* pMesh = item->polyhedron();

  if(!pMesh) return;

  // -----------------------------------
  // Create Mesh dialog
  // -----------------------------------
  QDialog dialog(mw);
  Ui::Meshing_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  // Connect checkboxes to spinboxes
  connect(ui.noApprox, SIGNAL(toggled(bool)),
          ui.approx,   SLOT(setEnabled(bool)));

  connect(ui.noFacetSizing, SIGNAL(toggled(bool)),
          ui.facetSizing,   SLOT(setEnabled(bool)));

  connect(ui.noAngle,    SIGNAL(toggled(bool)),
          ui.facetAngle, SLOT(setEnabled(bool)));

  connect(ui.noTetSizing, SIGNAL(toggled(bool)),
          ui.tetSizing,   SLOT(setEnabled(bool)));

  connect(ui.noTetShape, SIGNAL(toggled(bool)),
          ui.tetShape,   SLOT(setEnabled(bool)));

  // Set default parameters
  Scene_interface::Bbox bbox = item->bbox();
  ui.objectName->setText(item->name());
  ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                             .arg(bbox.width(),0,'g',3)
                             .arg(bbox.height(),0,'g',3)
                             .arg(bbox.depth(),0,'g',3) );

  double diag = bbox.diagonal_length();
  int decimals = 0;
  double sizing_default = get_approximate(diag * 0.05, 2, decimals);
  ui.facetSizing->setDecimals(-decimals+2);
  ui.facetSizing->setSingleStep(std::pow(10.,decimals));
  ui.facetSizing->setRange(diag * 10e-6, // min
                           diag); // max
  ui.facetSizing->setValue(sizing_default); // default value

  ui.tetSizing->setDecimals(-decimals+2);
  ui.tetSizing->setSingleStep(std::pow(10.,decimals));
  ui.tetSizing->setRange(diag * 10e-6, // min
                         diag); // max
  ui.tetSizing->setValue(sizing_default); // default value

  double approx_default = get_approximate(diag * 0.005, 2, decimals);
  ui.approx->setDecimals(-decimals+2);
  ui.approx->setSingleStep(std::pow(10.,decimals));
  ui.approx->setRange(diag * 10e-7, // min
                      diag); // max
  ui.approx->setValue(approx_default);

  // -----------------------------------
  // Get values
  // -----------------------------------
  int i = dialog.exec();
  if( i == QDialog::Rejected ) { return; }

  // 0 means parameter is not considered
  const double angle = !ui.noAngle->isChecked() ? 0 : ui.facetAngle->value();
  const double approx = !ui.noApprox->isChecked() ? 0 : ui.approx->value();
  const double facet_sizing = !ui.noFacetSizing->isChecked() ? 0 : ui.facetSizing->value();
  const double radius_edge = !ui.noTetShape->isChecked() ? 0 : ui.tetShape->value();
  const double tet_sizing = !ui.noTetSizing->isChecked() ? 0  : ui.tetSizing->value();
  const bool protect_features = ui.protect->isChecked();

  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_item* result_item = cgal_code_mesh_3(pMesh,
                                             item->name(),
                                             angle,
                                             facet_sizing,
                                             approx,
                                             tet_sizing,
                                             radius_edge,
                                             protect_features,
                                             scene);
  if(result_item) {
      result_item->setName(tr("%1 3d mesh (%2 %3 %4 %5)")
                           .arg(item->name())
                           .arg(angle)
                           .arg(facet_sizing)
                           .arg(tet_sizing)
                           .arg(approx));
      result_item->setColor(Qt::magenta);
      result_item->setRenderingMode(item->renderingMode());
      item->setVisible(false);
      scene->itemChanged(index);
      scene->addItem(result_item);
  }
  QApplication::restoreOverrideCursor();
}

#include "Polyhedron_demo_mesh_3_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
