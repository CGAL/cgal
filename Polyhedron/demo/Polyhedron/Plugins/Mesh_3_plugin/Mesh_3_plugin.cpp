#include "config.h"
#include "config_mesh_3.h"

#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include "Scene_c3t3_item.h"
#include <QInputDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <fstream>

#include "Scene_polyhedron_item.h"

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
#include "Scene_implicit_function_item.h"
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include "Scene_segmented_image_item.h"
#endif

#include "ui_Meshing_dialog.h"

// declare the CGAL function
CGAL::Three::Scene_item* cgal_code_mesh_3(const Polyhedron*,
                             QString filename,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing,
                             const double tet_shape,
                             const bool protect_features,
                             CGAL::Three::Scene_interface* scene);

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
Scene_item* cgal_code_mesh_3(const Implicit_function_interface* pfunction,
                             const double facet_angle,
                             const double facet_sizing,
                             const double facet_approx,
                             const double tet_sizing,
                             const double tet_shape,
                             CGAL::Three::Scene_interface* scene);
#endif

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
Scene_item* cgal_code_mesh_3(const Image* pImage,
                             const double facet_angle,
                             const double facet_sizing,
                             const double facet_approx,
                             const double tet_sizing,
                             const double tet_shape,
                             CGAL::Three::Scene_interface* scene);
#endif

using namespace CGAL::Three;
class Polyhedron_demo_mesh_3_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionMesh_3 = new QAction("Create a tetrahedral mesh", mw);
    if(actionMesh_3) {
      actionMesh_3->setProperty("subMenuName", "3D Mesh Generation");
      connect(actionMesh_3, SIGNAL(triggered()),
              this, SLOT(mesh_3()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMesh_3;
  }


  bool applicable(QAction*) const {
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  if(qobject_cast<Scene_implicit_function_item*>(scene->item(scene->mainSelectionIndex())))
    return true;
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  if( qobject_cast<Scene_segmented_image_item*>(scene->item(scene->mainSelectionIndex())))
    return true;
#endif
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
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* poly_item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  Scene_implicit_function_item* function_item =
    qobject_cast<Scene_implicit_function_item*>(scene->item(index));
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  Scene_segmented_image_item* image_item =
    qobject_cast<Scene_segmented_image_item*>(scene->item(index));
#endif

  Scene_item* item = NULL;
  bool features_protection_available = false;
  if (NULL != poly_item)
  {
    item = poly_item;
    features_protection_available = true;
  }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (NULL != function_item) { item = function_item; }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  else if (NULL != image_item)    { item = image_item; }
#endif

  if (NULL == item)
  {
    QMessageBox::warning(mw, tr(""),
      tr("Selected object can't be meshed"));
    return;
  }

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
  CGAL::Three::Scene_interface::Bbox bbox = item->bbox();
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

  ui.protect->setEnabled(features_protection_available);
  if (!features_protection_available)
    ui.protect->setChecked(false);

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
  Scene_item* temp_item = 0;
  if (NULL != poly_item)
  {
    Polyhedron* pMesh = poly_item->polyhedron();
    if (NULL == pMesh)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    temp_item = cgal_code_mesh_3(pMesh,
                                 item->name(),
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 radius_edge,
                                 protect_features,
                                 scene);
  }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (NULL != function_item)
  {
    const Implicit_function_interface* pFunction = function_item->function();
    if (NULL == pFunction)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    temp_item = cgal_code_mesh_3(pFunction,
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 radius_edge,
                                 scene);
  }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  else if (NULL != image_item)
  {
    const Image* pImage = image_item->image();
    if (NULL == pImage)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    temp_item = cgal_code_mesh_3(pImage,
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 radius_edge,
                                 scene);
  }
#endif

  Scene_c3t3_item *result_item = qobject_cast<Scene_c3t3_item*>(temp_item);
  if(result_item) {
      result_item->setName(tr("%1 3d mesh (%2 %3 %4 %5)")
                           .arg(item->name())
                           .arg(angle)
                           .arg(facet_sizing)
                           .arg(tet_sizing)
                           .arg(approx));
      result_item->setItemIsMulticolor(true);
      result_item->setRenderingMode(FlatPlusEdges);
      item->setVisible(false);
      result_item->set_data_item(item);
      scene->itemChanged(index);
      result_item->invalidate_buffers();
      int item_id = scene->addItem(result_item);
      scene->setSelectedItem(item_id);
  }
  QApplication::restoreOverrideCursor();
}

#include "Mesh_3_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
