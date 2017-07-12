#include "config.h"
#include "config_mesh_3.h"

#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "Messages_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include "Scene_c3t3_item.h"
#include <QInputDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <QDesktopServices>
#include <QUrl>
#include <fstream>

#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
#include "Scene_implicit_function_item.h"
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include "Scene_image_item.h"
#endif

#include "Meshing_thread.h"

#include "ui_Meshing_dialog.h"

using namespace CGAL::Three;

// Constants
const QColor default_mesh_color(45,169,70);

#include "Mesh_3_plugin_cgal_code.h" // declare functions `cgal_code_mesh_3`
#include "split_polylines.h"

class Mesh_3_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* msg_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionMesh_3 = new QAction("Create a Tetrahedral Mesh", mw);
    if(actionMesh_3) {
      actionMesh_3->setProperty("subMenuName", "Tetrahedral Mesh Generation");
      connect(actionMesh_3, SIGNAL(triggered()),
              this, SLOT(mesh_3_volume()));
    }

    actionMesh_3_surface = new QAction("Create a Surface Triangle Mesh", mw);
    if (actionMesh_3_surface){
      actionMesh_3_surface->setProperty("subMenuName", "Tetrahedral Mesh Generation");
      connect(actionMesh_3_surface, SIGNAL(triggered()),
              this, SLOT(mesh_3_surface()));
    }

    actionSplitPolylines = new QAction("Split polylines in a graph", mw);
    actionSplitPolylines->setProperty("subMenuName",
                                      "Tetrahedral Mesh Generation");
    connect(actionSplitPolylines, &QAction::triggered,
            this, &Mesh_3_plugin::splitPolylines);

    this->msg = msg_interface;
  }

  QList<QAction*> actions() const {
    return QList<QAction*>()
      << actionMesh_3
      << actionMesh_3_surface
      << actionSplitPolylines;
  }

  bool applicable(QAction* a) const {
    if(a == actionSplitPolylines) {
      return qobject_cast<Scene_polylines_item*>
        (scene->item(scene->mainSelectionIndex())) != 0;
    }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
    if(qobject_cast<Scene_implicit_function_item*>(scene->item(scene->mainSelectionIndex())) != NULL
      && a == actionMesh_3)
      return true;
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
    Q_FOREACH(int ind, scene->selectionIndices()){
      if( qobject_cast<Scene_image_item*>(scene->item(ind)))
        return true;
    }
#endif  
    Q_FOREACH(int ind, scene->selectionIndices()){
      Scene_polyhedron_item* poly_item
          = qobject_cast<Scene_polyhedron_item*>(scene->item(ind));
      Scene_surface_mesh_item* sm_item
          = qobject_cast<Scene_surface_mesh_item*>(scene->item(ind));
      if (NULL == poly_item)
      {
        if(NULL == sm_item)
          continue;
      }
      if (a == actionMesh_3)
      {
        if(poly_item)
          return poly_item->polyhedron()->is_closed();
        if(sm_item)
          return is_closed(*sm_item->polyhedron());
      }
      else
        return true;
    }
    return false;
  }

public Q_SLOTS:
  void mesh_3_volume();
  void mesh_3_surface();
  void mesh_3_surface_with_defaults() { mesh_3(true, true); }
  void splitPolylines();
  void meshing_done(Meshing_thread* t);
  void status_report(QString str);

private:
  void mesh_3(const bool surface_only, const bool use_defaults = false);
  void launch_thread(Meshing_thread* mesh_thread);
  void treat_result(Scene_item& source_item, Scene_c3t3_item& result_item) const;

private:
  QAction* actionMesh_3;
  QAction* actionMesh_3_surface;
  QAction* actionSplitPolylines;
  Messages_interface* msg;
  QMessageBox* message_box_;
  Scene_item* source_item_;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
}; // end class Mesh_3_plugin

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

void Mesh_3_plugin::splitPolylines() {
  Scene_item* main_item = scene->item(scene->mainSelectionIndex());
  Scene_polylines_item* polylines_item =
    qobject_cast<Scene_polylines_item*>(main_item);
  if(polylines_item == 0) return;

  Scene_polylines_item* new_item = new Scene_polylines_item;
  auto new_polylines = split_polylines(polylines_item->polylines);
  new_item->polylines =
    Polylines_container{new_polylines.begin(), new_polylines.end()};
  new_item->setName(tr("%1 (split)").arg(polylines_item->name()));
  scene->addItem(new_item);
}

void Mesh_3_plugin::mesh_3_surface()
{
  mesh_3(true);
}
void Mesh_3_plugin::mesh_3_volume()
{
  mesh_3(false);
}

void Mesh_3_plugin::mesh_3(const bool surface_only, const bool use_defaults)
{
  Scene_polyhedron_item* poly_item = NULL;
  Scene_polyhedron_item* bounding_poly_item = NULL;
  Scene_surface_mesh_item* sm_item = NULL;
  Scene_surface_mesh_item* bounding_sm_item = NULL;
  Scene_implicit_function_item* function_item = NULL;
  Scene_image_item* image_item = NULL;
  Scene_polylines_item* polylines_item = NULL;

  Q_FOREACH(int ind, scene->selectionIndices()) {

    if(poly_item == NULL)
    {
      poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(ind));
      if (poly_item != NULL
          && scene->selectionIndices().size() == 2
          && bounding_poly_item == NULL)
      {
        bounding_poly_item = qobject_cast<Scene_polyhedron_item*>(
            scene->item(scene->selectionIndices().back()));
        if (bounding_poly_item != NULL)
        {
          //if poly is bounding, and bounding_poly is a non-closed surface
          if (is_closed(*poly_item->polyhedron())
            && !is_closed(*bounding_poly_item->polyhedron()))
          {
            //todo : check poly_item is inside bounding_poly_item
            std::swap(poly_item, bounding_poly_item);
            //now bounding_poly_item is the bounding one
          }
        }
      }
    }
    if(sm_item == NULL)
    {
      sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(ind));
      if (sm_item != NULL
          && scene->selectionIndices().size() == 2
          && bounding_sm_item == NULL)
      {
        bounding_sm_item = qobject_cast<Scene_surface_mesh_item*>(
            scene->item(scene->selectionIndices().back()));
        if (bounding_sm_item != NULL)
        {
          if (is_closed(*sm_item->polyhedron())
            && !is_closed(*bounding_sm_item->polyhedron()))
          {
            //todo : check sm_item is inside bounding_sm_item
            std::swap(sm_item, bounding_sm_item);
            //now bounding_sm_item is the bounding one
          }
        }
      }
    }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
    if(function_item == NULL){
      function_item = qobject_cast<Scene_implicit_function_item*>(scene->item(ind));
    }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
    if(image_item == NULL){
      image_item = qobject_cast<Scene_image_item*>(scene->item(ind));
    }
#endif
    if(polylines_item == NULL){
      polylines_item = qobject_cast<Scene_polylines_item*>(scene->item(ind));
    }
  }
  Scene_item* item = NULL;
  bool features_protection_available = false;
  if(NULL != poly_item)
  {
    if (!poly_item->polyhedron()->is_pure_triangle())
    {
      QMessageBox::warning(mw, tr(""),
                           tr("Selected Scene_polyhedron_item is not triangulated."));
      return;
    }
    if (NULL != bounding_poly_item
      && !bounding_poly_item->polyhedron()->is_pure_triangle())
    {
      QMessageBox::warning(mw, tr(""),
        tr("Selected Scene_polyhedron_item is not triangulated."));
      return;
    }
    item = poly_item;
    features_protection_available = true;
  }
  if(NULL != sm_item)
  {
    if (!is_triangle_mesh(*sm_item->polyhedron()))
    {
      QMessageBox::warning(mw, tr(""),
                           tr("Selected Scene_surface_mesh__item is not triangulated."));
      return;
    }
    item = sm_item;
    features_protection_available = true;
  }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (NULL != function_item) { item = function_item; }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  else if (NULL != image_item)
  {
    item = image_item;
    features_protection_available = (polylines_item != NULL) || !image_item->isGray();

    bool fit_wrdtp = true;
    std::size_t img_wdim = image_item->image()->image()->wdim;
    WORD_KIND img_wordKind = image_item->image()->image()->wordKind;
    //check if the word type fits the hardcoded values in the plugin
    if(image_item->isGray())
    {
      if(img_wordKind != WK_FLOAT)
        fit_wrdtp = false;
      else
        if(img_wdim != 4)
          fit_wrdtp = false;
    }
    else
    {
      if(img_wordKind != WK_FIXED)
        fit_wrdtp = false;
      else
        if(img_wdim != 1)
          fit_wrdtp = false;
    }
    if(!fit_wrdtp)
    {
      QMessageBox::warning(mw, tr(""),
                           tr("Selected object can't be meshed because the image's word type is not supported by this plugin."));
      return;
    }
  }
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

  ui.advanced->setVisible(false);
  connect(ui.facetTopologyLabel,
          &QLabel::linkActivated,
          &QDesktopServices::openUrl);

  dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);
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

  connect(ui.protect, SIGNAL(toggled(bool)),
          ui.noEdgeSizing,   SLOT(setEnabled(bool)));

  connect(ui.protect, SIGNAL(toggled(bool)),
          ui.noEdgeSizing,   SLOT(setChecked(bool)));

  connect(ui.noEdgeSizing, SIGNAL(toggled(bool)),
          ui.edgeLabel,   SLOT(setEnabled(bool)));

  connect(ui.noEdgeSizing, SIGNAL(toggled(bool)),
          ui.edgeSizing,   SLOT(setEnabled(bool)));

  connect(ui.protect, SIGNAL(toggled(bool)),
          ui.sharpEdgesAngle, SLOT(setEnabled(bool)));

  connect(ui.protect, SIGNAL(toggled(bool)),
          ui.sharpEdgesAngleLabel, SLOT(setEnabled(bool)));

  connect(ui.protect, SIGNAL(toggled(bool)),
          ui.protectEdges, SLOT(setEnabled(bool)));

  // Set default parameters
  CGAL::Three::Scene_interface::Bbox bbox = item->bbox();
  ui.objectName->setText(item->name());
  ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                             .arg(bbox.xmax() - bbox.xmin(),0,'g',3)
                             .arg(bbox.ymax() - bbox.ymin(),0,'g',3)
                             .arg(bbox.zmax() - bbox.zmin(),0,'g',3) );

  double diag = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin()) + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin()) + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
  int decimals = 0;
  double sizing_default = get_approximate(diag * 0.05, 2, decimals);
  ui.facetSizing->setDecimals(-decimals+2);
  ui.facetSizing->setSingleStep(std::pow(10.,decimals));
  ui.facetSizing->setRange(diag * 10e-6, // min
                           diag); // max
  ui.facetSizing->setValue(sizing_default); // default value
  ui.edgeSizing->setValue(sizing_default);

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
  ui.protect->setChecked(features_protection_available);
  ui.protectEdges->setEnabled(features_protection_available);

  ui.initializationGroup->setVisible(image_item != NULL && !image_item->isGray());
  ui.grayImgGroup->setVisible(image_item != NULL && image_item->isGray());
  if (poly_item != NULL)
    ui.volumeGroup->setVisible(!surface_only && poly_item->polyhedron()->is_closed());
  else if (sm_item != NULL)
      ui.volumeGroup->setVisible(!surface_only && is_closed(*sm_item->polyhedron()));
  else
    ui.volumeGroup->setVisible(!surface_only);
  if ((poly_item == NULL && sm_item == NULL)|| polylines_item != NULL) {
    ui.sharpEdgesAngleLabel->setVisible(false);
    ui.sharpEdgesAngle->setVisible(false);

    ui.facetTopology->setEnabled(false);
    ui.facetTopology->setToolTip(tr("<b>Notice:</b> "
                                    "This option is only available with a"
                                    " polyhedron or a surface mesh, when features are detected"
                                    " automatically"));
  }
  ui.noEdgeSizing->setChecked(ui.protect->isChecked());
  ui.edgeLabel->setEnabled(ui.noEdgeSizing->isChecked());
  ui.edgeSizing->setEnabled(ui.noEdgeSizing->isChecked());

  if (features_protection_available)
  {
    if (NULL != poly_item || NULL != sm_item)
    {
      if (surface_only)
      {
        ui.protectEdges->addItem(QString("Sharp and Boundary edges"));
        ui.protectEdges->addItem(QString("Boundary edges only"));
      }
      else
        ui.protectEdges->addItem(QString("Sharp edges"));
    }
    else if(NULL != image_item)
    {
      if(polylines_item != NULL)
        ui.protectEdges->addItem(QString("Input polylines"));
      else
      {
        CGAL_assertion(!image_item->isGray());
        ui.protectEdges->addItem(QString("Polylines on cube"));
      }
    }
  }
  // -----------------------------------
  // Get values
  // -----------------------------------

  //reset cursor from the code for the scripts
  QApplication::restoreOverrideCursor();
  if(!use_defaults) {
    int i = dialog.exec();
    if( i == QDialog::Rejected ) { return; }
  }

  // 0 means parameter is not considered
  const double angle = !ui.noAngle->isChecked() ? 0 : ui.facetAngle->value();
  const double approx = !ui.noApprox->isChecked() ? 0 : ui.approx->value();
  const double facet_sizing = !ui.noFacetSizing->isChecked() ? 0 : ui.facetSizing->value();
  const double radius_edge = !ui.noTetShape->isChecked() ? 0 : ui.tetShape->value();
  const double tet_sizing = !ui.noTetSizing->isChecked() ? 0  : ui.tetSizing->value();
  const double edge_size = !ui.noEdgeSizing->isChecked() ? DBL_MAX : ui.edgeSizing->value();
  const bool protect_features = ui.protect->isChecked() && (ui.protectEdges->currentIndex() == 0);
  const bool protect_borders = ui.protect->isChecked() && (ui.protectEdges->currentIndex() == 1);
  const double sharp_edges_angle = ui.sharpEdgesAngle->value();
  const bool detect_connected_components = ui.detectComponents->isChecked();
  const int manifold =
    (ui.manifoldCheckBox->isChecked() ? 1 : 0)
    + (ui.facetTopology->isChecked() ? 2 : 0);
  const float iso_value = float(ui.iso_value_spinBox->value());
  const float value_outside = float(ui.value_outside_spinBox->value());
  const float inside_is_less =  float(ui.inside_is_less_checkBox->isChecked());


  QApplication::setOverrideCursor(Qt::WaitCursor);

  Meshing_thread* thread = NULL;

  // Polyhedron
  if ( NULL != poly_item )
  {
    Polyhedron* pMesh = poly_item->polyhedron();
    if (NULL == pMesh)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    Scene_polylines_item::Polylines_container plc;

    thread =    cgal_code_mesh_3(pMesh,
                                 (polylines_item == NULL)?plc:polylines_item->polylines,
                                 (bounding_poly_item == NULL)?NULL:bounding_poly_item->polyhedron(),
                                 item->name(),
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 edge_size,
                                 radius_edge,
                                 protect_features,
                                 protect_borders,//available only for poly_item
                                 sharp_edges_angle,
                                 manifold,
                                 surface_only,
                                 scene);
  }
  // Surface_mesh
  if ( NULL != sm_item )
  {
    SMesh* psMesh = sm_item->polyhedron();
    if (NULL == psMesh)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }
    typedef CGAL::Graph_with_descriptor_with_graph<SMesh> SMwgd;
    SMwgd *pMesh = new SMwgd(*psMesh);
    Scene_polylines_item::Polylines_container plc;
    SMwgd *pBMesh = (bounding_sm_item == NULL) ? NULL
                    : new SMwgd(*bounding_sm_item->polyhedron());

    thread =    cgal_code_mesh_3(pMesh,
                                 (polylines_item == NULL)?plc:polylines_item->polylines,
                                 pBMesh,
                                 item->name(),
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 edge_size,
                                 radius_edge,
                                 protect_features,
                                 protect_borders,
                                 sharp_edges_angle,
                                 manifold,
                                 surface_only,
                                 scene);
    delete pMesh;
  }
  // Image
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (NULL != function_item)
  {
    const Implicit_function_interface* pFunction = function_item->function();
    if (NULL == pFunction)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    thread =    cgal_code_mesh_3(pFunction,
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 edge_size,
                                 radius_edge,
                                 manifold,
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

    Scene_polylines_item::Polylines_container plc;

    thread =    cgal_code_mesh_3(pImage,
                                 (polylines_item == NULL)?plc:polylines_item->polylines,
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 edge_size,
                                 radius_edge,
                                 protect_features,
                                 manifold,
                                 scene,
                                 detect_connected_components,
                                 image_item->isGray(),
                                 iso_value,
                                 value_outside,
                                 inside_is_less);
  }
#endif

  if ( NULL == thread )
  {
    QMessageBox::critical(mw,tr(""),tr("ERROR: no thread created"));
    return;
  }

  // Launch thread
  source_item_ = item;
  launch_thread(thread);

  QApplication::restoreOverrideCursor();
}

void
Mesh_3_plugin::
launch_thread(Meshing_thread* mesh_thread)
{
  // -----------------------------------
  // Create message box with stop button
  // -----------------------------------
  message_box_ = new QMessageBox(QMessageBox::NoIcon,
                                 "Meshing",
                                 "Mesh generation in progress...",
                                 QMessageBox::Cancel,
                                 mw);

  message_box_->setDefaultButton(QMessageBox::Cancel);
  QAbstractButton* cancelButton = message_box_->button(QMessageBox::Cancel);
  cancelButton->setText(tr("Stop"));

  QObject::connect(cancelButton, SIGNAL(clicked()),
                   mesh_thread,  SLOT(stop()));

  message_box_->open();

  // -----------------------------------
  // Connect main thread to meshing thread
  // -----------------------------------
  QObject::connect(mesh_thread, SIGNAL(done(Meshing_thread*)),
                   this,        SLOT(meshing_done(Meshing_thread*)));

  QObject::connect(mesh_thread, SIGNAL(status_report(QString)),
                   this,        SLOT(status_report(QString)));

  // -----------------------------------
  // Launch mesher
  // -----------------------------------
  mesh_thread->start();
}


void
Mesh_3_plugin::
status_report(QString str)
{
  if ( NULL == message_box_ ) { return; }

  message_box_->setInformativeText(str);
}


void
Mesh_3_plugin::
meshing_done(Meshing_thread* thread)
{
  // Print message in console
  QString str = QString("Meshing of \"%1\" done in %2s<br>")
    .arg(source_item_->name())
    .arg(thread->time());

  Q_FOREACH( QString param, thread->parameters_log() )
  {
    str.append(QString("( %1 )<br>").arg(param));
  }

  Scene_c3t3_item* result_item = thread->item();
  const Scene_item::Bbox& bbox = result_item->bbox();
  str.append(QString("BBox (x,y,z): [ %1, %2 ], [ %3, %4 ], [ %5, %6 ], <br>")
    .arg(bbox.xmin())
    .arg(bbox.xmax())
    .arg(bbox.ymin())
    .arg(bbox.ymax())
    .arg(bbox.zmin())
    .arg(bbox.zmax()));

  msg->information(qPrintable(str));

  // Treat new c3t3 item
  treat_result(*source_item_, *result_item);

  // close message box
  message_box_->done(0);
  message_box_ = NULL;

  // free memory
  // TODO: maybe there is another way to do that
  delete thread;
}


void
Mesh_3_plugin::
treat_result(Scene_item& source_item,
             Scene_c3t3_item& result_item) const
{
  result_item.setName(tr("%1 [3D Mesh]").arg(source_item.name()));

  result_item.c3t3_changed();

  const Scene_item::Bbox& bbox = result_item.bbox();
  result_item.setPosition(float((bbox.xmin() + bbox.xmax())/2.f),
                          float((bbox.ymin() + bbox.ymax())/2.f),
                          float((bbox.zmin() + bbox.zmax())/2.f));

  result_item.setColor(default_mesh_color);
  result_item.setRenderingMode(source_item.renderingMode());
  result_item.set_data_item(&source_item);

  Q_FOREACH(int ind, scene->selectionIndices()) {
    scene->item(ind)->setVisible(false);
  }
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  scene->itemChanged(index);
  scene->setSelectedItem(-1);
  Scene_interface::Item_id new_item_id = scene->addItem(&result_item);
  scene->setSelectedItem(new_item_id);
}

#include "Mesh_3_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
