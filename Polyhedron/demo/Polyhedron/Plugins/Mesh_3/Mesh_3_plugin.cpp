#include "config.h"
#include "config_mesh_3.h"

#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Three.h>
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
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "mesh_3_plugin.json")

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
        (scene->item(scene->mainSelectionIndex())) != nullptr;
    }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
    if(qobject_cast<Scene_implicit_function_item*>
       (scene->item(scene->mainSelectionIndex())) != nullptr) {
      return true;
    }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
    if( qobject_cast<Scene_image_item*>
        (scene->item(scene->mainSelectionIndex())) != nullptr ) {
      return true;
    }
#endif
    for(int ind: scene->selectionIndices()){
      Scene_surface_mesh_item* sm_item
          = qobject_cast<Scene_surface_mesh_item*>(scene->item(ind));
      if(nullptr == sm_item)
        return false;
    }
    return true;
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
  void treat_result(Scene_item& source_item, Scene_c3t3_item* result_item) const;

private:
  QAction* actionMesh_3;
  QAction* actionMesh_3_surface;
  QAction* actionSplitPolylines;
  Messages_interface* msg;
  QMessageBox* message_box_;
  Scene_item* source_item_;
  QString source_item_name_;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;
  bool as_facegraph;
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
  QList<Scene_surface_mesh_item*> sm_items;
  Scene_surface_mesh_item* bounding_sm_item = nullptr;
  Scene_implicit_function_item* function_item = nullptr;
  Scene_image_item* image_item = nullptr;
  Scene_polylines_item* polylines_item = nullptr;

  for(int ind: scene->selectionIndices()) {
    Scene_surface_mesh_item* sm_item =
      qobject_cast<Scene_surface_mesh_item*>(scene->item(ind));
    if(sm_item) {
      sm_items.push_back(sm_item);
      if(is_closed(*sm_item->polyhedron())) {
        bounding_sm_item = sm_item;
      }
    }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
    else if(function_item == nullptr &&
            nullptr !=
            (function_item = qobject_cast<Scene_implicit_function_item*>(scene->item(ind))))
    {}
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
    else if(image_item == nullptr &&
            nullptr != (image_item = qobject_cast<Scene_image_item*>(scene->item(ind))))
    {}
#endif
    else if(polylines_item == nullptr &&
            nullptr != (polylines_item = qobject_cast<Scene_polylines_item*>(scene->item(ind))))
    {}
    else {
      QMessageBox::warning(mw, tr("Mesh_3 plugin"),
                           tr("Wrong selection of items"));
      return;
    }
  }
  Scene_item* item = nullptr;
  const bool more_than_one_item = sm_items.size() > 1;
  bool features_protection_available = false;
  if(!sm_items.empty())
  {
    for(auto sm_item : sm_items) {
      if(nullptr == sm_item->polyhedron()) {
        QApplication::restoreOverrideCursor();
        QMessageBox::critical(mw, tr("Mesh_3 plugin"),
                              tr("ERROR: no data in selected item %1").arg(sm_item->name()));
        return;
      }
      if (!is_triangle_mesh(*sm_item->polyhedron()))
      {
        QApplication::restoreOverrideCursor();
        QMessageBox::warning(mw, tr("Mesh_3 plugin"),
                             tr("Selected Scene_surface_mesh_item %1 is not triangulated.")
                             .arg(sm_item->name()));
        return;
      }
      if(sm_item->getNbIsolatedvertices() != 0)
      {
        QApplication::restoreOverrideCursor();
        QMessageBox::critical(mw, tr(""), tr("ERROR: there are isolated vertices in this mesh."));
        return;
      }
    }
    item = sm_items.front();
    features_protection_available = true;
  }
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (nullptr != function_item) { item = function_item; }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  else if (nullptr != image_item)
  {
    item = image_item;
    features_protection_available = true;

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

  if (nullptr == item)
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

  QString item_name = more_than_one_item ?
    QString("%1...").arg(item->name()) :
    item->name();

  // Set default parameters
  CGAL::Three::Scene_interface::Bbox bbox = item->bbox();
  if(more_than_one_item) {
    for(auto it: sm_items) {
      bbox = bbox + it->bbox();
    }
  }
  ui.objectName->setText(item_name);
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

  ui.facegraphCheckBox->setVisible(surface_only);
  ui.initializationGroup->setVisible(image_item != nullptr && !image_item->isGray());
  ui.grayImgGroup->setVisible(image_item != nullptr && image_item->isGray());
  if (!sm_items.empty())
      ui.volumeGroup->setVisible(!surface_only && nullptr != bounding_sm_item);
  else
    ui.volumeGroup->setVisible(!surface_only);
  if ((!sm_items.empty())|| polylines_item != nullptr) {
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
    if (!sm_items.empty())
    {
      if (surface_only)
      {
        ui.protectEdges->addItem(QString("Sharp and Boundary edges"));
        ui.protectEdges->addItem(QString("Boundary edges only"));
      }
      else
        ui.protectEdges->addItem(QString("Sharp edges"));
    }
    else if(nullptr != image_item)
    {
      if(polylines_item != nullptr)
        ui.protectEdges->addItem(QString("Input polylines"));
      else
      {
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
  as_facegraph = surface_only ? ui.facegraphCheckBox->isChecked() : false;

  Meshing_thread* thread = nullptr;
  if (!sm_items.empty())
  {
    QList<const SMesh*> polyhedrons;
    if(!surface_only) {
      sm_items.removeAll(bounding_sm_item);
    }
    std::transform(sm_items.begin(), sm_items.end(),
                   std::back_inserter(polyhedrons),
                   [](Scene_surface_mesh_item* item) {
                     return item->polyhedron();
                   });
    Scene_polylines_item::Polylines_container plc;
    SMesh *bounding_polyhedron =
      (bounding_sm_item == nullptr) ? nullptr : bounding_sm_item->polyhedron();

    thread =    cgal_code_mesh_3(polyhedrons,
                                 (polylines_item == nullptr)?plc:polylines_item->polylines,
                                 bounding_polyhedron,
                                 item_name,
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
  }
  // Image
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (nullptr != function_item)
  {
    const Implicit_function_interface* pFunction = function_item->function();
    if (nullptr == pFunction)
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
                                 surface_only,
                                 scene);
  }
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  else if (nullptr != image_item)
  {
    const Image* pImage = image_item->image();
    if (nullptr == pImage)
    {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    Scene_polylines_item::Polylines_container plc;

    thread =    cgal_code_mesh_3(pImage,
                                 (polylines_item == nullptr)?plc:polylines_item->polylines,
                                 angle,
                                 facet_sizing,
                                 approx,
                                 tet_sizing,
                                 edge_size,
                                 radius_edge,
                                 protect_features,
                                 manifold,
                                 surface_only,
                                 scene,
                                 detect_connected_components,
                                 image_item->isGray(),
                                 iso_value,
                                 value_outside,
                                 inside_is_less);
  }
#endif

  if ( nullptr == thread )
  {
    QMessageBox::critical(mw,tr(""),tr("ERROR: no thread created"));
    return;
  }

  // Launch thread
  source_item_ = item;
  source_item_name_ = item_name;
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

  QObject::connect(cancelButton, &QAbstractButton::clicked,
                   this, [mesh_thread](){
    mesh_thread->stop();
    mesh_thread->wait();
    QApplication::restoreOverrideCursor(); // restores cursor set in mesh_thread stop() function
  });

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
  if ( nullptr == message_box_ ) { return; }

  message_box_->setInformativeText(str);
}


void
Mesh_3_plugin::
meshing_done(Meshing_thread* thread)
{
  // Print message in console
  QString str = QString("Meshing of \"%1\" done in %2s<br>")
    .arg(source_item_name_)
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

  CGAL::Three::Three::information(qPrintable(str));

  // Treat new c3t3 item
  treat_result(*source_item_, result_item);

  // close message box
  message_box_->done(0);
  message_box_ = nullptr;

  // free memory
  // TODO: maybe there is another way to do that
  delete thread;
}


void
Mesh_3_plugin::
treat_result(Scene_item& source_item,
             Scene_c3t3_item* result_item) const
{
  if(!as_facegraph)
  {
    result_item->setName(tr("%1 [3D Mesh]").arg(source_item_name_));

    result_item->c3t3_changed();

    const Scene_item::Bbox& bbox = result_item->bbox();
    result_item->setPosition(float((bbox.xmin() + bbox.xmax())/2.f),
                            float((bbox.ymin() + bbox.ymax())/2.f),
                            float((bbox.zmin() + bbox.zmax())/2.f));

    result_item->setColor(default_mesh_color);
    result_item->setRenderingMode(source_item.renderingMode());
    result_item->set_data_item(&source_item);

    Q_FOREACH(int ind, scene->selectionIndices()) {
      scene->item(ind)->setVisible(false);
    }
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    scene->itemChanged(index);
    scene->setSelectedItem(-1);
    Scene_interface::Item_id new_item_id = scene->addItem(result_item);
    scene->setSelectedItem(new_item_id);
  }
  else
  {
    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item;
    CGAL::facets_in_complex_3_to_triangle_mesh(result_item->c3t3(), *new_item->face_graph());
    new_item->setName(tr("%1 [Remeshed]").arg(source_item_name_));
    Q_FOREACH(int ind, scene->selectionIndices()) {
      scene->item(ind)->setVisible(false);
    }
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    scene->itemChanged(index);
    scene->setSelectedItem(-1);
    Scene_interface::Item_id new_item_id = scene->addItem(new_item);
    new_item->invalidateOpenGLBuffers();
    new_item->redraw();
    scene->setSelectedItem(new_item_id);
    delete result_item;
  }
}

#include "Mesh_3_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
