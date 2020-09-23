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

#include <gsl/pointers>

// Small addition from GSL v2.0.0:
template <class T>
auto make_not_null(T&& t) {
    return gsl::not_null<std::remove_cv_t<std::remove_reference_t<T>>>{std::forward<T>(t)};
}

#include <boost/variant/variant.hpp>
#include <boost/optional/optional.hpp>
#include "Scene_polylines_item.h"

#ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
#include "Scene_implicit_function_item.h"
#endif
#ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
#include "Scene_image_item.h"
#include "Image_type.h"
#endif

#include "Meshing_thread.h"

#include "ui_Meshing_dialog.h"

using namespace CGAL::Three;

// Constants
const QColor default_mesh_color(45,169,70);

#include "Mesh_3_plugin_cgal_code.h" // declare functions `cgal_code_mesh_3`
#include "split_polylines.h"
#include <CGAL/Mesh_facet_topology.h>

class Mesh_3_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "mesh_3_plugin.json")

  Q_PROPERTY(double angle READ get_angle WRITE set_angle);
  Q_PROPERTY(double sharp_edges_angle_bound
             READ get_sharp_edges_angle_bound
             WRITE set_sharp_edges_angle_bound);
  Q_PROPERTY(double edges_sizing READ get_edges_sizing WRITE set_edges_sizing);
  Q_PROPERTY(double facets_sizing READ get_facets_sizing WRITE set_facets_sizing);
  Q_PROPERTY(double approx READ get_approx WRITE set_approx);
  Q_PROPERTY(double tets_sizing READ get_tets_sizing WRITE set_tets_sizing);
  Q_PROPERTY(double tets_shape READ get_tets_shape WRITE set_tets_shape);
  Q_PROPERTY(bool protect_features READ get_protect_features WRITE set_protect_features);
  Q_PROPERTY(bool protect_borders READ get_protect_borders WRITE set_protect_borders);
  Q_PROPERTY(bool manifold_criterion READ get_manifold_criterion WRITE set_manifold_criterion);

  typedef CGAL::Mesh_facet_topology Mesh_facet_topology;
  Q_ENUMS(Mesh_facet_topology)
  Q_PROPERTY(Mesh_facet_topology facet_topology
             READ get_facet_topology
             WRITE set_facet_topology)

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
    return !get_items_or_return_error_string();
  }

public Q_SLOTS:
  boost::optional<QString> get_items_or_return_error_string() const;
  void set_defaults();
  void mesh_3_volume();
  void mesh_3_surface();
  void mesh_3_surface_with_defaults() {
    mesh_3(Mesh_type::SURFACE_ONLY, Dialog_choice::NO_DIALOG);
  }
  void mesh_3_volume_with_defaults() {
    mesh_3(Mesh_type::VOLUME, Dialog_choice::NO_DIALOG);
  }
  void mesh_3(bool with_dialog) { // compatibility with old Qt Scripts
    return mesh_3(
        Mesh_type::VOLUME,
        with_dialog ? Dialog_choice::DIALOG : Dialog_choice::NO_DIALOG);
  }
  void splitPolylines();
  void meshing_done(Meshing_thread* t);
  void status_report(QString str);

public Q_SLOTS:
  void set_angle(const double v) { angle = v; };
  void set_sharp_edges_angle_bound(const double v) {
    sharp_edges_angle_bound = v;
  }
  void set_edges_sizing(const double v) { edges_sizing = v; };
  void set_facets_sizing(const double v) { facets_sizing = v; };
  void set_approx(const double v) { approx = v; };
  void set_tets_sizing(const double v) { tets_sizing = v; };
  void set_tets_shape(const double v) { tets_shape = v; };
  void set_manifold_criterion(const bool v) { manifold_criterion = v; }
  void set_facet_topology(const CGAL::Mesh_facet_topology v) {  facet_topology = v; }
  void set_protect_features(const bool v) { protect_features = v; };
  void set_protect_borders(const bool v) { protect_borders = v; };

  double get_angle() { return angle; };
  double get_sharp_edges_angle_bound() { return sharp_edges_angle_bound; }
  double get_edges_sizing() { return edges_sizing; };
  double get_facets_sizing() { return facets_sizing; };
  double get_approx() { return approx; };
  double get_tets_sizing() { return tets_sizing; };
  double get_tets_shape() { return tets_shape; };
  bool get_manifold_criterion() { return manifold_criterion; };
  CGAL::Mesh_facet_topology get_facet_topology() { return facet_topology; };
  bool get_protect_features() { return protect_features; };
  bool get_protect_borders() { return protect_borders; };


private:
  enum class Mesh_type : bool { VOLUME, SURFACE_ONLY };
  enum class Dialog_choice : bool { NO_DIALOG, DIALOG };
  void mesh_3(const Mesh_type mesh_type, const Dialog_choice dialog = Dialog_choice::DIALOG);
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

  double angle;
  double sharp_edges_angle_bound;
  int sizing_decimals;
  double approx;
  int approx_decimals;
  double edges_sizing;
  double facets_sizing;
  double tets_sizing;
  double tets_shape;
  bool manifold_criterion;
  CGAL::Mesh_facet_topology facet_topology;
  bool protect_features;
  bool protect_borders;

  struct Polyhedral_mesh_items {
    Polyhedral_mesh_items() noexcept
      : sm_items{}, bounding_sm_item(nullptr), polylines_item(nullptr) {}
    QList<gsl::not_null<Scene_surface_mesh_item*>> sm_items;
    Scene_surface_mesh_item* bounding_sm_item;
    Scene_polylines_item* polylines_item;
  };
  struct Image_mesh_items {
    Image_mesh_items(gsl::not_null<Scene_image_item*> ptr) : image_item(ptr) {}
    gsl::not_null<Scene_image_item*> image_item;
    Scene_polylines_item* polylines_item = nullptr;
  };
  struct Implicit_mesh_items {
    gsl::not_null<Scene_implicit_function_item*> function_item;
  };
  enum Item_types {
    POLYHEDRAL_MESH_ITEMS,
    IMAGE_MESH_ITEMS,
    IMPLICIT_MESH_ITEMS
  };
  mutable boost::optional<boost::variant<Polyhedral_mesh_items,
                                         Image_mesh_items,
                                         Implicit_mesh_items>>
      items;
  mutable bool features_protection_available = false;
  mutable Scene_item* item = nullptr;
  mutable CGAL::Three::Scene_interface::Bbox bbox = {};
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
  mesh_3(Mesh_type::SURFACE_ONLY);
}
void Mesh_3_plugin::mesh_3_volume()
{
  mesh_3(Mesh_type::VOLUME);
}

boost::optional<QString> Mesh_3_plugin::get_items_or_return_error_string() const {
  using boost::get;
  items = {};
  features_protection_available = false;
  item = nullptr;
  for (int ind : scene->selectionIndices()) {
    try {
      if (auto sm_item =
              qobject_cast<Scene_surface_mesh_item*>(scene->item(ind))) {
        if (!items) items = Polyhedral_mesh_items{};
        auto& poly_items = get<Polyhedral_mesh_items>(*items);
        auto& sm_items = poly_items.sm_items;
        sm_items.push_back(make_not_null(sm_item));
        if (is_closed(*sm_item->polyhedron())) {
          poly_items.bounding_sm_item = sm_item;
        }
      }
#  ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
      else if (auto function_item = qobject_cast<Scene_implicit_function_item*>(
                   scene->item(ind))) {
        if (!items)
          items = Implicit_mesh_items{make_not_null(function_item)};
        else
          return tr(
              "An implicit function cannot be mixed with other items type");
      }
#  endif
#  ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
      else if (auto image_item =
                   qobject_cast<Scene_image_item*>(scene->item(ind))) {
        if (!items)
          items = Image_mesh_items{make_not_null(image_item)};
        else
          return tr("An image items cannot be mixed with other items type");
      }
#  endif
      else if (auto polylines_item =
                   qobject_cast<Scene_polylines_item*>(scene->item(ind))) {
        if (!items) items = Polyhedral_mesh_items{};
        auto poly_items_ptr = get<Polyhedral_mesh_items>(&*items);
        if(poly_items_ptr) {
          if (poly_items_ptr->polylines_item) {
            return tr("Only one polyline item is accepted");
          } else {
            poly_items_ptr->polylines_item = polylines_item;
          }
        } else {
          auto image_items = get<Image_mesh_items>(*items);
          if (image_items.polylines_item) {
            return tr("Only one polyline item is accepted");
          } else {
            image_items.polylines_item = polylines_item;
          }
        }
      } else {
        return tr("Wrong selection of items");
      }
    } catch (const boost::bad_get&) { return tr("Wrong selection of items"); }
  } // end for loop on selected items
  if (!items) { return tr("Selected objects can't be meshed"); }
  item = nullptr;
  features_protection_available = false;
  if (auto poly_items = get<Polyhedral_mesh_items>(&*items)) {
    auto& sm_items = poly_items->sm_items;
    for (auto sm_item : sm_items) {
      if (nullptr == sm_item->polyhedron()) {
        return tr("ERROR: no data in selected item %1").arg(sm_item->name());
      }
      if (!is_triangle_mesh(*sm_item->polyhedron())) {
        return tr("Selected Scene_surface_mesh_item %1 is not triangulated.")
            .arg(sm_item->name());
      }
      if (sm_item->getNbIsolatedvertices() != 0) {
        return tr("ERROR: there are isolated vertices in this mesh.");
      }
    }
    if (!sm_items.empty()) item = sm_items.front();
    features_protection_available = true;
  }
#  ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  else if (auto implicit_mesh_items = get<Implicit_mesh_items>(&*items)) {
    item = implicit_mesh_items->function_item;
  }
#  endif
#  ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  else if (auto image_mesh_items = get<Image_mesh_items>(&*items)) {
    auto& image_item = image_mesh_items->image_item;
    item = image_item;
    features_protection_available = true;

    bool fit_wrdtp = true;
    std::size_t img_wdim = image_item->image()->image()->wdim;
    WORD_KIND img_wordKind = image_item->image()->image()->wordKind;
    // check if the word type fits the hardcoded values in the plugin
    if (image_item->isGray()) {
      if (img_wordKind != WK_FLOAT)
        fit_wrdtp = false;
      else if (img_wdim != 4)
        fit_wrdtp = false;
    } else {
      if (img_wordKind != WK_FIXED)
        fit_wrdtp = false;
      else if (img_wdim != 1)
        fit_wrdtp = false;
    }
    if (!fit_wrdtp) {
      return tr(
          "Selected object can't be meshed because the image's word type is "
          "not supported by this plugin.");
    }
  }
#  endif

  if(item) {
    bbox = item->bbox();
    if (auto poly_items = get<Polyhedral_mesh_items>(&*items)) {
      for (auto it : poly_items->sm_items) {
        bbox = bbox + it->bbox();
      }
      if (poly_items->polylines_item)
        bbox = bbox + poly_items->polylines_item->bbox();
    }
  }
  return {};
}

void Mesh_3_plugin::set_defaults() {
  auto error = get_items_or_return_error_string();
  if(error) return;
  double diag = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin()) + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin()) + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
  facets_sizing = get_approximate(diag * 0.05, 2, sizing_decimals);
  edges_sizing = facets_sizing;
  tets_sizing = facets_sizing;
  angle = 25.;
  sharp_edges_angle_bound = 60.;
  approx = get_approximate(diag * 0.005, 2, approx_decimals);
}

void Mesh_3_plugin::mesh_3(const Mesh_type mesh_type,
                           const Dialog_choice dialog_choice) {
  CGAL_assertion(static_cast<bool>(items));
  auto error_string = get_items_or_return_error_string();
  if (error_string) {
    QApplication::restoreOverrideCursor();
    QMessageBox::warning(mw, tr("Mesh_3 plugin"), *error_string);
    return;
  }
  using boost::get;
  const bool more_than_one_item =
      get<Polyhedral_mesh_items>(&*items) &&
      (get<Polyhedral_mesh_items>(&*items)->sm_items.size() > 1);

  Scene_image_item* image_item =
      get<Image_mesh_items>(&*items)
          ? get<Image_mesh_items>(&*items)->image_item.get()
          : nullptr;
  Scene_surface_mesh_item* bounding_sm_item =
      get<Polyhedral_mesh_items>(&*items)
          ? get<Polyhedral_mesh_items>(&*items)->bounding_sm_item
          : nullptr;
  Scene_polylines_item* polylines_item =
      get<Polyhedral_mesh_items>(&*items)
          ? get<Polyhedral_mesh_items>(&*items)->polylines_item
          : nullptr;
  Scene_implicit_function_item* function_item =
      get<Implicit_mesh_items>(&*items)
          ? get<Implicit_mesh_items>(&*items)->function_item.get()
          : nullptr;
  // -----------------------------------
  // Create Mesh dialog
  // -----------------------------------
  QDialog dialog(mw);
  Ui::Meshing_dialog ui;
  ui.setupUi(&dialog);

  ui.facetAngle->setRange(0.0, 30.0);
  ui.facetAngle->setValue(25.0);
  ui.edgeSizing->setMinimum(0.0);
  ui.sharpEdgesAngle->setMaximum(180);
  ui.iso_value_spinBox->setRange(-65536.0, 65536.0);
  ui.tetShape->setMinimum(1.0);

  ui.advanced->setVisible(false);
  connect(ui.facetTopologyLabel,
          &QLabel::linkActivated,
          &QDesktopServices::openUrl);

  dialog.setWindowFlags(Qt::Dialog | Qt::CustomizeWindowHint |
                        Qt::WindowCloseButtonHint);
  connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

  // Connect checkboxes to spinboxes
  connect(
      ui.noApprox, SIGNAL(toggled(bool)), ui.approx, SLOT(setEnabled(bool)));

  connect(ui.noFacetSizing,
          SIGNAL(toggled(bool)),
          ui.facetSizing,
          SLOT(setEnabled(bool)));

  connect(
      ui.noAngle, SIGNAL(toggled(bool)), ui.facetAngle, SLOT(setEnabled(bool)));

  connect(ui.noTetSizing,
          SIGNAL(toggled(bool)),
          ui.tetSizing,
          SLOT(setEnabled(bool)));

  connect(ui.noTetShape,
          SIGNAL(toggled(bool)),
          ui.tetShape,
          SLOT(setEnabled(bool)));

  connect(ui.protect,
          SIGNAL(toggled(bool)),
          ui.noEdgeSizing,
          SLOT(setEnabled(bool)));

  connect(ui.protect,
          SIGNAL(toggled(bool)),
          ui.noEdgeSizing,
          SLOT(setChecked(bool)));

  connect(ui.noEdgeSizing,
          SIGNAL(toggled(bool)),
          ui.edgeLabel,
          SLOT(setEnabled(bool)));

  connect(ui.noEdgeSizing,
          SIGNAL(toggled(bool)),
          ui.edgeSizing,
          SLOT(setEnabled(bool)));

  connect(ui.protect,
          SIGNAL(toggled(bool)),
          ui.sharpEdgesAngle,
          SLOT(setEnabled(bool)));

  connect(ui.protect,
          SIGNAL(toggled(bool)),
          ui.sharpEdgesAngleLabel,
          SLOT(setEnabled(bool)));

  connect(ui.protect,
          SIGNAL(toggled(bool)),
          ui.protectEdges,
          SLOT(setEnabled(bool)));

  QString item_name =
      more_than_one_item ? QString("%1...").arg(item->name()) : item->name();

  ui.objectName->setText(item_name);
  ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                             .arg(bbox.xmax() - bbox.xmin(),0,'g',3)
                             .arg(bbox.ymax() - bbox.ymin(),0,'g',3)
                             .arg(bbox.zmax() - bbox.zmin(),0,'g',3) );

  set_defaults();
  double diag = CGAL::sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin()) + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin()) + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
  ui.facetSizing->setRange(diag * 10e-6, // min
                           diag); // max
  ui.facetSizing->setValue(facets_sizing);
  ui.edgeSizing->setValue(edges_sizing);

  ui.tetSizing->setRange(diag * 10e-6, // min
                         diag);        // max
  ui.tetSizing->setValue(tets_sizing); // default value

  ui.approx->setRange(diag * 10e-7, // min
                      diag);        // max
  ui.approx->setValue(approx);

  ui.protect->setEnabled(features_protection_available);
  ui.protect->setChecked(features_protection_available);
  ui.protectEdges->setEnabled(features_protection_available);

  ui.facegraphCheckBox->setVisible(mesh_type == Mesh_type::SURFACE_ONLY);
  ui.initializationGroup->setVisible(image_item != nullptr &&
                                     !image_item->isGray());
  ui.grayImgGroup->setVisible(image_item != nullptr && image_item->isGray());
  if (items->which() == POLYHEDRAL_MESH_ITEMS)
    ui.volumeGroup->setVisible(mesh_type == Mesh_type::VOLUME &&
                               nullptr != bounding_sm_item);
  else
    ui.volumeGroup->setVisible(mesh_type == Mesh_type::VOLUME);
  ui.sharpEdgesAngle->setValue(sharp_edges_angle_bound);
  if (items->which() != POLYHEDRAL_MESH_ITEMS || polylines_item != nullptr) {
    ui.sharpEdgesAngleLabel->setVisible(false);
    ui.sharpEdgesAngle->setVisible(false);

    ui.facetTopology->setEnabled(false);
    ui.facetTopology->setToolTip(
        tr("<b>Notice:</b> "
           "This option is only available with a"
           " polyhedron or a surface mesh, when features are detected"
           " automatically"));
  }
  ui.noEdgeSizing->setChecked(ui.protect->isChecked());
  ui.edgeLabel->setEnabled(ui.noEdgeSizing->isChecked());
  ui.edgeSizing->setEnabled(ui.noEdgeSizing->isChecked());

  if (features_protection_available) {
    if (items->which() == POLYHEDRAL_MESH_ITEMS) {
      if (mesh_type == Mesh_type::SURFACE_ONLY) {
        ui.protectEdges->addItem(QString("Sharp and Boundary edges"));
        ui.protectEdges->addItem(QString("Boundary edges only"));
      } else
        ui.protectEdges->addItem(QString("Sharp edges"));
    } else if (items->which() == IMAGE_MESH_ITEMS) {
      if (polylines_item != nullptr)
        ui.protectEdges->addItem(QString("Input polylines"));
      else {
        ui.protectEdges->addItem(QString("Polylines on cube"));
      }
    }
  }
  // -----------------------------------
  // Get values
  // -----------------------------------

  // reset cursor from the code for the scripts
  QApplication::restoreOverrideCursor();
  if (dialog_choice == Dialog_choice::DIALOG) {
    int i = dialog.exec();
    if (i == QDialog::Rejected) { return; }
  }

  // 0 means parameter is not considered
  angle = !ui.noAngle->isChecked() ? 0 : ui.facetAngle->value();
  sharp_edges_angle_bound = ui.sharpEdgesAngle->value();
  std::cerr << "sharp_edges_angle_bound: " << sharp_edges_angle_bound << '\n';
  edges_sizing =
      !ui.noEdgeSizing->isChecked() ? DBL_MAX : ui.edgeSizing->value();
  facets_sizing = !ui.noFacetSizing->isChecked() ? 0 : ui.facetSizing->value();
  approx = !ui.noApprox->isChecked() ? 0 : ui.approx->value();
  tets_shape = !ui.noTetShape->isChecked() ? 0 : ui.tetShape->value();
  tets_sizing = !ui.noTetSizing->isChecked() ? 0 : ui.tetSizing->value();
  protect_features =
      ui.protect->isChecked() && (ui.protectEdges->currentIndex() == 0);
  protect_borders =
      ui.protect->isChecked() && (ui.protectEdges->currentIndex() == 1);
  const bool detect_connected_components = ui.detectComponents->isChecked();
  const int manifold = (ui.manifoldCheckBox->isChecked() ? 1 : 0) +
                       (ui.facetTopology->isChecked() ? 2 : 0);
  const float iso_value = float(ui.iso_value_spinBox->value());
  const float value_outside = float(ui.value_outside_spinBox->value());
  const bool inside_is_less = ui.inside_is_less_checkBox->isChecked();
  as_facegraph = (mesh_type == Mesh_type::SURFACE_ONLY)
                     ? ui.facegraphCheckBox->isChecked()
                     : false;

  Meshing_thread* thread = nullptr;
  switch (items->which()) {
  case POLYHEDRAL_MESH_ITEMS: {
    auto& poly_items = get<Polyhedral_mesh_items>(*items);
    auto& sm_items = poly_items.sm_items;
    const auto bounding_sm_item = poly_items.bounding_sm_item;
    const auto polylines_item = poly_items.polylines_item;
    QList<const SMesh*> polyhedrons;
    if(mesh_type != Mesh_type::SURFACE_ONLY) {
      sm_items.removeAll(make_not_null(bounding_sm_item));
    }
    std::transform(sm_items.begin(), sm_items.end(),
                   std::back_inserter(polyhedrons),
                   [](Scene_surface_mesh_item* item) {
                     return item->polyhedron();
                   });
    Scene_polylines_item::Polylines_container plc;
    SMesh* bounding_polyhedron = (bounding_sm_item == nullptr)
                                     ? nullptr
                                     : bounding_sm_item->polyhedron();

    thread = cgal_code_mesh_3(
        polyhedrons,
        (polylines_item == nullptr) ? plc : polylines_item->polylines,
        bounding_polyhedron,
        item_name,
        angle,
        facets_sizing,
        approx,
        tets_sizing,
        edges_sizing,
        tets_shape,
        protect_features,
        protect_borders,
        sharp_edges_angle_bound,
        manifold,
        mesh_type == Mesh_type::SURFACE_ONLY);
    break;
  }
  // Image
#  ifdef CGAL_MESH_3_DEMO_ACTIVATE_IMPLICIT_FUNCTIONS
  case IMPLICIT_MESH_ITEMS: {
    const Implicit_function_interface* pFunction = function_item->function();
    if (nullptr == pFunction) {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    thread = cgal_code_mesh_3(pFunction,
                              angle,
                              facets_sizing,
                              approx,
                              tets_sizing,
                              edges_sizing,
                              tets_shape,
                              manifold,
                              mesh_type == Mesh_type::SURFACE_ONLY);
    break;
  }
#  endif
#  ifdef CGAL_MESH_3_DEMO_ACTIVATE_SEGMENTED_IMAGES
  case IMAGE_MESH_ITEMS: {
    const Image* pImage = image_item->image();
    if (nullptr == pImage) {
      QMessageBox::critical(mw, tr(""), tr("ERROR: no data in selected item"));
      return;
    }

    Scene_polylines_item::Polylines_container plc;

    thread = cgal_code_mesh_3(
        pImage,
        (polylines_item == nullptr) ? plc : polylines_item->polylines,
        angle,
        facets_sizing,
        approx,
        tets_sizing,
        edges_sizing,
        tets_shape,
        protect_features,
        manifold,
        mesh_type == Mesh_type::SURFACE_ONLY,
        detect_connected_components,
        image_item->isGray(),
        iso_value,
        value_outside,
        inside_is_less);
    break;
  }
  default:
    CGAL::Three::Three::error(tr("Mesh_3 plugin"),
                              tr("This type of item is not handled!"));
    return;
  } // end switch
#  endif

  if (nullptr == thread) {
    QMessageBox::critical(mw, tr(""), tr("ERROR: no thread created"));
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
