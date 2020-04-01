#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Viewer_interface.h>
#include "Scene_surface_mesh_item.h"

#include "Scene_edit_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <QAction>
#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include <QShortcut>

#include "ui_Deform_mesh.h"

typedef Scene_surface_mesh_item Scene_facegraph_item;

using namespace CGAL::Three;
class Polyhedron_demo_edit_polyhedron_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  Polyhedron_demo_edit_polyhedron_plugin()
    : Polyhedron_demo_plugin_helper(), dock_widget(NULL)
  { }
  ~Polyhedron_demo_edit_polyhedron_plugin()
  {
    delete e_shortcut;
  }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  virtual void closure()
  {
    dock_widget->hide();
  }
  QList<QAction*> actions() const;
  bool applicable(QAction*) const;

public Q_SLOTS:
  void on_actionDeformation_triggered();
  /////// Dock window signal handlers //////
  // what they do is simply transmiting required 'action' to selected scene_edit_polyhedron_item object
  void on_AddCtrlVertPushButton_clicked();
  void on_PrevCtrlVertPushButton_clicked();
  void on_NextCtrlVertPushButton_clicked();
  void on_SelectAllVerticesPushButton_clicked();
  void on_DeleteCtrlVertPushButton_clicked();
  void on_ApplyAndClosePushButton_clicked();
  void on_DiscardChangesPushButton_clicked();
  void on_ClearROIPushButton_clicked();
  void on_ShowROICheckBox_stateChanged(int state);
  void on_ShowAsSphereCheckBox_stateChanged(int state);
  void on_ActivatePivotingCheckBox_stateChanged(int state);
  void on_ActivateFixedPlaneCheckBox_stateChanged(int state);
  void on_OverwritePushButton_clicked();
  void on_SaveROIPushButton_clicked();
  void on_ReadROIPushButton_clicked();
  void dock_widget_visibility_changed(bool visible);
  void on_Select_isolated_components_button_clicked();
  void on_Get_minimum_button_clicked();

  void on_BrushSpinBoxCtrlVert_changed(int);
  void on_BrushSpinBoxRoi_changed(int);
  void on_ROIRadioButton_toggled(bool);
  void on_importSelectionPushButton_clicked();
  void importSelection(Scene_polyhedron_selection_item*, Scene_edit_polyhedron_item*);
  void dispatchAction();
private:
  typedef CGAL::Three::Scene_interface::Item_id Item_id;
  std::vector<QColor> saved_color;
  Scene_edit_polyhedron_item* convert_to_edit_facegraph(Item_id, Scene_facegraph_item*);
  Scene_facegraph_item* convert_to_plain_facegraph(Item_id, Scene_edit_polyhedron_item*);
  void updateSelectionItems(Scene_facegraph_item* target);

  Ui::DeformMesh ui_widget;
  QDockWidget* dock_widget;

  QAction* actionDeformation;
  RenderingMode last_RM;
  QShortcut* e_shortcut;
}; // end Polyhedron_demo_edit_polyhedron_plugin

QList<QAction*> Polyhedron_demo_edit_polyhedron_plugin::actions() const {
  return QList<QAction*>() << actionDeformation;
}
bool Polyhedron_demo_edit_polyhedron_plugin::applicable(QAction*) const {
  Q_FOREACH(CGAL::Three::Scene_interface::Item_id i, scene->selectionIndices())
  {
    if(qobject_cast<Scene_facegraph_item*>(scene->item(i)))
      return true;
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item)
      return false;
    if (qobject_cast<Scene_facegraph_item*>(edit_item->sm_item()))
    {
      return true;
    }
  }
  return false;
}

void Polyhedron_demo_edit_polyhedron_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  mw = mainWindow;
  scene = scene_interface;
  actionDeformation = new QAction(
          "Surface Mesh Deformation"
        , mw);
  actionDeformation->setProperty("subMenuName", "Triangulated Surface Mesh Deformation");
  actionDeformation->setObjectName("actionDeformation");
  actionDeformation->setShortcutContext(Qt::ApplicationShortcut);
  autoConnectActions();
  e_shortcut = new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_E), mw);
  connect(e_shortcut, &QShortcut::activated, this, &Polyhedron_demo_edit_polyhedron_plugin::dispatchAction);
  connect(e_shortcut, &QShortcut::activatedAmbiguously, this, &Polyhedron_demo_edit_polyhedron_plugin::dispatchAction);

  // Connect Scene::newItem so that, if dock_widget is visible, convert
  // automatically polyhedron items to "edit polyhedron" items.

  ////////////////// Construct widget /////////////////////////////
  // First time, construct docking window
  dock_widget = new QDockWidget(
          "Surface Mesh Deformation"
        , mw);
  dock_widget->setVisible(false); // do not show at the beginning

  ui_widget.setupUi(dock_widget);
  dock_widget->setWindowTitle(tr(
                                  "Surface Mesh Deformation"
                                ));

  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget.AddCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_AddCtrlVertPushButton_clicked()));
  connect(ui_widget.PrevCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_PrevCtrlVertPushButton_clicked()));
  connect(ui_widget.NextCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_NextCtrlVertPushButton_clicked()));
  connect(ui_widget.SelectAllVerticesPushButton, SIGNAL(clicked()), this, SLOT(on_SelectAllVerticesPushButton_clicked()));
  connect(ui_widget.DeleteCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_DeleteCtrlVertPushButton_clicked()));
  connect(ui_widget.ApplyAndClosePushButton, SIGNAL(clicked()), this, SLOT(on_ApplyAndClosePushButton_clicked()));
  connect(ui_widget.DiscardChangesPushButton, SIGNAL(clicked()), this, SLOT(on_DiscardChangesPushButton_clicked()));
  connect(ui_widget.ClearROIPushButton, SIGNAL(clicked()), this, SLOT(on_ClearROIPushButton_clicked()));
  connect(ui_widget.ShowROICheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ShowROICheckBox_stateChanged(int)));
  connect(ui_widget.ShowAsSphereCheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ShowAsSphereCheckBox_stateChanged(int)));
  connect(ui_widget.ActivateFixedPlaneCheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ActivateFixedPlaneCheckBox_stateChanged(int)));
  connect(ui_widget.ActivatePivotingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ActivatePivotingCheckBox_stateChanged(int)));
  connect(ui_widget.OverwritePushButton, SIGNAL(clicked()), this, SLOT(on_OverwritePushButton_clicked()));
  connect(ui_widget.Select_isolated_components_button,  SIGNAL(clicked()), this, SLOT(on_Select_isolated_components_button_clicked()));
  connect(ui_widget.Get_minimum_button,  SIGNAL(clicked()), this, SLOT(on_Get_minimum_button_clicked()));

  connect(ui_widget.SaveROIPushButton, SIGNAL(clicked()), this, SLOT(on_SaveROIPushButton_clicked()));
  connect(ui_widget.ReadROIPushButton, SIGNAL(clicked()), this, SLOT(on_ReadROIPushButton_clicked()));
  connect(dock_widget, SIGNAL(visibilityChanged(bool)), this, SLOT(dock_widget_visibility_changed(bool)) );

  connect(ui_widget.BrushSpinBoxRoi, SIGNAL(valueChanged(int)), this, SLOT(on_BrushSpinBoxRoi_changed(int)));
  connect(ui_widget.BrushSpinBoxCtrlVert, SIGNAL(valueChanged(int)), this, SLOT(on_BrushSpinBoxCtrlVert_changed(int)));
  connect(ui_widget.ROIRadioButton, SIGNAL(toggled(bool)), this, SLOT(on_ROIRadioButton_toggled(bool)));
  connect(ui_widget.importSelectionPushButton, SIGNAL(clicked()), this, SLOT(on_importSelectionPushButton_clicked()));
  ///////////////////////////////////////////////////////////////////
}

void Polyhedron_demo_edit_polyhedron_plugin::on_actionDeformation_triggered()
{
  // dock widget should be constructed in init()
  if(dock_widget->isVisible()) { dock_widget->hide(); }
  else                         { dock_widget->show(); dock_widget->raise();}
}

/////// Dock window signal handlers //////
// what they do is simply transmitting required 'action' to selected scene_edit_polyhedron_item object
void Polyhedron_demo_edit_polyhedron_plugin::on_AddCtrlVertPushButton_clicked()
{

  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->create_ctrl_vertices_group();
}
void Polyhedron_demo_edit_polyhedron_plugin::on_PrevCtrlVertPushButton_clicked()
{
 int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->prev_ctrl_vertices_group();
  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_NextCtrlVertPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->next_ctrl_vertices_group();
  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_SelectAllVerticesPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->set_all_vertices_as_roi();
  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_DeleteCtrlVertPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->delete_ctrl_vertices_group();

  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ClearROIPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->clear_roi();
  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ApplyAndClosePushButton_clicked()
{
  dock_widget->setVisible(false);
}
void Polyhedron_demo_edit_polyhedron_plugin::on_DiscardChangesPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if (!edit_item) return;                             // the selected item is not of the right type

  edit_item->reset_deform_object();
  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item); //for redraw
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ShowROICheckBox_stateChanged(int /*state*/)
{
  for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }
    scene->itemChanged(edit_item);  // just for redraw
  }
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ShowAsSphereCheckBox_stateChanged(int state)
{
  for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }
    if(state == 0)
      edit_item->ShowAsSphere(false);
    else
      edit_item->ShowAsSphere(true);
    scene->itemChanged(edit_item);  // just for redraw
  }
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ActivatePivotingCheckBox_stateChanged(int state)
{
  for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    if(state == Qt::Checked) {
      edit_item->pivoting_begin();
    }
    else {
      edit_item->pivoting_end();
    }
    scene->itemChanged(edit_item);
  }
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ActivateFixedPlaneCheckBox_stateChanged(int)
{
    for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
    {
        Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
        if(!edit_item) { continue; }
        edit_item->update_frame_plane();
        edit_item->invalidateOpenGLBuffers();
        scene->itemChanged(edit_item);
    }
}
void Polyhedron_demo_edit_polyhedron_plugin::on_OverwritePushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->overwrite_deform_object();
}
void Polyhedron_demo_edit_polyhedron_plugin::on_Select_isolated_components_button_clicked() {
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  boost::optional<std::size_t> minimum =
    edit_item->select_isolated_components(ui_widget.Threshold_size_spin_box->value());
  if(minimum) {
    ui_widget.Threshold_size_spin_box->setValue((int) *minimum);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_Get_minimum_button_clicked() {
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  boost::optional<std::size_t> minimum = edit_item->get_minimum_isolated_component();
  if(minimum) {
    ui_widget.Threshold_size_spin_box->setValue((int) *minimum);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_SaveROIPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;

  QString fileName = QFileDialog::getSaveFileName(mw, "Save",
      "roi.txt", "Text (*.txt)");
  if(fileName.isNull()) { return; }

  edit_item->save_roi(fileName.toLocal8Bit().data());
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ReadROIPushButton_clicked()
{
  int item_id = scene->selectionIndices().front();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;

  QString fileName = QFileDialog::getOpenFileName(mw, "Read",
    "roi.txt", "Text (*.txt)");
  if(fileName.isNull()) { return; }

  edit_item->read_roi(fileName.toLocal8Bit().data());
  edit_item->invalidateOpenGLBuffers();
  scene->itemChanged(edit_item);
}


void Polyhedron_demo_edit_polyhedron_plugin::dock_widget_visibility_changed(bool visible)
{
  if(!visible)
  {
    for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
        i < end; ++i)
    {
      Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));

      if(edit_item) {
        edit_item->ShowAsSphere(false);
        Scene_facegraph_item* item = convert_to_plain_facegraph(i, edit_item);
        item->setRenderingMode(last_RM);
        updateSelectionItems(item);
        item->itemChanged();
      }
    }
  }
  else
  {
    ui_widget.ShowAsSphereCheckBox->setChecked(false);
    Scene_polyhedron_selection_item* selection_item = NULL;
    for(int i = 0; i<scene->numberOfEntries(); i++)
    {
      selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(i));
      if (selection_item)
        break;
      else
        selection_item = NULL;
    }
    Q_FOREACH(CGAL::Three::Scene_interface::Item_id i , scene->selectionIndices())
    {
      Scene_facegraph_item* poly_item = qobject_cast<Scene_facegraph_item*>(scene->item(i));
      if (poly_item &&
          CGAL::is_triangle_mesh(*poly_item->face_graph()))
      {
        bool is_valid = true;
        for(boost::graph_traits<Face_graph>::face_descriptor fd : faces(*poly_item->face_graph()))
        {
          if (CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(fd, *poly_item->face_graph()))
          {
            is_valid = false;
            break;
          }
        }
        if(!is_valid)
        {
          QMessageBox::warning(mw,
                               tr("Cannot edit degenerated facegraph_items"),
                               tr(" %1 has degenerated faces, therefore it is not editable.").arg(poly_item->name()));
          break;
        }
        last_RM = poly_item->renderingMode();
        if(!selection_item)
          convert_to_edit_facegraph(i, poly_item);
        else
          importSelection(selection_item, convert_to_edit_facegraph(i, poly_item));
      }
      else if(poly_item &&
              !CGAL::is_triangle_mesh(*poly_item->polyhedron()))
      {
        QMessageBox::warning(mw,
                             tr("Cannot edit non-triangle facegraph_items"),
                             tr(" %1 is not pure-triangle, therefore it is not editable.").arg(poly_item->name()));
      }
    }
  }
}


void Polyhedron_demo_edit_polyhedron_plugin::on_ROIRadioButton_toggled(bool value) {
  int k_ring = value ? ui_widget.BrushSpinBoxRoi->value() :
                       ui_widget.BrushSpinBoxCtrlVert->value();
  for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    edit_item->set_k_ring(k_ring);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_BrushSpinBoxCtrlVert_changed(int value) {
  if(ui_widget.ROIRadioButton->isChecked()) { return; }
  for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    edit_item->set_k_ring(value);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_BrushSpinBoxRoi_changed(int value) {
  if(!ui_widget.ROIRadioButton->isChecked()) { return; }
  for(CGAL::Three::Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    edit_item->set_k_ring(value);
  }
}

Scene_edit_polyhedron_item*
Polyhedron_demo_edit_polyhedron_plugin::convert_to_edit_facegraph(Item_id i,
                           Scene_facegraph_item* poly_item)
{
  QString poly_item_name = poly_item->name();
  Scene_edit_polyhedron_item* edit_poly = new Scene_edit_polyhedron_item(poly_item, &ui_widget, mw);

  edit_poly->setColor(poly_item->color());
  edit_poly->setName(QString("%1 (edit)").arg(poly_item->name()));
  edit_poly->setRenderingMode(Gouraud);
  poly_item->setName(poly_item_name); // Because it is changed when the
                                      // name of edit_poly is changed.
  int k_ring = ui_widget.ROIRadioButton->isChecked() ? ui_widget.BrushSpinBoxRoi->value() :
                                                       ui_widget.BrushSpinBoxCtrlVert->value();
  edit_poly->set_k_ring(k_ring);
  scene->setSelectedItem(-1);
  scene->replaceItem(i, edit_poly);
  scene->setSelectedItem(i);
  return edit_poly;
}

Scene_facegraph_item*
Polyhedron_demo_edit_polyhedron_plugin::convert_to_plain_facegraph(Item_id i,
                            Scene_edit_polyhedron_item* edit_item)
{
  Scene_facegraph_item* poly_item = edit_item->to_sm_item();
  scene->replaceItem(i, poly_item);
  delete edit_item;
  return poly_item;
}


void Polyhedron_demo_edit_polyhedron_plugin::on_importSelectionPushButton_clicked()
{

Scene_polyhedron_selection_item* selection_item = NULL;
Scene_edit_polyhedron_item* edit_item = NULL;
bool need_sel(true), need_edit(true);

// find selection_item and edit_item in selection
  Q_FOREACH(Item_id id, scene->selectionIndices())
  {
    if(need_sel)
    {
      Scene_polyhedron_selection_item* selection_test =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(id));

      if(selection_test)
      {
        selection_item = selection_test;
        need_sel = false;
      }
    }
    if(need_edit)
    {
      Scene_edit_polyhedron_item* edit_test =
          qobject_cast<Scene_edit_polyhedron_item*>(scene->item(id));
      if(edit_test)
      {
        edit_item = edit_test;
        need_edit = false;
      }
    }

    if(!need_sel && !need_edit)
      break;
  }
  if(!selection_item || !edit_item)
    return;
  importSelection(selection_item, edit_item);
}

void Polyhedron_demo_edit_polyhedron_plugin::importSelection(Scene_polyhedron_selection_item *selection_item, Scene_edit_polyhedron_item *edit_item)
{

  //converts the selection in selected points
  QVector<Scene_polyhedron_selection_item::fg_vertex_descriptor> sel_to_import;
  Q_FOREACH(Scene_polyhedron_selection_item::fg_vertex_descriptor vh, selection_item->selected_vertices)
    sel_to_import.push_back(vh);
  Q_FOREACH(Scene_polyhedron_selection_item::fg_edge_descriptor ed, selection_item->selected_edges)
  {
    Scene_polyhedron_selection_item::fg_vertex_descriptor vh = source(halfedge(ed, *selection_item->polyhedron()),*selection_item->polyhedron());
    if(!sel_to_import.contains(vh))
      sel_to_import.push_back(vh);

    vh = target(halfedge(ed, *selection_item->polyhedron()),*selection_item->polyhedron());
    if(!sel_to_import.contains(vh))
      sel_to_import.push_back(vh);
  }

  Q_FOREACH(Scene_polyhedron_selection_item::fg_face_descriptor fh, selection_item->selected_facets)
  {
    CGAL::Halfedge_around_face_circulator<Scene_facegraph_item::Face_graph> hafc(halfedge(fh, *selection_item->polyhedron()), *selection_item->polyhedron());
    CGAL::Halfedge_around_face_circulator<Scene_facegraph_item::Face_graph> end = hafc;
    CGAL_For_all(hafc, end)
    {
      if(!sel_to_import.contains(target(*hafc, *selection_item->polyhedron())))
        sel_to_import.push_back(target(*hafc, *selection_item->polyhedron()));
    }
  }

  //makes the selected points ROI
  Q_FOREACH(Scene_polyhedron_selection_item::fg_vertex_descriptor vh, sel_to_import)
  {
    edit_item->insert_roi_vertex(vh);
  }
  edit_item->invalidateOpenGLBuffers();
  if(selection_item->property("is_highlighting").toBool()){
    selection_item->setProperty("need_hl_restore", true);
    selection_item->set_highlighting(false);
  }
  selection_item->setVisible(false);
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
    v->update();
}

void Polyhedron_demo_edit_polyhedron_plugin::updateSelectionItems(Scene_facegraph_item* target)
{
  for(int i = 0; i<scene->numberOfEntries(); i++)
  {
    Scene_polyhedron_selection_item* sel_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(i));
    if(sel_item
       && sel_item->polyhedron() == target->polyhedron())
    {
      sel_item->invalidateOpenGLBuffers();
      if(!ui_widget.RemeshingCheckBox->isChecked()){
        sel_item->setVisible(true);
        if(sel_item->property("need_hl_restore").toBool()){
          sel_item->set_highlighting(true);
          sel_item->setProperty("need_hl_restore", false);
        }
      }
      else
        scene->erase(scene->item_id(sel_item));
    }
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::dispatchAction()
{
 if(applicable(actionDeformation))
   on_actionDeformation_triggered();
}
#include "Edit_polyhedron_plugin.moc"
