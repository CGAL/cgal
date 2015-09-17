#include "Polyhedron_demo_plugin_helper.h"
#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item.h"

#include <QAction>
#include <QMainWindow>
#include <QFileDialog>

#include <QGLViewer/qglviewer.h>

#include "ui_Deform_mesh.h"

class Polyhedron_demo_edit_polyhedron_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  Polyhedron_demo_edit_polyhedron_plugin() 
    : Polyhedron_demo_plugin_helper(), dock_widget(NULL)
  { }
  ~Polyhedron_demo_edit_polyhedron_plugin()
  { }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
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
  void on_ClearROIPushButton_clicked();
  void on_ShowROICheckBox_stateChanged(int state);
  void on_ShowAsSphereCheckBox_stateChanged(int state);  
  void on_ActivatePivotingCheckBox_stateChanged(int state);
  void on_OverwritePushButton_clicked();
  void on_SaveROIPushButton_clicked();
  void on_ReadROIPushButton_clicked();
  void dock_widget_visibility_changed(bool visible);
  void on_Select_isolated_components_button_clicked();
  void on_Get_minimum_button_clicked();

  void on_BrushSpinBoxCtrlVert_changed(int);
  void on_BrushSpinBoxRoi_changed(int);
  void on_ROIRadioButton_toggled(bool);
  void new_item_created(int item_id);

private:
  typedef Scene_interface::Item_id Item_id;

  Scene_edit_polyhedron_item* convert_to_edit_polyhedron(Item_id, Scene_polyhedron_item*);
  Scene_polyhedron_item* convert_to_plain_polyhedron(Item_id, Scene_edit_polyhedron_item*);

  Ui::DeformMesh ui_widget;
  QDockWidget* dock_widget;

  QAction* actionDeformation;
}; // end Polyhedron_demo_edit_polyhedron_plugin

QList<QAction*> Polyhedron_demo_edit_polyhedron_plugin::actions() const {
  return QList<QAction*>() << actionDeformation;
}
bool Polyhedron_demo_edit_polyhedron_plugin::applicable(QAction*) const { 
  Q_FOREACH(Scene_interface::Item_id i, scene->selectionIndices())
  {
    if(qobject_cast<Scene_polyhedron_item*>(scene->item(i)) 
        || qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i)))
      return true;
  }
  return false;
}

void Polyhedron_demo_edit_polyhedron_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
  mw = mainWindow;
  scene = scene_interface;

  actionDeformation = new QAction("Surface Mesh Deformation", mw);

  actionDeformation->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_E));
  connect(actionDeformation, SIGNAL(triggered()), this, SLOT(on_actionDeformation_triggered()));

  // Connect Scene::newItem so that, if dock_widget is visible, convert
  // automatically polyhedron items to "edit polyhedron" items.
  QObject* scene = dynamic_cast<QObject*>(scene_interface);
  if(scene) {
    connect(scene, SIGNAL(newItem(int)), this, SLOT(new_item_created(int)));
  } else {
    std::cerr << "ERROR " << __FILE__ << ":" << __LINE__ << " :"
              << " cannot convert scene_interface to scene!\n"; 
  }

  ////////////////// Construct widget /////////////////////////////
  // First time, construct docking window
  dock_widget = new QDockWidget("Mesh Deformation", mw);
  dock_widget->setVisible(false); // do not show at the beginning

  ui_widget.setupUi(dock_widget); 
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget.AddCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_AddCtrlVertPushButton_clicked()));
  connect(ui_widget.PrevCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_PrevCtrlVertPushButton_clicked()));
  connect(ui_widget.NextCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_NextCtrlVertPushButton_clicked()));
  connect(ui_widget.SelectAllVerticesPushButton, SIGNAL(clicked()), this, SLOT(on_SelectAllVerticesPushButton_clicked()));
  connect(ui_widget.DeleteCtrlVertPushButton, SIGNAL(clicked()), this, SLOT(on_DeleteCtrlVertPushButton_clicked()));
  connect(ui_widget.ApplyAndClosePushButton, SIGNAL(clicked()), this, SLOT(on_ApplyAndClosePushButton_clicked()));
  connect(ui_widget.ClearROIPushButton, SIGNAL(clicked()), this, SLOT(on_ClearROIPushButton_clicked()));
  connect(ui_widget.ShowROICheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ShowROICheckBox_stateChanged(int)));
  connect(ui_widget.ShowAsSphereCheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ShowAsSphereCheckBox_stateChanged(int)));  
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
  ///////////////////////////////////////////////////////////////////
}

void Polyhedron_demo_edit_polyhedron_plugin::on_actionDeformation_triggered()
{  
  // dock widget should be constructed in init()
  if(dock_widget->isVisible()) { dock_widget->hide(); }
  else                         { dock_widget->show(); }
}

/////// Dock window signal handlers //////
// what they do is simply transmitting required 'action' to selected scene_edit_polyhedron_item object
void Polyhedron_demo_edit_polyhedron_plugin::on_AddCtrlVertPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->create_ctrl_vertices_group();
}
void Polyhedron_demo_edit_polyhedron_plugin::on_PrevCtrlVertPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->prev_ctrl_vertices_group();
  edit_item->invalidate_buffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_NextCtrlVertPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->next_ctrl_vertices_group();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_SelectAllVerticesPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->set_all_vertices_as_roi();
  edit_item->invalidate_buffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_DeleteCtrlVertPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->delete_ctrl_vertices_group();
  edit_item->invalidate_buffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ClearROIPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->clear_roi();
  edit_item->invalidate_buffers();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ApplyAndClosePushButton_clicked()
{
  dock_widget->setVisible(false);
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ShowROICheckBox_stateChanged(int /*state*/)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }
    
    scene->itemChanged(edit_item);  // just for redraw   
  }  
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ShowAsSphereCheckBox_stateChanged(int /*state*/)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }
    
    scene->itemChanged(edit_item);  // just for redraw   
  }  
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ActivatePivotingCheckBox_stateChanged(int state)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
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
void Polyhedron_demo_edit_polyhedron_plugin::on_OverwritePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->overwrite_deform_object();
}
void Polyhedron_demo_edit_polyhedron_plugin::on_Select_isolated_components_button_clicked() {
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  boost::optional<std::size_t> minimum = 
    edit_item->select_isolated_components(ui_widget.Threshold_size_spin_box->value());
  if(minimum) {
    ui_widget.Threshold_size_spin_box->setValue((int) *minimum);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_Get_minimum_button_clicked() {
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  boost::optional<std::size_t> minimum = edit_item->get_minimum_isolated_component();
  if(minimum) {
    ui_widget.Threshold_size_spin_box->setValue((int) *minimum);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_SaveROIPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;  

  QString fileName = QFileDialog::getSaveFileName(mw, "Save", 
      "roi.txt", "Text (*.txt)");
  if(fileName.isNull()) { return; }

  edit_item->save_roi(fileName.toLocal8Bit().data());  
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ReadROIPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;  

  QString fileName = QFileDialog::getOpenFileName(mw, "Read", 
    "roi.txt", "Text (*.txt)");
  if(fileName.isNull()) { return; }

  edit_item->read_roi(fileName.toLocal8Bit().data());
  edit_item->invalidate_buffers();
  scene->itemChanged(edit_item); 
}
void Polyhedron_demo_edit_polyhedron_plugin::dock_widget_visibility_changed(bool visible)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(i));
    if (poly_item) 
    { poly_item->update_halfedge_indices(); }
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));

    if(visible && poly_item) {
      convert_to_edit_polyhedron(i, poly_item);
    } else if(!visible && edit_item) {
      convert_to_plain_polyhedron(i, edit_item);
    }
  }

  //QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  //if(visible)
  //{
  //  viewer->camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
  //}else
  //{
  //  viewer->camera()->setType(qglviewer::Camera::PERSPECTIVE);
  //}
}


void Polyhedron_demo_edit_polyhedron_plugin::on_ROIRadioButton_toggled(bool value) {
  int k_ring = value ? ui_widget.BrushSpinBoxRoi->value() : 
                       ui_widget.BrushSpinBoxCtrlVert->value();
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    edit_item->set_k_ring(k_ring);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_BrushSpinBoxCtrlVert_changed(int value) {
  if(ui_widget.ROIRadioButton->isChecked()) { return; }
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    edit_item->set_k_ring(value);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::on_BrushSpinBoxRoi_changed(int value) {
  if(!ui_widget.ROIRadioButton->isChecked()) { return; }
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }

    edit_item->set_k_ring(value);
  }
}

void Polyhedron_demo_edit_polyhedron_plugin::new_item_created(int item_id)
{
  if(dock_widget->isVisible()) {
    Scene_polyhedron_item* poly_item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
    if(poly_item) {
      convert_to_edit_polyhedron(item_id, poly_item);
    }
  }
}

Scene_edit_polyhedron_item* 
Polyhedron_demo_edit_polyhedron_plugin::convert_to_edit_polyhedron(Item_id i,
                           Scene_polyhedron_item* poly_item)
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
  scene->replaceItem(i, edit_poly);
  return edit_poly;
}

Scene_polyhedron_item* 
Polyhedron_demo_edit_polyhedron_plugin::convert_to_plain_polyhedron(Item_id i,
                            Scene_edit_polyhedron_item* edit_item) 
{
  Scene_polyhedron_item* poly_item = edit_item->to_polyhedron_item();
  scene->replaceItem(i, poly_item);
  delete edit_item;
  return poly_item;
}

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
