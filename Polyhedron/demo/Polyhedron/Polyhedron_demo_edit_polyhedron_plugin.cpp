
#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item.h"
#include <QAction>
#include <QKeySequence>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>
#include <QSettings>
#include <QFileDialog>

#include <QGLViewer/qglviewer.h>

#include <QDockWidget>
#include "ui_Deform_mesh.h"

// #include "MainWindow.h"

class Polyhedron_demo_edit_polyhedron_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  Polyhedron_demo_edit_polyhedron_plugin() 
    : Polyhedron_demo_plugin_helper(), ui_widget(NULL), dock_widget(NULL)
  { }
  ~Polyhedron_demo_edit_polyhedron_plugin()
  { }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const;
  bool applicable() const;

public slots:
  void on_actionDeformation_triggered();
  /////// Dock window signal handlers //////
  // what they do is simply transmiting required 'action' to selected scene_edit_polyhedron_item object
  void on_AddHandlePushButton_clicked();
  void on_PrevHandlePushButton_clicked();
  void on_NextHandlePushButton_clicked();
  void on_SelectAllVerticesPushButton_clicked();
  void on_DeleteHandlePushButton_clicked();  
  void on_ApplyAndClosePushButton_clicked();
  void on_ClearROIPushButton_clicked();
  void on_ShowROICheckBox_stateChanged(int state);
  void on_ShowAsSphereCheckBox_stateChanged(int state);  
  void on_ActivatePivotingCheckBox_stateChanged(int state);
  void on_OverwritePushButton_clicked();
  void on_SaveROIPushButton_clicked();
  void on_ReadROIPushButton_clicked();
  void dock_widget_visibility_changed(bool visible);
  ///////////////////////////////////////////
  void mesh_deformed(Scene_edit_polyhedron_item* edit_item);
  void mesh_repaint_needed(Scene_edit_polyhedron_item* edit_item);
  
  void item_destroyed();
  void new_item_created(int item_id);

private:
  typedef Scene_interface::Item_id Item_id;

  Scene_edit_polyhedron_item* convert_to_edit_polyhedron(Item_id, Scene_polyhedron_item*);
  Scene_polyhedron_item* convert_to_plain_polyhedron(Item_id, Scene_edit_polyhedron_item*);

  Ui::DeformMesh* ui_widget;
  QDockWidget* dock_widget;

  QAction* actionDeformation;
}; // end Polyhedron_demo_edit_polyhedron_plugin

QList<QAction*> Polyhedron_demo_edit_polyhedron_plugin::actions() const {
  return QList<QAction*>() << actionDeformation;
}
bool Polyhedron_demo_edit_polyhedron_plugin::applicable() const { 
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
  this->mw = mainWindow;
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
  ui_widget = new Ui::DeformMesh();

  ui_widget->setupUi(dock_widget); 
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget->AddHandlePushButton, SIGNAL(clicked()), this, SLOT(on_AddHandlePushButton_clicked()));
  connect(ui_widget->PrevHandlePushButton, SIGNAL(clicked()), this, SLOT(on_PrevHandlePushButton_clicked()));
  connect(ui_widget->NextHandlePushButton, SIGNAL(clicked()), this, SLOT(on_NextHandlePushButton_clicked()));
  connect(ui_widget->SelectAllVerticesPushButton, SIGNAL(clicked()), this, SLOT(on_SelectAllVerticesPushButton_clicked()));
  connect(ui_widget->DeleteHandlePushButton, SIGNAL(clicked()), this, SLOT(on_DeleteHandlePushButton_clicked()));
  connect(ui_widget->ApplyAndClosePushButton, SIGNAL(clicked()), this, SLOT(on_ApplyAndClosePushButton_clicked()));
  connect(ui_widget->ClearROIPushButton, SIGNAL(clicked()), this, SLOT(on_ClearROIPushButton_clicked()));
  connect(ui_widget->ShowROICheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ShowROICheckBox_stateChanged(int)));
  connect(ui_widget->ShowAsSphereCheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ShowAsSphereCheckBox_stateChanged(int)));  
  connect(ui_widget->ActivatePivotingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(on_ActivatePivotingCheckBox_stateChanged(int)));
  connect(ui_widget->OverwritePushButton, SIGNAL(clicked()), this, SLOT(on_OverwritePushButton_clicked()));
  
  connect(ui_widget->SaveROIPushButton, SIGNAL(clicked()), this, SLOT(on_SaveROIPushButton_clicked()));
  connect(ui_widget->ReadROIPushButton, SIGNAL(clicked()), this, SLOT(on_ReadROIPushButton_clicked()));
  connect(dock_widget, SIGNAL(visibilityChanged(bool)), this, SLOT(dock_widget_visibility_changed(bool)) );
  ///////////////////////////////////////////////////////////////////

  Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
}

void Polyhedron_demo_edit_polyhedron_plugin::on_actionDeformation_triggered()
{  
  // dock widget should be constructed in init()
  if(dock_widget != NULL)
  {
    dock_widget->show();
  }
}

/////// Dock window signal handlers //////
// what they do is simply transmiting required 'action' to selected scene_edit_polyhedron_item object
void Polyhedron_demo_edit_polyhedron_plugin::on_AddHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->create_handle_group();
}
void Polyhedron_demo_edit_polyhedron_plugin::on_PrevHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->prev_handle_group();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_NextHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->next_handle_group();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_SelectAllVerticesPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->set_all_vertices_as_roi();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_DeleteHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->delete_handle_group();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ClearROIPushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->clear_roi();
  scene->itemChanged(edit_item); // for repaint
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ApplyAndClosePushButton_clicked()
{
  dock_widget->setVisible(false);
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ShowROICheckBox_stateChanged(int state)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_edit_polyhedron_item* edit_item = qobject_cast<Scene_edit_polyhedron_item*>(scene->item(i));
    if(!edit_item) { continue; }
    
    scene->itemChanged(edit_item);  // just for redraw   
  }  
}
void Polyhedron_demo_edit_polyhedron_plugin::on_ShowAsSphereCheckBox_stateChanged(int state)
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
//////////////////////////////////
// slots get called by Scene_edit_polyhedron_item when mesh deformed or repaint needed
void Polyhedron_demo_edit_polyhedron_plugin::mesh_deformed(Scene_edit_polyhedron_item* edit_item)
{
  scene->itemChanged(edit_item); 
}
void Polyhedron_demo_edit_polyhedron_plugin::mesh_repaint_needed(Scene_edit_polyhedron_item* edit_item)
{
  scene->itemChanged(edit_item); 
}
//////////////////////////////////
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
  Scene_edit_polyhedron_item* edit_poly = new Scene_edit_polyhedron_item(poly_item, ui_widget);
  edit_poly->setColor(poly_item->color());
  edit_poly->setName(QString("%1 (edit)").arg(poly_item->name()));

  poly_item->setName(poly_item_name); // Because it is changed when the
                                      // name of edit_poly is changed.

  mw->installEventFilter(edit_poly); // filter mainwindows events for key(pressed/released)

  connect(edit_poly, SIGNAL(destroyed()), this, SLOT(item_destroyed()));
  connect(edit_poly, SIGNAL(mesh_deformed(Scene_edit_polyhedron_item*)), 
    this, SLOT(mesh_deformed(Scene_edit_polyhedron_item*)));
  connect(edit_poly, SIGNAL(mesh_repaint_needed(Scene_edit_polyhedron_item*)), 
    this, SLOT(mesh_repaint_needed(Scene_edit_polyhedron_item*)));

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


void Polyhedron_demo_edit_polyhedron_plugin::item_destroyed() { }

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin, Polyhedron_demo_edit_polyhedron_plugin)

#include "Polyhedron_demo_edit_polyhedron_plugin.moc"
