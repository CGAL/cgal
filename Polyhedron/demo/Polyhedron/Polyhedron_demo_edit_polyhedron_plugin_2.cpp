
#include "Polyhedron_demo_plugin_helper.h"

#include "Scene_polyhedron_item.h"
#include "Scene_edit_polyhedron_item_2.h"
#include <QAction>
#include <QKeySequence>
#include <QMainWindow>
#include <QMessageBox>
#include <QInputDialog>
#include <QSettings>

#include <QGLViewer/qglviewer.h>

#include <QDockWidget>
#include "ui_Deform_mesh_2.h"

class Polyhedron_demo_edit_polyhedron_plugin_2 : 
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  Polyhedron_demo_edit_polyhedron_plugin_2() 
    : Polyhedron_demo_plugin_helper(), size(0), edit_mode(false),
      dock_widget(NULL), ui_widget(NULL)
  {}

  ~Polyhedron_demo_edit_polyhedron_plugin_2();

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const;
  bool applicable() const;

  //! Applicable if any of the currently selected items is either a Polyhedron or an Edit polyhedron

public slots:
  void on_actionDeformation_triggered();
  void on_actionToggleEdit_triggered(bool);

  void on_AddHandlePushButton_clicked();
  void on_PrevHandlePushButton_clicked();
  void on_NextHandlePushButton_clicked();

  void item_destroyed();
  void new_item_created(int item_id);

private:
  typedef Scene_interface::Item_id Item_id;

  Scene_edit_polyhedron_item_2* 
  convert_to_edit_polyhedron(Item_id, Scene_polyhedron_item*);

  Scene_polyhedron_item* 
  convert_to_plain_polyhedron(Item_id, Scene_edit_polyhedron_item_2*);

  Ui::DeformMesh_2* ui_widget;
  QDockWidget* dock_widget;

  QAction* actionToggleEdit;
  QAction* actionDeformation;

  int size;
  bool edit_mode;
}; // end Polyhedron_demo_edit_polyhedron_plugin_2

Polyhedron_demo_edit_polyhedron_plugin_2::~Polyhedron_demo_edit_polyhedron_plugin_2()
{
  // IOY: note sure what it is doing but it constantly throw ex when I close the window

  //QSettings settings;
  //settings.beginGroup("Polyhedron edition");
  //settings.setValue("Deform_mesh widget area", 
  //                  this->mw->dockWidgetArea(widget));
  //settings.endGroup();
}

QList<QAction*> Polyhedron_demo_edit_polyhedron_plugin_2::actions() const {
  return QList<QAction*>() << actionDeformation;
  //return QList<QAction*>() << actionToggleEdit;
}

bool Polyhedron_demo_edit_polyhedron_plugin_2::applicable() const { 
  Q_FOREACH(Scene_interface::Item_id i, scene->selectionIndices())
  {
    if(qobject_cast<Scene_polyhedron_item*>(scene->item(i)) 
        || qobject_cast<Scene_edit_polyhedron_item_2*>(scene->item(i)))
      return true;
  }
  return false;
}

void Polyhedron_demo_edit_polyhedron_plugin_2::init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
  this->mw = mainWindow;
  actionDeformation = new QAction("Surface Mesh Deformation 2", mw);
  connect(actionDeformation, SIGNAL(triggered()), this, SLOT(on_actionDeformation_triggered()));
  
  actionToggleEdit = new QAction(tr("Toggle &edition of item(s)"), mainWindow);
  actionToggleEdit->setObjectName("actionToggleEdit");
  actionToggleEdit->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_E));
  actionToggleEdit->setCheckable(true);

  // Connect Scene::newItem so that, if edit_mode==true, convert
  // automatically polyhedron items to "edit polyhedron" items.
  QObject* scene = dynamic_cast<QObject*>(scene_interface);
  if(scene) {
    connect(scene, SIGNAL(newItem(int)),
            this, SLOT(new_item_created(int)));
  } else {
    std::cerr << "ERROR " << __FILE__ << ":" << __LINE__ << " :"
              << " cannot convert scene_interface to scene!\n"; 
  }

  ////////////////// Construct widget /////////////////////////////
  // First time, construct docking window
  dock_widget = new QDockWidget("Mesh Deformation", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  ui_widget = new Ui::DeformMesh_2();

  ui_widget->setupUi(dock_widget); 
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);
    
  // bind states of actionToggleEdit and editModeCb
  connect(actionToggleEdit, SIGNAL(triggered(bool)), ui_widget->editModeCb, SLOT(setChecked(bool)));
  connect(ui_widget->editModeCb, SIGNAL(clicked(bool)), actionToggleEdit, SLOT(setChecked(bool)));

  // make editModeCb actually trigger the slot
  connect(ui_widget->editModeCb, SIGNAL(clicked(bool)), this, SLOT(on_actionToggleEdit_triggered(bool)));

  connect(ui_widget->AddHandlePushButton, SIGNAL(clicked()), this, SLOT(on_AddHandlePushButton_clicked()));
  connect(ui_widget->PrevHandlePushButton, SIGNAL(clicked()), this, SLOT(on_PrevHandlePushButton_clicked()));
  connect(ui_widget->NextHandlePushButton, SIGNAL(clicked()), this, SLOT(on_NextHandlePushButton_clicked()));
  ///////////////////////////////////////////////////////////////////

  Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
}

void Polyhedron_demo_edit_polyhedron_plugin_2::on_actionDeformation_triggered()
{  
  // dock widget should be constructed in init()
  if(dock_widget != NULL)
  {
    dock_widget->show();
  }
}

void Polyhedron_demo_edit_polyhedron_plugin_2::on_AddHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item_2* edit_item = qobject_cast<Scene_edit_polyhedron_item_2*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->create_handle_group();
}

void Polyhedron_demo_edit_polyhedron_plugin_2::on_PrevHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item_2* edit_item = qobject_cast<Scene_edit_polyhedron_item_2*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->prev_handle_group();
  scene->itemChanged(edit_item); // for repaint
}

void Polyhedron_demo_edit_polyhedron_plugin_2::on_NextHandlePushButton_clicked()
{
  int item_id = scene->mainSelectionIndex();
  Scene_edit_polyhedron_item_2* edit_item = qobject_cast<Scene_edit_polyhedron_item_2*>(scene->item(item_id));
  if(!edit_item) return;                             // the selected item is not of the right type

  edit_item->next_handle_group();
  scene->itemChanged(edit_item); // for repaint
}

void Polyhedron_demo_edit_polyhedron_plugin_2::new_item_created(int item_id)
{
  if(edit_mode) {
    Scene_polyhedron_item* poly_item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
    if(poly_item) {
      convert_to_edit_polyhedron(item_id, poly_item);
    }
  }
}

Scene_edit_polyhedron_item_2* 
Polyhedron_demo_edit_polyhedron_plugin_2::convert_to_edit_polyhedron(Item_id i,
                           Scene_polyhedron_item* poly_item)
{
  QString poly_item_name = poly_item->name();
  Scene_edit_polyhedron_item_2* edit_poly = 
    new Scene_edit_polyhedron_item_2(poly_item);
  edit_poly->setColor(poly_item->color());
  edit_poly->setName(QString("%1 (edit)").arg(poly_item->name()));

  poly_item->setName(poly_item_name); // Because it is changed when the
                                      // name of edit_poly is changed.

  connect(edit_poly, SIGNAL(destroyed()), this, SLOT(item_destroyed()));

  edit_poly->ui_widget = ui_widget;

  scene->replaceItem(i, edit_poly);
  return edit_poly;
}

Scene_polyhedron_item* 
Polyhedron_demo_edit_polyhedron_plugin_2::convert_to_plain_polyhedron(Item_id i,
                            Scene_edit_polyhedron_item_2* edit_item) 
{
  Scene_polyhedron_item* poly_item = edit_item->to_polyhedron_item();
  scene->replaceItem(i, poly_item);
  delete edit_item;
  return poly_item;
}

void Polyhedron_demo_edit_polyhedron_plugin_2::on_actionToggleEdit_triggered(bool edit) {
  this->edit_mode = edit;
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    Scene_polyhedron_item* poly_item = 
      qobject_cast<Scene_polyhedron_item*>(scene->item(i));
    if (poly_item) poly_item->update_halfedge_indices();
    Scene_edit_polyhedron_item_2* edit_item = 
      qobject_cast<Scene_edit_polyhedron_item_2*>(scene->item(i));
    if(edit && poly_item) {
      convert_to_edit_polyhedron(i, poly_item);
    } else if(!edit && edit_item) {
      convert_to_plain_polyhedron(i, edit_item);
    }
  }
}

void Polyhedron_demo_edit_polyhedron_plugin_2::item_destroyed() {
  QObject* obj = sender(); // the item that is destroyed
  //Deform_map::iterator it = deform_map.find(obj);
  //if(it != deform_map.end()) {
  //  delete it->second.deform_mesh;  // TODO: uncomment that!
  //  deform_map.erase(it);
  //}
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_edit_polyhedron_plugin_2, Polyhedron_demo_edit_polyhedron_plugin_2)

#include "Polyhedron_demo_edit_polyhedron_plugin_2.moc"
