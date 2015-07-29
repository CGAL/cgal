#include <QtCore/qglobal.h>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include "Scene_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "ui_Point_set_selection_widget.h"
#include "Point_set_3.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include <map>


class Polyhedron_demo_point_set_selection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction*) const { 
      return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointSetSelection; }
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionPointSetSelection = new QAction(tr("Selection"), mw);
    connect(actionPointSetSelection, SIGNAL(triggered()), this, SLOT(selection_action()));

    dock_widget = new QDockWidget("Point Set Selection", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    add_dock_widget(dock_widget);

    connect(ui_widget.Selection_tool_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_tool_combo_box_changed(int)));
    connect(ui_widget.Selection_mode_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_mode_combo_box_changed(int)));
    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Invert_selection_button,  SIGNAL(clicked()), this, SLOT(on_Invert_selection_button_clicked()));
    connect(ui_widget.Erase_selected_points_button,  SIGNAL(clicked()), this, SLOT(on_Erase_selected_points_button_clicked()));
    connect(ui_widget.Create_point_set_item_button, SIGNAL(clicked()), this, SLOT(on_Create_point_set_item_button_clicked()));
  }

public Q_SLOTS:
  void selection_action() { 
    dock_widget->show();
    dock_widget->raise();
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item)
      return;
  }
  
  // Select all
  void on_Select_all_button_clicked() {
    Scene_points_with_normal_item* point_set_item = get_selected_item<Scene_points_with_normal_item>();
    if(!point_set_item)
      {
	print_message("Error: no point set selected!");
	return; 
      }

    point_set_item->selectAll();
  }
  
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }

    point_set_item->resetSelection();
  }

  void on_Erase_selected_points_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }

    point_set_item->deleteSelection();
  }

  void on_Invert_selection_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }

    point_set_item->invertSelection();
  }

  void on_Create_point_set_item_button_clicked() {
    Scene_points_with_normal_item* point_set_item
      = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!point_set_item) {
      print_message("Error: no point set selected!");
      return; 
    }
    if(point_set_item->isSelectionEmpty ()) {
      print_message("Error: there is no selected point in point set item!");
      return;
    }
    
    Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item();
    new_item->setName(QString("%1 (selected points)").arg(point_set_item->name()));
    new_item->set_has_normals (point_set_item->has_normals());
    new_item->setColor(point_set_item->color());
    new_item->setRenderingMode(point_set_item->renderingMode());
    new_item->setVisible(point_set_item->visible());

    typedef Point_set_3<Kernel> Point_set;
    for(typename Point_set::iterator it = point_set_item->point_set()->begin ();
	it != point_set_item->point_set()->end(); ++ it) {
      if (it->is_selected ())
	new_item->point_set()->push_back(*it);
    }
    new_item->resetSelection();
    new_item->changed();

    scene->addItem(new_item);
 }



private:
  Messages_interface* messages;
  QAction* actionPointSetSelection;

  QDockWidget* dock_widget;
  Ui::PointSetSelection ui_widget;
  Scene_points_with_normal_item* point_set_item;

}; // end Polyhedron_demo_point_set_selection_plugin

//Q_EXPORT_PLUGIN2(Polyhedron_demo_point_set_selection_plugin, Polyhedron_demo_point_set_selection_plugin)

#include "Polyhedron_demo_point_set_selection_plugin.moc"
