#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_interface.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Selection_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include <map>

class Polyhedron_demo_selection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
    Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable() const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
        || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex())); 
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionSelection; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m){
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionSelection = new QAction(tr("Selection"), mw);
    connect(actionSelection, SIGNAL(triggered()), this, SLOT(selection_action()));

    dock_widget = new QDockWidget("Selection", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Select_isolated_components_button,  SIGNAL(clicked()), this, SLOT(on_Select_isolated_components_button_clicked()));
    connect(ui_widget.Get_minimum_button,  SIGNAL(clicked()), this, SLOT(on_Get_minimum_button_clicked()));
    connect(ui_widget.Create_selection_item_button,  SIGNAL(clicked()), this, SLOT(on_Create_selection_item_button_clicked()));    
    connect(ui_widget.select_marked_edges_button,  SIGNAL(clicked()),
            this, SLOT(on_select_marked_edges_button_clicked()));

    QObject* scene = dynamic_cast<QObject*>(scene_interface);
    if(scene) { 
      connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)), this, SLOT(item_about_to_be_destroyed(Scene_item*)));
    } 
  }

  template<class SceneType>
  SceneType* get_selected_item() {
    int item_id = scene->mainSelectionIndex();
    SceneType* scene_item = qobject_cast<SceneType*>(scene->item(item_id));
    if(!scene_item) {
      // no selected SceneType, if there is only one in list use it, otherwise error
      int counter = 0;
      for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end && counter < 2; ++i) {
        if(SceneType* tmp = qobject_cast<SceneType*>(scene->item(i))) { 
          scene_item = tmp;
          counter++; 
        }
      }
      if(counter != 1) { return NULL; }
    }
    return scene_item;
  }

public slots:
  void selection_action() { 
    dock_widget->show();
  }
  // Select all
  void on_Select_all_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }
    selection_item->select_all();
  }
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }

    selection_item->clear();
  }
  // Isolated component related functions
  void on_Select_isolated_components_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }
    selection_item->select_isolated_components();
  }
  void on_Get_minimum_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }
    selection_item->get_minimum_isolated_component();
  }
  // Create selection item for selected polyhedron item
  void on_Create_selection_item_button_clicked() {
    Scene_polyhedron_item* poly_item = get_selected_item<Scene_polyhedron_item>();
    if(!poly_item) {
      print_message("Error: there is no selected polyhedron item!");
      return; 
    }

    QString poly_item_name = poly_item->name();
    Scene_polyhedron_selection_item* selection_poly =
      new Scene_polyhedron_selection_item(poly_item, &ui_widget);
    selection_item_map.insert(std::make_pair(poly_item, selection_poly));
    selection_poly->setName(QString("%1 (selection)").arg(poly_item->name()));
    selection_poly->setRenderingMode(Flat);

    mw->installEventFilter(selection_poly); // filter mainwindows events for key(pressed/released)
    scene->addItem(selection_poly);
  }
  void on_select_marked_edges_button_clicked()
  {
    for(Scene_interface::Item_id i = 0,
          end = scene->numberOfEntries();
        i < end; ++i)
    {
      Scene_polyhedron_selection_item* selection_item =
        qobject_cast<Scene_polyhedron_selection_item*>(scene->item(i));
      if(selection_item) {
        selection_item->select_marked_edges();
      }
    }
  }

  void item_about_to_be_destroyed(Scene_item* scene_item) {
    // if polyhedron item
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene_item);
    if(poly_item) {
      std::pair<Selection_item_map::iterator, Selection_item_map::iterator> res =
        selection_item_map.equal_range(poly_item);

      for(Selection_item_map::iterator begin = res.first; begin != res.second; ) {
        Scene_polyhedron_selection_item* selection_item = begin->second;
        selection_item_map.erase(begin++); // first erase from map, because scene->erase will cause a call to this function
        scene->erase( scene->item_id(selection_item) );
      }
    }
    // if polyhedron selection item
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene_item);
    if(selection_item) {
      Scene_polyhedron_item* poly_item = selection_item->polyhedron_item();
      std::pair<Selection_item_map::iterator, Selection_item_map::iterator> res =
        selection_item_map.equal_range(poly_item);
      for(Selection_item_map::iterator begin = res.first; begin != res.second; ++begin) {
        if(begin->second == selection_item) {
          selection_item_map.erase(begin); break;
        }
      }
    }
  }

private:
  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionSelection;

  QDockWidget* dock_widget;
  Ui::Selection ui_widget;
typedef std::multimap<Scene_polyhedron_item*, Scene_polyhedron_selection_item*> Selection_item_map;
  Selection_item_map selection_item_map;
}; // end Polyhedron_demo_selection_plugin

Q_EXPORT_PLUGIN2(Polyhedron_demo_selection_plugin, Polyhedron_demo_selection_plugin)

#include "Polyhedron_demo_selection_plugin.moc"
