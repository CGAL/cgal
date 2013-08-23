#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Scene_points_with_normal_item.h"

#include "Scene_interface.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "ui_Selection_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>

#include <map>

class Polyhedron_demo_selection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
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
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;

    actionSelection = new QAction(tr("Selection"), mw);
    connect(actionSelection, SIGNAL(triggered()), this, SLOT(selection_action()));

    dock_widget = new QDockWidget("Selection", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    add_dock_widget(dock_widget);

    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Select_isolated_components_button,  SIGNAL(clicked()), this, SLOT(on_Select_isolated_components_button_clicked()));
    connect(ui_widget.Get_minimum_button,  SIGNAL(clicked()), this, SLOT(on_Get_minimum_button_clicked()));
    connect(ui_widget.Create_selection_item_button,  SIGNAL(clicked()), this, SLOT(on_Create_selection_item_button_clicked()));    
    connect(ui_widget.Selection_type_combo_box, SIGNAL(currentIndexChanged(int)), 
            this, SLOT(on_Selection_type_combo_box_changed(int)));
    connect(ui_widget.Insertion_radio_button, SIGNAL(toggled(bool)), this, SLOT(on_Insertion_radio_button_toggled(bool)));
    connect(ui_widget.Brush_size_spin_box, SIGNAL(valueChanged(int)), this, SLOT(on_Brush_size_spin_box_changed(int)));
    connect(ui_widget.Create_point_set_item_button, SIGNAL(clicked()), this, SLOT(on_Create_point_set_item_button_clicked()));
    connect(ui_widget.Erase_selected_facets_button, SIGNAL(clicked()), this, SLOT(on_Erase_selected_facets_button_clicked()));
    connect(ui_widget.Expand_shrink_button, SIGNAL(clicked()), this, SLOT(on_Expand_shrink_button_clicked()));
    connect(ui_widget.Create_polyhedron_item_button, SIGNAL(clicked()), this, SLOT(on_Create_polyhedron_item_button_clicked()));
    QObject* scene = dynamic_cast<QObject*>(scene_interface);
    if(scene) { 
      connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)), this, SLOT(item_about_to_be_destroyed(Scene_item*)));
      connect(scene, SIGNAL(newItem(int)), this, SLOT(new_item_created(int)));
    } 
  }

public slots:
  void selection_action() { 
    dock_widget->show();
    dock_widget->raise();
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

    boost::optional<std::size_t> minimum = 
      selection_item->select_isolated_components(ui_widget.Threshold_size_spin_box->value());
    if(minimum) {
      ui_widget.Threshold_size_spin_box->setValue(*minimum);
    }
  }
  void on_Get_minimum_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }
    boost::optional<std::size_t> minimum = selection_item->get_minimum_isolated_component();
    if(minimum) {
      ui_widget.Threshold_size_spin_box->setValue(*minimum);
    }
  }
  // Create selection item for selected polyhedron item
  void on_Create_selection_item_button_clicked() {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    Scene_polyhedron_item* poly_item = get_selected_item<Scene_polyhedron_item>();
    if(!poly_item) {
      print_message("Error: there is no selected polyhedron item!");
      return; 
    }
    // all other arrangements (putting inside selection_item_map), setting names etc,
    // other params (e.g. k_ring) will be set inside new_item_created
    scene->addItem(new Scene_polyhedron_selection_item(poly_item, mw));
  }
  void on_Selection_type_combo_box_changed(int index) {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_active_handle_type(static_cast<Active_handle::Type>(index));
    }
  }
  void on_Insertion_radio_button_toggled(bool toggle){
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_is_insert(toggle);
    }
  }
  void on_Brush_size_spin_box_changed(int value) {
    for(Selection_item_map::iterator it = selection_item_map.begin(); it != selection_item_map.end(); ++it) {
      it->second->set_k_ring(value);
    }
  }

  void on_Create_point_set_item_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }
    if(selection_item->selected_vertices.empty()) {
      print_message("Error: there is no selected vertex in polyhedron selection item!");
      return;
    }
    Scene_points_with_normal_item* point_item = new Scene_points_with_normal_item();
    point_item->setName(QString("%1-points").arg(selection_item->name()));
    for(Scene_polyhedron_selection_item::Selection_set_vertex::iterator begin = selection_item->selected_vertices.begin(); 
       begin != selection_item->selected_vertices.end(); ++begin) {
       point_item->point_set()->push_back((*begin)->point());
    }
    
    scene->addItem(point_item);
  }
  void on_Erase_selected_facets_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }

    selection_item->erase_selected_facets();
  }
  void on_Create_polyhedron_item_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }

    Scene_polyhedron_item* poly_item = new Scene_polyhedron_item();
    if(selection_item->export_selected_facets_as_polyhedron(poly_item->polyhedron())) {
      poly_item->setName(QString("%1-facets").arg(selection_item->name()));
      poly_item->changed(); // for init()
      scene->addItem(poly_item);
    }
    else {
      delete poly_item;
      print_message("Error: polyhedron item is not created!");
    }
  }
  void on_Expand_shrink_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return; 
    }

    int steps = ui_widget.Expand_shrink_spin_box->value();
    selection_item->expand_or_shrink(steps);
  }
  // To handle empty selection items coming from loader
  void new_item_created(int item_id) {
    typedef Scene_polyhedron_selection_item::Active_handle Active_handle;
    Scene_polyhedron_selection_item* selection_item = 
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(item_id));
    if(!selection_item) { return; }

    Scene_polyhedron_item* poly_item = get_selected_item<Scene_polyhedron_item>();
    if(!poly_item) {
      CGAL_assertion(selection_item->polyhedron_item() == NULL); // which means it is coming from selection_io loader
      print_message("Error: please select corresponding polyhedron item from Geometric Objects list.");
      scene->erase(item_id);
      return;
    }

    if(selection_item->polyhedron_item() == NULL) { //coming from selection_io loader
      if(!selection_item->actual_load(poly_item, mw)) {
        print_message("Error: loading selection item is not successful!");
        scene->erase(item_id);
        return;
      }
    }
    // now set default params both for selection items coming from selection_io, or on_Create_selection_item_button_clicked
    Active_handle::Type aht = static_cast<Active_handle::Type>(ui_widget.Selection_type_combo_box->currentIndex());
    bool is_insert = ui_widget.Insertion_radio_button->isChecked();
    int k_ring = ui_widget.Brush_size_spin_box->value();

    selection_item->set_active_handle_type(aht);
    selection_item->set_is_insert(is_insert);
    selection_item->set_k_ring(k_ring);
    selection_item->setRenderingMode(Flat);
    if(selection_item->name() == "unamed") {
      selection_item->setName(tr("%1 (selection)").arg(poly_item->name()));
    }

    selection_item_map.insert(std::make_pair(poly_item, selection_item));
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
  Messages_interface* messages;
  QAction* actionSelection;

  QDockWidget* dock_widget;
  Ui::Selection ui_widget;
typedef std::multimap<Scene_polyhedron_item*, Scene_polyhedron_selection_item*> Selection_item_map;
  Selection_item_map selection_item_map;
}; // end Polyhedron_demo_selection_plugin

Q_EXPORT_PLUGIN2(Polyhedron_demo_selection_plugin, Polyhedron_demo_selection_plugin)

#include "Polyhedron_demo_selection_plugin.moc"
