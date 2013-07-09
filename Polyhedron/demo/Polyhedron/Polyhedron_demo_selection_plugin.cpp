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

    connect(dock_widget, SIGNAL(visibilityChanged(bool)), this, SLOT(dock_widget_visibility_changed(bool)) );
    connect(ui_widget.Select_all_button,  SIGNAL(clicked()), this, SLOT(on_Select_all_button_clicked()));
    connect(ui_widget.Clear_button,  SIGNAL(clicked()), this, SLOT(on_Clear_button_clicked()));
    connect(ui_widget.Show_selection_check_box, SIGNAL(stateChanged(int)), this, SLOT(on_Show_selection_check_box_stateChanged(int)));
    connect(ui_widget.Select_isolated_components_button,  SIGNAL(clicked()), this, SLOT(on_Select_isolated_components_button_clicked()));
    connect(ui_widget.Get_minimum_button,  SIGNAL(clicked()), this, SLOT(on_Get_minimum_button_clicked()));

    QObject* scene = dynamic_cast<QObject*>(scene_interface);
    if(scene) { connect(scene, SIGNAL(newItem(int)), this, SLOT(new_item_created(int))); } 
  }

  Scene_polyhedron_selection_item* get_selected_item() {
    int item_id = scene->mainSelectionIndex();
    Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(item_id));
    if(!selection_item) {
      print_message("Error: there is no selected polyhedron selection item!");
      return NULL;
    } 
    return selection_item;
  }

public slots:
  void selection_action() { 
    dock_widget->show();
  }
  // when dock is visible convert plain items to selection items
  void dock_widget_visibility_changed(bool visible)
  {
    for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
    {
      if(visible) { convert_to_selection_polyhedron(i); }
      else        { convert_to_plain_polyhedron(i);     }
    }
  }
  // converters between selection and plain polyhedron item
  Scene_polyhedron_selection_item* convert_to_selection_polyhedron(Scene_interface::Item_id i) {
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(i));
    if(!poly_item) { return NULL; }

    QString poly_item_name = poly_item->name();
    Scene_polyhedron_selection_item* selection_poly =
      new Scene_polyhedron_selection_item(poly_item, &ui_widget);
    selection_poly->setColor(poly_item->color());
    selection_poly->setName(QString("%1 (selection)").arg(poly_item->name()));
    selection_poly->setRenderingMode(poly_item->renderingMode());
    selection_poly->setVisible(poly_item->visible());

    mw->installEventFilter(selection_poly); // filter mainwindows events for key(pressed/released)
    scene->replaceItem(i, selection_poly);
    return selection_poly;
  }
  Scene_polyhedron_item* convert_to_plain_polyhedron(Scene_interface::Item_id i) {
    Scene_polyhedron_selection_item* selection_poly = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(i));
    if(!selection_poly) { return NULL; }

    Scene_polyhedron_item* poly_item = selection_poly->polyhedron_item();
    scene->replaceItem(i, poly_item);
    selection_poly->set_polyhedron_item(NULL); // do not let it delete item
    delete selection_poly;
    return poly_item;
  }
  // Select all
  void on_Select_all_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item();
    if(!selection_item) { return; }
    selection_item->select_all();
  }
  // Clear selection
  void on_Clear_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item();
    if(!selection_item) { return; }

    selection_item->clear();
  }
  // Isolated component related functions
  void on_Select_isolated_components_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item();
    if(!selection_item) { return; }
    selection_item->select_isolated_components();
  }
  void on_Get_minimum_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item();
    if(!selection_item) { return; }
    selection_item->get_minimum_isolated_component();
  }
  // redraw when show checkbox is changed
  void on_Show_selection_check_box_stateChanged(int /*state*/)
  {
    for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i) {
      Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(i));
      if(selection_item) 
      { scene->itemChanged(selection_item); } // just for redraw
    }  
  }
  // Convert new coming item to selection item
  void new_item_created(int item_id)
  {
    if(dock_widget->isVisible()) 
    { convert_to_selection_polyhedron(item_id); }
  }

private:
  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionSelection;

  QDockWidget* dock_widget;
  Ui::Selection ui_widget;

}; // end Polyhedron_demo_selection_plugin

Q_EXPORT_PLUGIN2(Polyhedron_demo_selection_plugin, Polyhedron_demo_selection_plugin)

#include "Polyhedron_demo_selection_plugin.moc"
