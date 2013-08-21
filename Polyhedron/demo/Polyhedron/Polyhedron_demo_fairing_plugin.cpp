#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Fairing_widget.h"
#include "Polyhedron_type.h"

#include <CGAL/Hole_filling.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QFileDialog>

#include <vector>
#include <algorithm>
#include <queue>

#include <boost/function_output_iterator.hpp>

class Polyhedron_demo_fairing_plugin :
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
  void print_message(QString message) { messages->information(message);}
  QList<QAction*> actions() const { return QList<QAction*>() << actionFairing; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionFairing = new QAction(tr("Fairing"), mw);
    connect(actionFairing, SIGNAL(triggered()), this, SLOT(fairing_action()));

    dock_widget = new QDockWidget("Fairing", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    connect(ui_widget.Fair_button,  SIGNAL(clicked()), this, SLOT(on_Fair_button_clicked()));  
    connect(ui_widget.Refine_button,  SIGNAL(clicked()), this, SLOT(on_Refine_button_clicked()));
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
  void fairing_action() {
    dock_widget->show();
  }

  void on_Fair_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item();
    if(!selection_item) { return; }

    if(selection_item->selected_vertices.empty()) {
      print_message("Error: please select a region of vertices!");
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int weight_index = ui_widget.weight_combo_box->currentIndex();

    if(weight_index == 1)
      CGAL::fair(*selection_item->polyhedron(), selection_item->selected_vertices.begin(),
        selection_item->selected_vertices.end(),
        CGAL::internal::Uniform_weight_fairing<Polyhedron>());
    if(weight_index == 0)
      CGAL::fair(*selection_item->polyhedron(), selection_item->selected_vertices.begin(),
        selection_item->selected_vertices.end(),
        CGAL::internal::Cotangent_weight_with_voronoi_area_fairing<Polyhedron>());
    selection_item->changed_with_poly_item();
    QApplication::restoreOverrideCursor();
  }

  void on_Refine_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item();
    if(!selection_item) { return; }

    if(selection_item->selected_facets.empty()) {
      print_message("Error: please select a region of facets!");
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    double alpha = ui_widget.Density_control_factor_spin_box->value();
    std::vector<Polyhedron::Facet_handle> new_facets;

    CGAL::refine(*selection_item->polyhedron(), selection_item->selected_facets.begin(),
      selection_item->selected_facets.end(), std::back_inserter(new_facets), Nop_out(), alpha);
    // add new facets to selection
    for(std::vector<Polyhedron::Facet_handle>::iterator it = new_facets.begin(); it != new_facets.end(); ++it) {
      selection_item->selected_facets.insert(*it);
    }
    selection_item->changed_with_poly_item();
    QApplication::restoreOverrideCursor();
  }

private:
  struct Nop_functor {
    template<class T>
    void operator()(const T & /*t*/) const {}
  };
  typedef boost::function_output_iterator<Nop_functor> Nop_out;

  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionFairing;

  QDockWidget* dock_widget;
  Ui::Fairing ui_widget;

}; // end Polyhedron_demo_fairing_plugin

Q_EXPORT_PLUGIN2(Polyhedron_demo_fairing_plugin, Polyhedron_demo_fairing_plugin)

#include "Polyhedron_demo_fairing_plugin.moc"
