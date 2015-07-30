#undef NDEBUG
#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "ui_Fairing_widget.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/refine.h>

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


class Polyhedron_demo_fairing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));  
  }
  void print_message(QString message) { messages->information(message);}
  QList<QAction*> actions() const { return QList<QAction*>() << actionFairing; }

  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionFairing = new QAction(tr("Fairing"), mw);
    connect(actionFairing, SIGNAL(triggered()), this, SLOT(fairing_action()));

    dock_widget = new QDockWidget("Fairing", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    add_dock_widget(dock_widget);

    connect(ui_widget.Fair_button,  SIGNAL(clicked()), this, SLOT(on_Fair_button_clicked()));  
    connect(ui_widget.Refine_button,  SIGNAL(clicked()), this, SLOT(on_Refine_button_clicked()));
  }

public Q_SLOTS:
  void fairing_action() {
    dock_widget->show();
    dock_widget->raise();
  }

  void on_Fair_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) { return; }

    if(selection_item->selected_vertices.empty()) {
      print_message("Error: please select a region of vertices!");
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int weight_index = ui_widget.Weight_combo_box->currentIndex();
    unsigned int continuity = ui_widget.Continuity_spin_box->value();

    if(weight_index == 1)
      CGAL::Polygon_mesh_processing::fair(*selection_item->polyhedron(),
        selection_item->selected_vertices,
        CGAL::Polygon_mesh_processing::parameters::weight_calculator(CGAL::internal::Uniform_weight_fairing<Polyhedron>(*selection_item->polyhedron())).
        fairing_continuity(continuity));
    if(weight_index == 0)
      CGAL::Polygon_mesh_processing::fair(*selection_item->polyhedron(),
        selection_item->selected_vertices,
        CGAL::Polygon_mesh_processing::parameters::fairing_continuity(continuity));
    selection_item->changed_with_poly_item();
    QApplication::restoreOverrideCursor();
  }

  void on_Refine_button_clicked() {
    Scene_polyhedron_selection_item* selection_item = get_selected_item<Scene_polyhedron_selection_item>();
    if(!selection_item) { return; }

    if(selection_item->selected_facets.empty()) {
      print_message("Error: please select a region of facets!");
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    double alpha = ui_widget.Density_control_factor_spin_box->value();
    std::vector<Polyhedron::Facet_handle> new_facets;

    CGAL::Polygon_mesh_processing::refine(*selection_item->polyhedron(),
      selection_item->selected_facets,
      std::back_inserter(new_facets),
      CGAL::Emptyset_iterator(),
      CGAL::Polygon_mesh_processing::parameters::density_control_factor(alpha));
    // add new facets to selection
    for(std::vector<Polyhedron::Facet_handle>::iterator it = new_facets.begin(); it != new_facets.end(); ++it) {
      selection_item->selected_facets.insert(*it);
    }
    selection_item->changed_with_poly_item();
    QApplication::restoreOverrideCursor();
  }

private:
  Messages_interface* messages;
  QAction* actionFairing;

  QDockWidget* dock_widget;
  Ui::Fairing ui_widget;

}; // end Polyhedron_demo_fairing_plugin

// Q_EXPORT_PLUGIN2(Polyhedron_demo_fairing_plugin, Polyhedron_demo_fairing_plugin)

#include "Polyhedron_demo_fairing_plugin.moc"
