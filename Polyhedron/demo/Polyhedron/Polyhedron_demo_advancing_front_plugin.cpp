#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include <Scene_polygon_soup_item.h>


#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include "ui_Polyhedron_demo_advancing_front_plugin.h"

// Poisson reconstruction method:
// Reconstructs a surface mesh from a point set and writes facet indices into polygon soup.
void advancing_front_reconstruct(const Point_set& points,
                                 double sm_perimeter,
                                 double sm_area,
                                 Scene_polygon_soup_item*);

class Polyhedron_demo_advancing_front_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionAdvancingFrontReconstruction;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionAdvancingFrontReconstruction = new QAction(tr("Advancing Front reconstruction"), mainWindow);
    actionAdvancingFrontReconstruction->setObjectName("actionAdvancingFrontReconstruction");
    
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  //! Applicate for Point_sets with normals.
  bool applicable() const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionAdvancingFrontReconstruction;
  }

public slots:
  void on_actionAdvancingFrontReconstruction_triggered();
}; // end class Polyhedron_demo_advancing_front_plugin


class Polyhedron_demo_advancing_front_plugin_dialog : public QDialog, private Ui::AdvancingFrontDialog
{
  Q_OBJECT
  public:
    Polyhedron_demo_advancing_front_plugin_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
      
    }

    double trianglePerimeter() const { return m_inputPerimeter->value(); }
    double triangleArea() const { return m_inputArea->value(); }
};

void Polyhedron_demo_advancing_front_plugin::on_actionAdvancingFrontReconstruction_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(point_set_item)
  {
    // Gets point set
    Point_set* points = point_set_item->point_set();
    if(!points) return;

    // Gets options
    Polyhedron_demo_advancing_front_plugin_dialog dialog;
    if(!dialog.exec())
      return;
    const double sm_perimeter     = dialog.trianglePerimeter();
    const double sm_area    = dialog.triangleArea();


    QApplication::setOverrideCursor(Qt::WaitCursor);

    // Add polyhedron to scene
    Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item();
    
    for(Point_set::iterator it = points->begin(); it!= points->end(); ++it){
      new_item->new_vertex(it->x(), it->y(), it->z());
    }

    // Reconstruct point set as a polyhedron
    advancing_front_reconstruct(*points, sm_perimeter, sm_area, new_item);


    new_item->setName(tr("%1 Advancing Front (%2 %3)")
                      .arg(point_set_item->name())
                      .arg(sm_perimeter)
                      .arg(sm_area));
    new_item->setColor(Qt::lightGray);
    scene->addItem(new_item);
    
    // Hide point set
    point_set_item->setVisible(false);
    scene->itemChanged(index);
    

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_advancing_front_plugin, Polyhedron_demo_advancing_front_plugin)

#include "Polyhedron_demo_advancing_front_plugin.moc"
