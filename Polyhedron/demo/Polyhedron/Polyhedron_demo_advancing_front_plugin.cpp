#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include <Scene_polyhedron_item.h>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include "ui_Polyhedron_demo_advancing_front_plugin.h"

struct Perimeter {

  double bound;

  Perimeter(double bound)
    : bound(bound)
  {}

  bool operator()(const Kernel::Point_3& p, const Kernel::Point_3& q, const Kernel::Point_3& r) const
  {
    if(bound == 0){
      return false;
    }
    double d  = sqrt(squared_distance(p,q));
    if(d>bound) return true;
    d += sqrt(squared_distance(p,r)) ;
    if(d>bound) return true;
    d+= sqrt(squared_distance(q,r));
    return d>bound;
  }
};


class Polyhedron_demo_advancing_front_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionAdvancingFrontReconstruction;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionAdvancingFrontReconstruction = new QAction(tr("Advancing Front reconstruction"), mainWindow);
    actionAdvancingFrontReconstruction->setObjectName("actionAdvancingFrontReconstruction");
    
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionAdvancingFrontReconstruction;
  }

public Q_SLOTS:
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


    QApplication::setOverrideCursor(Qt::WaitCursor);

    // Add polyhedron to scene
    
    // Reconstruct point set as a polyhedron
    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(Polyhedron());
    Polyhedron& P = * const_cast<Polyhedron*>(new_item->polyhedron());
    Perimeter filter(sm_perimeter);
      CGAL::advancing_front_surface_reconstruction((points)->begin(), points->end(), P, filter);


    new_item->setName(tr("%1 Advancing Front (%2)")
                      .arg(point_set_item->name())
                      .arg(sm_perimeter));
    new_item->setColor(Qt::lightGray);
    scene->addItem(new_item);
    
    // Hide point set
    point_set_item->setVisible(false);
    scene->itemChanged(index);
    

    QApplication::restoreOverrideCursor();
  
  }
}

#include "Polyhedron_demo_advancing_front_plugin.moc"
