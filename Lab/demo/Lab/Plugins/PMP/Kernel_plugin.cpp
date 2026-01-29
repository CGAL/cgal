#include <fstream>
#include <QApplication>
#include <QMainWindow>
#include <QMessageBox>
#include <QAction>
#include <QStringList>

#include "Scene_surface_mesh_item.h"
#include "SMesh_type.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/kernel.h>

#include "Kernel_type.h"
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;
typedef Kernel::FT FT;
using namespace CGAL::Three;
class CGAL_Lab_mesh_kernel_plugin :
  public QObject,
  public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

public:
  QList<QAction*> actions() const {
    return _actions;
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

   void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
   {
     scene = scene_interface;
     mw = mainWindow;
     QAction* actionKernel = new QAction(tr("Compute Mesh Kernel"), mainWindow);
     actionKernel->setProperty("subMenuName","Polygon Mesh Processing");
     connect(actionKernel, SIGNAL(triggered()),
             this, SLOT(on_actionKernel_triggered()));
     _actions << actionKernel;
   }
public Q_SLOTS:
  void on_actionKernel_triggered();
private:
  QList<QAction*> _actions;
  CGAL::Three::Scene_interface* scene;
  QMainWindow* mw;

}; // end CGAL_Lab_mesh_kernel_plugin


void CGAL_Lab_mesh_kernel_plugin::on_actionKernel_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* sm_item =
    qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(sm_item)
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    SMesh* sMesh = sm_item->polyhedron();
    SMesh* kernel = new SMesh;

    CGAL::Polygon_mesh_processing::kernel(*sMesh, *kernel);

    if(is_empty(*kernel)){
      QApplication::restoreOverrideCursor();
      std::cout << "done (empty kernel)" << std::endl;
      QMessageBox::information(mw, tr("Empty Kernel"),
                               tr("The kernel of the polyhedron \"%1\" is empty.").
                               arg(sm_item->name()));
      return;
    }

    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(kernel);
    new_item->setName(tr("%1 (kernel)").arg(sm_item->name()));
    new_item->setColor(Qt::magenta);
    new_item->setRenderingMode(sm_item->renderingMode());
    sm_item->setRenderingMode(Wireframe);

    scene->addItem(new_item);
    scene->itemChanged(sm_item);

    QApplication::restoreOverrideCursor();
  }
}

#include "Kernel_plugin.moc"
