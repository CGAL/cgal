#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include "Scene_polyhedron_item.h"
#include "Scene_segmented_image_item.h"
#include "Image_type.h"
#include <QInputDialog>
#include <QFileDialog>
#include <fstream>

// declare the CGAL function
Scene_item* cgal_code_mesh_3(const Polyhedron*,
                             QString filename,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing);

Scene_item* cgal_code_mesh_3(const Image*,
                             QString filename,
                             const double angle,
                             const double sizing,
                             const double approx,
                             const double tets_sizing);

class Polyhedron_demo_mesh_3_plugin : 
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);
public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionMesh_3 = new QAction("Create a tetrahedral mesh", mw);
    if(actionMesh_3) {
      connect(actionMesh_3, SIGNAL(triggered()),
              this, SLOT(mesh_3()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMesh_3;
  }
public slots:
  void mesh_3();

private:
  QAction* actionMesh_3;
}; // end class Polyhedron_demo_mesh_3_plugin

void Polyhedron_demo_mesh_3_plugin::mesh_3()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_item* item = 0;
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Scene_segmented_image_item* image_item = 
    qobject_cast<Scene_segmented_image_item*>(scene->item(index));

  if(poly_item)
  {
    item = poly_item; 
    Polyhedron* pMesh = poly_item->polyhedron();
    if(!pMesh) return;
  }
  else if(image_item)
  {
    item = image_item;
    const Image* image = image_item->image();
    if(!image) return;
  }

  if(item) {
    // TODO: get parameters using ONE dialog box
    // sizing and approximation parameters should be expressed as ratio of 
    // scene bbox diagonal.

    double diag = scene->len_diagonal();

    bool ok;
    const double angle = 
      QInputDialog::getDouble(mw, tr("Min triangle angle"),
                              tr("Angle:"),
                              25., // default value
                              1., // min
                              30., // max
                              2, // decimals
                              &ok);
    if(!ok) return;

    const double sizing = 
      QInputDialog::getDouble(mw, "Facets sizing",
      "Facets size bound:",
      diag * 0.05, // default value
      diag * 10e-6, // min
      diag, // max
      4, // decimals
      &ok);
    if(!ok) return;

    const double tets_sizing = 
      QInputDialog::getDouble(mw, "Tetrahedra sizing",
      "Tetrahedra size bound:",
      diag * 0.05, // default value
      diag * 10e-6, // min
      diag, // max
      4, // decimals
      &ok);
    if(!ok) return;

    const double approx = 
      QInputDialog::getDouble(mw, "Approximation error",
      "Error:",
      diag * 0.005, // default value
      diag * 10e-7, // min
      diag, // max
      6, // decimals
      &ok);
    if(!ok) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    Scene_item* result_item = 0;
    if(poly_item) {
      Polyhedron* pMesh = poly_item->polyhedron();
      
      if(!pMesh) return;

      result_item = cgal_code_mesh_3(pMesh,
                                     item->name(),
                                     angle,
                                     sizing,
                                     approx,
                                     tets_sizing);
    } else if(image_item) {
      const Image* pImage = image_item->image();

      if(!pImage) return;

      result_item = cgal_code_mesh_3(pImage,
                                     item->name(),
                                     angle,
                                     sizing,
                                     approx,
                                     tets_sizing);
    }
    if(result_item) {
      result_item->setName(tr("%1 3d mesh (%2 %3 %4 %5)")
                           .arg(item->name())
                           .arg(angle)
                           .arg(sizing)
                           .arg(tets_sizing)
                           .arg(approx));
      result_item->setColor(Qt::magenta);
      result_item->setRenderingMode(item->renderingMode());
      item->setVisible(false);
      scene->itemChanged(index);
      scene->addItem(result_item);
    }
    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_mesh_3_plugin, Polyhedron_demo_mesh_3_plugin);

#include "Polyhedron_demo_mesh_3_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
