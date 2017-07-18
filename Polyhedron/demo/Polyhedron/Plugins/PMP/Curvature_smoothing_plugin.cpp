
#include <QtCore/qglobal.h>
#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>



using namespace CGAL::Three;
class Polyhedron_demo_smothing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")


public:
    void init(QMainWindow* mainWindow,
              Scene_interface* scene_interface,
              Messages_interface*)
    {
      this->scene = scene_interface;
      this->mw = mainWindow;

      actionCurvatureSmoothing_ = new QAction("Mean curvature flow smoothing", mw);
      actionCurvatureSmoothing_->setProperty("subMenuName", "Polygon Mesh Processing");
      if (actionCurvatureSmoothing_) {
        connect(actionCurvatureSmoothing_, SIGNAL(triggered()),
          this, SLOT(smooth()));
      }
    }

    QList<QAction*> actions() const {
        return QList<QAction*>() << actionCurvatureSmoothing_;
      }

    bool applicable(QAction*) const
    {
      const Scene_interface::Item_id index = scene->mainSelectionIndex();
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
      else
        return false;
    }



public Q_SLOTS:
    void smooth()
    {
        const Scene_interface::Item_id index = scene->mainSelectionIndex();
        Scene_polyhedron_item* poly_item =
          qobject_cast<Scene_polyhedron_item*>(scene->item(index));

        // wait cursor
        QApplication::setOverrideCursor(Qt::WaitCursor);
        QTime time;
        time.start();

        std::cout << "Smoothing..." << std::endl;

        if(poly_item)
        {
            Polyhedron& pmesh = *poly_item->polyhedron();
            CGAL::Polygon_mesh_processing::curvature_flow(pmesh);

            // to fix
            if(poly_item)
            {
                poly_item->invalidateOpenGLBuffers();
                Q_EMIT poly_item->itemChanged();
            }
            else
            {
                std::cerr<<"selection_item!"<<std::endl;
            }
            //

        }
        else
        {
            std::cerr << "No smoothing!" << std::endl;
        }

        std::cout << " ok (" << time.elapsed() << " ms)" << std::endl;

        // default cursor
         QApplication::restoreOverrideCursor();
    }



private:
    Scene_interface *scene;
    QMainWindow* mw;
    QAction* actionCurvatureSmoothing_;



};



#include "Curvature_smoothing_plugin.moc"





