
#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>

#include <CGAL/Polygon_mesh_processing/smoothing.h>
#include <CGAL/Polygon_mesh_processing/random_pertubation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Smoothing_plugin.h"


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

      actionSmoothing_ = new QAction("Smoothing", mw);
      actionSmoothing_->setProperty("subMenuName", "Polygon Mesh Processing");
      if (actionSmoothing_) {
        connect(actionSmoothing_, SIGNAL(triggered()),
          this, SLOT(smooth()));
      }
    }

    QList<QAction*> actions() const {
        return QList<QAction*>() << actionSmoothing_;
      }

    bool applicable(QAction*) const
    {
      const Scene_interface::Item_id index = scene->mainSelectionIndex();
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
      else if (qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index)))
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
        Scene_polyhedron_selection_item* selection_item =
          qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

        //Create dialog box
        QDialog dialog(mw);
        Ui::Smoothing_dialog ui
          = smooth_dialog(&dialog, poly_item, selection_item);

        //Get values
        int i = dialog.exec();
        if (i == QDialog::Rejected)
        {
          std::cout << "Smoothing aborted" << std::endl;
          return;
        }

        // wait cursor
        QApplication::setOverrideCursor(Qt::WaitCursor);
        QTime time;
        time.start();

        std::cout << "Smoothing..." << std::endl;

        if(poly_item || selection_item)
        {
            Polyhedron& pmesh = *poly_item->polyhedron();
            CGAL::Polygon_mesh_processing::compatible_remeshing(pmesh, faces(pmesh), edges(pmesh),
                                                                CGAL::Polygon_mesh_processing::parameters::all_default());

        }
        else
        {
            std::cerr << "No smooth!" << std::endl;
        }

        std::cout << " ok (" << time.elapsed() << " ms)" << std::endl;

        // default cursor
         QApplication::restoreOverrideCursor();
    }




    Ui::Smoothing_dialog smooth_dialog(QDialog* dialog, Scene_polyhedron_item* poly_item, Scene_polyhedron_selection_item* selection_item)
    {

        Ui::Smoothing_dialog ui;
        ui.setupUi(dialog);
        connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
        connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

        //Set default parameters
        Scene_interface::Bbox bbox = poly_item != NULL ? poly_item->bbox()
          : (selection_item != NULL ? selection_item->bbox()
          : scene->bbox());
        ui.objectName->setText(poly_item != NULL ? poly_item->name()
          : (selection_item != NULL ? selection_item->name()
          : QString("Remeshing parameters")));

        ui.objectNameSize->setText(
          tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
          .arg(bbox.xmax() - bbox.xmin(), 0, 'g', 3)
          .arg(bbox.ymax() - bbox.ymin(), 0, 'g', 3)
          .arg(bbox.zmax() - bbox.zmin(), 0, 'g', 3));

        double diago_length = CGAL::sqrt((bbox.xmax() - bbox.xmin())*(bbox.xmax() - bbox.xmin())
          + (bbox.ymax() - bbox.ymin())*(bbox.ymax() - bbox.ymin())
          + (bbox.zmax() - bbox.zmin())*(bbox.zmax() - bbox.zmin()));
        double log = std::log10(diago_length);
        unsigned int nb_decimals = (log > 0) ? 5 : (std::ceil(-log) + 3);


        return ui;
    }


private:
    Scene_interface *scene;
    QMainWindow* mw;
    QAction* actionSmoothing_;



};



#include "Smoothing_plugin.moc"


















