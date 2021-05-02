#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QStringList>
#include <QInputDialog>
#include <QtPlugin>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_surface_mesh_item.h"

#include "Polyhedron_demo_detect_sharp_corners.h"

typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef CGAL::Kernel_traits<Scene_surface_mesh_item::Face_graph::Point>::Kernel Kernel;

typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
typedef boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;


using namespace CGAL::Three;

class Polyhedron_demo_detect_sharp_corners_plugin :
        public QObject,
        public Polyhedron_demo_plugin_interface {
Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(
          IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "detect_sharp_corners_plugin.json")

public:
    void init(QMainWindow *mainWindow, Scene_interface *scene_interface, Messages_interface *) {
        this->scene = scene_interface;
        this->mw = mainWindow;
        actionSharpCorners = new QAction("Detect Sharp Features", mw);
        actionSharpCorners->setObjectName("detectSharpFeaturesAction");
        if (actionSharpCorners) {
            connect(actionSharpCorners, SIGNAL(triggered()),
                    this, SLOT(detectSharpCornersWithInputDialog()));
        }
    }

    bool applicable(QAction *) const {
        Q_FOREACH(int index, scene->selectionIndices()) {
                Scene_facegraph_item *item =
                        qobject_cast<Scene_facegraph_item *>(scene->item(index));
                if (item) return true;
            }
        return false;
    }

    QList<QAction *> actions() const {
        return QList<QAction *>() << actionSharpCorners;
    }

public Q_SLOTS:

    void detectSharpCorners(bool input_dialog = false, double angle = 60);

    void detectSharpCornersWithInputDialog();

protected:
    Kernel::Vector_3 facet_normal(face_descriptor f);

    bool is_sharp(halfedge_descriptor he);

private:
    QAction *actionSharpCorners;
    CGAL::Three::Scene_interface *scene;
    QMainWindow *mw;
}; // end Polyhedron_demo_detect_sharp_corners_plugin

void Polyhedron_demo_detect_sharp_corners_plugin::detectSharpCornersWithInputDialog() {
    detectSharpCorners(true);
}

namespace PMP = CGAL::Polygon_mesh_processing;

void Polyhedron_demo_detect_sharp_corners_plugin::detectSharpCorners(bool input_dialog,
                                                                     double angle) {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    typedef std::pair<int, FaceGraph *> Poly_tuple;

    // Get selected items
    QList<Poly_tuple> polyhedrons;
    Q_FOREACH(int index, scene->selectionIndices()) {
            Scene_facegraph_item *item =
                    qobject_cast<Scene_facegraph_item *>(scene->item(index));
            if (!item)
                return;

            FaceGraph *pMesh = item->polyhedron();
            if (!pMesh)
                return;
            item->show_feature_edges(true);
            polyhedrons << std::make_pair(index, pMesh);
        }

    QApplication::restoreOverrideCursor();
    if (input_dialog) {
        bool ok = true;
        angle = QInputDialog::getDouble(NULL,
                                        tr("Sharp corners max angle"),
                                        tr("Angle in degrees between 0 and 180:"),
                                        angle, // value
                                        0.,          // min
                                        180., // max
                                        2,          // decimals
                                        &ok);
        if (!ok) return;
    }
    // Detect edges
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    std::size_t first_patch = 1;
    Q_FOREACH(Poly_tuple tuple, polyhedrons) {
            Scene_facegraph_item *item =
                    qobject_cast<Scene_facegraph_item *>(scene->item(tuple.first));
            FaceGraph *pMesh = tuple.second;
            if (!pMesh)
                continue;

            boost::property_map<FaceGraph, CGAL::vertex_is_feature_t>::type vif
                    = get(CGAL::vertex_is_feature, *pMesh);


            PMP::detect_sharp_corners(pMesh,angle,vif);
            std::size_t sharp_corners_counter;
            vertex_descriptor v;

                if (get(vif,v))
                    ++sharp_corners_counter;

            //update item
            item->setItemIsMulticolor(true);
            item->computeItemColorVectorAutomatically(true);
            item->invalidateOpenGLBuffers();

            // update scene
            scene->itemChanged(tuple.first);
        }

    // default cursor
    QApplication::restoreOverrideCursor();
}

#include "Detect_sharp_corners_plugin.moc"



