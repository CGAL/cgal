#include "ui_Approx_decomposition_widget.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include "Scene.h"
#include "Color_map.h"

//#define CGAL_APPROX_DECOMPOSITION_VERBOSE
#include <CGAL/approx_decomposition.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/property_map.h>
#include <CGAL/Random.h>

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>
#include <QDebug>
#include <QObject>
#include <QDockWidget>

using namespace CGAL::Three;

typedef CGAL::Real_timer Timer;

class Polyhedron_demo_approx_decomposition_plugin : public QObject, public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
    QList<QAction*> actions() const
    {
        return QList<QAction*>() << m_decomposition_action;
    }

    bool applicable(QAction*) const
    {
        return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
               qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    }

    void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
    {
        this->scene = scene_interface;
        this->mw = mainWindow;

        // actions
        m_decomposition_action = new QAction("Approximate Decomposition", mw);
        m_decomposition_action->setProperty("subMenuName", "Triangulated Surface Mesh Segmentation");
        m_decomposition_action->setObjectName("decomposition_action");

        autoConnectActions();

        // ui
        m_dock_widget = new QDockWidget("Approximate decomposition parameters", mw);
        m_dock_widget->setVisible(false); // do not show at the beginning
        m_ui_widget.setupUi(m_dock_widget);
        this->mw->addDockWidget(Qt::LeftDockWidgetArea, m_dock_widget);

        // signal-slot bindings
        connect(m_ui_widget.decompose_button, SIGNAL(clicked()), this, SLOT(on_decompose_button_clicked()));
    }

    virtual void closure()
    {
        m_dock_widget->hide();
    }

//    template<class FacegraphItem>
//    void apply_Decomposition_button_clicked(FacegraphItem* item);

public Q_SLOTS:
    void on_decomposition_action_triggered()
    {
        m_dock_widget->show();
    }

    void on_decompose_button_clicked()
    {
        std::cout << "Decomposing..." << std::endl;
        
        CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
        
        Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
        if (item)
        {
            decompose(item);
            return;
        }
        
        Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
        if (sm_item)
        {
            decompose(sm_item);
            return;
        }
    }

private:
    QAction* m_decomposition_action;
    QDockWidget* m_dock_widget;
    Ui::decomposition_widget m_ui_widget;

    template <class FacegraphItem>
    void decompose(FacegraphItem* item)
    {
        typedef typename FacegraphItem::Face_graph Facegraph;
        
        QApplication::setOverrideCursor(Qt::WaitCursor);

        // parameters
        std::size_t min_number_of_clusters = m_ui_widget.min_clusters_spin_box->value();
        double concavity_threshold = m_ui_widget.concavity_threshold_spin_box->value();
        bool extract_segmentation = m_ui_widget.segmentation_check_box->isChecked();
        bool use_concavity_colors = m_ui_widget.concavity_colors_check_box->isChecked();

        FacegraphItem* segmentation_item = item;

        // create new item and use it for segmentation if the flag is set
        FacegraphItem* new_item = new FacegraphItem(*item->face_graph());
        new_item->setFlatPlusEdgesMode();
        segmentation_item = new_item;
        
        // decompose mesh
        Facegraph& mesh = *segmentation_item->face_graph();

        typedef typename boost::graph_traits<Facegraph>::face_descriptor face_descriptor;
        typedef std::map<face_descriptor, std::size_t> Clusters_id_map;
        Clusters_id_map clusters_map;
        typedef boost::associative_property_map<Clusters_id_map> Clusters_id_pmap;
        Clusters_id_pmap clusters_pmap(clusters_map);

        Timer timer;

#ifndef CGAL_LINKED_WITH_TBB
        std::cout << "Running sequentially. For performance reasons it's recommended to run in parralel using TBB." << std::endl;
#endif

        timer.start();
        std::size_t clusters_num = CGAL::convex_decomposition<Facegraph, Clusters_id_pmap, CGAL::Parallel_tag>(mesh, clusters_pmap, concavity_threshold, min_number_of_clusters);
        timer.stop();
        
        std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
        std::cout << "Number of clusters: " << clusters_num << std::endl;
      
        // extract segmentation if the flag is set
        if (extract_segmentation)
        {
            // colorize segmentation
            typedef typename boost::property_map<Facegraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
            Patch_id_pmap patch_pmap = get(CGAL::face_patch_id_t<int>(), mesh);

            std::vector<QColor>& colors = segmentation_item->color_vector();
            colors.clear();

            BOOST_FOREACH(face_descriptor face, faces(mesh))
            {
                put(patch_pmap, face, static_cast<int>(clusters_map[face]));
            }
            for (std::size_t i = 0; i < clusters_num; ++i)
            {
                QColor color(CGAL::get_default_random().get_int(41, 255),
                             CGAL::get_default_random().get_int(41, 255),
                             CGAL::get_default_random().get_int(41, 255));
                colors.push_back(color);
            }
            segmentation_item->setItemIsMulticolor(true);
            segmentation_item->setName(tr("%1-segmentation-[%2,%3]").arg(clusters_num).arg(concavity_threshold).arg(min_number_of_clusters));

            // add to the scene
            scene->addItem(segmentation_item);
            
            // refresh item
            segmentation_item->invalidateOpenGLBuffers();
            scene->itemChanged(scene->item_id(segmentation_item));
        }

        // extract decomposition


        // setup default view        
        item->setVisible(false);
        segmentation_item->setVisible(false);
        scene->setSelectedItem(scene->item_id(segmentation_item));

        QApplication::restoreOverrideCursor();
    }
};

#include "Approx_decomposition_plugin.moc"
