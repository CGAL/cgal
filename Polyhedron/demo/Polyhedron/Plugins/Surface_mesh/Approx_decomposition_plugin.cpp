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
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/iterator.h>
#include <boost/unordered_map.hpp>

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
        m_decomposition_action = new QAction("Approximate Convex Decomposition", mw);
        m_decomposition_action->setProperty("subMenuName", "Triangulated Surface Mesh Segmentation");
        m_decomposition_action->setObjectName("decomposition_action");

        autoConnectActions();

        // colors
        init_gradient_colors();

        // ui
        m_dock_widget = new QDockWidget("Approximate convex decomposition widget", mw);
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

    std::vector<QColor> m_gradient_colors;

    template <class FacegraphItem>
    void decompose(FacegraphItem* item)
    {
        typedef typename FacegraphItem::Face_graph Facegraph;
        typedef typename boost::graph_traits<Facegraph>::face_descriptor face_descriptor;
        typedef typename boost::graph_traits<Facegraph>::vertex_descriptor vertex_descriptor;
        
        typedef typename boost::property_map<Facegraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
        
        QApplication::setOverrideCursor(Qt::WaitCursor);

        // parameters
        std::size_t min_number_of_clusters = m_ui_widget.min_clusters_spin_box->value();
        double concavity_threshold = m_ui_widget.concavity_threshold_spin_box->value();
        bool extract_segmentation = m_ui_widget.segmentation_check_box->isChecked();
        bool use_concavity_colors = m_ui_widget.concavity_colors_check_box->isChecked();

        // create new item and use it for segmentation if the flag is set
        FacegraphItem* segmentation_item = new FacegraphItem(*item->face_graph());
        segmentation_item->setFlatPlusEdgesMode();
        
        // decompose mesh
        Facegraph& mesh = *segmentation_item->face_graph();

        typedef std::map<face_descriptor, std::size_t> Clusters_id_map;
        Clusters_id_map clusters_map;
        typedef boost::associative_property_map<Clusters_id_map> Clusters_id_pmap;
        Clusters_id_pmap clusters_pmap(clusters_map);

        Timer timer;

#ifndef CGAL_LINKED_WITH_TBB
        std::cout << "Running sequentially. For performance reasons it's recommended to run in parralel using TBB." << std::endl;
        typedef CGAL::Sequential_tag Concurrency_tag;
#else
        typedef CGAL::Parallel_tag Concurrency_tag;
#endif

        timer.start();
        std::size_t clusters_num = CGAL::approximate_convex_decomposition<Concurrency_tag>(mesh, clusters_pmap, concavity_threshold, min_number_of_clusters);
        timer.stop();
        
        std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
        std::cout << "Number of clusters: " << clusters_num << std::endl;
     
        std::vector<QColor> colors;
        generate_colors(colors, clusters_num);

        // extract segmentation if the flag is set
        if (extract_segmentation)
        {
            // colorize segmentation
            Patch_id_pmap patch_pmap = get(CGAL::face_patch_id_t<int>(), mesh);

            BOOST_FOREACH(face_descriptor face, faces(mesh))
            {
                put(patch_pmap, face, static_cast<int>(clusters_map[face]));
            }

            set_color_read_only(segmentation_item);
            segmentation_item->color_vector() = colors;
            segmentation_item->setItemIsMulticolor(true);
            segmentation_item->setName(tr("%1-segmentation-[%2,%3]").arg(clusters_num).arg(concavity_threshold).arg(min_number_of_clusters));

            // add to the scene
            scene->addItem(segmentation_item);
            
            // refresh item
            segmentation_item->invalidateOpenGLBuffers();
            scene->itemChanged(scene->item_id(segmentation_item));
        }

        // extract decomposition
        FacegraphItem* decomposition_item = new FacegraphItem();
        {
            decomposition_item->setFlatPlusEdgesMode();
            Facegraph& decomposition_mesh = *decomposition_item->face_graph();

            typedef CGAL::Face_filtered_graph<Facegraph> Filtered_graph;

            // add convex hulls
            for (std::size_t i = 0; i < clusters_num; ++i)
            {
                // construct convex hull of the i-th cluster
                Filtered_graph filtered_mesh(mesh, i, clusters_pmap);
                Facegraph cluster;
                CGAL::copy_face_graph(filtered_mesh, cluster);

                Facegraph conv_hull;
                std::vector<Point_3> pts;

                if (num_vertices(cluster) > 3)
                {
                    BOOST_FOREACH(vertex_descriptor vert, vertices(cluster))
                    {
                        pts.push_back(get(CGAL::vertex_point, cluster)[vert]);
                    }

                    CGAL::convex_hull_3(pts.begin(), pts.end(), conv_hull);
                }
                else
                {
                    conv_hull = cluster;
                }

                boost::unordered_map<face_descriptor, face_descriptor> f2f;
                typedef std::pair<face_descriptor, face_descriptor> Faces_pair;

                // add the convex hull
                CGAL::copy_face_graph(conv_hull, decomposition_mesh, CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(), std::inserter(f2f, f2f.end()));
                
                // assign patch id to the convex hull
                Patch_id_pmap patch_pmap = get(CGAL::face_patch_id_t<int>(), decomposition_mesh);

                BOOST_FOREACH(Faces_pair faces_pair, f2f)
                {
                    put(patch_pmap, faces_pair.second, static_cast<int>(i));
                }
            }
          
            set_color_read_only(decomposition_item);
            if (use_concavity_colors)
            {

                std::vector<QColor>& colors = decomposition_item->color_vector();
                for (std::size_t i = 0; i < clusters_num; ++i)
                {
                    double concavity = CGAL::concavity_value<Concurrency_tag>(mesh, clusters_pmap, i);
                    int step = std::min(255, int(concavity / concavity_threshold * 255));
                    colors.push_back(m_gradient_colors[step]);
                } 
            }
            else
            {
                decomposition_item->color_vector() = colors;
            }
            decomposition_item->setItemIsMulticolor(true);
            decomposition_item->setName(tr("%1-decomposition-[%2,%3]").arg(clusters_num).arg(concavity_threshold).arg(min_number_of_clusters));

            // add to the scene
            scene->addItem(decomposition_item);
                
            // refresh item
            decomposition_item->invalidateOpenGLBuffers();
            scene->itemChanged(scene->item_id(decomposition_item));
        }

        // setup default view        
        item->setVisible(false);
        segmentation_item->setVisible(false);
        scene->setSelectedItem(scene->item_id(item));

        if (!extract_segmentation)
        {
            delete segmentation_item;
        }

        QApplication::restoreOverrideCursor();
    }

    void generate_colors(std::vector<QColor>& colors, std::size_t cnt)
    {
        colors.clear();

        const std::size_t custom_cnt = 14;
        colors.push_back(QColor(173, 35, 35));
        colors.push_back(QColor(87, 87, 87));
        colors.push_back(QColor(42, 75, 215));
        colors.push_back(QColor(29, 105, 20));
        colors.push_back(QColor(129, 74, 25));
        colors.push_back(QColor(129, 38, 192));
        colors.push_back(QColor(160, 160, 160));
        colors.push_back(QColor(129, 197, 122));
        colors.push_back(QColor(157, 175, 255));
        colors.push_back(QColor(41, 208, 208));
        colors.push_back(QColor(255, 146, 51));
        colors.push_back(QColor(255, 238, 51));
        colors.push_back(QColor(233, 222, 187));
        colors.push_back(QColor(255, 205, 243));

        colors.resize(std::min(cnt, custom_cnt));

        for (std::size_t i = custom_cnt; i < cnt; ++i)
        {
            QColor color(CGAL::get_default_random().get_int(41, 255),
                         CGAL::get_default_random().get_int(41, 255),
                         CGAL::get_default_random().get_int(41, 255));
            colors.push_back(color);
        }
        
    }

    void init_gradient_colors()
    {
        m_gradient_colors.resize(256);
        int r = 0, g = 0, b = 255;
        for (int i = 0; i <= 255; ++i)
        {
            if (i > 128 && i <= 192) { r = static_cast<int>( ((i - 128) / (192.0 - 128)) * 255 ); }
            if (i > 0 && i <= 98)    { g = static_cast<int>( ((i) / (98.0)) * 255 ); }
            if (i > 191 && i <=255)  { g = 255 - static_cast<int>( ((i - 191) / (255.0 - 191)) * 255 ); }
            if (i > 64 && i <= 127)  { b = 255 - static_cast<int>( ((i - 64) / (127.0 - 64)) * 255 ); }
            m_gradient_colors[i] = QColor(r, g, b);
        }
    }

    void set_color_read_only(Scene_polyhedron_item* poly)
    {
        poly->set_color_vector_read_only(true);
    }

    void set_color_read_only(Scene_surface_mesh_item*)
    {}
};

#include "Approx_decomposition_plugin.moc"
