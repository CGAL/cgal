#include "ui_Approximate_convex_segmentation_widget.h"

#include "Scene_surface_mesh_item.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#include "Scene.h"
#include "Color_map.h"

//#define CGAL_APPROXIMATE_CONVEX_SEGMENTATION_VERBOSE
#include <CGAL/approximate_convex_segmentation.h>
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
#include <cmath>

using namespace CGAL::Three;

typedef CGAL::Real_timer Timer;

#ifndef CGAL_LINKED_WITH_TBB
    typedef CGAL::Sequential_tag Concurrency_tag;
#else
    typedef CGAL::Parallel_tag Concurrency_tag;
#endif

class Polyhedron_demo_approximate_convex_segmentation_plugin : public QObject, public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  QList<QAction*> actions() const
  {
    return QList<QAction*>() << m_segmentation_action;
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
    m_segmentation_action = new QAction("Approximate Convex Segmentation", mw);
    m_segmentation_action->setProperty("subMenuName", "Triangulated Surface Mesh Segmentation");
    m_segmentation_action->setObjectName("approximate_convex_segmentation_action");

    autoConnectActions();

    // colors
    init_gradient_colors();

    // ui
    m_segmentation_widget = new QDockWidget("Approximate convex segmentation widget", mw);
    m_segmentation_widget->setVisible(false); // do not show at the beginning
    m_segmentation_ui.setupUi(m_segmentation_widget);
    this->mw->addDockWidget(Qt::LeftDockWidgetArea, m_segmentation_widget);

    // signal-slot bindings
    connect(m_segmentation_ui.segmentize_button, SIGNAL(clicked()), this, SLOT(on_segmentize_button_clicked()));
    connect(m_segmentation_ui.concavity_values_button, SIGNAL(clicked()), this, SLOT(on_concavity_values_button_clicked()));
  }

  virtual void closure()
  {
    m_segmentation_widget->hide();
  }

public Q_SLOTS:
  void on_approximate_convex_segmentation_action_triggered()
  {
    m_segmentation_widget->show();
  }

  void on_segmentize_button_clicked()
  {
    std::cout << "Segmentizing..." << std::endl;
    
    CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
    
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if (item)
    {
      segmentize(item);
      return;
    }
    
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
    if (sm_item)
    {
      segmentize(sm_item);
      return;
    }
  }

  void on_concavity_values_button_clicked()
  {
    std::cout << "Computing concavity values..." << std::endl;

    CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
    
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(scene->item(index));
    if (item)
    {
      compute_concavity_values(item);
      return;
    }
    
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(index));
    if (sm_item)
    {
      compute_concavity_values(sm_item);
      return;
    }
  }

private:
  QAction* m_segmentation_action;
  QDockWidget* m_segmentation_widget;
  Ui::approx_segmentation m_segmentation_ui;

  std::vector<QColor> m_gradient_colors;

  template <class FacegraphItem>
  void segmentize(FacegraphItem* item)
  {
    typedef typename FacegraphItem::Face_graph Facegraph;
    typedef typename boost::graph_traits<Facegraph>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<Facegraph>::vertex_descriptor vertex_descriptor;
    
    typedef typename boost::property_map<Facegraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // parameters
    std::size_t min_number_of_segments = m_segmentation_ui.min_segments_spin_box->value();
    double concavity_threshold = m_segmentation_ui.concavity_threshold_spin_box->value();
    bool extract_convex_hulls = m_segmentation_ui.convex_hulls_check_box->isChecked();
    bool use_concavity_colors = m_segmentation_ui.concavity_colors_check_box->isChecked();
    bool use_shortest_method = m_segmentation_ui.shortest_method_check_box->isChecked();
    bool enable_postprocessing_segments = m_segmentation_ui.postprocess_segments_check_box->isChecked();
    double small_segment_threshold = m_segmentation_ui.small_segment_threshold_spin_box->value();

    // create a new item and use it for segmentation
    FacegraphItem* segmentation_item = new FacegraphItem(*item->face_graph());
    segmentation_item->setFlatPlusEdgesMode();
    
    // segmenation mesh
    Facegraph& segmentation_mesh = *segmentation_item->face_graph();

    typedef std::map<face_descriptor, std::size_t> Segments_id_map;
    Segments_id_map segments_map;
    typedef boost::associative_property_map<Segments_id_map> Segments_id_pmap;
    Segments_id_pmap segments_pmap(segments_map);

    boost::vector_property_map<Facegraph> convex_hulls_pmap;

    Timer timer;

#ifndef CGAL_LINKED_WITH_TBB
    std::cout << "Running sequentially. For performance reasons it's recommended to run in parralel using TBB." << std::endl;
#endif

    timer.start();
    std::size_t segments_num =
      CGAL::approximate_convex_segmentation<Concurrency_tag>(segmentation_mesh,
                                                             segments_pmap,
                                                             concavity_threshold,
                                                             CGAL::parameters::min_number_of_segments(min_number_of_segments).
                                                             segments_convex_hulls(convex_hulls_pmap).
                                                             use_closest_point(use_shortest_method).
                                                             postprocess_segments(enable_postprocessing_segments).
                                                             small_segment_threshold(small_segment_threshold));
    timer.stop();
    
    std::cout << "Elapsed time: " << timer.time() << " seconds" << std::endl;
    std::cout << "Number of segments: " << segments_num << std::endl;
  
    std::vector<QColor> colors;
    generate_colors(colors, segments_num);

    // extract segmentation
    {
      // colorize segmentation
      Patch_id_pmap patch_pmap = get(CGAL::face_patch_id_t<int>(), segmentation_mesh);

      BOOST_FOREACH(face_descriptor face, faces(segmentation_mesh))
      {
        put(patch_pmap, face, static_cast<int>(segments_map[face]));
      }

      set_color_read_only(segmentation_item);
      segmentation_item->color_vector() = colors;
      segmentation_item->setItemIsMulticolor(true);
      segmentation_item->setName(tr("%1-segmentation-[%2,%3]").arg(segments_num).arg(concavity_threshold).arg(min_number_of_segments));

      // add to the scene
      scene->addItem(segmentation_item);
      segmentation_item->setFlatPlusEdgesMode();
      
      // refresh item
      segmentation_item->invalidateOpenGLBuffers();
      scene->itemChanged(scene->item_id(segmentation_item));
    }

    // extract convex hulls if the corresponding flag is set 
    if (extract_convex_hulls)
    {
      FacegraphItem* convex_hulls_item = new FacegraphItem();
      
      convex_hulls_item->setFlatPlusEdgesMode();
      Facegraph& convex_hulls_mesh = *convex_hulls_item->face_graph();

      typedef CGAL::Face_filtered_graph<Facegraph> Filtered_graph;

      // add convex hulls
      for (std::size_t i = 0; i < segments_num; ++i)
      {
        Facegraph& convex_hull = convex_hulls_pmap[i];

        boost::unordered_map<face_descriptor, face_descriptor> f2f;
        typedef std::pair<face_descriptor, face_descriptor> Faces_pair;

        // add the convex hull
        CGAL::copy_face_graph(convex_hull, convex_hulls_mesh, CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(), std::inserter(f2f, f2f.end()));
        
        // assign patch id to the convex hull
        Patch_id_pmap patch_pmap = get(CGAL::face_patch_id_t<int>(), convex_hulls_mesh);

        BOOST_FOREACH(Faces_pair faces_pair, f2f)
        {
          put(patch_pmap, faces_pair.second, static_cast<int>(i));
        }
      }
     
      set_color_read_only(convex_hulls_item);
      if (use_concavity_colors)
      {
        std::vector<QColor>& colors = convex_hulls_item->color_vector();
        for (std::size_t i = 0; i < segments_num; ++i)
        {
          double concavity = CGAL::concavity_values<Concurrency_tag>(segmentation_mesh, segments_pmap, i);
          int step = std::min(255, int(concavity / concavity_threshold * 255));
          colors.push_back(m_gradient_colors[step]);
        } 
      }
      else
      {
        convex_hulls_item->color_vector() = colors;
      }
      convex_hulls_item->setItemIsMulticolor(true);
      convex_hulls_item->setName(tr("%1-convex-hulls-[%2,%3]").arg(segments_num).arg(concavity_threshold).arg(min_number_of_segments));

      // add to the scene
      scene->addItem(convex_hulls_item);
      convex_hulls_item->setVisible(false);
      convex_hulls_item->setFlatPlusEdgesMode();
        
      // refresh item
      convex_hulls_item->invalidateOpenGLBuffers();
      scene->itemChanged(scene->item_id(convex_hulls_item));
    }

    // setup default view
    item->setVisible(false);
    scene->setSelectedItem(scene->item_id(item));

    QApplication::restoreOverrideCursor();
  }

  template <class FacegraphItem>
  void compute_concavity_values(FacegraphItem* item)
  {
    typedef typename FacegraphItem::Face_graph Facegraph;
    typedef typename boost::graph_traits<Facegraph>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<Facegraph>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<Facegraph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<Facegraph, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    // parameters
    bool use_shortest_method = m_segmentation_ui.shortest_method_check_box->isChecked();

    // create new item for concavity values
    FacegraphItem* concavity_values_item = new FacegraphItem(*item->face_graph());
      
    concavity_values_item->setFlatPlusEdgesMode();
    Facegraph& concavity_values_mesh = *concavity_values_item->face_graph();

    //for each input vertex compute its distance to the convex hull
    typedef std::map<vertex_descriptor, double> Vertex_double_map;
    Vertex_double_map distances_map;
    boost::associative_property_map<Vertex_double_map> distances_pmap(distances_map);

    double concavity_value = CGAL::concavity_values<Concurrency_tag>(concavity_values_mesh, distances_pmap, CGAL::parameters::use_closest_point(use_shortest_method));
    std::cout << "Concavity value of the mesh: " << concavity_value << std::endl;

    // assign patch id to each face and colorize it
    Patch_id_pmap patch_pmap = get(CGAL::face_patch_id_t<int>(), concavity_values_mesh);

    int face_id = 0;
    double max_distance = 0.;
    std::vector<double> distances;
    BOOST_FOREACH(face_descriptor face, faces(concavity_values_mesh))
    {
      put(patch_pmap, face, static_cast<int>(face_id++));

      double dist = 0;
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(face, concavity_values_mesh), concavity_values_mesh))
      {   
        dist += distances_map[target(hd, concavity_values_mesh)];
      }
      dist /= 3.;

      max_distance = std::max(max_distance, dist);
      distances.push_back(dist);
    }
      
    set_color_read_only(concavity_values_item);
    std::vector<QColor>& colors = concavity_values_item->color_vector();
    colors.clear();
    
    for (int i = 0; i < face_id; ++i)
    {
      int step = distances[i] > 0 ? std::min(255, int(std::pow(distances[i] / max_distance, 0.3) * 255)) : 0;
      colors.push_back(m_gradient_colors[step]);
    }
   
    concavity_values_item->setItemIsMulticolor(true);
    concavity_values_item->setName(tr("concavity-values"));

    // add to the scene
    scene->addItem(concavity_values_item);
    item->setVisible(false);
    concavity_values_item->setVisible(true);
    concavity_values_item->setFlatMode();
      
    // refresh item
    concavity_values_item->invalidateOpenGLBuffers();
    scene->itemChanged(scene->item_id(concavity_values_item));

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

#include "Approximate_convex_segmentation_plugin.moc"
