//#define CGAL_PMP_REMESHING_VERBOSE
//#define CGAL_PMP_REMESHING_DEBUG
//#define CGAL_PMP_REMESHING_VERY_VERBOSE

#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QDialog>
#include <QtPlugin>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>

#ifdef CGAL_LINKED_WITH_TBB
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#endif

#include "ui_Isotropic_remeshing_dialog.h"

using namespace CGAL::Three;
class Polyhedron_demo_isotropic_remeshing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionIsotropicRemeshing_ = new QAction("Isotropic Remeshing", mw);
    actionIsotropicRemeshing_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionIsotropicRemeshing_) {
      connect(actionIsotropicRemeshing_, SIGNAL(triggered()),
        this, SLOT(isotropic_remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionIsotropicRemeshing_;
  }

  bool applicable(QAction*) const
  {
    if (scene->selectionIndices().size() == 1)
    {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    }

    Q_FOREACH(int index, scene->selectionIndices())
    {
      //if one polyhedron is found in the selection, it's fine
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
    }
    return false;
  }

public Q_SLOTS:
  void isotropic_remeshing()
  {
    if (scene->selectionIndices().size() > 1)
    {
      isotropic_remeshing_of_several_polyhedra();
      return;
    }
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if (poly_item || selection_item)
    {
      // Create dialog box
      QDialog dialog(mw);
      Ui::Isotropic_remeshing_dialog ui
        = remeshing_dialog(&dialog, poly_item, selection_item);

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }
      bool edges_only = ui.splitEdgesOnly_checkbox->isChecked();
      double target_length = ui.edgeLength_dspinbox->value();
      unsigned int nb_iter = ui.nbIterations_spinbox->value();
      bool protect = ui.protect_checkbox->isChecked();
      bool smooth_features = ui.smooth1D_checkbox->isChecked();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      time.start();

      typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
      typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
      typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

      const Polyhedron& pmesh = (poly_item != NULL)
        ? *poly_item->polyhedron()
        : *selection_item->polyhedron();
      boost::property_map<Polyhedron, CGAL::face_index_t>::type fim
        = get(CGAL::face_index, pmesh);
      unsigned int id = 0;
      BOOST_FOREACH(face_descriptor f, faces(pmesh))
      {
        put(fim, f, id++);
      }

      if (selection_item)
      {
        if (edges_only)
        {
          std::vector<edge_descriptor> edges;
          BOOST_FOREACH(edge_descriptor e, selection_item->selected_edges)
          {
            if (selection_item->selected_facets.find(face(halfedge(e, pmesh), pmesh))
                 != selection_item->selected_facets.end()
             || selection_item->selected_facets.find(face(opposite(halfedge(e, pmesh), pmesh), pmesh))
                 != selection_item->selected_facets.end())
              edges.push_back(e);
          }
          BOOST_FOREACH(face_descriptor f, selection_item->selected_facets)
          {
            BOOST_FOREACH(halfedge_descriptor he, halfedges_around_face(halfedge(f, pmesh), pmesh))
            {
              if (selection_item->selected_facets.find(face(opposite(he, pmesh), pmesh))
                  == selection_item->selected_facets.end())
              edges.push_back(edge(he, pmesh));
            }
          }
          if (!edges.empty())
            CGAL::Polygon_mesh_processing::split_long_edges(
              edges
              , target_length
              , *selection_item->polyhedron()
              , PMP::parameters::geom_traits(Kernel())
              .edge_is_constrained_map(selection_item->constrained_edges_pmap()));
          else
            std::cout << "No selected or boundary edges to be split" << std::endl;
        }
        else
        {
         CGAL::Polygon_mesh_processing::isotropic_remeshing(
           selection_item->selected_facets
         , target_length
         , *selection_item->polyhedron()
         , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
         .protect_constraints(protect)
         .edge_is_constrained_map(selection_item->constrained_edges_pmap())
         .smooth_along_features(smooth_features)
         .vertex_is_constrained_map(selection_item->constrained_vertices_pmap()));
        }
        selection_item->poly_item_changed();
        selection_item->clear<face_descriptor>();
        selection_item->changed_with_poly_item();
      }
      else if (poly_item)
      {
        if (edges_only)
        {
          std::vector<halfedge_descriptor> border;
          CGAL::Polygon_mesh_processing::border_halfedges(
            faces(*poly_item->polyhedron()),
            pmesh,
            std::back_inserter(border));
          std::vector<edge_descriptor> border_edges;
          BOOST_FOREACH(halfedge_descriptor h, border)
            border_edges.push_back(edge(h, pmesh));

          if (!border_edges.empty())
            CGAL::Polygon_mesh_processing::split_long_edges(
              border_edges
              , target_length
              , *poly_item->polyhedron()
              , PMP::parameters::geom_traits(Kernel()));
          else
            std::cout << "No border to be split" << std::endl;
        }
        else
        {
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
           faces(*poly_item->polyhedron())
         , target_length
         , *poly_item->polyhedron()
         , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
         .protect_constraints(protect)
         .smooth_along_features(smooth_features));
        }
        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();
      }
      else{
        std::cout << "Can't remesh that type of thing" << std::endl;
      }
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      // default cursor
      QApplication::restoreOverrideCursor();
    }
  }

  void isotropic_remeshing_of_several_polyhedra()
  {
    // Remeshing parameters
    bool edges_only = false;
    double target_length = 0.;
    unsigned int nb_iter = 1;
    bool protect = false;
    bool smooth_features = true;

    std::vector<Scene_polyhedron_item*> selection;
    BOOST_FOREACH(int index, scene->selectionIndices())
    {
      Scene_polyhedron_item* poly_item =
        qobject_cast<Scene_polyhedron_item*>(scene->item(index));

      if (poly_item == NULL)
      {
        std::cout << scene->item(index)->name().data()
          << " is not a Polyhedron, remeshing skipped\n";
        continue;
      }
      else
      {
        selection.push_back(poly_item);

        if (target_length == 0.)//parameters have not been set yet
        {
        QDialog dialog(mw);
        Ui::Isotropic_remeshing_dialog ui = remeshing_dialog(&dialog, poly_item);
        ui.objectName->setText(QString::number(scene->selectionIndices().size())
          .append(QString(" items to be remeshed")));
        int i = dialog.exec();
        if (i == QDialog::Rejected)
        {
          std::cout << "Remeshing aborted" << std::endl;
          return;
        }

        edges_only = ui.splitEdgesOnly_checkbox->isChecked();
        target_length = ui.edgeLength_dspinbox->value();
        nb_iter = ui.nbIterations_spinbox->value();
        protect = ui.protect_checkbox->isChecked();
        smooth_features = ui.smooth1D_checkbox->isChecked();
        }
      }
    }

    if(target_length == 0.)//parameters have not been set
    {                      // i.e. no item is a polyhedron
      std::cout << "Remeshing aborted" << std::endl;
      return;
    }

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);
    int total_time = 0;

#ifdef CGAL_LINKED_WITH_TBB
    QTime time;
    time.start();

      tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, selection.size()),
        Remesh_polyhedron_item_for_parallel_for<Remesh_polyhedron_item>(
          selection, edges_only, target_length, nb_iter, protect, smooth_features));

    total_time = time.elapsed();
#else
    Remesh_polyhedron_item remesher(edges_only,
      target_length, nb_iter, protect, smooth_features);
    BOOST_FOREACH(Scene_polyhedron_item* poly_item, selection)
    {
      QTime time;
      time.start();

        remesher(poly_item);

      total_time += time.elapsed();
      std::cout << "Remeshing of " << poly_item->name().data()
                << " done in " << time.elapsed() << " ms" << std::endl;
    }
#endif
    std::cout << "Remeshing of all selected items done in "
      << total_time << " ms" << std::endl;

    BOOST_FOREACH(Scene_polyhedron_item* poly_item, selection)
    {
      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();
    }
    
    // default cursor
    QApplication::restoreOverrideCursor();
  }

private:
  struct Remesh_polyhedron_item
  {
    typedef boost::graph_traits<Polyhedron>::edge_descriptor     edge_descriptor;
    typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
    typedef boost::graph_traits<Polyhedron>::face_descriptor     face_descriptor;

    bool edges_only_;
    double target_length_;
    unsigned int nb_iter_;
    bool protect_;
    bool smooth_features_;

  protected:
    void remesh(Scene_polyhedron_item* poly_item) const
    {
      //fill face_index property map
      boost::property_map<Polyhedron, CGAL::face_index_t>::type fim
        = get(CGAL::face_index, *poly_item->polyhedron());
      unsigned int id = 0;
      BOOST_FOREACH(face_descriptor f, faces(*poly_item->polyhedron()))
      { put(fim, f, id++); }

      if (edges_only_)
      {
        std::vector<halfedge_descriptor> border;
        CGAL::Polygon_mesh_processing::border_halfedges(
          faces(*poly_item->polyhedron())
          , *poly_item->polyhedron()
          , std::back_inserter(border));
        std::vector<edge_descriptor> border_edges;
        BOOST_FOREACH(halfedge_descriptor h, border)
          border_edges.push_back(edge(h, *poly_item->polyhedron()));

        CGAL::Polygon_mesh_processing::split_long_edges(
            border_edges
          , target_length_
          , *poly_item->polyhedron());
      }
      else
      {
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
            faces(*poly_item->polyhedron())
          , target_length_
          , *poly_item->polyhedron()
          , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter_)
          .protect_constraints(protect_)
          .smooth_along_features(smooth_features_));
      }
    }

  public:
    Remesh_polyhedron_item(
      const bool edges_only,
      const double target_length,
      const unsigned int nb_iter,
      const bool protect,
      const bool smooth_features)
      : edges_only_(edges_only)
      , target_length_(target_length)
      , nb_iter_(nb_iter)
      , protect_(protect)
      , smooth_features_(smooth_features)
    {}

    Remesh_polyhedron_item(const Remesh_polyhedron_item& remesh)
      : edges_only_(remesh.edges_only_)
      , target_length_(remesh.target_length_)
      , nb_iter_(remesh.nb_iter_)
      , protect_(remesh.protect_)
      , smooth_features_(remesh.smooth_features_)
    {}

    void operator()(Scene_polyhedron_item* poly_item) const
    {
      remesh(poly_item);
    }
  };

#ifdef CGAL_LINKED_WITH_TBB
  template<typename RemeshFunctor>
  struct Remesh_polyhedron_item_for_parallel_for
    : RemeshFunctor
  {
    const std::vector<Scene_polyhedron_item*>& selection_;

  public:
    // Constructor
    Remesh_polyhedron_item_for_parallel_for(
      const std::vector<Scene_polyhedron_item*>& selection,
      const bool edges_only,
      const double target_length,
      const unsigned int nb_iter,
      const bool protect,
      const bool smooth_features)
      : RemeshFunctor(edges_only, target_length, nb_iter, protect, smooth_features)
      , selection_(selection)
    {
      ;
    }

    // Constructor
    Remesh_polyhedron_item_for_parallel_for(
      const Remesh_polyhedron_item_for_parallel_for &remesh)
      : RemeshFunctor(remesh)
      , selection_(remesh.selection_)
    {}

    // operator()
    void operator()(const tbb::blocked_range<size_t>& r) const
    {
      for (size_t i = r.begin(); i != r.end(); ++i)
        RemeshFunctor::remesh(selection_[i]);
    }
  };
#endif

  Ui::Isotropic_remeshing_dialog
  remeshing_dialog(QDialog* dialog,
                   Scene_polyhedron_item* poly_item,
                   Scene_polyhedron_selection_item* selection_item = NULL)
  {
    Ui::Isotropic_remeshing_dialog ui;
    ui.setupUi(dialog);
    connect(ui.buttonBox, SIGNAL(accepted()), dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()), dialog, SLOT(reject()));

    //connect checkbox to spinbox
    connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
            ui.nbIterations_spinbox, SLOT(setDisabled(bool)));
    connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
            ui.protect_checkbox, SLOT(setDisabled(bool)));
    connect(ui.protect_checkbox, SIGNAL(toggled(bool)),
            ui.smooth1D_checkbox, SLOT(setDisabled(bool)));
    connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
            ui.smooth1D_checkbox, SLOT(setDisabled(bool)));

    //Set default parameters
    Scene_interface::Bbox bbox = poly_item != NULL ? poly_item->bbox()
      : (selection_item != NULL ? selection_item->bbox()
        : scene->bbox());
    ui.objectName->setText(poly_item != NULL ? poly_item->name()
      : (selection_item != NULL ? selection_item->name()
        : QString("Remeshing parameters")));

    ui.objectNameSize->setText(
      tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
      .arg(bbox.width(), 0, 'g', 3)
      .arg(bbox.height(), 0, 'g', 3)
      .arg(bbox.depth(), 0, 'g', 3));

    double diago_length = bbox.diagonal_length();
    ui.edgeLength_dspinbox->setDecimals(3);
    ui.edgeLength_dspinbox->setSingleStep(0.001);
    ui.edgeLength_dspinbox->setRange(1e-6 * diago_length, //min
      2.   * diago_length);//max
    ui.edgeLength_dspinbox->setValue(0.05 * diago_length);

    std::ostringstream oss;
    oss << "Diagonal length of the Bbox of the selection to remesh is ";
    oss << diago_length << "." << std::endl;
    oss << "Default is 5% of it" << std::endl;
    ui.edgeLength_dspinbox->setToolTip(QString::fromStdString(oss.str()));

    ui.nbIterations_spinbox->setSingleStep(1);
    ui.nbIterations_spinbox->setRange(1/*min*/, 1000/*max*/);
    ui.nbIterations_spinbox->setValue(1);

    ui.protect_checkbox->setChecked(false);
    ui.smooth1D_checkbox->setChecked(true);

    return ui;
  }


private:
  QAction* actionIsotropicRemeshing_;

}; // end Polyhedron_demo_isotropic_remeshing_plugin

//Q_EXPORT_PLUGIN2(Polyhedron_demo_isotropic_remeshing_plugin,
//                 Polyhedron_demo_isotropic_remeshing_plugin)

#include "Isotropic_remeshing_plugin.moc"
