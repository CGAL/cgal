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

    actionIsotropicRemeshing_ = new QAction("Isotropic remeshing", mw);
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
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void isotropic_remeshing()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if (poly_item || selection_item)
    {
      // Create dialog box
      QDialog dialog(mw);
      Ui::Isotropic_remeshing_dialog ui;
      ui.setupUi(&dialog);
      connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
      connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

      //connect checkbox to spinbox
      connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
              ui.nbIterations_spinbox, SLOT(setDisabled(bool)));
      connect(ui.splitEdgesOnly_checkbox, SIGNAL(toggled(bool)),
              ui.protect_checkbox, SLOT(setDisabled(bool)));

      //Set default parameters
      bool p_ = (poly_item != NULL);
      Scene_interface::Bbox bbox = p_ ? poly_item->bbox() : selection_item->bbox();
      ui.objectName->setText(p_ ? poly_item->name() : selection_item->name());
      ui.objectNameSize->setText(
        tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
        .arg(bbox.width(),  0, 'g', 3)
        .arg(bbox.height(), 0, 'g', 3)
        .arg(bbox.depth(),  0, 'g', 3));

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

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }
      bool edges_only = ui.splitEdgesOnly_checkbox->isChecked();
      double target_length = ui.edgeLength_dspinbox->value();
      int nb_iter = ui.nbIterations_spinbox->value();
      bool protect = ui.protect_checkbox->isChecked();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      time.start();

      typedef boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
      typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
      typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

      const Polyhedron& pmesh = p_ ? *poly_item->polyhedron()
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
        std::vector<edge_descriptor> updated_selected_edges;
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
          CGAL::Polygon_mesh_processing::split_long_edges(
            *selection_item->polyhedron()
            , edges
            , target_length
            , std::back_inserter(updated_selected_edges)
            , PMP::parameters::geom_traits(Kernel()));
        }
        else
        {
        std::vector<bool> selected(
          selection_item->polyhedron()->size_of_halfedges()/2,
          false);

        if (selection_item->selected_edges.empty())
          CGAL::Polygon_mesh_processing::isotropic_remeshing(
          *selection_item->polyhedron()
          , selection_item->selected_facets
          , target_length
          , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
          .protect_constraints(protect));
        else
          CGAL::Polygon_mesh_processing::isotropic_remeshing(
         *selection_item->polyhedron()
         , selection_item->selected_facets
         , target_length
         , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
         .protect_constraints(protect)
         .edge_is_constrained_map(selection_item->selected_edges_pmap(selected))
         );

        }
        selection_item->poly_item_changed();
        selection_item->clear_all();
        selection_item->selected_edges.insert(updated_selected_edges.begin(),
                                              updated_selected_edges.end());
        selection_item->changed_with_poly_item();
      }
      else if (poly_item)
      {
        if (edges_only)
        {
          std::vector<halfedge_descriptor> border;
          CGAL::Polygon_mesh_processing::border_halfedges(
            faces(*poly_item->polyhedron()),
            std::back_inserter(border),
            pmesh);
          std::vector<edge_descriptor> border_edges;
          BOOST_FOREACH(halfedge_descriptor h, border)
            border_edges.push_back(edge(h, pmesh));

          CGAL::Polygon_mesh_processing::split_long_edges(*poly_item->polyhedron()
                                                        , border_edges
                                                        , target_length);
        }
        else
        {
        CGAL::Polygon_mesh_processing::isotropic_remeshing(
         *poly_item->polyhedron()
         , faces(*poly_item->polyhedron())
         , target_length
         , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
         .protect_constraints(protect));
        }
        poly_item->invalidate_buffers();
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

private:
  QAction* actionIsotropicRemeshing_;

}; // end Polyhedron_demo_isotropic_remeshing_plugin

//Q_EXPORT_PLUGIN2(Polyhedron_demo_isotropic_remeshing_plugin,
//                 Polyhedron_demo_isotropic_remeshing_plugin)

#include "Isotropic_remeshing_plugin.moc"
