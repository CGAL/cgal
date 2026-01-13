#include <QtCore/qglobal.h>

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>


#include "Scene_surface_mesh_item.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/utility.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <QElapsedTimer>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QDialog>
#include <QtPlugin>
#include <QMessageBox>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>
#include <cmath>
#include <unordered_set>

#include "ui_Remesh_planar_patches_dialog.h"


namespace PMP = CGAL::Polygon_mesh_processing;

template <class SM>
struct Mesh_map
{
  typedef SM value_type;
  typedef SM& reference;
  typedef Scene_surface_mesh_item* key_type;
  typedef boost::lvalue_property_map_tag category;

  SM& operator[](Scene_surface_mesh_item* poly_item) const
  {
    return *poly_item->polyhedron();
  }
};

using namespace CGAL::Three;
class CGAL_Lab_remesh_planar_patches_plugin :
  public QObject,
  public CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0" FILE "remesh_planar_patches_plugin.json")

  typedef Scene_surface_mesh_item::Face_graph Mesh;
  typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

  typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionRemeshPlanarPatches_ = new QAction("Remesh Planar Patches", mw);
    actionRemeshPlanarPatches_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionRemeshPlanarPatches_) {
      connect(actionRemeshPlanarPatches_, SIGNAL(triggered()),
        this, SLOT(remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRemeshPlanarPatches_;
  }

  bool applicable(QAction*) const
  {
    if (scene->selectionIndices().size() == 1)
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

    for(int index : scene->selectionIndices())
    {
      if (qobject_cast<Scene_surface_mesh_item*>(scene->item(index)))
        return true;
    }
    return false;
  }

public Q_SLOTS:
  void remeshing()
  {
    if (scene->selectionIndices().size() > 1)
    {
      // Create dialog box
      QDialog dialog(mw);
      ui.setupUi(&dialog);
      connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
      connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

      ui.create_new_item_checkbox->setEnabled(false);
      ui.use_region_growing_checkbox->setEnabled(false);

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }

      bool do_not_triangulate_faces = ui.notriangulation_checkbox->isChecked();
      double cos_threshold = ui.cos_dspinbox->value();

      std::vector<Scene_surface_mesh_item*> meshes;
      for(int index : scene->selectionIndices())
      {
        Scene_surface_mesh_item* poly_item =
          qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

        if (poly_item != nullptr)
          meshes.push_back(poly_item);
      }

      PMP::decimate_meshes_with_common_interfaces(meshes, -cos_threshold, Mesh_map<Mesh>(), do_not_triangulate_faces);

      for (Scene_surface_mesh_item* poly_item : meshes)
      {
        poly_item->invalidateOpenGLBuffers();
        Q_EMIT poly_item->itemChanged();
      }

      return;
    }

    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_surface_mesh_item* poly_item =
      qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

    if (poly_item)
    {
      // Create dialog box
      QDialog dialog(mw);
      ui.setupUi(&dialog);
      connect(ui.buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
      connect(ui.buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

      // Get values
      int i = dialog.exec();
      if (i == QDialog::Rejected)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }

      bool do_not_triangulate_faces = ui.notriangulation_checkbox->isChecked();
      bool create_new_item = ui.create_new_item_checkbox->isChecked();
      double cos_threshold = ui.cos_dspinbox->value();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QElapsedTimer time;
      time.start();

      Mesh& pmesh = *poly_item->polyhedron();

      if (!CGAL::is_triangle_mesh(pmesh))
      {
        QApplication::restoreOverrideCursor();
        if (QMessageBox::Ok ==
            QMessageBox::question(mw, tr("Error - Triangulate Faces?"),
              tr("The input mesh is not a triangulated surface mesh.\n"
                 "Do you wish to triangulate faces first, or cancel remeshing ?"),
              (QMessageBox::Ok | QMessageBox::Cancel), QMessageBox::Ok))
        {
          QApplication::setOverrideCursor(Qt::WaitCursor);
          PMP::triangulate_faces(pmesh);
        }
        else
        {
          return;
        }
      }

      if (ui.use_region_growing_checkbox->isChecked())
      {
        typedef boost::property_map<Mesh, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
        Patch_id_pmap in_fpmap = get(CGAL::face_patch_id_t<int>(), pmesh);
        std::vector<std::size_t> corner_id_map(num_vertices(pmesh), -1);
        std::vector<bool> ecm(num_edges(pmesh), false);
        boost::vector_property_map<CGAL::Epick::Vector_3> normal_map;
        std::size_t nb_regions =
          PMP::region_growing_of_planes_on_faces(pmesh,
                                                 in_fpmap,
                                                 CGAL::parameters::cosine_of_maximum_angle(cos_threshold).
                                                                   region_primitive_map(normal_map).
                                                                   maximum_distance(ui.dist_dspinbox->value()).
                                                                   postprocess_regions(ui.postprocess_regions_checkbox->isChecked()));
        std::size_t nb_corners =
          PMP::detect_corners_of_regions(pmesh,
                                         in_fpmap,
                                         nb_regions,
                                         CGAL::make_random_access_property_map(corner_id_map),
                                         CGAL::parameters::cosine_of_maximum_angle(cos_threshold).
                                                           maximum_distance(ui.dist_dspinbox->value()).
                                                           edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));
        if (create_new_item)
        {
          Scene_surface_mesh_item* new_item=new Scene_surface_mesh_item();
          Mesh& out = *new_item->polyhedron();

          Patch_id_pmap out_fpmap = get(CGAL::face_patch_id_t<int>(), out);
//TODO: use the return type
          PMP::remesh_almost_planar_patches(pmesh,
                                            out,
                                            nb_regions, nb_corners,
                                            in_fpmap,
                                            CGAL::make_random_access_property_map(corner_id_map),
                                            CGAL::make_random_access_property_map(ecm),
                                            CGAL::parameters::patch_normal_map(normal_map),
                                            CGAL::parameters::do_not_triangulate_faces(do_not_triangulate_faces).face_patch_map(out_fpmap));


          new_item->setName(tr("%1_remeshed").arg(poly_item->name()));
          scene->setSelectedItem( scene->addItem(new_item) );

          poly_item->setItemIsMulticolor(true);
          poly_item->computeItemColorVectorAutomatically(true);
          poly_item->invalidateOpenGLBuffers();
          Q_EMIT poly_item->itemChanged();
          new_item->setItemIsMulticolor(true);
          new_item->computeItemColorVectorAutomatically(false);
          new_item->color_vector()=poly_item->color_vector(); // colors are not deterministic
          new_item->invalidateOpenGLBuffers();
          Q_EMIT new_item->itemChanged();
        }
        else
        {
          PMP::remesh_almost_planar_patches(pmesh,
                                            pmesh,
                                            nb_regions, nb_corners,
                                            in_fpmap,
                                            CGAL::make_random_access_property_map(corner_id_map),
                                            CGAL::make_random_access_property_map(ecm),
                                            CGAL::parameters::patch_normal_map(normal_map),
                                            CGAL::parameters::visitor([](Mesh& pmesh){pmesh.clear_without_removing_property_maps ();})
                                                             .do_not_triangulate_faces(do_not_triangulate_faces));
          pmesh.remove_property_map<Mesh::Face_index, int>(in_fpmap);
          poly_item->invalidateOpenGLBuffers();

          Q_EMIT poly_item->itemChanged();
        }
      }
      else
      {
        if (create_new_item)
        {
          typedef boost::property_map<Mesh, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
          Scene_surface_mesh_item* new_item=new Scene_surface_mesh_item();
          Mesh& out = *new_item->polyhedron();


          Patch_id_pmap in_fpmap = get(CGAL::face_patch_id_t<int>(), pmesh);
          Patch_id_pmap out_fpmap = get(CGAL::face_patch_id_t<int>(), out);

          PMP::remesh_planar_patches(pmesh,
                                     out,
                                     CGAL::parameters::cosine_of_maximum_angle(cos_threshold)
                                          .face_patch_map(in_fpmap),
                                     CGAL::parameters::face_patch_map(out_fpmap)
                                                      .do_not_triangulate_faces(do_not_triangulate_faces));


          new_item->setName(tr("%1_remeshed").arg(poly_item->name()));
          scene->setSelectedItem( scene->addItem(new_item) );

          poly_item->setItemIsMulticolor(true);
          poly_item->computeItemColorVectorAutomatically(true);
          poly_item->invalidateOpenGLBuffers();
          Q_EMIT poly_item->itemChanged();
          new_item->setItemIsMulticolor(true);
          new_item->computeItemColorVectorAutomatically(false);
          new_item->color_vector()=poly_item->color_vector(); // colors are not deterministic
          new_item->invalidateOpenGLBuffers();
          Q_EMIT new_item->itemChanged();

        }
        else
        {
          PMP::remesh_planar_patches(pmesh,
                                     pmesh,
                                     CGAL::parameters::cosine_of_maximum_angle(cos_threshold),
                                     CGAL::parameters::visitor([](Mesh& pmesh){pmesh.clear_without_removing_property_maps ();})
                                                      .do_not_triangulate_faces(do_not_triangulate_faces));

          poly_item->invalidateOpenGLBuffers();

          Q_EMIT poly_item->itemChanged();
        }
      }
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      // default cursor
      QApplication::restoreOverrideCursor();
    }
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;

  QAction* actionRemeshPlanarPatches_;
  Ui::Remesh_planar_patches_dialog ui;
}; // end CGAL_Lab_remesh_planar_patches_plugin

#include "Remesh_planar_patches_plugin.moc"
