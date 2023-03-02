#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>


#include "Scene_surface_mesh_item.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/utility.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

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




using namespace CGAL::Three;
class Polyhedron_demo_remesh_planar_patches_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "isotropic_remeshing_plugin.json")

  typedef Scene_surface_mesh_item Scene_surface_mesh_item;
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

    Q_FOREACH(int index, scene->selectionIndices())
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

      struct Mesh_map
      {
        typedef Mesh value_type;
        typedef Mesh& reference;
        typedef Scene_surface_mesh_item* key_type;
        typedef boost::lvalue_property_map_tag category;

        Mesh& operator[](Scene_surface_mesh_item* poly_item) const
        {
          return *poly_item->polyhedron();
        }
      };

      CGAL::Polygon_mesh_processing::decimate_meshes_with_common_interfaces(meshes, cos_threshold, Mesh_map(), do_not_triangulate_faces);

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
          CGAL::Polygon_mesh_processing::triangulate_faces(pmesh);
        }
        else
        {
          return;
        }
      }

      if (create_new_item)
      {
        typedef boost::property_map<Mesh, CGAL::face_patch_id_t<int> >::type Patch_id_pmap;
        Scene_surface_mesh_item* new_item=new Scene_surface_mesh_item();
        Mesh& out = *new_item->polyhedron();


        Patch_id_pmap in_fpmap = get(CGAL::face_patch_id_t<int>(), pmesh);
        Patch_id_pmap out_fpmap = get(CGAL::face_patch_id_t<int>(), out);

        CGAL::Polygon_mesh_processing::remesh_planar_patches(pmesh,
                                                             out,
                                                             CGAL::parameters::cosine_of_maxium_angle(cos_threshold)
                                                                  .face_patch_map(in_fpmap)
                                                                  .do_not_triangulate_faces(do_not_triangulate_faces),
                                                             CGAL::parameters::face_patch_map(out_fpmap));


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
        CGAL::Polygon_mesh_processing::remesh_planar_patches(pmesh,
                                                             pmesh,
                                                             CGAL::parameters::cosine_of_maxium_angle(cos_threshold)
                                                                              .do_not_triangulate_faces(do_not_triangulate_faces),
                                                             CGAL::parameters::visitor([](Mesh& pmesh){pmesh.clear_without_removing_property_maps ();}));

        poly_item->invalidateOpenGLBuffers();

        Q_EMIT poly_item->itemChanged();
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
}; // end Polyhedron_demo_remesh_planar_patches_plugin

#include "Remesh_planar_patches_plugin.moc"
