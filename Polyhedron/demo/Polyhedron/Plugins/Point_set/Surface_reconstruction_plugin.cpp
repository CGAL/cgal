#include <QTime>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QObject>

#include <fstream>

#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "SMesh_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/compute_average_spacing.h>

#include "ui_Surface_reconstruction_plugin.h"
#include "CGAL/Kernel_traits.h"

typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

SMesh* advancing_front (const Point_set& points,
                        double longest_edge,
                        double radius_ratio_bound,
                        double beta,
                        bool structuring);

SMesh* poisson_reconstruct (Point_set& points,
                            Kernel::FT sm_angle, // Min triangle angle (degrees).
                            Kernel::FT sm_radius, // Max triangle size w.r.t. point set average spacing.
                            Kernel::FT sm_distance, // Approximation error w.r.t. point set average spacing.
                            bool conjugate_gradient,
                            bool use_two_passes,
                            bool do_not_fill_holes);

void scale_space (const Point_set& points,
                  std::vector<Scene_polygon_soup_item*>& items,
                  bool jet_smoother, 
                  unsigned int iterations,
                  unsigned int neighbors, unsigned int fitting, unsigned int monge,
                  unsigned int neighborhood_size, unsigned int samples,
                  bool advancing_front_mesher,
                  bool generate_smooth,
                  double longest_edge, double radius_ratio_bound, double beta_angle,
                  bool separate_shells, bool force_manifold);

class Polyhedron_demo_surface_reconstruction_plugin_dialog : public QDialog, private Ui::SurfaceReconstructionDialog
{
  Q_OBJECT
public:
  Polyhedron_demo_surface_reconstruction_plugin_dialog(QWidget* /*parent*/ = 0)
  {
    setupUi(this);
  }

  unsigned int method () const
  {
    return tabWidget->currentIndex();
  }

  void disable_poisson()
  {
    tabWidget->setTabEnabled(1, false);
    tabWidget->setTabText(1, QString("Poisson (requires oriented normals)"));
  }

  void enable_structuring()
  {
    m_use_structuring->setEnabled(true);
    m_use_structuring->setText(QString("Use Point Set Structuring"));
  }
  
  // Advancing front
  double longest_edge () const { return m_longestEdge->value (); }
  double radius_ratio_bound () const { return m_radiusRatioBound->value (); }
  double beta_angle () const { return m_betaAngle->value (); }
  bool structuring() const { return m_use_structuring->isChecked(); }

  // Scale Space
  bool scalespace_js() const { return m_scalespace_jet->isChecked(); }
  unsigned int iterations () const { return m_iterations->value (); }
  unsigned int neighbors () const { return m_neighbors->value(); }
  unsigned int fitting () const { return m_fitting->value(); }
  unsigned int monge () const { return m_monge->value(); }
  unsigned int neighborhood_size () const { return m_neighborhood_size->value (); }
  unsigned int samples () const { return m_samples->value (); }
  bool scalespace_af() const { return m_scalespace_af->isChecked(); }
  bool generate_smoothed () const { return m_genSmooth->isChecked (); }
  double longest_edge_2 () const { return m_longestEdge_2->value (); }
  double radius_ratio_bound_2 () const { return m_radiusRatioBound_2->value (); }
  double beta_angle_2 () const { return m_betaAngle_2->value (); }
  bool separate_shells () const { return m_genShells->isChecked (); }
  bool force_manifold () const { return m_forceManifold->isChecked (); }

  // Poisson
  double angle () const { return m_inputAngle->value (); }
  double radius () const { return m_inputRadius->value (); }
  double distance () const { return m_inputDistance->value (); }
  bool conjugate_gradient() const { return m_conjugate_gradient->isChecked(); }
  bool two_passes () const { return m_inputTwoPasses->isChecked (); }
  bool do_not_fill_holes () const { return m_doNotFillHoles->isChecked (); }
  
};

class Polyhedron_demo_surface_reconstruction_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  QAction* actionSurfaceReconstruction;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionSurfaceReconstruction = new QAction(tr("Surface Reconstruction"), mainWindow);
    actionSurfaceReconstruction->setObjectName("actionSurfaceReconstruction");
    autoConnectActions();

  }

  void advancing_front_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void scale_space_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  void poisson_reconstruction (const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog);
  
  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSurfaceReconstruction;
  }

private:

public Q_SLOTS:
  void on_actionSurfaceReconstruction_triggered();
}; // end class Polyhedron_surface_reconstruction_plugin


void Polyhedron_demo_surface_reconstruction_plugin::on_actionSurfaceReconstruction_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      //generate the dialog box to set the options
      Polyhedron_demo_surface_reconstruction_plugin_dialog dialog;
      dialog.setWindowFlags(Qt::Dialog|Qt::CustomizeWindowHint|Qt::WindowCloseButtonHint);

      if (!pts_item->point_set()->has_normal_map())
        dialog.disable_poisson();
      if (pts_item->point_set()->has_property_map<int> ("shape"))
        dialog.enable_structuring();
      
      if(!dialog.exec())
        return;

      unsigned int method = dialog.method ();
      switch (method)
        {
        case 0:
          advancing_front_reconstruction (dialog);
          break;
        case 1:
          poisson_reconstruction (dialog);
          break;
        case 2:
          scale_space_reconstruction (dialog);
          break;
        default:
          std::cerr << "Error: unkown method." << std::endl;
          return;
        }
      

    }
}

void Polyhedron_demo_surface_reconstruction_plugin::advancing_front_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      std::cerr << "Advancing front reconstruction... ";
      
      // Reconstruct point set as a polyhedron
      SMesh* mesh = advancing_front (*points,
                                     dialog.longest_edge(),
                                     dialog.radius_ratio_bound(),
                                     CGAL_PI * dialog.beta_angle() / 180.,
                                     dialog.structuring());
      if (mesh)
      {
        // Add polyhedron to scene
        Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(mesh);
        new_item->setName(tr("%1 (advancing front)").arg(pts_item->name()));
        new_item->setColor(Qt::lightGray);
        new_item->invalidateOpenGLBuffers();
        scene->addItem(new_item);

        // Hide point set
        pts_item->setVisible(false);
        scene->itemChanged(index);
      }

      QApplication::restoreOverrideCursor();
    }
}


void Polyhedron_demo_surface_reconstruction_plugin::scale_space_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      std::cout << "Scale scape surface reconstruction...";
      
      std::vector<Scene_polygon_soup_item*> reco_items;

      scale_space (*points, reco_items,
                   dialog.scalespace_js(),
                   dialog.iterations(),
                   dialog.neighbors(), dialog.fitting(), dialog.monge(),
                   dialog.neighborhood_size (), dialog.samples(),
                   dialog.scalespace_af(),
                   dialog.generate_smoothed (),
                   dialog.longest_edge_2(), dialog.radius_ratio_bound_2(),
                   CGAL_PI * dialog.beta_angle_2 () / 180.,
                   dialog.separate_shells (), dialog.force_manifold ());

      for (std::size_t i = 0; i < reco_items.size (); ++ i)
        {
          if (!(dialog.scalespace_af()))
          {
            if (dialog.force_manifold () && i > reco_items.size () - 3)
            {
              if (dialog.generate_smoothed () && i % 2)
                reco_items[i]->setName(tr("%1 (scale space smooth garbage)").arg(scene->item(index)->name()));
              else
                reco_items[i]->setName(tr("%1 (scale space garbage)").arg(scene->item(index)->name()));
            }
            else
            {
              if (dialog.generate_smoothed ())
                {
                  if (i % 2)
                    reco_items[i]->setName(tr("%1 (scale space smooth shell %2)").arg(scene->item(index)->name()).arg((i+1)/2));
                  else
                    reco_items[i]->setName(tr("%1 (scale space shell %2)").arg(scene->item(index)->name()).arg(((i+1)/2)+1));
                }
              else
                reco_items[i]->setName(tr("%1 (scale space shell %2)").arg(scene->item(index)->name()).arg(i+1));
            }
          }
          else
          {
            if (dialog.generate_smoothed ())
            {
              if (i % 2)
                reco_items[i]->setName(tr("%1 (scale space smooth)").arg(scene->item(index)->name()));
              else
                reco_items[i]->setName(tr("%1 (scale space)").arg(scene->item(index)->name()));
            }
            else
              reco_items[i]->setName(tr("%1 (scale space)").arg(scene->item(index)->name()));
          }
          scene->addItem (reco_items[i]);
        }

      QApplication::restoreOverrideCursor();
    }
}


void Polyhedron_demo_surface_reconstruction_plugin::poisson_reconstruction
(const Polyhedron_demo_surface_reconstruction_plugin_dialog& dialog)
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(point_set_item)
    {
      // Gets point set
      Point_set* points = point_set_item->point_set();
      if(!points) return;

      const double sm_angle     = dialog.angle ();
      const double sm_radius    = dialog.radius ();
      const double sm_distance  = dialog.distance ();
      bool conjugate_gradient = dialog.conjugate_gradient();
      bool use_two_passes = dialog.two_passes();
      bool do_not_fill_holes = dialog.do_not_fill_holes();

      QApplication::setOverrideCursor(Qt::WaitCursor);

      // Reconstruct point set as a polyhedron
      SMesh* mesh = poisson_reconstruct (*points, sm_angle, sm_radius, sm_distance,
                                         conjugate_gradient, use_two_passes,
                                         do_not_fill_holes);
      if (mesh)
      {
        // Add polyhedron to scene
        Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(mesh);
        new_item->setName(tr("%1 Poisson").arg(point_set_item->name()));
        new_item->setColor(Qt::lightGray);
        scene->addItem(new_item);

        // Hide point set
        point_set_item->setVisible(false);
        scene->itemChanged(index);
      }

      QApplication::restoreOverrideCursor();
    }
}

#include "Surface_reconstruction_plugin.moc"
