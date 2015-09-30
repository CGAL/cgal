#include <QTime>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QObject>

#include "Scene_polygon_soup_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "ui_Polyhedron_demo_meta_reconstruction_plugin.h"

namespace MetaReconstruction
{
  
  unsigned int scale_of_anisotropy (const Point_set& points)
  {
    unsigned int scale = 6;

    for (; scale < 400; scale = static_cast<unsigned int>(scale * 1.5))
      {
	

      }

    
    return 6;
  }

  
  unsigned int scale_of_noise (const Point_set& points, unsigned int scale_min = 6)
  {

    return 6;
  }

  void simplify_point_set (const Point_set& points, unsigned int scale)
  {

  }

  void smooth_point_set (const Point_set& points, unsigned int scale)
  {

  }
}

class Polyhedron_demo_meta_reconstruction_plugin_dialog : public QDialog, private Ui::MetaReconstructionOptionsDialog
{
  Q_OBJECT
public:
  Polyhedron_demo_meta_reconstruction_plugin_dialog(QWidget* /*parent*/ = 0)
  {
    setupUi(this);
  }

  bool boundaries() const { return m_boundaries->isChecked(); }
  bool interpolate() const { return m_interpolate->isChecked(); }
};

#include <CGAL/Scale_space_surface_reconstruction_3.h>

class Polyhedron_demo_meta_reconstruction_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionScaleSpaceReconstruction;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionScaleSpaceReconstruction = new QAction(tr("Meta surface reconstruction"), mainWindow);
    actionScaleSpaceReconstruction->setObjectName("actionMetaReconstruction");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionScaleSpaceReconstruction;
  }

public Q_SLOTS:
  void on_actionMetaReconstruction_triggered();
}; // end class Polyhedron_meta_reconstruction_plugin


void Polyhedron_demo_meta_reconstruction_plugin::on_actionMetaReconstruction_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
    {
      // Gets point set
      Point_set* points = pts_item->point_set();

      //generate the dialog box to set the options
      Polyhedron_demo_meta_reconstruction_plugin_dialog dialog;
      if(!dialog.exec())
	return;

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      std::cout << "Meta surface reconstruction with the following requirements:" << std::endl
		<< (dialog.boundaries() ? " * Output shape has boundaries" : " * Output shape is closed") << std::endl
		<< (dialog.interpolate() ? " * Output shape passes through input points"
		    : " * Output shape approximates input points") << std::endl;


      std::cout << "Analysis of input point set:" << std::endl;
      time.start();

      // TODO: analyse if point set is isotropic
      unsigned int aniso_scale = MetaReconstruction::scale_of_anisotropy (*points);
      bool isotropic = (aniso_scale > 6);

      // TODO: analyse if point set is noisy
      unsigned int noise_scale = MetaReconstruction::scale_of_noise (*points, aniso_scale);
      bool noisy = (noise_scale > 6);
      
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      
      if (!(dialog.interpolate()) && (noisy || isotropic))
	{
	  points = new Point_set();
	  std::copy (pts_item->point_set()->begin(), pts_item->point_set()->end(),
		     std::back_inserter (*points));
	  
	  std::cout << "Preprocessing:" << std::endl;
	  time.restart();

	  if (isotropic)
	    {
	      // TODO: simplify point set
	      MetaReconstruction::simplify_point_set (*points, aniso_scale);
	    }
	  if (noisy)
	    {
	      // TODO: smooth point set
	      MetaReconstruction::smooth_point_set (*points, noise_scale);
	    }

	  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	}
      

      if (dialog.interpolate())
	{
	  if (noisy)
	    { 
	      std::cout << "Scale space reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      std::cout << "Advancing front reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }

	}
      else
	{
	  if (dialog.boundaries())
	    {
	      std::cout << "Scale space reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	  else
	    {
	      std::cout << "Poisson reconstruction:" << std::endl;
	      time.restart();
	      // TODO

	      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
	    }
	}

      if (!(dialog.interpolate()) && (noisy || isotropic))
	delete points;

      // default cursor
      QApplication::restoreOverrideCursor();
    }
}

#include "Polyhedron_demo_meta_reconstruction_plugin.moc"
