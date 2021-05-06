#include "config.h"
#include "Kernel_type.h"
#include "Scene_polylines_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/Polygon_mesh_processing/polyline_snapping.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

#include <boost/function_output_iterator.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;
using CGAL::cpp11::get;

struct Get_snapping_points_map
{
  typedef CGAL::cpp11::tuple<std::size_t, std::size_t,
                             std::size_t, std::size_t,
                             Kernel::FT, Kernel::FT> OutputType;

  typedef std::map<std::pair<std::size_t, std::size_t>, Kernel::Point_3> Map_pts;

  std::vector<std::vector<Kernel::Point_3> >& polylines;
  boost::shared_ptr<Map_pts> map_pts;
    
  Get_snapping_points_map (std::vector<std::vector<Kernel::Point_3> >& polylines)
    : polylines (polylines), map_pts (new Map_pts()) { }

  void operator() (const OutputType& s)
  {
    Kernel::Point_3 pa
      = CGAL::barycenter (polylines[get<0>(s)][get<1>(s)],     (1. - get<4>(s)),
                          polylines[get<0>(s)][get<1>(s) + 1],       get<4>(s));
    Kernel::Point_3 pb
      = CGAL::barycenter (polylines[get<2>(s)][get<3>(s)],     (1. - get<5>(s)),
                          polylines[get<2>(s)][get<3>(s) + 1],       get<5>(s));

    Kernel::Point_3 p = CGAL::midpoint (pa, pb);

    map_pts->insert (std::make_pair (std::make_pair (get<0>(s), get<1>(s)), p));
    map_pts->insert (std::make_pair (std::make_pair (get<2>(s), get<3>(s)), p));
  }
};


using namespace CGAL::Three;
class Polyhedron_demo_snap_polylines_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionSnapPolylines;

public:
  
  void init (QMainWindow* mainWindow,
             CGAL::Three::Scene_interface* scene_interface,
             Messages_interface*)
  {
    scene = scene_interface;
    mw = mainWindow;
    actionSnapPolylines = new QAction(tr("Snap Polylines"), mainWindow);
    actionSnapPolylines->setObjectName("actionSnapPolylines");
    actionSnapPolylines->setProperty("subMenuName","Operations on Polylines");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionSnapPolylines;
  }

  bool applicable(QAction*) const {
    Q_FOREACH(int index, scene->selectionIndices())
      if (!qobject_cast<Scene_polylines_item*>(scene->item(index)))
        return false;
    return true;
  }

public Q_SLOTS:
  
  void on_actionSnapPolylines_triggered();

};

void Polyhedron_demo_snap_polylines_plugin::on_actionSnapPolylines_triggered()
{
  bool ok;

  const double tolerance =
    QInputDialog::getDouble((QWidget*)mw,
                            tr("Snap Polylines"), // dialog title
                            tr("Tolerance:"), // field label
                            0.005 * scene->len_diagonal(), // default value
                            0.000001, // min
                            100000., // max
                            6, // decimals
                            &ok);
  if(!ok) 
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);
  QApplication::processEvents();

  std::vector<std::vector<Kernel::Point_3> > polylines;
  std::ostringstream oss;
  oss << "Snapped ";
    
  Q_FOREACH(int index, scene->selectionIndices())
  {
    Scene_polylines_item* item =
      qobject_cast<Scene_polylines_item*>(scene->item(index));
    if(item)
    {
      if (polylines.empty())
        oss << "[";
      else
        oss << " ; ";
      
      oss << item->name().toStdString();
      std::copy (item->polylines.begin(), item->polylines.end(),
                 std::back_inserter (polylines));
    }
  }
  oss << "]";

  std::cerr << polylines.size() << " input polyline(s)" << std::endl;
  if (!polylines.empty())
  {
    Get_snapping_points_map snapping_pts(polylines);

    PMP::polyline_snapping (polylines, tolerance,
                            boost::make_function_output_iterator
                            (snapping_pts),
                            Kernel());

    std::cerr << snapping_pts.map_pts->size() / 2
              << " output snapping point(s)" << std::endl;
      
    Scene_polylines_item* new_item
      = new Scene_polylines_item;

    for (std::size_t i = 0; i < polylines.size(); ++ i)
    {
      new_item->polylines.push_back (std::vector<Kernel::Point_3>());
      std::vector<Kernel::Point_3>& poly = new_item->polylines.back();

      for (std::size_t j = 0; j < polylines[i].size(); ++ j)
      {
        poly.push_back (polylines[i][j]);
        Get_snapping_points_map::Map_pts::iterator
          it = snapping_pts.map_pts->find (std::make_pair (i, j));
        if (it == snapping_pts.map_pts->end())
          continue;
        poly.push_back (it->second);
      }
    }
      
    new_item->setName (QString("%1").arg(oss.str().c_str()));
    scene->addItem(new_item);
  }

  QApplication::restoreOverrideCursor();
}



#include "Snap_polylines_plugin.moc"
