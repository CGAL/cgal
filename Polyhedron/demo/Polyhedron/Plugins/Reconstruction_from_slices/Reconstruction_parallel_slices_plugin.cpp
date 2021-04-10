#include <QTime>
#include <QApplication>
#include <QAction>
#include <QStringList>
#include <QMainWindow>
#include <QtPlugin>
#include <QMessageBox>


#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"
#include "Polyhedron_type.h"

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Reconstruction_from_parallel_slices_3.h>
#include <CGAL/Reconstruction_from_parallel_slices_3/contour_providers.h>


using namespace CGAL::Three;
class Polyhedron_demo_reconstruction_parallel_slices_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  QList<QAction*> _actions;
  Scene_interface* scene;
  QMainWindow* mw;
  int detect_constant_coordinate(Scene_polylines_item* polylines_item);
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    QAction* actionReconstructionParallelSlices = new QAction(tr("Reconstruction from parallel slices"), mainWindow);
    actionReconstructionParallelSlices->setObjectName("actionReconstructionParallelSlices");
    connect(actionReconstructionParallelSlices, SIGNAL(triggered()), this, SLOT(on_actionReconstructionParallelSlices_triggered()));
    _actions << actionReconstructionParallelSlices;
    scene = scene_interface;
    mw=mainWindow;
  }

  QList<QAction*> actions()const {return _actions;}

  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionReconstructionParallelSlices";
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionReconstructionParallelSlices_triggered();

}; // end Polyhedron_demo_reconstruction_parallel_slices_plugin


int Polyhedron_demo_reconstruction_parallel_slices_plugin::detect_constant_coordinate(Scene_polylines_item* polylines_item)
{
  typedef Scene_polylines_item::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline;
  typedef std::list<Polyline> Polylines_container;
  Polylines_container& polylines=polylines_item->polylines;
  
  Polylines_container::size_type nb_polylines=polylines.size();
  
  switch (nb_polylines)
  {
    case 0:
    case 1:
      QMessageBox::warning(mw, tr("Error"),
                     tr("The selected polylines item must contains at least two contours."));
      return -1;
    default:
    {
      Polylines_container::iterator it_poly=polylines.begin();
      std::vector<int> const_coords;
      const_coords.reserve(3);
      const_coords.push_back(0);
      const_coords.push_back(1);
      const_coords.push_back(2);
      
      
      while(it_poly!=polylines.end())
      {
        Polyline& polyline = *it_poly++;
        if ( polyline.size()==1 ) continue;
        
        Point_3 pt=polyline[0];
        
        for (Polyline::iterator it_pt=CGAL::cpp11::next( polyline.begin() );it_pt!=polyline.end();++it_pt)
        {
          std::vector<int> to_keep;
          std::vector<int>::size_type nc=const_coords.size();
          for (std::vector<int>::size_type i=0;i<nc;++i)
            if ( pt[ const_coords[i] ] == (*it_pt)[ const_coords[i] ] ) to_keep.push_back( const_coords[i] );
          switch( to_keep.size() )
          {
            case 0:
              QMessageBox::warning(mw, tr("Error"),
                             tr("The selected polylines item are not contained in parallel slices."));
              return -1;
            case 1:
              return to_keep[0];
            default:
              const_coords=to_keep;
          }
        }
      }
      QMessageBox::warning(mw, tr("Error"),
                     tr("The selected polylines item are not contained in parallel slices."));
      return -1;
    }
  }
}

struct Visitor_update{
  Scene_interface* m_scene;
  int m_item_index;
  
  Visitor_update(Scene_interface* scene,int item_index):m_scene(scene),m_item_index(item_index){}
  void one_layer_is_finished() { Q_EMIT m_scene->itemChanged(m_item_index); }
};

void Polyhedron_demo_reconstruction_parallel_slices_plugin::on_actionReconstructionParallelSlices_triggered()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polygon_as_vector_of_Point_3_in_axis_aligned_planes<K,Scene_polylines_item::Polylines_container> Contour_reader;
  typedef CGAL::Incremental_slice_writer_into_polyhedron<Polyhedron,K,Visitor_update> Slice_writer;
  
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polylines_item* polylines_item = 
    qobject_cast<Scene_polylines_item*>(scene->item(index));
  
  if(polylines_item)
  {
    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::cout << "Analysing input...";
    int cst_coord = detect_constant_coordinate(polylines_item);
    if ( cst_coord == -1 ) return;
    
    switch (cst_coord)
    {
      case 0:
        std::cout << "  Constant coordinate is X"<< std::endl;
        break;
      case 1:
        std::cout << "  Constant coordinate is Y"<< std::endl;
        break;
      case 2:
        std::cout << "  Constant coordinate is Z"<< std::endl;
        break;
      default:
        std::cout << "  Invalid input."<< std::endl;
        break;
        return;
    }
    
    
    QTime time;
    time.start();
    std::cout << "Reconstructing ...";

    // add a new polyhedron
    Polyhedron *pReconst = new Polyhedron;

    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pReconst);
    new_item->setName(tr("%1 (reconstruction)").arg(polylines_item->name()));
    int item_index=scene->addItem(new_item);    
    
    Visitor_update vis_update(scene,item_index);
    //slice writer writing in a polyhedron
    Slice_writer writer(*pReconst,vis_update);

    Contour_reader reader(polylines_item->polylines, cst_coord);
    CGAL::Reconstruction_from_parallel_slices_3<Slice_writer> reconstruction;
    reconstruction.run(reader,writer,cst_coord);

    new_item->invalidateOpenGLBuffers();
    scene->itemChanged(item_index);

    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;



    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

#include "Reconstruction_parallel_slices_plugin.moc"
