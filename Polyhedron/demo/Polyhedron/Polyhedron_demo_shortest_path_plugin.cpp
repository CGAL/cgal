#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_type.h"
#include "Scene.h"

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QAction>
#include <QDebug>
#include <QObject>
#include <QDockWidget>
//#include <QtConcurrentRun>
#include <map>
#include <algorithm>
#include <vector>

class Polyhedron_demo_shortest_path_plugin : 
    public QObject,
    public Polyhedron_demo_plugin_helper
{
    Q_OBJECT
        Q_INTERFACES(Polyhedron_demo_plugin_interface)
private:

    typedef boost::property_map<Polyhedron, boost::vertex_index_t>::type VertexIndexMap;
    typedef boost::property_map<Polyhedron, CGAL::halfedge_index_t>::type HalfedgeIndexMap;
    typedef boost::property_map<Polyhedron, CGAL::face_index_t>::type FaceIndexMap;
    typedef boost::property_map<Polyhedron, CGAL::vertex_point_t>::type VertexPointMap;
    
    typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron> Polyhedron_shortest_path_traits;
    typedef CGAL::Polyhedron_shortest_path<Polyhedron_shortest_path_traits /*, VertexIndexMap, HalfedgeIndexMap, FaceIndexMap, VertexPointMap*/> Polyhedron_shortest_path;
    
    struct ShortestPathsPointsVisitor
    {
      typedef std::vector<typename Polyhedron_shortest_path::Point_3> Container; 
      Container& m_container;

      ShortestPathsPointsVisitor(Container& container)
        : m_container(container)
      {
      }
      
      void point(const Polyhedron_shortest_path::Point_3& point)
      {
        std::cout << point << std::endl;
        m_container.push_back(point);
      }
    };
    
    
    
    typedef std::map<Scene_polyhedron_item*, Polyhedron_shortest_path* > Shortest_paths_map;
public:

    QList<QAction*> actions() const 
    {
        return QList<QAction*>() << actionComputeShortestPath;
    }

    bool applicable() const 
    {
      return 
        // Eventually this should be applicable to the polyhedron itself? I just want to get it working on the simplest possible case first
        qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    }
    
    void init(QMainWindow* mainWindow, Scene_interface* scene_interface) 
    {
        this->scene = scene_interface;
        this->mw = mainWindow;
        actionComputeShortestPath = new QAction("Compute Shortest Path", this->mw);
        
        connect(actionComputeShortestPath, SIGNAL(triggered()),this, SLOT(on_actionComputeShortestPath_triggered()));

        // This is for later
        /*if( Scene* scene = dynamic_cast<Scene*>(scene_interface) ) {
            connect(scene, SIGNAL(itemAboutToBeDestroyed(Scene_item*)), this, SLOT(itemAboutToBeDestroyed(Scene_item*)));
        }*/
    }
    
    void check_and_set_ids(Polyhedron* polyhedron);
    
    public slots:
        void on_actionComputeShortestPath_triggered();

private:
    QAction* actionComputeShortestPath;
    //QDockWidget* dock_widget;
    //Ui::Todo ui_widget;
};

void Polyhedron_demo_shortest_path_plugin::on_actionComputeShortestPath_triggered()
{ 
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  std::cout << "Selection Index: " << index << std::endl;
  
  Scene_polylines_item* lines_item = 
    qobject_cast<Scene_polylines_item*>(scene->item(index));
  
  Scene_polyhedron_selection_item* selection_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));
  
  std::cout << "Has selection_item : " << (selection_item ? "Yes" : "No") << std::endl;
  
  if (selection_item)
  {
    std::cout << "Has polyhedron : " << (selection_item->polyhedron() ? "Yes" : "No") << std::endl;
    
    if (selection_item->polyhedron())
    {
      check_and_set_ids(selection_item->polyhedron());
    
      if (selection_item->selected_vertices.begin() != selection_item->selected_vertices.end())
      {
        bool foundOrigin = false;
        Polyhedron::Vertex_iterator origin;
        bool foundDestination = false;
        Polyhedron::Vertex_iterator destination;
      
        std::cout << "Selected Vertices: " << std::endl;
        for(Scene_polyhedron_selection_item::Selection_set_vertex::iterator it = selection_item->selected_vertices.begin(); it != selection_item->selected_vertices.end(); ++it) 
        {
          if (!foundOrigin)
          {
            foundOrigin = true;
            origin = *it;
            std::cout << "Origin = " << (*it)->id() << std::endl;
          }
          else if (!foundDestination)
          {
            foundDestination = true;
            destination = *it;
            std::cout << "Destination = " << (*it)->id() << std::endl;
          }
          else
          {
            std::cout << "Leftovers: " << (*it)->id() << std::endl;
          }
        }
        std::cout << std::endl;
        
        if (foundOrigin && foundDestination)
        {
          typedef typename boost::graph_traits<Polyhedron> GraphTraits;

          Polyhedron_shortest_path shortestPath(*selection_item->polyhedron());//, 
            //CGAL::get(boost::vertex_index, *selection_item->polyhedron()), 
            //CGAL::get(CGAL::halfedge_index, *selection_item->polyhedron()),
            //CGAL::get(CGAL::face_index, *selection_item->polyhedron()),
            //CGAL::get(CGAL::vertex_point, *selection_item->polyhedron()));
          
          std::vector<typename Polyhedron_shortest_path::vertex_descriptor> vertices;
          
          typename Polyhedron_shortest_path::vertex_iterator begin, end;
          boost::tie(begin, end) = boost::vertices(*selection_item->polyhedron());
          
          for (typename Polyhedron_shortest_path::vertex_iterator curr = begin; curr != end; ++curr)
          {
            vertices.push_back(*curr);
          }
          
          std::cout << "Computing shortest paths between " << origin->id() << " and " << destination->id() << std::endl;
          
          shortestPath.compute_shortest_paths(vertices[origin->id()]);
          
          Scene_polylines_item* polylines = new Scene_polylines_item();
          
          polylines->polylines.push_back(Scene_polylines_item::Polyline());
          
          std::cout << "Origin: " << shortestPath.get_vertex_location(origin) << std::endl;
          std::cout << "Destination: " << shortestPath.get_vertex_location(destination) << std::endl;

          std::cout << "Past" << std::endl;
          
          ShortestPathsPointsVisitor visitor(polylines->polylines.back());
          
          shortestPath.shortest_path_points(vertices[destination->id()], visitor);
          
          for (Scene_polylines_item::Polyline::iterator it = polylines->polylines.back().begin(); it != polylines->polylines.back().end(); ++it)
          {
            std::cout << *it << std::endl;
          }
          
          scene->addItem(polylines);
          
          std::cout << "Done" << std::endl;
        }
      }
      else
      {
        std::cout << "No selected vertices." << std::endl;
      }
    }
  }
}


void Polyhedron_demo_shortest_path_plugin::check_and_set_ids(Polyhedron* polyhedron)
{
  Polyhedron::Vertex_iterator testVertex1 = polyhedron->vertices_begin();
  Polyhedron::Vertex_iterator testVertex2 = ++polyhedron->vertices_begin();
  
  if(testVertex1->id() == testVertex2->id())
  {
    std::size_t vertexId = 0;
    for(Polyhedron::Vertex_iterator currentVertex = polyhedron->vertices_begin();
        currentVertex != polyhedron->vertices_end(); ++currentVertex, ++vertexId)
    {
        currentVertex->id() = vertexId;
    }
  }
  
  Polyhedron::Halfedge_iterator testHalfedge1 = polyhedron->halfedges_begin();
  Polyhedron::Halfedge_iterator testHalfedge2 = ++polyhedron->halfedges_begin();
  
  if (testHalfedge1->id() == testHalfedge2->id())
  {
    std::size_t halfedgeId = 0;
    for(Polyhedron::Halfedge_iterator currentHalfedge = polyhedron->halfedges_begin();
        currentHalfedge != polyhedron->halfedges_end(); ++currentHalfedge, ++halfedgeId)
    {
        currentHalfedge->id() = halfedgeId;
    }
  }
  
  Polyhedron::Facet_iterator testFacet1 = polyhedron->facets_begin();
  Polyhedron::Facet_iterator testFacet2 = ++polyhedron->facets_begin();
  
  if (testFacet1->id() == testFacet2->id())
  {
    std::size_t facetId = 0;
    for(Polyhedron::Facet_iterator currentFacet = polyhedron->facets_begin();
        currentFacet != polyhedron->facets_end(); ++currentFacet, ++facetId)
    {
        currentFacet->id() = facetId;
    }
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_shortest_path_plugin, Polyhedron_demo_shortest_path_plugin)

#include "Polyhedron_demo_shortest_path_plugin.moc"
