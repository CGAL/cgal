// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>
#include <CGAL/Polyhedron_shortest_path/Internal/misc_functions.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#define UNUSED(X) (void)sizeof(X)

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
typedef Traits::Point_3 Point_3;
typedef Traits::Vector_3 Vector_3;
typedef Traits::FT FT;
typedef Traits::Barycentric_coordinate Barycentric_coordinate;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::halfedge_iterator halfedge_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

typedef boost::property_map<typename Traits::Polyhedron, CGAL::vertex_external_index_t>::type VertexExternalIndexMap;
typedef boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type HalfedgeExternalIndexMap;
typedef boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type FaceExternalIndexMap;
typedef boost::property_map<typename Traits::Polyhedron, CGAL::vertex_point_t>::type VertexPointMap;

enum Sequence_item_type
{
  SEQUENCE_ITEM_VERTEX,
  SEQUENCE_ITEM_EDGE,
  SEQUENCE_ITEM_FACE,
};

template <class Traits>
struct Sequence_item
{
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef typename Traits::FT FT;
  
  Sequence_item_type type;
  size_t index;
  Barycentric_coordinate faceAlpha;
  FT edgeAlpha;
  vertex_descriptor vertex;
  halfedge_descriptor halfedge;
  face_descriptor face;
};

Point_3 get_item_location(const Sequence_item<Traits>& item, Polyhedron_shortest_path& helper)
{
  if (item.type == SEQUENCE_ITEM_VERTEX)
  {
    return helper.get_vertex_location(item.vertex);
  }
  else if (item.type == SEQUENCE_ITEM_EDGE)
  {
    return helper.get_edge_location(item.halfedge, item.edgeAlpha);
  }
  else
  {
    return helper.get_face_location(item.face, item.faceAlpha);
  }
}

template <class Traits, 
  class VIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_external_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type,
  class VPM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_point_t>::type>
struct Edge_sequence_collector
{
  typedef typename Traits::Polyhedron Polyhedron;
  typedef typename Traits::FT FT;
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef VIM VertexIndexMap;
  typedef HIM HalfedgeIndexMap;
  typedef FIM FaceIndexMap;
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;

  VertexIndexMap m_vertexIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  FaceIndexMap m_faceIndexMap;
  
  std::vector<Sequence_item<Traits> > m_sequence;
  
  Edge_sequence_collector(Polyhedron& p)
    : m_vertexIndexMap(CGAL::get(boost::vertex_external_index, p))
    , m_halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, p))
    , m_faceIndexMap(CGAL::get(CGAL::face_external_index, p))
  {
  }

  Edge_sequence_collector(VertexIndexMap& vertexIndexMap, HalfedgeIndexMap& halfedgeIndexMap, FaceIndexMap& faceIndexMap)
    : m_vertexIndexMap(vertexIndexMap)
    , m_halfedgeIndexMap(halfedgeIndexMap)
    , m_faceIndexMap(faceIndexMap)
  {
  }
  
  void edge(halfedge_descriptor he, FT alpha)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_EDGE;
    item.index = m_halfedgeIndexMap[he];
    item.edgeAlpha = alpha;
    item.halfedge = he;
    m_sequence.push_back(item);
  }
  
  void vertex(vertex_descriptor v)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_VERTEX;
    item.index = m_vertexIndexMap[v];
    item.vertex = v;
    m_sequence.push_back(item);
  }
  
  void face(face_descriptor f, Barycentric_coordinate alpha)
  {
    Sequence_item<Traits> item;
    item.type = SEQUENCE_ITEM_FACE;
    item.index = m_faceIndexMap[f];
    item.faceAlpha = alpha;
    item.face = f;
    m_sequence.push_back(item);
  }
};

class Edge_sequence_comparator
{
private:
  std::vector<halfedge_descriptor> m_edges;
  Polyhedron_3& m_polyhedron;

  const Sequence_item<Traits>& second_last(Edge_sequence_collector<Traits>* list) const
  {
    return list->m_sequence[list->m_sequence.size() - 2];
  }
  
  size_t edge_index(halfedge_descriptor e) const
  {
    for (size_t i = 0; i < m_edges.size(); ++i)
    {
      if (m_edges[i] == CGAL::opposite(e, m_polyhedron))
      {
        return i;
      }
    }
    
    assert(false && "Error, edge was not found");
  }
  
  size_t vertex_index(vertex_descriptor v) const
  {
    for (size_t i = 0; i < m_edges.size(); ++i)
    {
      if (CGAL::source(m_edges[i], m_polyhedron) == v)
      {
        return i;
      }
    }
    
    assert(false && "Error, vertex was not found");
  }
  
  size_t get_index(const Sequence_item<Traits>& item) const
  {
    if (item.type == SEQUENCE_ITEM_VERTEX)
    {
      return vertex_index(item.vertex);
    }
    else if (item.type == SEQUENCE_ITEM_EDGE)
    {
      return edge_index(item.halfedge);
    }
    else
    {
      assert(false && "Error, face item found");
    }
  }
  
public:
  Edge_sequence_comparator(Polyhedron_3& polyhedron, face_descriptor locationFace)
    : m_polyhedron(polyhedron)
  {
    m_edges.push_back(CGAL::halfedge(locationFace, m_polyhedron));
    m_edges.push_back(CGAL::next(m_edges[0], m_polyhedron));
    m_edges.push_back(CGAL::next(m_edges[1], m_polyhedron));
  }
  
  Edge_sequence_comparator(Polyhedron_3& polyhedron, vertex_descriptor locationVertex)
    : m_polyhedron(polyhedron)
  {
    halfedge_descriptor startEdge = CGAL::halfedge(locationVertex, m_polyhedron);
    
    halfedge_descriptor currentEdge = startEdge;
    
    do
    {
      currentEdge = CGAL::opposite(currentEdge, m_polyhedron);
      m_edges.push_back(CGAL::next(currentEdge, m_polyhedron));
      currentEdge = CGAL::next(m_edges.back(), m_polyhedron);
    }
    while(currentEdge != startEdge);
  }
  
  Edge_sequence_comparator(Polyhedron_3& polyhedron, halfedge_descriptor locationEdge)
    : m_polyhedron(polyhedron)
  {
    m_edges.push_back(CGAL::next(locationEdge, m_polyhedron));
    m_edges.push_back(CGAL::next(m_edges[0], m_polyhedron));
    m_edges.push_back(CGAL::next(CGAL::opposite(locationEdge, m_polyhedron), m_polyhedron));
    m_edges.push_back(CGAL::next(m_edges[2], m_polyhedron));
  }

  bool operator() (Edge_sequence_collector<Traits>* lhs, Edge_sequence_collector<Traits>* rhs) const
  {
    const Sequence_item<Traits>& left = second_last(lhs);
    const Sequence_item<Traits>& right = second_last(rhs);
    
    size_t leftIndex = get_index(left);
    size_t rightIndex = get_index(right);
    
    if (leftIndex < rightIndex)
    {
      return true;
    }
    else if (rightIndex < leftIndex)
    {
      return false;
    }
    else if (left.type == SEQUENCE_ITEM_VERTEX)
    {
      return true;
    }
    else if (right.type == SEQUENCE_ITEM_VERTEX)
    {
      return false;
    }
    else
    {
      return left.edgeAlpha < right.edgeAlpha;
    }
  }
};

void run_isolines_no_id(Traits& traits, Polyhedron_3& polyhedron, Polyhedron_shortest_path& shortestPaths, Polyhedron_shortest_path::Face_location_pair faceLocation, const Edge_sequence_comparator& edgeComparator, FT delta, const std::string& outPathsFilename, const std::string& outIsolinesFilename)
{
  Traits::Compute_squared_distance_3 compute_squared_distance_3(traits.compute_squared_distance_3_object());
  shortestPaths.compute_shortest_paths(faceLocation.first, faceLocation.second);

  vertex_iterator verticesCurrent, verticesEnd;
  
  std::ofstream outPaths(outPathsFilename.c_str());
  
  VertexExternalIndexMap vertexIndexMap(CGAL::get(boost::vertex_external_index, polyhedron));
  HalfedgeExternalIndexMap halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, polyhedron));
  FaceExternalIndexMap faceIndexMap(CGAL::get(CGAL::face_external_index, polyhedron));

  std::vector<Edge_sequence_collector<Traits>*> edgeSequences;
  
  FT maxDistance(0.0);
  
  for (boost::tie(verticesCurrent, verticesEnd) = boost::vertices(polyhedron); verticesCurrent != verticesEnd; ++verticesCurrent)
  {
    Edge_sequence_collector<Traits>* collector = new Edge_sequence_collector<Traits>(vertexIndexMap, halfedgeIndexMap, faceIndexMap);
    
    shortestPaths.shortest_path_sequence(*verticesCurrent, *collector);
    
    maxDistance = std::max(maxDistance, shortestPaths.shortest_distance_to_vertex(*verticesCurrent));
    
    if (collector->m_sequence.size() >= 2)
    {
      edgeSequences.push_back(collector);
    }
    
    outPaths << collector->m_sequence.size() + 1;
    
    outPaths << " " << shortestPaths.get_vertex_location(*verticesCurrent);
    
    for (size_t i = 0; i < collector->m_sequence.size(); ++i)
    {
      outPaths << " " << get_item_location(collector->m_sequence[i], shortestPaths);
    }
    
    outPaths << std::endl;
  }
  
  outPaths.close();

  std::sort(edgeSequences.begin(), edgeSequences.end(), edgeComparator);

    size_t numPaths = edgeSequences.size();
  
  std::vector<Point_3> currentLocations(numPaths);
  std::vector<size_t> currentItems(numPaths);
  std::vector<bool> currentlyCompleted(numPaths);
  size_t numPathsCompleted = 0;
  
  Point_3 startPoint = shortestPaths.get_face_location(faceLocation.first, faceLocation.second);
  
  for (size_t i = 0; i < edgeSequences.size(); ++i)
  {
    currentLocations[i] = startPoint;
    currentItems[i] = edgeSequences[i]->m_sequence.size() - 2;
    currentlyCompleted[i] = false;
  }
  
  FT currentSampleDistance = FT(0.0);
  
  std::ofstream outIsoLines(outIsolinesFilename.c_str());
  
  while (numPathsCompleted < numPaths - 1)
  {
    std::vector<Point_3> points;
    
    currentSampleDistance += delta;
    
    // advance each path to the next sample point
    for (size_t i = 0; i < numPaths; ++i)
    {
      FT remainingDistance = delta;
      
      while (!currentlyCompleted[i] && remainingDistance > FT(0.0))
      {
        Point_3 nextLocation = get_item_location(edgeSequences[i]->m_sequence[currentItems[i]], shortestPaths);
        FT distance = CGAL::sqrt(compute_squared_distance_3(nextLocation, currentLocations[i]));
        
        //std::cout << distance << std::endl;
        
        if (distance < remainingDistance)
        {
          remainingDistance -= distance;

          if (currentItems[i] == 0)
          {
            currentlyCompleted[i] = true;
            ++numPathsCompleted;
          }
          else
          {
            --currentItems[i];
            currentLocations[i] = nextLocation;
          }
        }
        else
        {
          Vector_3 direction = nextLocation - currentLocations[i];
          Vector_3 scaled = direction / distance;
          
          currentLocations[i] = currentLocations[i] + scaled * remainingDistance;
          remainingDistance = FT(0.0);
        }
      }
    }
    
    //std::cout << "Current Distance: " << currentSampleDistance << std::endl;
    //std::cout << "Num Paths Cleared: " << numPathsCompleted << std::endl;
    
    if (numPathsCompleted >= numPaths - 1)
    {
      break;
    }
    
    size_t firstPath = 0;
    
    // connect pairs of adjacent paths
    while (currentlyCompleted[firstPath])
    {
      ++firstPath;
    }
    
    size_t currentPath = firstPath;
    
    while (currentPath < numPaths)
    {
      size_t nextPath = currentPath + 1;
      
      while (currentlyCompleted[nextPath % numPaths])
      {
        ++nextPath;
      }
      
      points.push_back(currentLocations[currentPath]);
      
      currentPath = nextPath;
    }
    
    outIsoLines << points.size();
    
    for (size_t i = 0; i < points.size(); ++i)
    {
      outIsoLines << " " << points[i];
    }
    
    outIsoLines << " " << points[0] << std::endl;
  }
  
  outIsoLines.close();

  for (size_t i = 0; i < edgeSequences.size(); ++i)
  {
    delete edgeSequences[i];
  }
}


int main(int argc, char** argv)
{
  namespace po = boost::program_options ;
  po::options_description options;

  options.add_options()
    ("help,h", "Display help message")
    ("polyhedron,p", po::value<std::string>(), "Polyhedron input file")
    ("vertex,v", po::value<size_t>(), "Randomization seed value")
    ("face,f", po::value<size_t>(), "Number of trials to run")
    ("location,l", po::value<std::string>(), "Number of source points per trial")
    ("delta,d", po::value<double>()->default_value(0.1), "Number of queries to run per trial")
    ("output-paths,o", po::value<std::string>()->default_value("shortestpaths.cgal"), "Filename to write shortest paths output")
    ("output-isolines,i", po::value<std::string>()->default_value("isolines.cgal"), "Filename to write isolines")
    ;
    
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, options), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << options << std::endl;
  }
  else if (vm.count("polyhedron"))
  {
    bool useVertex = false;
    bool useFace = false;
    
    if (vm.count("vertex"))
    {
      useVertex = true;
    }
    else if (vm.count("face"))
    {
      if (vm.count("location"))
      {
        useFace = true;
      }
      else
      {
        std::cerr << "Must specify face location with face." << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Must specify either a face or a vertex." << std::endl;
      return EXIT_FAILURE;
    }
  
    if (useVertex || useFace)
    {
       FT delta(vm["delta"].as<double>());
      
      if (delta < FT(0.0))
      {
        std::cerr << "Must specify non-negative isoline delta." << std::endl;
        return EXIT_FAILURE;
      }
    
      Polyhedron_shortest_path::Face_location_pair faceLocation;
      vertex_descriptor targetVertex;
    
      Polyhedron_3 polyhedron;
      std::ifstream inFile(vm["polyhedron"].as<std::string>().c_str());
      inFile >> polyhedron;
      inFile.close();
      
      Traits traits;
      Polyhedron_shortest_path shortestPaths(polyhedron, traits);

      if (useVertex)
      {
        size_t vertexIndex = vm["vertex"].as<size_t>();
        if (vertexIndex >= boost::num_vertices(polyhedron))
        {
          std::cerr << "Vertex index out of range (#Vertices = " << boost::num_vertices(polyhedron) << ")." << std::endl;
          return EXIT_FAILURE;
        }

        vertex_iterator currentVertex, endVertex;
        boost::tie(currentVertex, endVertex) = boost::vertices(polyhedron);
        
        while (vertexIndex > 0)
        {
          ++currentVertex;
          --vertexIndex;
        }
        
        faceLocation = shortestPaths.get_vertex_as_face_location(*currentVertex);
        targetVertex = *currentVertex;
      }
      else
      {
        size_t faceIndex = vm["face"].as<size_t>();
        if (faceIndex >= CGAL::num_faces(polyhedron))
        {
          std::cerr << "Face index out of range (#Face = " << CGAL::num_faces(polyhedron) << ")." << std::endl;
          return EXIT_FAILURE;
        }
        
        face_iterator currentFace, endFace;
        boost::tie(currentFace, endFace) = CGAL::faces(polyhedron);
        while (faceIndex > 0)
        {
          ++currentFace;
          --faceIndex;
        }
        
        std::istringstream istr(vm["location"].as<std::string>().c_str());
        
        Barycentric_coordinate location;
        
        istr >> location;
        
        FT coordinateSum = location[0] + location[1] + location[2];
        
        if (CGAL::abs(coordinateSum - FT(1.0)) > FT(0.00001))
        {
          std::cerr << "Invalid face location (must be a valid, normalized, barycentric coordinate)." << std::endl;
          return EXIT_FAILURE;
        }
        
        faceLocation = Polyhedron_shortest_path::Face_location_pair(*currentFace, location);
      }
      
      run_isolines_no_id(traits, polyhedron, shortestPaths, faceLocation, useVertex ? Edge_sequence_comparator(polyhedron, targetVertex) : Edge_sequence_comparator(polyhedron, faceLocation.first), delta, vm["output-paths"].as<std::string>(), vm["output-isolines"].as<std::string>());
    }
  }
  else
  {
    std::cout << "Please specify a polyhedron to use.  Use option --help for more details." << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
