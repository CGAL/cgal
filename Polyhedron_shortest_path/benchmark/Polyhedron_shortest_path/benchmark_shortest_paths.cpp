// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Random.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>

#define UNUSED(X) (void)sizeof(X)

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
typedef Traits::Barycentric_coordinate Barycentric_coordinate;
typedef Traits::Construct_barycentric_coordinate Construct_barycentric_coordinate;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::halfedge_iterator halfedge_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;


typedef boost::timer::nanosecond_type nanosecond_type;

template <class IntType>
class SampledValue
{
private:
  IntType m_minimum;
  IntType m_maximum;
  IntType m_sum;
  IntType m_numSamples;
  
public:
  SampledValue()
  {
    reset();
  }
  
  void reset()
  {
    m_minimum = IntType(0);
    m_maximum = IntType(0);
    m_sum = IntType(0);
    m_numSamples = IntType(0);
  }

  IntType minimum() const
  {
    return m_minimum;
  }
  
  IntType maximum() const
  {
    return m_maximum;
  }
  
  IntType sum() const
  {
    return m_sum;
  }
  
  IntType average_integer() const
  {
    return m_sum / m_numSamples;
  }
  
  IntType average_remainder() const
  {
    return m_sum % m_numSamples;
  }
  
  double average_float() const
  {
    return double(average_integer()) + (double(average_remainder()) / double(m_numSamples));
  }
  
  void add_sample(IntType sample)
  {
    if (m_numSamples == 0)
    {
      m_minimum = sample;
      m_maximum = sample;
    }
    else
    {
      m_minimum = std::min(m_minimum, sample);
      m_maximum = std::max(m_maximum, sample);
    }
    m_sum += sample;
    ++m_numSamples;
  }
};

struct Benchmark_data
{
  size_t numVertices;
  size_t numFaces;
  size_t numEdges;

  SampledValue<boost::timer::nanosecond_type> constructionTime;
  SampledValue<boost::timer::nanosecond_type> queryTime;
  
#if !defined(NDEBUG)
  SampledValue<size_t> peakMemoryUsage;
#endif

  void reset()
  {
    numVertices = 0;
    numEdges = 0;
    numFaces = 0;
    constructionTime.reset();
    queryTime.reset();
#if !defined(NDEBUG)
    peakMemoryUsage.reset();
#endif
  }
};

void print_results(std::ostream& stream, const std::string& filename, const Benchmark_data& outData)
{
  stream << "Filename | " << filename << std::endl;
  stream << "Num Vertices | " << outData.numVertices << std::endl;
  stream << "Num Edges | " << outData.numEdges << std::endl;
  stream << "Num Faces | " << outData.numFaces << std::endl;
  
  stream << "Construction  | " << outData.constructionTime.average_float() / 1.0e9 << " | " << double(outData.constructionTime.minimum()) / 1.0e9 << " | " << double(outData.constructionTime.maximum()) / 1.0e9 << " |" << std::endl;
  stream << "Query         | " << outData.queryTime.average_float() / 1.0e9 << " | " << double(outData.queryTime.minimum()) / 1.0e9 << " | " << double(outData.queryTime.maximum()) / 1.0e9 << " |" << std::endl;
#if !defined(NDEBUG)
  stream << "Memory (Peak) | " << outData.peakMemoryUsage.average_float() / 1.0e6 << " | " << double(outData.peakMemoryUsage.minimum()) / 1.0e6 << " | " << double(outData.peakMemoryUsage.maximum()) / 1.0e6 << " |" << std::endl;
#endif
}

Barycentric_coordinate random_coordinate(CGAL::Random& rand)
{
  Construct_barycentric_coordinate construct_barycentric_coordinate;
  Traits::FT u = rand.uniform_01<Traits::FT>();
  Traits::FT v = rand.uniform_real(Traits::FT(0.0), Traits::FT(1.0) - u);
  return construct_barycentric_coordinate(u, v, Traits::FT(1.0) - u - v);
}

void run_benchmarks_no_id(CGAL::Random& rand, size_t numTrials, size_t numSources, size_t numQueries, Polyhedron_3& polyhedron, Benchmark_data& outData)
{
  boost::timer::cpu_timer timer;
  outData.reset();
  
  outData.numVertices = boost::num_vertices(polyhedron);
  outData.numEdges = boost::num_edges(polyhedron);
  outData.numFaces = CGAL::num_faces(polyhedron);

  face_iterator startFace, endFace;
  boost::tie(startFace, endFace) = CGAL::faces(polyhedron);
  
  std::vector<face_descriptor> allFaces;
  
  for (face_iterator currentFace = startFace; currentFace != endFace; ++currentFace)
  {
    allFaces.push_back(*currentFace);
  }
  
  Traits traits;
  Polyhedron_shortest_path shortestPaths(polyhedron, traits);
  
  for (size_t i = 0; i < numTrials; ++i)
  {
    std::vector<Polyhedron_shortest_path::Face_location> sourcePoints;
    
    while (sourcePoints.size() < numSources)
    {
      face_descriptor sourceFace = allFaces[rand.get_int(0, allFaces.size())];
      Barycentric_coordinate sourceLocation = random_coordinate(rand);
      sourcePoints.push_back(Polyhedron_shortest_path::Face_location(sourceFace, sourceLocation));
    }

    timer.start();
    shortestPaths.construct_sequence_tree(sourcePoints.begin(), sourcePoints.end());
    timer.stop();
    
    boost::timer::cpu_times elapsed = timer.elapsed();
    
    outData.constructionTime.add_sample(elapsed.wall);
    
#if !defined(NDEBUG)
    outData.peakMemoryUsage.add_sample(std::max(shortestPaths.peak_memory_usage(), shortestPaths.current_memory_usage()));
#endif
    
    for (size_t j = 0; j < numQueries; ++j)
    {
      face_descriptor sourceFace = allFaces[rand.get_int(0, allFaces.size())];
      Barycentric_coordinate sourceLocation = random_coordinate(rand);
      
      timer.start();
      Traits::FT distance = shortestPaths.shortest_distance_to_source_points(sourceFace, sourceLocation).first;
      timer.stop();
      UNUSED(distance);
      
      boost::timer::cpu_times elapsed = timer.elapsed();
      
      outData.queryTime.add_sample(elapsed.wall);
    }
  }
}

int main(int argc, char* argv[])
{
  namespace po = boost::program_options;
  po::options_description options;

  options.add_options()
    ("help,h", "Display help message")
    ("polyhedron,p", po::value<std::string>(), "Polyhedron input file")
    ("randomseed,r", po::value<size_t>()->default_value(0), "Randomization seed value")
    ("trials,t", po::value<size_t>()->default_value(10), "Number of trials to run")
    ("numpoints,n", po::value<size_t>()->default_value(10), "Number of source points per trial")
    ("queries,q", po::value<size_t>()->default_value(10), "Number of queries to run per trial")
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
    Polyhedron_3 polyhedron;
    std::ifstream inFile(vm["polyhedron"].as<std::string>().c_str());
    inFile >> polyhedron;
    inFile.close();
    
    Benchmark_data results;
    
    CGAL::Random rand(vm["randomseed"].as<size_t>());
    
    run_benchmarks_no_id(
      rand,
      vm["trials"].as<size_t>(),
      vm["numpoints"].as<size_t>(),
      vm["queries"].as<size_t>(),
      polyhedron,
      results);
      
    print_results(std::cout, vm["polyhedron"].as<std::string>(), results);
  }
  else
  {
    std::cerr << "Please specify a polyhedron to use.  Use option --help for more details." << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}