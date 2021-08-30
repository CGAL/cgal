#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>

#include <CGAL/Random.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path_traits.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path.h>

#define UNUSED(X) (void)sizeof(X)

namespace po = boost::program_options;

typedef boost::timer::nanosecond_type nanosecond_type;

enum Kernel_type
{
  KERNEL_UNKNOWN,
  KERNEL_IPICK,
  KERNEL_EPICK,
  KERNEL_EPECK,
};

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
      m_minimum = (std::min)(m_minimum, sample);
      m_maximum = (std::max)(m_maximum, sample);
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
  stream << "Query         | " << 1.0 / (outData.queryTime.average_float() / 1.0e9) << " | " << double(outData.queryTime.minimum()) / 1.0e9 << " | " << double(outData.queryTime.maximum()) / 1.0e9 << " |" << std::endl;
#if !defined(NDEBUG)
  stream << "Memory (Peak) | " << outData.peakMemoryUsage.average_float() / 1.0e6 << " | " << double(outData.peakMemoryUsage.minimum()) / 1.0e6 << " | " << double(outData.peakMemoryUsage.maximum()) / 1.0e6 << " |" << std::endl;
#endif
}

template <class Traits>
typename Traits::Barycentric_coordinates random_coordinates(CGAL::Random& rand)
{
  typedef typename Traits::FT FT;
  typename Traits::Construct_barycentric_coordinates construct_barycentric_coordinates;
  FT u = rand.uniform_real(FT(0), FT(1));
  FT v = rand.uniform_real(FT(0), FT(1) - u);
  return construct_barycentric_coordinates(u, v, FT(1) - u - v);
}

template <class Kernel>
void run_benchmarks(CGAL::Random& rand, size_t numTrials, size_t numSources, size_t numQueries, const std::string& polyhedronFile, Benchmark_data& outData)
{
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
  typedef typename Traits::Barycentric_coordinates Barycentric_coordinates;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
  typedef typename Surface_mesh_shortest_path::Face_location Face_location;
  typedef typename Surface_mesh_shortest_path::FT FT;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::face_iterator face_iterator;

  std::ifstream inFile(polyhedronFile.c_str());

  if (!inFile)
  {
    throw std::runtime_error(std::string("Model file \"") + polyhedronFile + std::string("\" does not exist."));
  }

  Polyhedron_3 polyhedron;

  inFile >> polyhedron;
  inFile.close();

  CGAL::set_halfedgeds_items_id(polyhedron);

  boost::timer::cpu_timer timer;
  outData.reset();

  outData.numVertices = num_vertices(polyhedron);
  outData.numEdges = num_edges(polyhedron);
  outData.numFaces = num_faces(polyhedron);

  face_iterator startFace, endFace;
  boost::tie(startFace, endFace) = faces(polyhedron);

  std::vector<face_descriptor> allFaces;

  for (face_iterator currentFace = startFace; currentFace != endFace; ++currentFace)
  {
    allFaces.push_back(*currentFace);
  }

  Traits traits;
  Surface_mesh_shortest_path shortestPaths(polyhedron, traits);

  for (size_t i = 0; i < numTrials; ++i)
  {
    std::vector<Face_location> sourcePoints;

    while (sourcePoints.size() < numSources)
    {
      face_descriptor sourceFace = allFaces[rand.get_int(0, allFaces.size())];
      Barycentric_coordinates sourceLocation = random_coordinate<Traits>(rand);
      sourcePoints.push_back(Face_location(sourceFace, sourceLocation));
    }

    shortestPaths.clear();
    shortestPaths.add_source_points(sourcePoints.begin(), sourcePoints.end());

    timer.start();
    shortestPaths.build_sequence_tree();
    timer.stop();

    boost::timer::cpu_times elapsed = timer.elapsed();

    outData.constructionTime.add_sample(elapsed.wall);

#if !defined(NDEBUG)
    outData.peakMemoryUsage.add_sample((std::max)(shortestPaths.peak_memory_usage(), shortestPaths.current_memory_usage()));
#endif

    for (size_t j = 0; j < numQueries; ++j)
    {
      face_descriptor sourceFace = allFaces[rand.get_int(0, allFaces.size())];
      Barycentric_coordinates sourceLocation = random_coordinate<Traits>(rand);

      timer.start();
      FT distance = shortestPaths.shortest_distance_to_source_points(sourceFace, sourceLocation).first;
      timer.stop();
      UNUSED(distance);

      boost::timer::cpu_times elapsed = timer.elapsed();

      outData.queryTime.add_sample(elapsed.wall);
    }
  }
}

Kernel_type parse_kernel_type(const std::string& s)
{
  if (boost::iequals(s, "ipick"))
  {
    return KERNEL_IPICK;
  }
  else if (boost::iequals(s, "epick"))
  {
    return KERNEL_EPICK;
  }
  else if (boost::iequals(s, "epeck"))
  {
    return KERNEL_EPECK;
  }
  else
  {
    return KERNEL_UNKNOWN;
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
    ("trials,t", po::value<size_t>()->default_value(20), "Number of trials to run")
    ("numpoints,n", po::value<size_t>()->default_value(1), "Number of source points per trial")
    ("queries,q", po::value<size_t>()->default_value(100), "Number of queries to run per trial")
    ("kernel,k", po::value<std::string>()->default_value("epick"), "Kernel to use.  One of \'ipick\', \'epick\', \'epeck\'")
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
    Benchmark_data results;

    CGAL::Random rand(vm["randomseed"].as<size_t>());
    size_t numTrials = vm["trials"].as<size_t>();
    size_t numPoints = vm["numpoints"].as<size_t>();
    size_t numQueries = vm["queries"].as<size_t>();
    std::string polyhedronFile = vm["polyhedron"].as<std::string>();

    try
    {
      switch (parse_kernel_type(vm["kernel"].as<std::string>()))
      {
        case KERNEL_IPICK:
          run_benchmarks<CGAL::Simple_cartesian<double> >(rand, numTrials, numPoints, numQueries, polyhedronFile, results);
          break;
        case KERNEL_EPICK:
          run_benchmarks<CGAL::Exact_predicates_inexact_constructions_kernel>(rand, numTrials, numPoints, numQueries, polyhedronFile, results);
          break;
        case KERNEL_EPECK:
          run_benchmarks<CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt>(rand, numTrials, numPoints, numQueries, polyhedronFile, results);
          break;
        default:
          std::cerr << "Invalid kernel type: \"" << vm["kernel"].as<std::string>() << "\"" << std::endl;
          return EXIT_FAILURE;
      }
    }
    catch (std::exception& e)
    {
      std::cerr << "Runtime error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }

    print_results(std::cout, polyhedronFile, results);
  }
  else
  {
    std::cerr << "Please specify a polyhedron to use.  Use option --help for more details." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
