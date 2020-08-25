#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Classification/Point_set_neighborhood.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Fuzzy_sphere.h>

#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>

#include <CGAL/hierarchy_simplify_point_set.h>

#include <boost/filesystem.hpp>

#include <CGAL/OpenGR/gret_sdp.h>

#include <vector>
#include <string>
#include <numeric>


// For computations 3D space
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Point_range = std::vector<Point_3>;
using Pwn = std::pair<Point_3, Vector_3>;
using Pwn_range = std::vector<Pwn>;
using Point_map = CGAL::First_of_pair_property_map<Pwn>;
using Normal_map = CGAL::Second_of_pair_property_map<Pwn>;

// KD tree in N dimension
using Search_traits_base = CGAL::Search_traits_3<Kernel>;
using Point_3_map = typename CGAL::Pointer_property_map<Point_3>::type;;
using Search_traits = CGAL::Search_traits_adapter<std::size_t, Point_3_map, Search_traits_base>;
using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
using Tree = typename Knn::Tree;
using Splitter = typename Knn::Splitter;
using Distance = typename Knn::Distance;
using Tree_ptr = std::unique_ptr<Tree>;
using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_traits>;

// convenience definitions for correspondences
using Correspondence = std::pair<size_t, size_t>;
using CorrespondenceRange = std::vector<Correspondence>;

// extracts point clouds and transformations from a stanford config file
template <typename PointMap, typename PointRange, typename Transformation>
void extractPCAndTrFromStandfordConfFile(
        const std::string &confFilePath,
        std::vector<Transformation>& transforms,
        std::vector<PointRange>& point_clouds,
        const PointMap& point_map
        ){
    using namespace boost;
    using namespace std;

    typedef Eigen::Transform<double, 3, Eigen::Affine> Transform;

    std::vector<string> files;
    
    if(!(filesystem::exists(confFilePath) && filesystem::is_regular_file(confFilePath)))
      throw std::runtime_error("Config file does not exist or is no regular file.");

    // extract the working directory for the configuration path
    const std::string workingDir = filesystem::path(confFilePath).parent_path().native();    
    if(!filesystem::exists(confFilePath))
      throw std::runtime_error("Directory \"" + workingDir + "\" does not exist.");

    // read the configuration file and call the matching process
    std::string line;
    std::ifstream confFile;
    confFile.open(confFilePath);
    if(!confFile.is_open())
      throw std::runtime_error("Could not open config file.");

    while ( getline (confFile,line) )
    {
        std::istringstream iss (line);
        std::vector<string> tokens{istream_iterator<string>{iss},
                              istream_iterator<string>{}};

        // here we know that the tokens are:
        // [0]: keyword, must be bmesh
        // [1]: 3D object filename
        // [2-4]: target translation with previous object
        // [5-8]: target quaternion with previous object

        if (tokens.size() == 9){
            if (tokens[0].compare("bmesh") == 0){
                std::string inputfile = filesystem::path(confFilePath).parent_path().string()+string("/")+tokens[1];
                if(!(filesystem::exists(inputfile) && filesystem::is_regular_file(inputfile)))
                  throw std::runtime_error("File \"" + inputfile + "\" does not exist or is no regular file.");

                // build the Eigen rotation matrix from the rotation and translation stored in the files
                Eigen::Matrix<double, 3, 1> tr (
                            std::atof(tokens[2].c_str()),
                            std::atof(tokens[3].c_str()),
                            std::atof(tokens[4].c_str()));

                Eigen::Quaternion<double> quat(
                            std::atof(tokens[8].c_str()), // eigen starts by w
                            std::atof(tokens[5].c_str()),
                            std::atof(tokens[6].c_str()),
                            std::atof(tokens[7].c_str()));

                quat.normalize();

                Transform transform (Transform::Identity());
                transform.rotate(quat);
                transform.translate(-tr);

                transforms.emplace_back(
                transform(0,0), transform(0,1), transform(0,2), transform(0,3),
                transform(1,0), transform(1,1), transform(1,2), transform(1,3),
                transform(2,0), transform(2,1), transform(2,2), transform(2,3)
                );

                files.push_back(inputfile);
            }
        }
    }
    confFile.close();

    std::ifstream pc_file;
    int num_point_clouds = files.size();
    point_clouds.resize(num_point_clouds);
    for(int i = 0; i < num_point_clouds; i++){
        const string& file = files[i];
        pc_file.open(file);
        if(!pc_file ||
            !CGAL::read_ply_points(pc_file, std::back_inserter(point_clouds[i]), CGAL::parameters::point_map(point_map)))
        {
          std::cerr << "Error: cannot read file " << file << std::endl;
          throw std::exception();
        } 
        pc_file.close();
    }
}

void constructKdTree(Point_range& point_range, Tree_ptr& tree, Distance& distance){
  Point_3_map point_map = CGAL::make_property_map(point_range);
  tree.reset(
    new Tree(
    boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(point_range.size()),
    Splitter(), Search_traits(point_map)
    ));
  tree->build();
  distance = Distance(point_map);
}

template <typename PointMap, typename TransformationRange>
void computeCorrespondences(const std::vector<Pwn_range>& point_clouds, const PointMap& point_map, const TransformationRange& transformations, std::vector<CorrespondenceRange>& correspondences, 
                            int sampling_number = 200, const double max_dist = 0.00001){
  int num_point_clouds = point_clouds.size();

  Point_range merged_point_cloud;
  std::vector<Point_range> transformed_point_clouds(num_point_clouds);
  // construct transformed point clouds
  for (size_t i = 0; i < num_point_clouds; i++){
    for (size_t j = 0; j < point_clouds[i].size(); j++){
      Point_3 point = get(point_map, point_clouds[i][j]).transform(transformations[i].inverse());
      transformed_point_clouds[i].push_back(point);
      merged_point_cloud.push_back(point);
    }
  }

  sampling_number = merged_point_cloud.size() < sampling_number? merged_point_cloud.size() : sampling_number;
  merged_point_cloud.erase(CGAL::hierarchy_simplify_point_set(merged_point_cloud, 
                            CGAL::parameters::size((double)merged_point_cloud.size()/(double)sampling_number)
                            .maximum_variation(0.33)),
                            merged_point_cloud.end());
  
  int num_global_coordinates = merged_point_cloud.size();

  // construct a KD tree in N dimensions 
  std::vector<Tree_ptr> trees(num_point_clouds);
  std::vector<Distance> distances(num_point_clouds);
  for (size_t i = 0; i < num_point_clouds; i++)
    constructKdTree(transformed_point_clouds[i], trees[i], distances[i]);

  // construct correspondences
  for (size_t i = 0; i < merged_point_cloud.size(); i++) {
    CorrespondenceRange correspondence;
    const Point_3& query_point = merged_point_cloud[i];
    for (size_t j = 0; j < num_point_clouds; j++) {
      Knn knn(*trees[j], query_point, 1, 0, true, distances[j]);
      double dist = knn.begin()->second;
      if(dist < max_dist){
        std::size_t nn = knn.begin()->first;
        correspondence.emplace_back(j, nn);
      } 
    }
    if(!correspondence.empty())
      correspondences.push_back(std::move(correspondence));
  }       
}

int main (int argc, char** argv)
{
  const char* config_fname = (argc>1)?argv[1]:"data/bun.conf";
  
  std::vector<Pwn_range> point_clouds;
  Point_map point_map;
  std::vector<Kernel::Aff_transformation_3> ground_truth_transformations;

  extractPCAndTrFromStandfordConfFile(config_fname, ground_truth_transformations, point_clouds, point_map);
  int num_point_clouds = point_clouds.size();

  std::vector<CorrespondenceRange> correspondences;
  computeCorrespondences(point_clouds, point_map, ground_truth_transformations, correspondences);

  // EITHER call the registration method GRET-SDP from OpenGR to get the transfromations for registration
  std::vector<Kernel::Aff_transformation_3> transformations(num_point_clouds);
  CGAL::OpenGR::compute_registration_transformations(point_clouds, correspondences, CGAL::parameters::point_map(point_map).normal_map(Normal_map()), transformations);

  // OR call  the registration method GRET-SDP from OpenGR and apply transformations directly 
  int num_total_points = std::accumulate(point_clouds.begin(), point_clouds.end(), 0, [](int sum, const auto& pc){ return sum + pc.size(); });
  Pwn_range registered_point_cloud(num_total_points);
  CGAL::OpenGR::register_point_clouds(point_clouds, correspondences, CGAL::parameters::point_map(point_map).normal_map(Normal_map()), registered_point_cloud);

  std::ofstream out("registered_point_clouds.ply");
  if (!out ||
    !CGAL::write_ply_points(
      out, registered_point_cloud,
      CGAL::parameters::point_map(point_map)))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}