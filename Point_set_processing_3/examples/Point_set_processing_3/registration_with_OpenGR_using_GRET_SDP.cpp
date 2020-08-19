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

// For computations 3D space
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Point_range = std::vector<Point_3>;
// Point with global coordinate index
using Indexed_Point = std::pair<Point_3, int>;
using Point_map = CGAL::First_of_pair_property_map<Indexed_Point>;
using Index_map = CGAL::Second_of_pair_property_map<Indexed_Point>;

// KD tree in N dimension
using Search_traits_base = CGAL::Search_traits_3<Kernel>;
using Point_3_map = typename CGAL::Pointer_property_map<Point_3>::type;
using Search_traits = CGAL::Search_traits_adapter<std::size_t, Point_3_map, Search_traits_base>;
using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
using Tree = typename Knn::Tree;
using Splitter = typename Knn::Splitter;
using Distance = typename Knn::Distance;
using Tree_ptr = std::unique_ptr<Tree>;
using Fuzzy_sphere = CGAL::Fuzzy_sphere<Search_traits>;

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
    
    //VERIFY (filesystem::exists(confFilePath) && filesystem::is_regular_file(confFilePath));

    // extract the working directory for the configuration path
    const std::string workingDir = filesystem::path(confFilePath).parent_path().native();
    //VERIFY (filesystem::exists(workingDir));

    // read the configuration file and call the matching process
    std::string line;
    std::ifstream confFile;
    confFile.open(confFilePath);
    //VERIFY (confFile.is_open());

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
                //VERIFY(filesystem::exists(inputfile) && filesystem::is_regular_file(inputfile));

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

void constructKdTree(Point_range& feature_range, Tree_ptr& tree, Distance& distance){
  Point_3_map point_d_map = CGAL::make_property_map(feature_range);
  tree.reset(
    new Tree(
    boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(feature_range.size()),
    Splitter(), Search_traits(point_d_map)
    ));
  tree->build();
  distance = Distance(point_d_map);
}

template <typename PointRange, typename TransformationRange, typename IndexedPointRange>
int computeCorrespondences(const std::vector<PointRange>& point_clouds, const TransformationRange& transformations, std::vector<IndexedPointRange>& patches, 
                            int sampling_number = 200, const double max_dist = 0.00001){
  using PointType = typename PointRange::value_type;
  int num_point_clouds = point_clouds.size();

  Point_range merged_point_cloud;
  std::vector<Point_range> transformed_point_clouds(num_point_clouds);
  // construct transformed patches
  for (size_t i = 0; i < num_point_clouds; i++){
    for (size_t j = 0; j < point_clouds[i].size(); j++){
      const PointType point = point_clouds[i][j].transform(transformations[i].inverse());
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
    Point_3& query_point = merged_point_cloud[i];
    for (size_t j = 0; j < patches.size(); j++) {
      Knn knn(*trees[j], query_point, 1, 0, true, distances[j]);
      double dist = knn.begin()->second;
      if(dist < max_dist){
        std::size_t nn = knn.begin()->first;
        patches[j].emplace_back(point_clouds[j][nn] , i);
      } 
    }
  }       

  return num_global_coordinates;
}

template <typename PointRange, typename TransformationRange>
void transformAndMergePointSets(const std::vector<PointRange>& point_clouds, const TransformationRange& transformations, PointRange& registered_point_cloud){
  for (size_t i = 0; i < point_clouds.size(); i++)
    for (size_t j = 0; j < point_clouds[i].size(); j++)
      registered_point_cloud.push_back(point_clouds[i][j].transform(transformations[i]));
}


int main (int argc, char** argv)
{
  const char* config_fname = (argc>1)?argv[1]:"gret-sdp-data/bun.conf";

  std::vector<Point_range> point_clouds;
  CGAL::Identity_property_map<Point_3> point_map;
  std::vector<Kernel::Aff_transformation_3> ground_truth_transformations;

  extractPCAndTrFromStandfordConfFile(config_fname, ground_truth_transformations, point_clouds, point_map);
  int num_point_clouds = point_clouds.size();

  std::vector<std::vector<Indexed_Point>> patches(num_point_clouds);
  int num_global_coordinates = computeCorrespondences(point_clouds, ground_truth_transformations, patches);

  CGAL::OpenGR::GRET_SDP<Kernel> matcher;
  matcher.registerPatches(patches, num_global_coordinates, CGAL::parameters::point_map(Point_map())
                                                .vertex_index_map(Index_map()));

  std::vector<Kernel::Aff_transformation_3> computed_transformations;
  matcher.getTransformations(computed_transformations);

  Point_range registered_point_cloud;
  transformAndMergePointSets(point_clouds, computed_transformations, registered_point_cloud);

  std::ofstream out("registered_point_clouds.ply");
  if (!out ||
    !CGAL::write_ply_points(
      out, registered_point_cloud,
      CGAL::parameters::point_map(CGAL::Identity_property_map<Point_3>())))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}