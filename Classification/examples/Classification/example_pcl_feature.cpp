#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/shot.h>
#include <pcl/common/io.h>
#include <pcl/features/don.h>
#include <pcl/filters/conditional_removal.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

typedef CGAL::Simple_cartesian<double>                                          Kernel;
typedef Kernel::Point_3                                                         Point;
typedef Kernel::Iso_cuboid_3                                                    Iso_cuboid_3;
typedef std::vector<Point>                                                      Point_range;
typedef CGAL::Identity_property_map<Point>                                      Pmap;

namespace Classification = CGAL::Classification;

typedef Classification::ETHZ::Random_forest_classifier                          Classifier;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef CGAL::Point_set_3<Point>                                                Point_set;
typedef Point_set::Point_map                                                    Point_map;
typedef Point_set::Property_map<int>                                            Imap;
typedef Point_set::Property_map<unsigned char>                                  UCmap;
typedef Classification::Point_set_neighborhood<Kernel, Point_set, Point_map>    Neighborhood;

typedef pcl::SHOT352                                                            pcl_SHOT;

///////////////////////////////////////////////////////////////////
//! [Feature]

// User-defined feature that communicates the given dimension of
// the SHOT feature descriptor computed using PCL.
class cgal_SHOT_dim : public CGAL::Classification::Feature_base
{
  std::vector<float> m_values;
public:
  cgal_SHOT_dim(std::size_t i, const pcl::PointCloud<pcl_SHOT>::Ptr shotFeatures)
  {
    this->set_name("SHOT[" + std::to_string(i) + "]");
    m_values.reserve(shotFeatures->size());
    for (std::size_t j = 0; j < shotFeatures->size(); ++j) {
      m_values.push_back(shotFeatures->at(j).descriptor[i]);
    }
  }

  float value(std::size_t pt_index)
  {
    return m_values.at(pt_index);
  }
};

// User-defined feature that communicates the given dimension of
// the DifferenceOfNormalsEstimation descriptor computed using PCL.
class cgal_DON_dim : public CGAL::Classification::Feature_base
{
    std::vector<float> m_values;
public:
    cgal_DON_dim(std::size_t i, const pcl::PointCloud<pcl::Normal>::Ptr DoNFeatures)
    {
        this->set_name("DoN[" + std::to_string(i) + "]");
        m_values.reserve(DoNFeatures->size());
        for (std::size_t j = 0; j < DoNFeatures->size(); ++j) {
            m_values.push_back(DoNFeatures->at(j).normal[i]);
        }
    }

    float value(std::size_t pt_index)
    {
        return m_values.at(pt_index);
    }
};

class cgal_SHOT : public CGAL::Classification::Feature_base_multi_dim
{
    pcl::PointCloud<pcl_SHOT>::Ptr m_features;
public:
    cgal_SHOT(const pcl::PointCloud<pcl_SHOT>::Ptr shotFeatures) : m_features(shotFeatures)
    {
        this->set_name("SHOT");
    }

    float value(std::size_t pt_index, std::size_t dim_index)
    {
        return m_features->at(pt_index).descriptor[dim_index];
    }
};


//! [Feature]
///////////////////////////////////////////////////////////////////

int main (int argc, char** argv)
{
  const std::string filename(argc > 1 ? argv[1] : CGAL::data_file_path("points_3/b9_training.ply"));

  std::cerr << "Reading input" << std::endl;
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in)
  {
    std::cerr << "Error: cannot read \"" << filename << "\"." << std::endl;
  }

  Point_set pts;
  in >> pts;

  Imap label_map;
  bool lm_found = false;
  std::tie(label_map, lm_found) = pts.property_map<int>("label");
  if (!lm_found)
  {
    std::cerr << "Error: \"label\" property not found in input file." << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Bbox_3 bbox = CGAL::bbox_3(pts.points().begin(), pts.points().end());
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  std::cout << "diag_length = " << diag_length << std::endl;
  const double scale = diag_length * 0.075f;

  std::cerr << "Converting point cloud" << std::endl;
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
  cloud->width = pts.number_of_points();
  cloud->height = 1;
  cloud->is_dense = false;
  cloud->resize(cloud->width * cloud->height);
  for (std::size_t i = 0; i < pts.number_of_points(); ++i) {
    cloud->at(i) = pcl::PointXYZ(pts.point(i).x(), pts.point(i).y(), pts.point(i).z());
  }
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
  normalEstimation.setInputCloud(cloud);

  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  normalEstimation.setSearchMethod(tree);

  std::cerr << "Computing normals" << std::endl;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud< pcl::Normal>);
  normalEstimation.setRadiusSearch(scale);
  normalEstimation.compute(*normals);

  std::cerr << "Computing features" << std::endl;
  pcl::SHOTEstimation<pcl::PointXYZ, pcl::Normal, pcl_SHOT> shotEstimation;
  shotEstimation.setInputCloud(cloud);
  shotEstimation.setInputNormals(normals);

  shotEstimation.setSearchMethod(tree);
  pcl::PointCloud<pcl_SHOT>::Ptr shotFeatures(new pcl::PointCloud<pcl_SHOT>);
  shotEstimation.setKSearch(0);
  shotEstimation.setRadiusSearch(scale * 2.0f);

  shotEstimation.compute(*shotFeatures);


  Feature_set features;
  constexpr std::size_t dims = (std::size_t)pcl_SHOT::descriptorSize();
  features.add_multidimensional_feature2<cgal_SHOT>(dims, shotFeatures);

  //return EXIT_SUCCESS;

  pcl::PointCloud<pcl::Normal>::Ptr normals_large_scale(new pcl::PointCloud<pcl::Normal>);
  normalEstimation.setRadiusSearch(scale * 2.0f);
  normalEstimation.compute(*normals_large_scale);

  pcl::DifferenceOfNormalsEstimation<pcl::PointXYZ, pcl::Normal, pcl::Normal> don;
  don.setInputCloud(cloud);
  don.setNormalScaleSmall(normals);
  don.setNormalScaleLarge(normals_large_scale);

  pcl::PointCloud<pcl::Normal>::Ptr doncloud(new pcl::PointCloud<pcl::Normal>);
  doncloud->width = cloud->width;
  doncloud->height = cloud->height;
  doncloud->is_dense = cloud->is_dense;
  doncloud->resize(doncloud->width * doncloud->height);
 
  don.initCompute();
  don.computeFeature(*doncloud);


  std::vector<Feature_handle> shot_feature_handles = features.add_multidimensional_feature<cgal_SHOT_dim>(dims, shotFeatures);

  features.add_multidimensional_feature<cgal_DON_dim>(3, doncloud);

  std::cout << "features.size() = " << features.size() << std::endl;
  features.remove(shot_feature_handles);
  std::cout << "features.size() = " << features.size() << std::endl;
  for (auto& feature : features) {
      std::cout << feature->name() << std::endl;
  }
  shot_feature_handles = features.add_multidimensional_feature<cgal_SHOT_dim>(dims, shotFeatures);
  std::cout << "features.size() = " << features.size() << std::endl;
  features.remove(shot_feature_handles);
  std::vector<Feature_handle> don_feature_handles = features.add_multidimensional_feature<cgal_DON_dim>(3, doncloud);
  std::cout << "features.size() = " << features.size() << std::endl;
  features.remove(don_feature_handles);
  std::cout << "features.size() = " << features.size() << std::endl;
  for (auto& feature : features) {
      std::cout << feature->name() << std::endl;
  }

  Neighborhood neighborhood(pts, pts.point_map());

  Label_set labels;
  Label_handle ground = labels.add("ground");
  Label_handle vegetation = labels.add("vegetation");
  Label_handle roof = labels.add("roof");

  Classifier classifier (labels, features);

  std::cerr << "Training classifier" << std::endl;
  classifier.train (pts.range(label_map));

  std::cerr << "Classifying" << std::endl;
  std::vector<int> label_indices(pts.number_of_points(), -1);
  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (pts, pts.point_map(), labels, classifier,
     neighborhood.k_neighbor_query(12),
     0.5, 1, label_indices);

  UCmap red = pts.add_property_map<unsigned char>("red", 0).first;
  UCmap green = pts.add_property_map<unsigned char>("green", 0).first;
  UCmap blue = pts.add_property_map<unsigned char>("blue", 0).first;
  for (std::size_t i = 0; i < label_indices.size(); ++i)
  {
    Label_handle label = labels[label_indices[i]];
    const CGAL::IO::Color& color = label->color();
    red[i] = color.red();
    green[i] = color.green();
    blue[i] = color.blue();
  }

  CGAL::IO::write_PLY("classification.ply", pts, CGAL::parameters::use_binary_mode(true));

  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
