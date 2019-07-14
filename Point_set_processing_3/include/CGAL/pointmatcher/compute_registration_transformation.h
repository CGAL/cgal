// TODO: Copyright info

#ifndef CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H
#define CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Point_set_3.h>

#include <pointmatcher/PointMatcher.h>

namespace CGAL {

namespace pointmatcher {

namespace internal {

template<typename Scalar,
         typename PointRange,
         typename PointMap,
         typename VectorMap,
         typename PM_matrix>
void
copy_cgal_points_to_pm_matrix
(const PointRange& prange, PointMap point_map, VectorMap vector_map, PM_matrix& pm_points, PM_matrix& pm_normals)
{
  int idx = 0;
  for(const auto& p : prange)
  {
    // position
    const auto& pos = get(point_map, p);
    pm_points(0, idx) = pos.x();
    pm_points(1, idx) = pos.y();
    pm_points(2, idx) = pos.z();
    pm_points(3, idx) = Scalar(1.);
    
    // normal
    const auto& normal = get (vector_map, p);
    pm_normals(0, idx) = normal.x();
    pm_normals(1, idx) = normal.y();
    pm_normals(2, idx) = normal.z();
    
    ++idx;
  }
}

} // end of namespace internal

template <class Kernel,
          class PointRange1,
          class PointRange2,
          class PointMap1,
          class PointMap2,
          class VectorMap1,
          class VectorMap2>
typename Kernel::Aff_transformation_3
compute_registration_transformation(const PointRange1& range1, const PointRange2& range2,
                                    PointMap1 point_map1, PointMap2 point_map2,
                                    VectorMap1 vector_map1, VectorMap2 vector_map2)
{
  using Scalar    = typename Kernel::FT;
  
  using PM                  = PointMatcher<Scalar>;
  using PM_cloud            = typename PM::DataPoints;
  using PM_matrix           = typename PM::Matrix;
  using PM_labels           = typename PM_cloud::Labels;
  using PM_transform        = typename PM::Transformation;
  using PM_transform_params = typename PM::TransformationParameters;
  
  // ref_points: 1, points: 2
  std::size_t nb_ref_points = range1.size();
  std::size_t nb_points     = range2.size();
  
  PM_matrix ref_points_pos_matrix    = PM_matrix (4, nb_ref_points);
  PM_matrix ref_points_normal_matrix = PM_matrix (3, nb_ref_points);
  PM_matrix points_pos_matrix    = PM_matrix (4, nb_points);
  PM_matrix points_normal_matrix = PM_matrix (3, nb_points);
  
  // convert cgal points to pointmatcher points
  internal::copy_cgal_points_to_pm_matrix<Scalar>(range1,
                                                  point_map1,
                                                  vector_map1,
                                                  ref_points_pos_matrix, // out
                                                  ref_points_normal_matrix); // out

  internal::copy_cgal_points_to_pm_matrix<Scalar>(range2,
                                                  point_map2,
                                                  vector_map2,
                                                  points_pos_matrix, // out
                                                  points_normal_matrix); // out
                                        
  auto construct_PM_cloud = [](const PM_matrix& positions, const PM_matrix& normals) -> PM_cloud
  {
    PM_cloud cloud;
    
    cloud.addFeature("x", positions.row(0));
    cloud.addFeature("y", positions.row(1));
    cloud.addFeature("z", positions.row(2));
    cloud.addFeature("pad", positions.row(3));
    cloud.addDescriptor("normals", normals);
    
    return cloud;
  };

  PM_cloud ref_cloud = construct_PM_cloud(ref_points_pos_matrix, ref_points_normal_matrix);
  PM_cloud cloud     = construct_PM_cloud(points_pos_matrix,     points_normal_matrix);

  typename PM::ICP icp;
  
  // TODO: Make it configurable? Use named parameters to have possible settings?
  icp.setDefault();
  
  PM_transform_params transform_params = PM_transform_params::Identity(4,4);
  try 
  {
		const PM_transform_params prior = transform_params;
    std::cerr << "Cloud points nb: " << cloud.getNbPoints() << std::endl;
    std::cerr << "Ref Cloud points nb: " << ref_cloud.getNbPoints() << std::endl;
		transform_params = icp(cloud, ref_cloud, prior);
    std::cerr << transform_params << std::endl;
	}
	catch (typename PM::ConvergenceError& error)
	{
		std::cerr << "ERROR PM::ICP failed to converge: " << std::endl;
		std::cerr << "   " << error.what() << std::endl;
    // TODO: What to do?
		//continue;
	}

	// Rigid transformation
	std::shared_ptr<PM_transform> transform = PM::get().REG(Transformation).create("RigidTransformation");
  transform_params = transform->correctParameters(transform_params);

  typename Kernel::Aff_transformation_3 cgal_transform
    (transform_params(0,0), transform_params(0,1), transform_params(0,2), transform_params(0,3),
     transform_params(1,0), transform_params(1,1), transform_params(1,2), transform_params(1,3),
     transform_params(2,0), transform_params(2,1), transform_params(2,2), transform_params(2,3));
  
  return cgal_transform;
}

} } // end of namespace CGAL::pointmatcher

#endif // CGAL_POINTMATCHER_COMPUTE_REGISTRATION_TRANSFORMATION_H
