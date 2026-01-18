
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
#include <vtkNIFTIImageReader.h>
#include <vtkImageReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkDemandDrivenPipeline.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/IO/read_vtk_image_data.h>

#include <boost/lexical_cast.hpp>
#include <boost/functional.hpp>

#include <filesystem>

typedef short Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

template <bool less_than_iso>
class Is_interior {
  double iso;
public:
  Is_interior(double iso): iso(iso) {}

  template <typename T>
  int operator()(T v) const {
    return int( (v < iso) == less_than_iso);
  }
};

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
  // Loads image

  // Usage: mesh_3D_gray_vtk_image <nii file or dicom directory> iso_level=1  facet_size=1  facet_distance=0.1  cell_size=1 interior_less_than_isolevel=true

  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/squircle.nii");

  vtkSmartPointer<vtkImageData> vtk_image = nullptr;
  Image_word_type iso = (argc>2)? boost::lexical_cast<Image_word_type>(argv[2]): 1;
  double fs = (argc>3)? boost::lexical_cast<double>(argv[3]): 1;
  double fd = (argc>4)? boost::lexical_cast<double>(argv[4]): 0.1;
  double cs = (argc>5)? boost::lexical_cast<double>(argv[5]): 1;
  bool less = (argc>6)? boost::lexical_cast<bool>(argv[6]): true;

  fs::path path(fname);

  if(fs::is_regular_file(path)){
    std::cout << "regular file" << std::endl;
    if (path.has_extension()){
      fs::path stem = path.stem();
      if ((path.extension() == ".nii") || (stem.has_extension() && (stem.extension() == ".nii") && (path.extension() == ".gz"))) {
        auto reader = vtkSmartPointer<vtkNIFTIImageReader>::New();
        reader->SetFileName(fname.c_str());
        reader->Update();
        vtk_image = reader->GetOutput();
        vtk_image->Print(std::cerr);
      }
    }
  }
  else if (fs::is_directory(path)) {
    auto dicom_reader = vtkSmartPointer<vtkDICOMImageReader>::New();
    dicom_reader->SetDirectoryName(argv[1]);

    vtkDemandDrivenPipeline* executive =
      vtkDemandDrivenPipeline::SafeDownCast(dicom_reader->GetExecutive());
    if (executive)
      {
        executive->SetReleaseDataFlag(0, 0); // where 0 is the port index
      }

    auto smoother = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smoother->SetStandardDeviations(1., 1., 1.);
    smoother->SetInputConnection(dicom_reader->GetOutputPort());
    smoother->Update();
    vtk_image = smoother->GetOutput();
    vtk_image->Print(std::cerr);
  }
  if(vtk_image == nullptr){
    std::cout << "No image loaded" << std::endl;
    return 0;
  }
  CGAL::Image_3 image = CGAL::IO::read_vtk_image_data(vtk_image);
  if(image.image() == nullptr){
    std::cerr << "could not create a CGAL::Image_3 from the vtk image\n";
    return 0;
  }
  /// [Domain creation]
  namespace params = CGAL::parameters;

  Mesh_domain domain = less ? Mesh_domain::create_gray_image_mesh_domain(image,
                                params::image_values_to_subdomain_indices(Is_interior<true>(iso)).
                                value_outside(iso+1))
                            : Mesh_domain::create_gray_image_mesh_domain(image,
                                params::image_values_to_subdomain_indices(Is_interior<false>(iso)).
                                value_outside(iso+1));
  /// [Domain creation]

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_size(fs).facet_distance(fd).
                         cell_radius_edge_ratio(3).cell_size(cs));

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  return 0;
}
