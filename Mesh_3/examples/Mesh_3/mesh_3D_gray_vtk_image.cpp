
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
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
#include <CGAL/read_vtk_image_data.h>

#include <boost/lexical_cast.hpp>
#include <boost/functional.hpp>

typedef short Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

class Less {
  double iso;
public:
  Less(double iso): iso(iso) {}

  template <typename T>
  int operator()(T v) const {
    return int(v < iso);
  }
};

int main(int argc, char* argv[])
{
  // Loads image
  if(argc == 1){
    std::cerr << "Usage:  " << argv[0] << " <directory with dicom data> iso_level=1  facet_size=1  facet_distance=0.1  cell_size=1\n";
    return 0;
  }

  Image_word_type iso = (argc>2)? boost::lexical_cast<Image_word_type>(argv[2]): 1;
  double fs = (argc>3)? boost::lexical_cast<double>(argv[3]): 1;
  double fd = (argc>4)? boost::lexical_cast<double>(argv[4]): 0.1;
  double cs = (argc>5)? boost::lexical_cast<double>(argv[5]): 1;

  vtkDICOMImageReader*dicom_reader = vtkDICOMImageReader::New();
  dicom_reader->SetDirectoryName(argv[1]);

  vtkDemandDrivenPipeline*executive =
    vtkDemandDrivenPipeline::SafeDownCast(dicom_reader->GetExecutive());
  if (executive)
    {
      executive->SetReleaseDataFlag(0, 0); // where 0 is the port index
    }

  vtkImageGaussianSmooth* smoother = vtkImageGaussianSmooth::New();
  smoother->SetStandardDeviations(1., 1., 1.);
  smoother->SetInputConnection(dicom_reader->GetOutputPort());
  smoother->Update();
  vtkImageData* vtk_image = smoother->GetOutput();
  vtk_image->Print(std::cerr);

  CGAL::Image_3 image = CGAL::read_vtk_image_data(vtk_image);
  if(image.image() == 0){
    std::cerr << "could not create a CGAL::Image_3 from the vtk image\n";
    return 0;
  }
  /// [Domain creation]
  // To avoid verbose function and named parameters call
  using namespace CGAL::parameters;

  Mesh_domain domain = Mesh_domain::create_gray_image_mesh_domain
    (image,
     image_values_to_subdomain_indices = Less(iso),
     value_outside = 0);
  /// [Domain creation]

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=fs, facet_distance=fd,
                         cell_radius_edge_ratio=3, cell_size=cs);

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}
