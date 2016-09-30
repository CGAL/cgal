#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Planimetric_grid.h>
#include <CGAL/Data_classification/Attribute.h>
#include <CGAL/Data_classification/Attribute_vertical_dispersion.h>
#include <CGAL/Data_classification/Attribute_elevation.h>
#include <CGAL/Data_classification/Attribute_verticality.h>
#include <CGAL/Data_classification/Attribute_distance_to_plane.h>
#include <CGAL/Data_classification/Attributes_eigen.h>
#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point>::iterator Iterator;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Point_set_classification<Kernel, Iterator, Pmap> Classification;

typedef CGAL::Data_classification::Planimetric_grid<Kernel, Iterator, Pmap>      Planimetric_grid;
typedef CGAL::Data_classification::Neighborhood<Kernel, Iterator, Pmap>          Neighborhood;
typedef CGAL::Data_classification::Local_eigen_analysis<Kernel, Iterator, Pmap>  Local_eigen_analysis;
typedef CGAL::Data_classification::Attribute Attribute;
typedef CGAL::Data_classification::Attribute_vertical_dispersion<Kernel, Iterator, Pmap> Dispersion;
typedef CGAL::Data_classification::Attribute_elevation<Kernel, Iterator, Pmap>           Elevation;
typedef CGAL::Data_classification::Attribute_verticality<Kernel, Iterator, Pmap>         Verticality;
typedef CGAL::Data_classification::Attribute_distance_to_plane<Kernel, Iterator, Pmap>   Distance_to_plane;
typedef CGAL::Data_classification::Attribute_linearity<Kernel, Iterator, Pmap>           Linearity;
typedef CGAL::Data_classification::Attribute_planarity<Kernel, Iterator, Pmap>           Planarity;
typedef CGAL::Data_classification::Attribute_sphericity<Kernel, Iterator, Pmap>          Sphericity;
typedef CGAL::Data_classification::Attribute_omnivariance<Kernel, Iterator, Pmap>        Omnivariance;
typedef CGAL::Data_classification::Attribute_anisotropy<Kernel, Iterator, Pmap>          Anisotropy;
typedef CGAL::Data_classification::Attribute_eigentropy<Kernel, Iterator, Pmap>          Eigentropy;
typedef CGAL::Data_classification::Attribute_sum_eigenvalues<Kernel, Iterator, Pmap>     Sum_eigen;
typedef CGAL::Data_classification::Attribute_surface_variation<Kernel, Iterator, Pmap>   Surface_variation;



class My_ply_interpreter
{
  std::vector<Point>& points;
  std::vector<int>& labels;
    
public:
  My_ply_interpreter (std::vector<Point>& points,
                      std::vector<int>& labels)
    : points (points), labels (labels)
  { }

  // Init and test if input file contains the right properties
  bool is_applicable (CGAL::Ply_reader& reader)
  {
    return reader.does_tag_exist<double> ("x")
      && reader.does_tag_exist<double> ("y")
      && reader.does_tag_exist<double> ("z")
      && reader.does_tag_exist<int> ("label");
  }

  // Describes how to process one line (= one point object)
  void process_line (CGAL::Ply_reader& reader)
  {
    double x = 0., y = 0., z = 0.;
    int l = 0;

    reader.assign (x, "x");
    reader.assign (y, "y");
    reader.assign (z, "z");
    reader.assign (l, "label");

    points.push_back (Point (x, y, z));
    labels.push_back(l);
  }

};


int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/b9.ply");
  std::ifstream in (filename.c_str());
  std::vector<Point> pts;
  std::vector<int> labels;

  std::cerr << "Reading input" << std::endl;
  My_ply_interpreter interpreter (pts, labels);
  if (!in
      || !(CGAL::read_ply_custom_points (in, interpreter, Kernel())))
    {
      std::cerr << "Error: cannot read " << filename << std::endl;
      return EXIT_FAILURE;
    }

  double grid_resolution = 1.;
  double radius_neighbors = 5.;
  double radius_dtm = 30;

  Classification psc (pts.begin (), pts.end(), Pmap());
  
  std::cerr << "Computing useful structures" << std::endl;
  Iso_cuboid_3 bbox = CGAL::bounding_box (pts.begin(), pts.end());
  Neighborhood neighborhood (pts.begin(), pts.end(), Pmap());  
  std::cerr << "Computing attributes" << std::endl;

  for (std::size_t i = 0; i < 4; ++ i)
    {
      Planimetric_grid grid (pts.begin(), pts.end(), Pmap(), bbox, grid_resolution);

      Local_eigen_analysis eigen (pts.begin(), pts.end(), Pmap(), neighborhood, radius_neighbors);
      
      psc.add_segmentation_attribute (new Dispersion (pts.begin(), pts.end(), Pmap(), grid,
                                                      grid_resolution,
                                                      radius_neighbors));
      psc.add_segmentation_attribute (new Elevation (pts.begin(), pts.end(), Pmap(), bbox, grid,
                                                     grid_resolution,
                                                     radius_neighbors,
                                                     radius_dtm));
      psc.add_segmentation_attribute (new Verticality (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Distance_to_plane (pts.begin(), pts.end(), Pmap(), eigen));
      psc.add_segmentation_attribute (new Linearity (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Planarity (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Sphericity (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Omnivariance (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Anisotropy (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Eigentropy (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Sum_eigen (pts.begin(), pts.end(), eigen));
      psc.add_segmentation_attribute (new Surface_variation (pts.begin(), pts.end(), eigen));
      grid_resolution *= 2;
      radius_neighbors *= 2;
      radius_dtm *= 2;
    }
  
  // Add types to PSC
  CGAL::Data_classification::Type ground ("ground");
  CGAL::Data_classification::Type vege ("vegetation");
  CGAL::Data_classification::Type roof ("roof");
  CGAL::Data_classification::Type facade ("facade");

  for (std::size_t i = 0; i < labels.size(); ++ i)
    {
      switch (labels[i])
        {
        case 0:
          vege.training_set().push_back (i);
          break;
        case 1:
          ground.training_set().push_back (i);
          break;
        case 2:
          roof.training_set().push_back (i);
          break;
        case 3:
          facade.training_set().push_back (i);
          break;
        default:
          break;
        }
    }

  psc.add_classification_type (&vege);
  psc.add_classification_type (&ground);
  psc.add_classification_type (&roof);
  psc.add_classification_type (&facade);

  //  psc.set_multiplicative(true);
  std::cerr << "Training" << std::endl;
  psc.training();

  psc.run_with_graphcut (neighborhood, 0.5);
  //psc.run_quick();
  
  // Save the output in a colored PLY format

  std::ofstream f ("classification.ply");
  f.precision(18);
  f << "ply" << std::endl
    << "format ascii 1.0" << std::endl
    << "element vertex " << pts.size() << std::endl
    << "property double x" << std::endl
    << "property double y" << std::endl
    << "property double z" << std::endl
    << "property uchar red" << std::endl
    << "property uchar green" << std::endl
    << "property uchar blue" << std::endl
    << "property int label" << std::endl
    << "end_header" << std::endl;
  
  for (std::size_t i = 0; i < pts.size(); ++ i)
    {
      f << pts[i] << " ";
      
      CGAL::Data_classification::Type* type = psc.classification_type_of (i);
      if (type == &vege)
        f << "0 255 27 0" << std::endl;
      else if (type == &ground)
        f << "245 180 0 1" << std::endl;
      else if (type == &roof)
        f << "255 0 170 2" << std::endl;
      else if (type == &facade)
        f << "100 0 255 3" << std::endl;
      else
        {
          f << "0 0 0 -1" << std::endl;
          //          std::cerr << "Error: unknown classification type" << std::endl;
        }
    }

  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
