#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Planimetric_grid.h>
#include <CGAL/Data_classification/Segmentation_attribute_vertical_dispersion.h>
#include <CGAL/Data_classification/Segmentation_attribute_elevation.h>
#include <CGAL/Data_classification/Segmentation_attribute_verticality.h>
#include <CGAL/Data_classification/Segmentation_attribute_distance_to_plane.h>
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
typedef CGAL::Segmentation_attribute Attribute;
typedef CGAL::Segmentation_attribute_vertical_dispersion<Kernel, Iterator, Pmap> Dispersion;
typedef CGAL::Segmentation_attribute_elevation<Kernel, Iterator, Pmap>           Elevation;
typedef CGAL::Segmentation_attribute_verticality<Kernel, Iterator, Pmap>         Verticality;
typedef CGAL::Segmentation_attribute_distance_to_plane<Kernel, Iterator, Pmap>   Distance_to_plane;


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


struct Linearity : Attribute
{
  std::vector<double> attrib;

  Linearity (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back ((ev[2] - ev[1]) / (ev[2] + 1e-10));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "linearity"; }
};

struct Planarity : Attribute
{
  std::vector<double> attrib;

  Planarity (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back ((ev[1] - ev[0]) / (ev[2] + 1e-10));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "planarity"; }
};

struct Sphericity : Attribute
{
  std::vector<double> attrib;

  Sphericity (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back (ev[0] / (ev[2] + 1e-10));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "sphericity"; }
};

struct Omnivariance : Attribute
{
  std::vector<double> attrib;

  Omnivariance (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back (std::pow (std::fabs(ev[0] * ev[1] * ev[2]), 0.333333333));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "omnivariance"; }
};

struct Anisotropy : Attribute
{
  std::vector<double> attrib;

  Anisotropy (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back ((ev[2] - ev[0]) / (ev[2] + 1e-10));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "anisotropy"; }
};

struct Eigentropy : Attribute
{
  std::vector<double> attrib;

  Eigentropy (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back (- ev[0] * std::log(std::fabs(ev[0]) + 1e-10)
                          - ev[1] * std::log(std::fabs(ev[1]) + 1e-10)
                          - ev[2] * std::log(std::fabs(ev[2]) + 1e-10));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "eigentropy"; }
};

struct SumEigen : Attribute
{
  std::vector<double> attrib;

  SumEigen (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back (ev[0] + ev[1] + ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "sum_eigen"; }
};

struct SurfaceVariation : Attribute
{
  std::vector<double> attrib;

  SurfaceVariation (Iterator begin, Iterator end, Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        attrib.push_back (ev[0] / (ev[0] + ev[1] + ev[2] + 1e-10));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "surface_variation"; }
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

  double grid_resolution = 1.544;
  double radius_neighbors = 7.781;
  double radius_dtm = 60;

  std::cerr << "Computing useful structures" << std::endl;
  Iso_cuboid_3 bbox = CGAL::bounding_box (pts.begin(), pts.end());
  Planimetric_grid grid (pts.begin(), pts.end(), Pmap(), bbox, grid_resolution);
  Neighborhood neighborhood (pts.begin(), pts.end(), Pmap());
  Local_eigen_analysis eigen (pts.begin(), pts.end(), Pmap(), neighborhood, radius_neighbors);
  
  std::cerr << "Computing attributes" << std::endl;
  Dispersion disp (pts.begin(), pts.end(), Pmap(), grid,
                   grid_resolution,
                   radius_neighbors);
  Elevation elev (pts.begin(), pts.end(), Pmap(), bbox, grid,
                  grid_resolution,
                  radius_neighbors,
                  radius_dtm);
  Verticality verti (pts.begin(), pts.end(), eigen);
  Distance_to_plane d2p (pts.begin(), pts.end(), Pmap(), eigen);

  Classification psc (pts.begin (), pts.end(), Pmap());
  
  psc.add_segmentation_attribute (&disp);
  psc.add_segmentation_attribute (&elev);
  psc.add_segmentation_attribute (&verti);
  psc.add_segmentation_attribute (&d2p);
  psc.add_segmentation_attribute (new Linearity (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new Planarity (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new Sphericity (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new Omnivariance (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new Anisotropy (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new Eigentropy (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new SumEigen (pts.begin(), pts.end(), eigen));
  psc.add_segmentation_attribute (new SurfaceVariation (pts.begin(), pts.end(), eigen));

  // Add types to PSC
  CGAL::Classification_type ground ("ground");
  CGAL::Classification_type vege ("vegetation");
  CGAL::Classification_type roof ("roof");
  CGAL::Classification_type facade ("facade");

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
  
  std::cerr << "Training" << std::endl;
  psc.training(neighborhood);

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
    << "property int label" << std::endl
    << "end_header" << std::endl;
  
  for (std::size_t i = 0; i < pts.size(); ++ i)
    {
      f << pts[i] << " ";
      
      CGAL::Classification_type* type = psc.classification_type_of (i);
      if (type == &vege)
        f << "0" << std::endl;
      else if (type == &ground)
        f << "1" << std::endl;
      else if (type == &roof)
        f << "2" << std::endl;
      else if (type == &facade)
        f << "3" << std::endl;
      else
        {
          f << "-1" << std::endl;
          //          std::cerr << "Error: unknown classification type" << std::endl;
        }
    }

  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
