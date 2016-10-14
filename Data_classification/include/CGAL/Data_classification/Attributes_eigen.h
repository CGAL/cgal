#ifndef CGAL_DATA_CLASSIFICATION_ATTRIBUTES_EIGEN_H
#define CGAL_DATA_CLASSIFICATION_ATTRIBUTES_EIGEN_H

#include <vector>

namespace CGAL {

namespace Data_classification {
  
template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_linearity : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_linearity (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back ((ev[2] - ev[1]) / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "linearity"; }
};

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_planarity : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_planarity (RandomAccessIterator begin,
                       RandomAccessIterator end,
                       Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back ((ev[1] - ev[0]) / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "planarity"; }
};

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_sphericity : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_sphericity (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back (ev[0] / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "sphericity"; }
};

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_omnivariance : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_omnivariance (RandomAccessIterator begin,
                          RandomAccessIterator end,
                          Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
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

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_anisotropy : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_anisotropy (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back ((ev[2] - ev[0]) / ev[2]);
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "anisotropy"; }
};

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_eigentropy : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_eigentropy (RandomAccessIterator begin,
                        RandomAccessIterator end,
                        Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);
        if (ev[0] < 1e-15
            || ev[1] < 1e-15
            || ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back (- ev[0] * std::log(ev[0])
                            - ev[1] * std::log(ev[1])
                            - ev[2] * std::log(ev[2]));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "eigentropy"; }
};

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_sum_eigenvalues : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_sum_eigenvalues (RandomAccessIterator begin,
                             RandomAccessIterator end,
                             Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      attrib.push_back (eigen.sum_eigenvalues(i));

    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "sum_eigen"; }
};

template <typename Kernel, typename RandomAccessIterator, typename PointPMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Attribute_surface_variation : public Attribute
{
  typedef Data_classification::Local_eigen_analysis<Kernel, RandomAccessIterator,
                                                    PointPMap, DiagonalizeTraits> Local_eigen_analysis;
  std::vector<double> attrib;
public:
  Attribute_surface_variation (RandomAccessIterator begin,
                               RandomAccessIterator end,
                               Local_eigen_analysis& eigen)
  {
    std::size_t size = (std::size_t)(end - begin);
    attrib.reserve (size);
    for (std::size_t i = 0; i < size; ++ i)
      {
        const typename Local_eigen_analysis::Eigenvalues& ev = eigen.eigenvalue(i);

        if (ev[0] + ev[1] + ev[2] < 1e-15)
          attrib.push_back (0.);
        else
          attrib.push_back (ev[0] / (ev[0] + ev[1] + ev[2]));
      }
    this->compute_mean_max (attrib, mean, this->max);
  }
  virtual double value (std::size_t pt_index)
  {
    return attrib[pt_index];
  }
  virtual std::string id() { return "surface_variation"; }
};


} // namespace Data_classification

} // namespace CGAL

#endif // CGAL_DATA_CLASSIFICATION_ATTRIBUTES_EIGEN_H
