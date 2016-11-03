#ifndef CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
#define CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H

#include <vector>

#include <CGAL/Classification/Image.h>

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassification

    \brief Class that precomputes a 2D planimetric grid used for
    digital terrain modeling.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
  */

template <typename Kernel, typename RandomAccessIterator, typename PointMap>
class Planimetric_grid
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;
  
  typedef Image<std::vector<std::size_t> > Image_indices;
  typedef Image<bool> Image_bool;

  Image_indices m_grid;
  std::vector<std::size_t> m_x;
  std::vector<std::size_t> m_y;
  
public:

  /// \cond SKIP_IN_MANUAL
  Planimetric_grid () { }
  /// \endcond

  /*
    \brief Constructs a planimetric grid based on the input range.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param bbox The bounding box of the input range
    \param grid_resolution Resolution of the planimetric grid
  */
  Planimetric_grid (const RandomAccessIterator& begin,
                    const RandomAccessIterator& end,
                    PointMap point_map,
                    const Iso_cuboid_3& bbox,
                    const FT grid_resolution)
  {

    std::size_t size = (std::size_t)(end - begin);
    std::size_t width = (std::size_t)((bbox.xmax() - bbox.xmin()) / grid_resolution) + 1;
    std::size_t height = (std::size_t)((bbox.ymax() - bbox.ymin()) / grid_resolution) + 1;
    m_grid = Image_indices (width, height);

    for (std::size_t i = 0; i < size; ++ i)
      {
        const Point_3& p = get(point_map, begin[i]);
        m_x.push_back ((std::size_t)((p.x() - bbox.xmin()) / grid_resolution));
        m_y.push_back ((std::size_t)((p.y() - bbox.ymin()) / grid_resolution));
        m_grid(m_x.back(), m_y.back()).push_back (i);
      }

  }

  std::size_t width() const { return m_grid.width(); }
  std::size_t height() const { return m_grid.height(); }

  /*!
    \brief Returns the indices of points lying in the given indexed cell.
  */
  const std::vector<std::size_t>& indices(std::size_t x, std::size_t y) const { return m_grid(x,y); }
  /*!
    \brief Returns `true` if the indexed cell is to be used for classification.
  */
  bool mask(std::size_t x, std::size_t y) const { return (!(m_grid(x,y).empty())); }

  /*!
    \brief Returns the `x` coordinate of the indexed point in the grid.
  */
  std::size_t x(std::size_t index) const { return m_x[index]; }
  /*!
    \brief Returns the `y` coordinate of the indexed point in the grid.
  */
  std::size_t y(std::size_t index) const { return m_y[index]; }
};
  

}
  
}


#endif // CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
