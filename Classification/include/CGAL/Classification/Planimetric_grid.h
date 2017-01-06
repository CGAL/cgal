#ifndef CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
#define CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H

#include <vector>

#include <CGAL/Classification/Image.h>

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationDataStructures

    \brief Class that precomputes a 2D planimetric grid.

    The grid is composed of squared cells with a user-defined size,
    each cell containing the list of indices of the points whose
    projection along the Z-axis lies within this cell. The mapping
    from each point to the cell it lies in is also stored.

    \tparam Kernel model of \cgal Kernel.
    \tparam Range range of items, model of `ConstRange`. Its iterator type
    is `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `Range` and value type
    is `Point_3<Kernel>`.
  */

template <typename Kernel, typename Range, typename PointMap>
class Planimetric_grid
{
public:
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;

private:
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

    \param input input range.
    \param point_map property map to access the input points.
    \param bbox bounding box of the input range.
    \param grid_resolution resolution of the planimetric grid.
  */
  Planimetric_grid (const Range& input,
                    PointMap point_map,
                    const Iso_cuboid_3& bbox,
                    const double grid_resolution)
  {
    std::size_t width = (std::size_t)((bbox.xmax() - bbox.xmin()) / grid_resolution) + 1;
    std::size_t height = (std::size_t)((bbox.ymax() - bbox.ymin()) / grid_resolution) + 1;
    m_grid = Image_indices (width, height);

    for (std::size_t i = 0; i < input.size(); ++ i)
      {
        const Point_3& p = get(point_map, *(input.begin()+i));
        m_x.push_back ((std::size_t)((p.x() - bbox.xmin()) / grid_resolution));
        m_y.push_back ((std::size_t)((p.y() - bbox.ymin()) / grid_resolution));
        m_grid(m_x.back(), m_y.back()).push_back (i);
      }

  }

  /*!
    \brief Returns the number of cells along the X-axis.
  */
  std::size_t width() const { return m_grid.width(); }
  /*!
    \brief Returns the number of cells along the Y-axis.
  */
  std::size_t height() const { return m_grid.height(); }

  /*!
    \brief Returns the indices of the points lying in the cell at position `(x,y)`.
  */
  const std::vector<std::size_t>& indices(std::size_t x, std::size_t y) const { return m_grid(x,y); }
  
  /*!
    \brief Returns `false` if the cell at position `(x,y)` is empty, `true` otherwise.
  */
  bool mask(std::size_t x, std::size_t y) const { return (!(m_grid(x,y).empty())); }

  /*!
    \brief Returns the `x` grid coordinate of the point at position `index`.
  */
  std::size_t x(std::size_t index) const { return m_x[index]; }
  /*!
    \brief Returns the `y` grid coordinate of the point at position `index`.
  */
  std::size_t y(std::size_t index) const { return m_y[index]; }
};
  

}
  
}


#endif // CGAL_CLASSIFICATION_PLANIMETRIC_GRID_H
