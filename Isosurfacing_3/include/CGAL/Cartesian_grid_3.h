#ifndef CGAL_CARTESIAN_GRID_3_H
#define CGAL_CARTESIAN_GRID_3_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Image_3.h>

#include <type_traits>

namespace CGAL {

// TODO: not sure if anything other than float works
template <class Traits>
class Cartesian_grid_3 {
public:
    typedef typename Traits::FT FT;

public:
    Cartesian_grid_3(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim, const Bbox_3 &bbox) {
        create_image(xdim, ydim, zdim, bbox);
    }

    Cartesian_grid_3(const Image_3 &image) : grid(image) {}

    FT value(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return grid.value(x, y, z);
    }

    FT &value(const std::size_t x, const std::size_t y, const std::size_t z) {
        FT *data = (FT *)grid.image()->data;
        return data[(z * ydim() + y) * xdim() + x];
    }

    std::size_t xdim() const {
        return grid.xdim();
    }
    std::size_t ydim() const {
        return grid.ydim();
    }
    std::size_t zdim() const {
        return grid.zdim();
    }

    // TODO: better return types
    double voxel_x() const {
        return grid.vx();
    }
    double voxel_y() const {
        return grid.vy();
    }
    double voxel_z() const {
        return grid.vz();
    }

    float offset_x() const {
        return grid.tx();
    }
    float offset_y() const {
        return grid.ty();
    }
    float offset_z() const {
        return grid.tz();
    }

private:
    void create_image(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim, const Bbox_3 &bbox);

private:
    Image_3 grid;
};

template <typename T>
void Cartesian_grid_3<T>::create_image(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim,
                                       const Bbox_3 &bbox) {

    WORD_KIND wordkind;
    if (std::is_floating_point<FT>::value)
        wordkind = WK_FLOAT;
    else
        wordkind = WK_FIXED;

    SIGN sign;
    if (std::is_signed<FT>::value)
        sign = SGN_SIGNED;
    else
        sign = SGN_UNSIGNED;

    const double vx = bbox.x_span() / xdim;
    const double vy = bbox.y_span() / ydim;
    const double vz = bbox.z_span() / zdim;

    _image *im = _createImage(xdim, ydim, zdim,
                              1,           // vectorial dimension
                              vx, vy, vz,  // voxel size
                              sizeof(FT),  // image word size in bytes
                              wordkind,    // image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
                              sign);       // image word sign

    if (im == nullptr || im->data == nullptr) {
        throw std::bad_alloc();  // TODO: idk?
    }

    im->tx = bbox.xmin();
    im->ty = bbox.ymin();
    im->tz = bbox.zmin();

    grid = Image_3(im, Image_3::OWN_THE_DATA);
}

}  // end namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_3_H