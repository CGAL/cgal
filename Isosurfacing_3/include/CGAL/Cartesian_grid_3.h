// Copyright (c) 2022 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julian Stahl

#ifndef CGAL_CARTESIAN_GRID_3_H
#define CGAL_CARTESIAN_GRID_3_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Image_3.h>

#include <array>
#include <type_traits>
#include <vector>

namespace CGAL {

template <class GeomTraits>
class Cartesian_grid_3 {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Vector_3 Vector;

public:
    Cartesian_grid_3(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim, const Bbox_3& bbox)
        : sizes{xdim, ydim, zdim}, bbox(bbox) {

        values.resize(xdim * ydim * zdim);
        gradients.resize(xdim * ydim * zdim);

        const FT d_x = bbox.x_span() / (xdim - 1);
        const FT d_y = bbox.y_span() / (ydim - 1);
        const FT d_z = bbox.z_span() / (zdim - 1);
        spacing = Vector(d_x, d_y, d_z);
    }

    Cartesian_grid_3(const Image_3& image) {
        from_image(image);
    }

    FT operator()(const std::array<std::size_t, 3>& idx) const {
        return values[linear_index(idx[0], idx[1], idx[2])];
    }

    FT value(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return values[linear_index(x, y, z)];
    }

    FT& value(const std::size_t x, const std::size_t y, const std::size_t z) {
        return values[linear_index(x, y, z)];
    }

    Vector gradient(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return gradients[linear_index(x, y, z)];
    }

    Vector& gradient(const std::size_t x, const std::size_t y, const std::size_t z) {
        return gradients[linear_index(x, y, z)];
    }

    std::size_t xdim() const {
        return sizes[0];
    }
    std::size_t ydim() const {
        return sizes[1];
    }
    std::size_t zdim() const {
        return sizes[2];
    }

    const Bbox_3& get_bbox() const {
        return bbox;
    }

    const Vector& get_spacing() const {
        return spacing;
    }

private:
    std::size_t linear_index(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return (z * ydim() + y) * xdim() + x;
    }

    void from_image(const Image_3& image);
    Image_3 to_image() const;

private:
    std::vector<FT> values;
    std::vector<Vector> gradients;

    std::array<std::size_t, 3> sizes;

    Bbox_3 bbox;
    Vector spacing;
};

template <typename GeomTraits>
void Cartesian_grid_3<GeomTraits>::from_image(const Image_3& image) {
    const FT max_x = image.tx() + (image.xdim() - 1) * image.vx();
    const FT max_y = image.ty() + (image.ydim() - 1) * image.vy();
    const FT max_z = image.tz() + (image.zdim() - 1) * image.vz();
    bbox = Bbox_3(image.tx(), image.ty(), image.tz(), max_x, max_y, max_z);

    spacing = Vector(image.vx(), image.vy(), image.vz());

    sizes[0] = image.xdim();
    sizes[1] = image.ydim();
    sizes[2] = image.zdim();

    values.resize(xdim() * ydim() * zdim());
    gradients.resize(xdim() * ydim() * zdim());

    for (std::size_t x = 0; x < sizes[0]; x++) {
        for (std::size_t y = 0; y < sizes[1]; y++) {
            for (std::size_t z = 0; z < sizes[2]; z++) {

                value(x, y, z) = image.value(x, y, z);
            }
        }
    }
}

template <typename GeomTraits>
Image_3 Cartesian_grid_3<GeomTraits>::to_image() const {

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

    const double vx = bbox.x_span() / (xdim() - 1);
    const double vy = bbox.y_span() / (ydim() - 1);
    const double vz = bbox.z_span() / (zdim() - 1);

    _image* im = _createImage(xdim(), ydim(), zdim(),
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

    FT* data = (FT*)im->data;
    for (std::size_t x = 0; x < xdim(); x++) {
        for (std::size_t y = 0; y < ydim(); y++) {
            for (std::size_t z = 0; z < zdim(); z++) {

                data[(z * ydim() + y) * xdim() + x] = value(x, y, z);
            }
        }
    }

    return Image_3(im, Image_3::OWN_THE_DATA);
}

}  // end namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_3_H
