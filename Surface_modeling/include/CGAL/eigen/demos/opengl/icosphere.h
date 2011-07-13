// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_ICOSPHERE_H
#define EIGEN_ICOSPHERE_H

#include <Eigen/Core>
#include <vector>

class IcoSphere
{
  public:
    IcoSphere(unsigned int levels=1);
    const std::vector<Eigen::Vector3f>& vertices() const { return mVertices; }
    const std::vector<int>& indices(int level) const;
    void draw(int level);
  protected:
    void _subdivide();
    std::vector<Eigen::Vector3f> mVertices;
    std::vector<std::vector<int>*> mIndices;
    std::vector<int> mListIds;
};

#endif // EIGEN_ICOSPHERE_H
