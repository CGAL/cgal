// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. Eigen itself is part of the KDE project.
//
// Copyright (C) 2010 Benoit Jacob <jacob.benoit.1@gmail.com>
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

#include "main.h"

void test_eigen2_bug_132() {
  int size = 100;
  MatrixXd A(size, size);
  VectorXd b(size), c(size);
  {
    VectorXd y = A.transpose() * (b-c);  // bug 132: infinite recursion in coeffRef
    VectorXd z = (b-c).transpose() * A;  // bug 132: infinite recursion in coeffRef
  }

  // the following ones weren't failing, but let's include them for completeness:
  {
    VectorXd y = A * (b-c);
    VectorXd z = (b-c).transpose() * A.transpose();
  }
}
