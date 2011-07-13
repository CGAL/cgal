// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
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

#include "binary_library.h"
#include "stdio.h"

void demo_MatrixXd()
{
  struct C_MatrixXd *matrix1, *matrix2, *result;
  printf("*** demo_MatrixXd ***\n");
  
  matrix1 = MatrixXd_new(3, 3);
  MatrixXd_set_zero(matrix1);
  MatrixXd_set_coeff(matrix1, 0, 1, 2.5);
  MatrixXd_set_coeff(matrix1, 1, 0, 1.4);
  printf("Here is matrix1:\n");
  MatrixXd_print(matrix1);

  matrix2 = MatrixXd_new(3, 3);
  MatrixXd_multiply(matrix1, matrix1, matrix2);
  printf("Here is matrix1*matrix1:\n");
  MatrixXd_print(matrix2);

  MatrixXd_delete(matrix1);
  MatrixXd_delete(matrix2);
}

// this helper function takes a plain C array and prints it in one line
void print_array(double *array, int n)
{
  struct C_Map_MatrixXd *m = Map_MatrixXd_new(array, 1, n);
  Map_MatrixXd_print(m);
  Map_MatrixXd_delete(m);
}

void demo_Map_MatrixXd()
{
  struct C_Map_MatrixXd *map;
  double array[5];
  int i;
  printf("*** demo_Map_MatrixXd ***\n");
  
  for(i = 0; i < 5; ++i) array[i] = i;
  printf("Initially, the array is:\n");
  print_array(array, 5);
  
  map = Map_MatrixXd_new(array, 5, 1);
  Map_MatrixXd_add(map, map, map);
  Map_MatrixXd_delete(map);

  printf("Now the array is:\n");
  print_array(array, 5);
}

int main()
{
  demo_MatrixXd();
  demo_Map_MatrixXd();
}
