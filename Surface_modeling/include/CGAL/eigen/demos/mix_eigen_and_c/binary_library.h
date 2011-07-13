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

// This is a pure C header, no C++ here.
// The functions declared here will be implemented in C++ but
// we don't have to know, because thanks to the extern "C" syntax,
// they will be compiled to C object code.

#ifdef __cplusplus
extern "C"
{
#endif

  // just dummy empty structs to give different pointer types,
  // instead of using void* which would be type unsafe
  struct C_MatrixXd {};
  struct C_Map_MatrixXd {};

  // the C_MatrixXd class, wraps some of the functionality
  // of Eigen::MatrixXd.
  struct C_MatrixXd* MatrixXd_new(int rows, int cols);
  void    MatrixXd_delete     (struct C_MatrixXd *m);
  double* MatrixXd_data       (struct C_MatrixXd *m);
  void    MatrixXd_set_zero   (struct C_MatrixXd *m);
  void    MatrixXd_resize     (struct C_MatrixXd *m, int rows, int cols);
  void    MatrixXd_copy       (struct C_MatrixXd *dst,
                               const struct C_MatrixXd *src);
  void    MatrixXd_copy_map   (struct C_MatrixXd *dst,
                               const struct C_Map_MatrixXd *src);  
  void    MatrixXd_set_coeff  (struct C_MatrixXd *m,
                               int i, int j, double coeff);
  double  MatrixXd_get_coeff  (const struct C_MatrixXd *m,
                               int i, int j);
  void    MatrixXd_print      (const struct C_MatrixXd *m);
  void    MatrixXd_add        (const struct C_MatrixXd *m1,
                               const struct C_MatrixXd *m2,
                               struct C_MatrixXd *result);  
  void    MatrixXd_multiply   (const struct C_MatrixXd *m1,
                               const struct C_MatrixXd *m2,
                               struct C_MatrixXd *result);
  
  // the C_Map_MatrixXd class, wraps some of the functionality
  // of Eigen::Map<MatrixXd>
  struct C_Map_MatrixXd* Map_MatrixXd_new(double *array, int rows, int cols);
  void   Map_MatrixXd_delete     (struct C_Map_MatrixXd *m);
  void   Map_MatrixXd_set_zero   (struct C_Map_MatrixXd *m);
  void   Map_MatrixXd_copy       (struct C_Map_MatrixXd *dst,
                                  const struct C_Map_MatrixXd *src);
  void   Map_MatrixXd_copy_matrix(struct C_Map_MatrixXd *dst,
                                  const struct C_MatrixXd *src);  
  void   Map_MatrixXd_set_coeff  (struct C_Map_MatrixXd *m,
                                  int i, int j, double coeff);
  double Map_MatrixXd_get_coeff  (const struct C_Map_MatrixXd *m,
                                  int i, int j);
  void   Map_MatrixXd_print      (const struct C_Map_MatrixXd *m);
  void   Map_MatrixXd_add        (const struct C_Map_MatrixXd *m1,
                                  const struct C_Map_MatrixXd *m2,
                                  struct C_Map_MatrixXd *result);  
  void   Map_MatrixXd_multiply   (const struct C_Map_MatrixXd *m1,
                                  const struct C_Map_MatrixXd *m2,
                                  struct C_Map_MatrixXd *result);

#ifdef __cplusplus
} // end extern "C"
#endif