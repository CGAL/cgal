// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of LCC-Lab.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_JACOBIAN_H
#define LCC_JACOBIAN_H

#include <array>
#include <utility>
#include <CGAL/MP_Float.h>
#include <CGAL/Kernel/global_functions.h>

// #include <tools/Orientation.h>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
auto normalize(T & V)
{
  auto const slen=V.squared_length();
  auto const d=CGAL::approximate_sqrt(slen);
  V=V/d;
}
////////////////////////////////////////////////////////////////////////////////
template<typename Vector>
double scale_jacobian(Vector v1, Vector v2, Vector v3)
{
  double res;
  normalize(v1);
  normalize(v2);
  normalize(v3);
  res=CGAL::determinant(v1, v2, v3);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
/*       4----7
 *      /|   /|
 *     5----6 |
 *     | 3--|-2
 *     |/   |/
 *     0----1
 */
template<typename LCC>
double scale_jacobian_of_hexa(LCC& lcc,
                              std::array<typename LCC::Vertex_attribute_descriptor, 8>&
                              vertices)
{
//  // Compute scaled jacobi of a triangulated hexa given 8 vertices
//  if(!is_positive_tetra(lcc.point_of_vertex_attribute(vertices[0]),
//                        lcc.point_of_vertex_attribute(vertices[1]),
//                        lcc.point_of_vertex_attribute(vertices[4]),
//                        lcc.point_of_vertex_attribute(vertices[5])))
//  {
//    std::swap(vertices[1], vertices[3]);
//    std::swap(vertices[4], vertices[6]);
//  }

  using Vector=typename LCC::Vector;
  Vector v0(lcc.point_of_vertex_attribute(vertices[0]),
      lcc.point_of_vertex_attribute(vertices[1]));
  Vector v1(lcc.point_of_vertex_attribute(vertices[0]),
      lcc.point_of_vertex_attribute(vertices[3]));
  Vector v2(lcc.point_of_vertex_attribute(vertices[0]),
      lcc.point_of_vertex_attribute(vertices[5]));
  double res=scale_jacobian(v0, v1, v2);

  v0=Vector(lcc.point_of_vertex_attribute(vertices[3]),
      lcc.point_of_vertex_attribute(vertices[0]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[3]),
      lcc.point_of_vertex_attribute(vertices[2]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[3]),
      lcc.point_of_vertex_attribute(vertices[4]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

  v0=Vector(lcc.point_of_vertex_attribute(vertices[2]),
      lcc.point_of_vertex_attribute(vertices[3]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[2]),
      lcc.point_of_vertex_attribute(vertices[1]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[2]),
      lcc.point_of_vertex_attribute(vertices[7]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

  v0=Vector(lcc.point_of_vertex_attribute(vertices[1]),
      lcc.point_of_vertex_attribute(vertices[2]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[1]),
      lcc.point_of_vertex_attribute(vertices[0]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[1]),
      lcc.point_of_vertex_attribute(vertices[6]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

  v0=Vector(lcc.point_of_vertex_attribute(vertices[4]),
      lcc.point_of_vertex_attribute(vertices[7]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[4]),
      lcc.point_of_vertex_attribute(vertices[5]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[4]),
      lcc.point_of_vertex_attribute(vertices[3]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

  v0=Vector(lcc.point_of_vertex_attribute(vertices[5]),
      lcc.point_of_vertex_attribute(vertices[4]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[5]),
      lcc.point_of_vertex_attribute(vertices[6]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[5]),
      lcc.point_of_vertex_attribute(vertices[0]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

  v0=Vector(lcc.point_of_vertex_attribute(vertices[6]),
      lcc.point_of_vertex_attribute(vertices[5]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[6]),
      lcc.point_of_vertex_attribute(vertices[7]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[6]),
      lcc.point_of_vertex_attribute(vertices[1]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

  v0=Vector(lcc.point_of_vertex_attribute(vertices[7]),
      lcc.point_of_vertex_attribute(vertices[6]));
  v1=Vector(lcc.point_of_vertex_attribute(vertices[7]),
      lcc.point_of_vertex_attribute(vertices[4]));
  v2=Vector(lcc.point_of_vertex_attribute(vertices[7]),
      lcc.point_of_vertex_attribute(vertices[2]));
  res=std::min(res, scale_jacobian(v0, v1, v2));

 /* for(int i=0; i<8; ++i)
  { std::cout<<lcc.point_of_vertex_attribute(vertices[i])<<"  "; }
  // std::cout<<" -> "<<res<<std::endl;
  // CGAL::draw(temp);
  std::cout<<"VOL: "<<volume_of_hexa
             (lcc.point_of_vertex_attribute(vertices[0]),
      lcc.point_of_vertex_attribute(vertices[1]),
      lcc.point_of_vertex_attribute(vertices[2]),
      lcc.point_of_vertex_attribute(vertices[3]),
      lcc.point_of_vertex_attribute(vertices[5]),
      lcc.point_of_vertex_attribute(vertices[6]),
      lcc.point_of_vertex_attribute(vertices[7]),
      lcc.point_of_vertex_attribute(vertices[4]))<<" JABOBI: "<<res<<std::endl;

  if(volume_of_hexa
     (lcc.point_of_vertex_attribute(vertices[0]),
lcc.point_of_vertex_attribute(vertices[1]),
lcc.point_of_vertex_attribute(vertices[2]),
lcc.point_of_vertex_attribute(vertices[3]),
lcc.point_of_vertex_attribute(vertices[5]),
lcc.point_of_vertex_attribute(vertices[6]),
lcc.point_of_vertex_attribute(vertices[7]),
lcc.point_of_vertex_attribute(vertices[4]))==1)
  {
  LCC3 temp;
  temp.make_hexahedron(lcc.point_of_vertex_attribute(vertices[0]),
                       lcc.point_of_vertex_attribute(vertices[1]),
                       lcc.point_of_vertex_attribute(vertices[2]),
                       lcc.point_of_vertex_attribute(vertices[3]),
                       lcc.point_of_vertex_attribute(vertices[4]),
                       lcc.point_of_vertex_attribute(vertices[5]),
                       lcc.point_of_vertex_attribute(vertices[6]),
                       lcc.point_of_vertex_attribute(vertices[7]));
  CGAL::draw(temp);
  } */

  return res;
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double jacobian_tri_hex_for_dart_v1(LCC& lcc, typename LCC::Dart_descriptor d)
{
  // Get the vertices of virtual hexahedra, considering virtually merging
  // face(d) and face(beta2(d)), and face(beta1,2(d)) and face(beta1,2,1,2(d)).
  /*       4----7
   *      /|   /|
   *     5----6 |
   *     | 3--|-2
   *     |/   |/
   *     0----1
   */
  std::array<typename LCC::Vertex_attribute_descriptor, 8> vertices;
  vertices[0]=lcc.vertex_attribute(lcc.beta(d, 0));
  vertices[1]=lcc.vertex_attribute(d);
  vertices[2]=lcc.vertex_attribute(lcc.beta(d, 2, 0));
  vertices[3]=lcc.vertex_attribute(lcc.beta(d, 1));

  vertices[4]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 1, 2, 0));
  vertices[5]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0));
  vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0));
  if(vertices[6]==vertices[5])
  { vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0, 2, 0)); }
  vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0));
  if(vertices[7]==vertices[6])
  { vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0, 2, 0)); }

  return scale_jacobian_of_hexa(lcc, vertices);
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double jacobian_tri_hex_for_dart_v2(LCC& lcc, typename LCC::Dart_descriptor d)
{
  // Get the vertices of virtual hexahedra, considering virtually merging
  // face(d) and face(beta2(d)), and face(beta1,2(d)) and face(beta1,2,0,2(d)).
  /*       4----7
   *      /|   /|
   *     5----6 |
   *     | 3--|-2
   *     |/   |/
   *     0----1
   */
  std::array<typename LCC::Vertex_attribute_descriptor, 8> vertices;
  vertices[0]=lcc.vertex_attribute(lcc.beta(d, 0));
  vertices[1]=lcc.vertex_attribute(d);
  vertices[2]=lcc.vertex_attribute(lcc.beta(d, 2, 0));
  vertices[3]=lcc.vertex_attribute(lcc.beta(d, 1));

  vertices[4]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0));
  vertices[5]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0, 2, 0));
  vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0));
  if(vertices[6]==vertices[5])
  { vertices[6]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 0, 2, 0)); }
  vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0));
  if(vertices[7]==vertices[6])
  { vertices[7]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 2, 0, 2, 0)); }

  return scale_jacobian_of_hexa(lcc, vertices);
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double jacobian_of_triangulated_hexahedron(LCC& lcc, typename LCC::Dart_descriptor d)
{
  return std::max({jacobian_tri_hex_for_dart_v1(lcc, d),
                   jacobian_tri_hex_for_dart_v1(lcc, lcc.beta(d, 1)),
                   jacobian_tri_hex_for_dart_v1(lcc, lcc.beta(d, 0)),
                   jacobian_tri_hex_for_dart_v2(lcc, d),
                   jacobian_tri_hex_for_dart_v2(lcc, lcc.beta(d, 1)),
                   jacobian_tri_hex_for_dart_v2(lcc, lcc.beta(d, 0))});
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
double jacobian_of_hexahedron(LCC& lcc, typename LCC::Dart_descriptor d)
{
  assert(lcc.is_volume_combinatorial_hexahedron(d));

  std::array<typename LCC::Vertex_attribute_descriptor, 8> vertices;
  vertices[0]=lcc.vertex_attribute(d);
  vertices[1]=lcc.vertex_attribute(lcc.beta(d, 0));
  vertices[2]=lcc.vertex_attribute(lcc.beta(d, 0, 2, 1, 1));
  vertices[3]=lcc.vertex_attribute(lcc.beta(d, 2, 1, 1));

  vertices[4]=lcc.vertex_attribute(lcc.beta(d, 2, 0));
  vertices[5]=lcc.vertex_attribute(lcc.beta(d, 2));
  vertices[6]=lcc.vertex_attribute(lcc.beta(d, 1, 1));
  vertices[7]=lcc.vertex_attribute(lcc.beta(d, 1, 2, 0));
  
  return scale_jacobian_of_hexa(lcc, vertices);
}
////////////////////////////////////////////////////////////////////////////////
template<typename LCC>
void display_all_jacobians(LCC& lcc)
{
  /// Special function that work only for triangulated hexahedra
  /// The code commented out allow to triangulate the face of hex if you wand
  
  // Triangulate all squares
/*   std::vector<Dart_descriptor> faces_to_subdivide;
  faces_to_subdivide.reserve(lcc.number_of_darts()/4);
  for(auto itd=lcc.one_dart_per_cell<2>().begin(),
        itdend=lcc.one_dart_per_cell<2>().end(); itd!=itdend; ++itd)
  { faces_to_subdivide.push_back(itd); }

  for(Dart_descriptor itd: faces_to_subdivide)
  { lcc.insert_cell_1_in_cell_2(itd, lcc.beta<1,1>(itd)); }
*/

  // Display jacobian of triangulated hexahedra
  for(auto it=lcc.template one_dart_per_cell<3>().begin(),
      itend=lcc.template one_dart_per_cell<3>().end(); it!=itend; ++it)
  { std::cout<<jacobian_of_triangulated_hexahedron(lcc, it)<<std::endl; }
}
////////////////////////////////////////////////////////////////////////////////
/// Display the jacobian of all hexahedra, print other for other volumes.
template<typename LCC>
void display_all_jacobians_of_hex(LCC& lcc, const std::string& other="X")
{
  // Display jacobian of hexahedra
  for(auto it=lcc.template one_dart_per_cell<3>().begin(),
      itend=lcc.template one_dart_per_cell<3>().end(); it!=itend; ++it)
  {
    if(lcc.is_volume_combinatorial_hexahedron(it))
    { std::cout<<jacobian_of_hexahedron(lcc, it)<<std::endl; }
    else
    { std::cout<<other<<std::endl; }
  }
}

////////////////////////////////////////////////////////////////////////////////
#endif // LCC_JACOBIAN_H
