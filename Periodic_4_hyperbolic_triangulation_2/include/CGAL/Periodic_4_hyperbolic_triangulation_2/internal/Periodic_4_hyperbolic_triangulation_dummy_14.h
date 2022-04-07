// Copyright (c) 2016-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Iordan Iordanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/Aff_transformation_2.h>

namespace CGAL {

template < class GT, class TDS >
inline std::vector<typename Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>::Vertex_handle >
Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>::
insert_dummy_points(bool rational)
{
  typedef typename GT::FT FT;

  typedef Aff_transformation_2<GT> Transformation;
  Transformation rotate(ROTATION, FT(1)/CGAL::sqrt(FT(2)), FT(1)/CGAL::sqrt(FT(2)));
  Transformation rot8(ROTATION, -CGAL::sqrt(FT(2)-CGAL::sqrt(FT(2)))/FT(2), CGAL::sqrt(FT(2)+CGAL::sqrt(FT(2)))/FT(2));
  Transformation opposite(ROTATION, FT(0), FT(-1)); // sin(theta), cos(theta)

  int fcount = 32;    // Faces count
  int vcount = 14;    // Vertices count

  dummy_points.clear();

  if(rational)
  {
    // Push back the origin
    dummy_points.push_back(Dummy_point(FT(0), FT(0)));

    // Compute rational approximation of internal points
    dummy_points.push_back(Dummy_point( FT(1)/FT(2),  FT(-4)/FT(19))); // 0
    dummy_points.push_back(Dummy_point( FT(1)/FT(2),  FT(4)/FT(19))); // 1
    dummy_points.push_back(Dummy_point( FT(4)/FT(19), FT(1)/FT(2) )); // 2
    dummy_points.push_back(Dummy_point( FT(-4)/FT(19), FT(1)/FT(2) )); // 3
    dummy_points.push_back(Dummy_point( FT(-1)/FT(2),  FT(4)/FT(19))); // 4
    dummy_points.push_back(Dummy_point( FT(-1)/FT(2),  FT(-4)/FT(19))); // 5
    dummy_points.push_back(Dummy_point( FT(-4)/FT(19), FT(-1)/FT(2) )); // 6
    dummy_points.push_back(Dummy_point( FT(4)/FT(19), FT(-1)/FT(2) )); // 7

    // Compute rational approximations of the midpoints
    dummy_points.push_back(Dummy_point(FT(-9)/FT(14),  FT(0)   ));      // 4
    dummy_points.push_back(Dummy_point(FT(-5)/FT(11),  FT(-5)/FT(11))); // 5
    dummy_points.push_back(Dummy_point(FT(0),          FT(-9)/FT(14))); // 6
    dummy_points.push_back(Dummy_point(FT(5)/FT(11),   FT(-5)/FT(11))); // 7

    // The vertex v_0
    dummy_points.push_back(Dummy_point(FT(97)/FT(125), FT(-26)/FT(81) )); // 0
  }
  else
  { // Algebraic dummy points
    FT F0 = FT(0);
    FT F1 = FT(1);
    FT F2 = FT(2);

    dummy_points.push_back(Dummy_point(FT(0), FT(0)));   // origin

    FT i1(CGAL::sqrt(FT(2) - sqrt(FT(2))));
    FT i2(FT(2)*i1);
    FT i3(FT(2)*CGAL::sqrt(FT(2)) - i2);
    FT i4(-1 + i3);
    Point a0(FT(CGAL::sqrt(i4)), FT(0));          // internal point
    a0 = rot8(a0);
    dummy_points.push_back(Dummy_point(a0));
    for(int k=1; k<8; ++k)
    {
      a0 = rotate(a0);
      dummy_points.push_back(Dummy_point(a0));
    }

    dummy_points.push_back(Dummy_point(FT(-CGAL::sqrt(CGAL::sqrt(F2) - F1)), F0)); //  μ(s_0)
    dummy_points.push_back(Dummy_point(FT(-CGAL::sqrt((CGAL::sqrt(F2) - F1) / F2)),
                                       FT(-CGAL::sqrt((CGAL::sqrt(F2) - F1) / F2))));     //  μ(s_1)
    dummy_points.push_back(Dummy_point(F0, FT(-CGAL::sqrt(CGAL::sqrt(F2) - F1)))); //  μ(s_2)
    dummy_points.push_back(Dummy_point(FT(CGAL::sqrt((CGAL::sqrt(F2) - F1) / F2)),
                                       -FT(CGAL::sqrt((CGAL::sqrt(F2) - F1) / F2)))); //  μ(s_3)
    dummy_points.push_back(Dummy_point(FT(CGAL::sqrt(CGAL::sqrt(F2) + F1)/ F2),
                                       - FT(CGAL::sqrt(CGAL::sqrt(F2) - F1)/ F2))); //  v_0
  }

  Hyperbolic_translation off[32][3];
  for(int i=0; i<32; ++i)
    for(int j=0; j<3; ++j)
      off[i][j] = Hyperbolic_translation();

  int tri[32][3];
  for(int i=0; i<4; ++i)
  {
    Hyperbolic_translation vo1, vo2;
    if(i == 0) { vo1 = Hyperbolic_translation(1, 6, 3);    vo2 = Hyperbolic_translation();         }
    if(i == 1) { vo1 = Hyperbolic_translation(1, 4);       vo2 = Hyperbolic_translation(1, 6, 3);  }
    if(i == 2) { vo1 = Hyperbolic_translation(3);          vo2 = Hyperbolic_translation(1, 4);     }
    if(i == 3) { vo1 = Hyperbolic_translation(4, 1, 6, 3); vo2 = Hyperbolic_translation(3);        }

    int i0 = 0;
    int fn = i;
    int i1 = (fn+1);
    int i2 = (fn+2);
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    //----------------------
    fn = i + 4;
    i1 = ((i+4) % 8) + 1;
    i2 = ((i+5) % 8) + 1;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    //----------------------
    fn = i + 8;
    i0 = i+9;
    i1 = i+2;
    i2 = i+1;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    off[fn][0] = Hyperbolic_translation(i);
    //----------------------
    fn = i + 12;
    i0 = i+9;
    i1 = ((i+5)%8)+1;
    i2 = ((i+4)%8)+1;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    //off[fn][0] = Hyperbolic_translation(i+4);
    //----------------------
    fn = 2*(i + 8);
    i0 = i+1;
    i1 = ((i+5)%8)+1;
    i2 = i+9;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    off[fn][1] = Hyperbolic_translation(i);
    off[fn][2] = Hyperbolic_translation(i);
    //----------------------
    fn = 2*(i + 8) + 1;
    i0 = i+2;
    i1 = i+9;
    i2 = i+5;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    off[fn][1] = Hyperbolic_translation(i);
    off[fn][2] = Hyperbolic_translation(i);
    //----------------------
    fn = 2*(i + 12);
    i0 = ((i+5)%8)+1;
    i1 = i+1;
    i2 = 13;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    off[fn][0] = Hyperbolic_translation(i);
    off[fn][2] = vo2;
    //----------------------
    fn = 2*(i + 12) + 1;
    i0 = i+5;
    i1 = 13;
    i2 = i+2;
    tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
    off[fn][0] = Hyperbolic_translation(i);
    off[fn][1] = vo1;
  }

  Vertex_handle vertices[14];
  for(int i=0; i<vcount; ++i)
  {
    vertices[i] = tds().create_vertex();
    vertices[i]->set_point(dummy_points[i]());
    dummy_points[i].set_vertex(vertices[i]);
  }

  Face_handle faces[32];
  for(int i=0; i<fcount; ++i)
  {
    faces[i] = tds().create_face(vertices[tri[i][0]], vertices[tri[i][1]], vertices[tri[i][2]]);
    for(int j=0; j<3; ++j)
      faces[i]->set_translation(j, off[i][j]);
  }

  tds().set_dimension(2);

  for(int i=0; i<4; ++i)
  {
    int fn = i;
    int n0 = i+8;
    int n1 = i+1;
    tds().set_adjacency(faces[fn], 0, faces[n0], 0);
    tds().set_adjacency(faces[fn], 1, faces[n1], 2);

    fn = i+4;
    n0 = i+12;
    n1 = (i+5) % 8;
    tds().set_adjacency(faces[fn], 0, faces[n0], 0);
    tds().set_adjacency(faces[fn], 1, faces[n1], 2);

    fn = i+8;
    n0 = 2*(i+8);
    n1 = 2*(i+8)+1;
    tds().set_adjacency(faces[fn], 1, faces[n0], 1);
    tds().set_adjacency(faces[fn], 2, faces[n1], 2);

    fn = 2*(i+8);
    n0 = (i+12);
    n1 = 2*n0;
    tds().set_adjacency(faces[fn], 0, faces[n0], 2);
    tds().set_adjacency(faces[fn], 2, faces[n1], 2);

    fn = 2*(i+8) + 1;
    n0 = (i+12);
    n1 = 2*n0 + 1;
    tds().set_adjacency(faces[fn], 0, faces[n0], 1);
    tds().set_adjacency(faces[fn], 1, faces[n1], 1);


    if(i < 3)
    {
      fn = 2*(i+12);
      n0 = 24 + ((i+1)*2 + 1);
      tds().set_adjacency(faces[fn], 1, faces[n0], 2);

      fn = 2*(i+12)+1;
      n0 = 25 + (i*2 + 1);
      tds().set_adjacency(faces[fn], 0, faces[n0], 0);
    }
  }

  int fn = 24;
  int n0 = 30;
  tds().set_adjacency(faces[fn], 0, faces[n0], 1);

  fn = 25;
  n0 = 31;
  tds().set_adjacency(faces[fn], 2, faces[n0], 0);

  vertices[0]->set_face(faces[0]);
  // vertices 1-8
  for(int i = 1; i < 9; ++i)
    vertices[i]->set_face(faces[i+7]);

  // vertices 9-12
  for(int i = 9; i < 13; ++i)
    vertices[i]->set_face(faces[i+3]);

  vertices[13]->set_face(faces[30]);

  std::vector<Vertex_handle> ret(vcount);
  for(int i=0; i<vcount; ++i)
    ret.push_back(vertices[i]);

  for(Face_iterator fit = tds().faces_begin(); fit != tds().faces_end(); ++fit)
    this->make_canonical(fit);

  CGAL_triangulation_assertion(is_valid(true));

  return ret;
}


} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H

