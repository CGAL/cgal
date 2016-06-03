// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/candidate-packages/Triangulation_2/include/CGAL/Delaunay_triangulation_2.h $
// $Id: Delaunay_triangulation_2.h 57509 2010-07-15 09:14:09Z sloriot $
//
//
// Author(s)     : Iordan Iordanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H

#include <CGAL/Aff_transformation_2.h>

namespace CGAL {

    template<class Point, class Face_handle>
    Point barycenter(Face_handle fh) {
        Point p0 = fh->offset(0).apply(fh->vertex(0)->point());
        Point p1 = fh->offset(1).apply(fh->vertex(1)->point());
        Point p2 = fh->offset(2).apply(fh->vertex(2)->point());
        return Point( (p0.x() + p1.x() + p2.x())/3, (p0.y() + p1.y() + p2.y())/3 );
    }
    
    template<class Point, class Face_handle>
    Point midpoint(Face_handle fh, int i, int j) {
        Point p0 = fh->offset(i).apply(fh->vertex(i)->point());
        Point p1 = fh->offset(j).apply(fh->vertex(j)->point());
        return Point( (p0.x() + p1.x())/2, (p0.y() + p1.y())/2 );
    }

    template<class Point, class Face_handle>
    Point vertex(Face_handle fh, int i) {
        Point p0 = fh->offset(i).apply(fh->vertex(i)->point());
        return Point( p0.x(), p0.y() );
    }

    template < class GT, class TDS >
    inline std::vector<typename Periodic_4_hyperbolic_triangulation_2<GT, TDS>::Vertex_handle >
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
    insert_dummy_points() {

        typedef typename GT::FT FT;

        typedef Aff_transformation_2<GT> Transformation; 
        Transformation rotate(ROTATION, FT(1)/CGAL::sqrt(FT(2)), FT(1)/CGAL::sqrt(FT(2))); 


        int fcount = 32;    // Faces count
        int vcount = 14;    // Vertices count

        FT F0 = FT(0);
        FT F1 = FT(1);
        FT F2 = FT(2);

        std::vector<typename GT::Point_2> pts;
        pts.push_back(Point(FT( CGAL::sqrt(CGAL::sqrt(F2) + F1)/ F2 ), - FT( CGAL::sqrt(CGAL::sqrt(F2) - F1)/ F2) ));           //  v_0
        pts.push_back(Point(FT( CGAL::sqrt(CGAL::sqrt(F2) - F1)), F0));                                                         //  μ(s_0)
        pts.push_back(Point(FT( CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2) ), FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)) ));      //  μ(s_1)
        pts.push_back(Point(F0, FT( CGAL::sqrt(CGAL::sqrt(F2) - F1))));                                                         //  μ(s_2)
        pts.push_back(Point(-FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)), FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)) ));       //  μ(s_3)

        FT num_x = FT( F2 + CGAL::sqrt(F2) - CGAL::sqrt(FT(6)) );
        FT den_x = FT( FT(4)*CGAL::sqrt(CGAL::sqrt(F2) - F1) );
        FT xi = num_x/den_x;
        FT yi = FT( (CGAL::sqrt(F2) - F2*CGAL::sqrt(FT(3)) + CGAL::sqrt(FT(6))) / (FT(4)*CGAL::sqrt(CGAL::sqrt(F2) - F1)) );

        Point a0(xi, yi);

        /* Consistency reasons (has to be corresponding to vertex 5 of the figure) */
        a0 = rotate(a0);
        a0 = rotate(a0);
        a0 = rotate(a0);

        pts.push_back(a0);
        for (int k = 1; k < 8; k++) {
            a0 = rotate(a0);
            pts.push_back(a0);
        }

        pts.push_back(Point( FT(0), FT(0) ));

        
        int tri[32][3] = { 
                {  0,   1,   9},    //  0
                {  1,  10,   9},    //  1
                {  0,  10,   1},    //  2
                {  0,   2,  10},    //  3
                {  0,  11,   2},    //  4
                {  0,   3,  11},    //  5
                {  0,  12,   3},    //  6
                {  0,   4,  12},    //  7
                {  0,   5,   4},    //  8
                {  0,   1,   5},    //  9
                {  0,   6,   1},    // 10
                {  0,   2,   6},    // 11
                {  0,   7,   2},    // 12
                {  0,   3,   7},    // 13
                {  0,   8,   3},    // 14
                {  0,   4,   8},    // 15
                {  0,   9,   4},    // 16
                {  2,  11,  10},    // 17
                {  3,  12,  11},    // 18
                {  4,   5,  12},    // 19
                {  1,   6,   5},    // 20
                {  2,   7,   6},    // 21
                {  3,   8,   7},    // 22
                {  4,   9,   8},    // 23
                {  9,  10,  13},    // 24
                { 10,  11,  13},    // 25
                { 11,  12,  13},    // 26
                {  5,  13,  12},    // 27
                {  5,   6,  13},    // 28
                {  6,   7,  13},    // 29
                {  7,   8,  13},    // 30
                {  8,   9,  13},    // 31
        };


        Offset off[32][3] = {
            { Offset(),         Offset(),   Offset() },     //  0
            { Offset(),         Offset(),   Offset() },     //  1
            { Offset(1,6,3),    Offset(),   Offset() },     //  2
            { Offset(1,6,3),    Offset(),   Offset() },     //  3
            { Offset(1,4),      Offset(),   Offset() },     //  4
            { Offset(1,4),      Offset(),   Offset() },     //  5
            { Offset(3),        Offset(),   Offset() },     //  6
            { Offset(3),        Offset(),   Offset() },     //  7
            { Offset(3,6,1,4),  Offset(),   Offset() },     //  8
            { Offset(3,6,1,4),  Offset(4),  Offset() },     //  9
            { Offset(4),        Offset(),   Offset(4)},     // 10
            { Offset(4),        Offset(5),  Offset() },     // 11
            { Offset(6,3),      Offset(),   Offset(5)},     // 12
            { Offset(6,3),      Offset(6),  Offset() },     // 13
            { Offset(6,1,4),    Offset(),   Offset(6)},     // 14
            { Offset(6,1,4),    Offset(7),  Offset() },     // 15
            { Offset(),         Offset(),   Offset(7)},     // 16
            { Offset(),         Offset(),   Offset() },     // 17
            { Offset(),         Offset(),   Offset() },     // 18
            { Offset(),         Offset(),   Offset() },     // 19
            { Offset(4),        Offset(),   Offset() },     // 20
            { Offset(5),        Offset(),   Offset() },     // 21
            { Offset(6),        Offset(),   Offset() },     // 22
            { Offset(7),        Offset(),   Offset() },     // 23
            { Offset(),         Offset(),   Offset() },     // 24
            { Offset(),         Offset(),   Offset() },     // 25
            { Offset(),         Offset(),   Offset() },     // 26
            { Offset(),         Offset(),   Offset() },     // 27
            { Offset(),         Offset(),   Offset() },     // 28
            { Offset(),         Offset(),   Offset() },     // 29
            { Offset(),         Offset(),   Offset() },     // 30
            { Offset(),         Offset(),   Offset() },     // 31
        };
        
        Offset noff[32][3] = {
            { Offset(),         Offset(),   Offset(0)},     //  0
            { Offset(),         Offset(),   Offset() },     //  1
            { Offset(),         Offset(0),  Offset() },     //  2
            { Offset(),         Offset(),   Offset(1)},     //  3
            { Offset(),         Offset(1),  Offset() },     //  4
            { Offset(),         Offset(),   Offset(2)},     //  5
            { Offset(),         Offset(2),  Offset() },     //  6
            { Offset(),         Offset(),   Offset(3)},     //  7
            { Offset(),         Offset(3),  Offset() },     //  8
            { Offset(),         Offset(),   Offset(4)},     //  9
            { Offset(),         Offset(4),  Offset() },     // 10
            { Offset(),         Offset(),   Offset(5)},     // 11
            { Offset(),         Offset(5),  Offset() },     // 12
            { Offset(),         Offset(),   Offset(6)},     // 13
            { Offset(),         Offset(6),  Offset() },     // 14
            { Offset(),         Offset(),   Offset(7)},     // 15
            { Offset(),         Offset(7),  Offset() },     // 16
            { Offset(),         Offset(),   Offset() },     // 17
            { Offset(),         Offset(),   Offset() },     // 18
            { Offset(),         Offset(),   Offset() },     // 19
            { Offset(),         Offset(),   Offset() },     // 20
            { Offset(),         Offset(),   Offset() },     // 21
            { Offset(),         Offset(),   Offset() },     // 22
            { Offset(),         Offset(),   Offset() },     // 23
            { Offset(),         Offset(),   Offset() },     // 24
            { Offset(),         Offset(),   Offset() },     // 25
            { Offset(),         Offset(),   Offset() },     // 26
            { Offset(),         Offset(),   Offset() },     // 27
            { Offset(),         Offset(),   Offset() },     // 28
            { Offset(),         Offset(),   Offset() },     // 29
            { Offset(),         Offset(),   Offset() },     // 30
            { Offset(),         Offset(),   Offset() },     // 31
        };

        int nbr[32][3] = {
                {  1,  16,  10 },   //  0
                { 24,   0,   2 },   //  1
                {  1,   9,   3 },   //  2
                { 17,   2,  12 },   //  3
                { 17,  11,   5 },   //  4
                { 18,   4,  14 },   //  5
                { 18,  13,   7 },   //  6 
                { 19,   6,  16 },   //  7
                { 19,  15,   9 },   //  8
                { 20,   8,   2 },   //  9
                { 20,   0,  11 },   // 10
                { 21,  10,   4 },   // 11
                { 21,   3,  13 },   // 12
                { 22,  12,   6 },   // 13
                { 22,   5,  15 },   // 14
                { 23,  14,   8 },   // 15
                { 23,   7,   0 },   // 16
                { 25,   3,   4 },   // 17
                { 26,   5,   6 },   // 18
                { 27,   7,   8 },   // 19
                { 28,   9,  10 },   // 20
                { 29,  11,  12 },   // 21
                { 30,  13,  14 },   // 22
                { 31,  15,  16 },   // 23
                { 25,  31,   1 },   // 24
                { 26,  24,  17 },   // 25
                { 27,  25,  18 },   // 26
                { 26,  19,  28 },   // 27
                { 29,  27,  20 },   // 28
                { 30,  28,  21 },   // 29
                { 31,  29,  22 },   // 30
                { 24,  30,  23 },   // 31
        };


        Vertex_handle vertices[14];
        for (int i = 0; i < vcount; i++) {
            vertices[i] = _tds.create_vertex();
            vertices[i]->set_point(pts[i]);
            vertices[i]->set_idx(i);
        }


        Face_handle faces[32];
        for (int i = 0; i < fcount; i++) {
            int x, y, z;
            faces[i] = _tds.create_face(vertices[tri[i][0]], vertices[tri[i][1]], vertices[tri[i][2]]);
            faces[i]->set_number(i);
        }


        // This needs to be done AFTER the faces are created
        for (int i = 0; i < fcount; i++) {
            faces[i]->set_neighbors(            faces[nbr[i][0]],   faces[nbr[i][1]],   faces[nbr[i][2]] );
            faces[i]->set_offsets(              off[i][0],          off[i][1],          off[i][2]        );
            faces[i]->set_neighbor_face_offsets(noff[i][0],         noff[i][1],         noff[i][2]       );
        }


        vertices[ 0]->set_face(faces[ 0]);
        vertices[ 1]->set_face(faces[ 2]);
        vertices[ 2]->set_face(faces[ 4]);
        vertices[ 3]->set_face(faces[ 6]);
        vertices[ 4]->set_face(faces[ 8]);
        vertices[ 5]->set_face(faces[20]);
        vertices[ 6]->set_face(faces[21]);
        vertices[ 7]->set_face(faces[22]);
        vertices[ 8]->set_face(faces[23]);
        vertices[ 9]->set_face(faces[ 1]);
        vertices[10]->set_face(faces[17]);
        vertices[11]->set_face(faces[18]);
        vertices[12]->set_face(faces[19]);
        vertices[13]->set_face(faces[28]);

        _tds.set_dimension(2);

        std::vector<Vertex_handle> ret(vcount);
        for (int i = 0; i < vcount; i++) {
            ret.push_back(vertices[i]);
        }


        return ret;
        }
        
    } // namespace CGAL
    
#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H
    
