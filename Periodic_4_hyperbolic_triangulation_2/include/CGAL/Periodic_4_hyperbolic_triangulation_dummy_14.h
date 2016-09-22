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
        Transformation rot8(ROTATION, -CGAL::sqrt(FT(2)-CGAL::sqrt(FT(2)))/FT(2), CGAL::sqrt(FT(2)+CGAL::sqrt(FT(2)))/FT(2));
        Transformation opposite(ROTATION, FT(0), FT(-1)); // sin(theta), cos(theta)

        int fcount = 32;    // Faces count
        int vcount = 14;    // Vertices count

        FT F0 = FT(0);
        FT F1 = FT(1);
        FT F2 = FT(2);

        std::vector<typename GT::Point_2> pts;

        pts.push_back(Point( FT(0), FT(0) ));   // origin -> pt 0

        FT i1(CGAL::sqrt(FT(2) - sqrt(FT(2))));
        FT i2(FT(2)*i1);
        FT i3(FT(2)*CGAL::sqrt(FT(2)) - i2);
        FT i4(-1 + i3);
        Point a0(FT( CGAL::sqrt( i4 )), FT(0));          // internal point
        
        a0 = rot8(a0);
        pts.push_back(a0);
        for (int k = 1; k < 8; k++) {
            a0 = rotate(a0);
            pts.push_back(a0);   // pt 1 - 8
        }


        pts.push_back(Point(-FT( CGAL::sqrt(CGAL::sqrt(F2) - F1)), F0));                                                          //  μ(s_0)  pt 9
        pts.push_back(Point(-FT( CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2) ), -FT(CGAL::sqrt((CGAL::sqrt(F2) - F1) / F2)) ));       //  μ(s_1)  pt 10
        pts.push_back(Point(F0,                                             -FT(CGAL::sqrt(CGAL::sqrt(F2) - F1))));               //  μ(s_2)  pt 11
        pts.push_back(Point(FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)),    -FT(CGAL::sqrt((CGAL::sqrt(F2) - F1) / F2)) ));       //  μ(s_3)  pt 12
        pts.push_back(Point(FT( CGAL::sqrt(CGAL::sqrt(F2) + F1)/ F2 ),      -FT(CGAL::sqrt(CGAL::sqrt(F2) - F1)/ F2) ));          //  v_0     pt 13
        
        // for (int i = 0; i < 4; i++) {
        //     pts[12-i] = opposite(pts[12-i]);
        // }

        Offset off[32][3]; 
        for (int i = 0; i < 32; i++) {
            for (int j = 0; j < 3; j++) {
                off[i][j] = Offset();
            }
        }

        int tri[32][3];
        for (int i = 0; i < 4; i++) {
            Offset vo1, vo2;
            if (i == 0) { vo1 = Offset(1, 6, 3);    vo2 = Offset();         }
            if (i == 1) { vo1 = Offset(1, 4);       vo2 = Offset(1, 6, 3);  }
            if (i == 2) { vo1 = Offset(3);          vo2 = Offset(1, 4);     }
            if (i == 3) { vo1 = Offset(4, 1, 6, 3); vo2 = Offset(3);        }

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
            off[fn][0] = Offset(i);
            //----------------------
            fn = i + 12;
            i0 = i+9;
            i1 = ((i+5)%8)+1;
            i2 = ((i+4)%8)+1;
            tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
            //off[fn][0] = Offset(i+4);
            //----------------------
            fn = 2*(i + 8);
            i0 = i+1;
            i1 = ((i+5)%8)+1;
            i2 = i+9;
            tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
            off[fn][1] = Offset(i);
            off[fn][2] = Offset(i);
            //----------------------
            fn = 2*(i + 8) + 1;
            i0 = i+2;
            i1 = i+9;
            i2 = i+5;
            tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
            off[fn][1] = Offset(i);
            off[fn][2] = Offset(i);
            //----------------------
            fn = 2*(i + 12);
            i0 = ((i+5)%8)+1;
            i1 = i+1;
            i2 = 13;
            tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
            off[fn][0] = Offset(i);
            off[fn][2] = vo2;
            //----------------------
            fn = 2*(i + 12) + 1;
            i0 = i+5;
            i1 = 13;
            i2 = i+2;
            tri[fn][0] = i0; tri[fn][1] = i1; tri[fn][2] = i2;
            off[fn][0] = Offset(i);
            off[fn][1] = vo1;
        }



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
            faces[i]->set_offsets(      off[i][0],           off[i][1],           off[i][2]);
            faces[i]->set_number(i);
        }

        _tds.set_dimension(2);

        for (int i = 0; i < 4; i++) {
            int fn = i;
            int n0 = i+8;
            int n1 = i+1;
            _tds.set_adjacency( faces[fn], 0, faces[n0], 0 );
            _tds.set_adjacency( faces[fn], 1, faces[n1], 2 );
            
            fn = i+4;
            n0 = i+12;
            n1 = (i+5) % 8;
            _tds.set_adjacency( faces[fn], 0, faces[n0], 0 );
            _tds.set_adjacency( faces[fn], 1, faces[n1], 2 );
           
            fn = i+8;
            n0 = 2*(i+8);
            n1 = 2*(i+8)+1;
            _tds.set_adjacency( faces[fn], 1, faces[n0], 1 );
            _tds.set_adjacency( faces[fn], 2, faces[n1], 2 );
           
            fn = 2*(i+8);
            n0 = (i+12);
            n1 = 2*n0;
            _tds.set_adjacency( faces[fn], 0, faces[n0], 2 );
            _tds.set_adjacency( faces[fn], 2, faces[n1], 2 );
           
            fn = 2*(i+8) + 1;
            n0 = (i+12);
            n1 = 2*n0 + 1;
            _tds.set_adjacency( faces[fn], 0, faces[n0], 1 );
            _tds.set_adjacency( faces[fn], 1, faces[n1], 1 );
           

            if (i < 3) {
                fn = 2*(i+12);
                n0 = 24 + ((i+1)*2 + 1);
                _tds.set_adjacency( faces[fn], 1, faces[n0], 2 );

                fn = 2*(i+12)+1;
                n0 = 25 + (i*2 + 1);
                _tds.set_adjacency( faces[fn], 0, faces[n0], 0 );
            }
        }

        int fn = 24;
        int n0 = 30;
        _tds.set_adjacency( faces[fn], 0, faces[n0], 1 );

        fn = 25;
        n0 = 31;
        _tds.set_adjacency( faces[fn], 2, faces[n0], 0 );


        vertices[0]->set_face(faces[0]);
        // vertices 1-8
        for (int i = 1; i < 9; i++) {
            vertices[i]->set_face(faces[i+7]);
        }
        // vertices 9-12
        for (int i = 9; i < 13; i++) {
            vertices[i]->set_face(faces[i+3]);
        }  
        vertices[13]->set_face(faces[30]);

        std::vector<Vertex_handle> ret(vcount);
        for (int i = 0; i < vcount; i++) {
            ret.push_back(vertices[i]);
        }


        for (int i = 0; i < 32; i++) {
            cout << "Face " << i << ": vertices " << tri[i][0] << ", " << tri[i][1] << ", " << tri[i][2] << "    offsets " << off[i][0] << ", " << off[i][1] << ", " << off[i][2] << endl;
        }


        CGAL_triangulation_assertion(is_valid(true));

        return ret;
        }


    } // namespace CGAL
    
#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_DUMMY_14_H
    
