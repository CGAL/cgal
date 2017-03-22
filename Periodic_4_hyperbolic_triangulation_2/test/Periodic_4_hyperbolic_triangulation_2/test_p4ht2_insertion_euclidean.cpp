// Copyright (c) 2016-2017 INRIA Nancy Grand-Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Iordan Iordanov <iordan.iordanov@loria.fr>

#include <CGAL/basic.h>
#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

#include <CGAL/Timer.h>


typedef double                                                                  NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>                                  Triangulation;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;

typedef CGAL::Creator_uniform_2<double, Point>                                  Creator;

using namespace std;

int main(int argc, char** argv) {    

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number_of_points_to_insert] [optional: number_of_iterations]" << endl;
        return -1;
    }

    int N = atoi(argv[1]);
    int iters = 1;
    if (argc == 3) {
        iters = atoi(argv[2]);
    }

    cout << "---- for best results, make sure that you have compiled me in Release mode ----" << endl;
    
    double extime = 0.0;

    for (int exec = 1; exec <= iters; exec++) {
        std::vector<Point> pts;
        pts.reserve(N);
        CGAL::Random_points_in_disc_2<Point, Creator> g(1.0);
        CGAL::cpp11::copy_n( g, N, std::back_inserter(pts));

        cout << "iteration " << exec << ": inserting into triangulation (rational dummy points)... "; cout.flush();
        Triangulation tr;
        //tr.insert_dummy_points(true);  

        CGAL::Timer tt;
        tt.start();
        tr.insert(pts.begin(), pts.end());
        //for (int j = 0; j < pts.size(); j++) {
        //    tr.insert(pts[j]);
        //}
        tt.stop();
        cout << "DONE! (# of vertices = " << tr.number_of_vertices() << ", time = " << tt.time() << " secs)" << endl;
        extime += tt.time();
    }

    double avgtime = extime / (double)iters;
    cout << "---------------------------------------" << endl;
    cout << "Average execution time over " << iters << " iterations: " << avgtime << " secs" << endl << endl;

    return 0;
}
