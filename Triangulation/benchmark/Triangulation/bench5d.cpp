//  Switch off the static filters for 4, 5, and 6 points.
//  #define CGAL_NO_STATIC_FILTER_456

#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

typedef CGAL::Epick_d< CGAL::Dimension_tag<5> >  K;
typedef CGAL::Triangulation_vertex<K>                                                     Vertex;
typedef CGAL::Triangulation_ds_full_cell<void,CGAL::TDS_full_cell_mirror_storage_policy>  DS_full_cell;
typedef CGAL::Triangulation_full_cell<K,CGAL::No_full_cell_data, DS_full_cell>            Full_cell;
typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<5>,
                                            Vertex,
                                            Full_cell>                                    TDS;
typedef CGAL::Triangulation<K,TDS>                                                        Triangulation;


int main()
{
    const int D = 5;   // we work in Euclidean 5-space

    std::vector<Triangulation::Point> points;
    std::ifstream in("points_5.txt");
    Triangulation::Point p;
    int d;
    in >> d;
    assert(d == D);
    int n;
    in >> n;
    points.reserve(n);
    while (in >> p) {
        points.push_back(p);
    }

    CGAL::Timer timer;
    timer.start();
    Triangulation t(D);                      // create triangulation
    assert(t.empty());

    //std::array<double, 5> ar = { 0, 0, 0, 0, 0 };
    //Triangulation::Point orig(ar.begin(), ar.end());
    //t.insert(orig);  // insert origin
    t.insert(points.begin(), points.end());  // compute triangulation

    std::cout << "Time taken to insert "<< points.size() << " points: " << timer.time() << " seconds." << std::endl;
    timer.reset();
    assert( t.is_valid() );
    // - - - - - - - - - - - - - - - - - - - - - - - - STEP 2
    typedef Triangulation::Face Face;
    typedef std::vector<Face> Faces;
    Faces edges;
    std::back_insert_iterator<Faces> out(edges);
    t.tds().incident_faces(t.infinite_vertex(), 1, out);
    // collect faces of dimension 1 (edges) incident to the infinite vertex
    std::cout << "There are " << edges.size()
              << " vertices on the convex hull." << std::endl;
    std::cout << "Time taken to find edges: " << timer.time() << " seconds." << std::endl;
     return 0;
}