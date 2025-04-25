
#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

typedef CGAL::Epick_d< CGAL::Dimension_tag<6> >  K;
typedef CGAL::Triangulation_vertex<K, std::ptrdiff_t>                                     Vertex;
typedef CGAL::Triangulation_ds_full_cell<void,CGAL::TDS_full_cell_mirror_storage_policy>  DS_full_cell;
typedef CGAL::Triangulation_full_cell<K,CGAL::No_full_cell_data, DS_full_cell>            Full_cell;
typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<6>,
                                            Vertex,
                                            Full_cell>                                    TDS;
typedef CGAL::Triangulation<K,TDS>                                                        Triangulation;


int main(int argc, char* argv[])
{
    const int D = 6;   // we work in Euclidean 6-space

    std::vector<Triangulation::Point> points;
    std::ifstream in(argv[1]);
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


    t.insert_and_index(points.begin(), points.end());  // compute triangulation

    std::cout << "Time taken to insert "<< points.size() << " points: " << timer.time() << " seconds." << std::endl;
    timer.reset();
    assert( t.is_valid() );

    std::vector<Triangulation::Full_cell_handle> infinite_cells;
    t.tds().incident_full_cells(t.infinite_vertex(), std::back_inserter(infinite_cells));

    std::set<std::pair<std::ptrdiff_t, std::ptrdiff_t>> edges;

    for(auto ch : infinite_cells) {
        std::vector<Triangulation::Vertex_handle> vertices;
        std::cout << "cell" << std::endl;
        vertices.reserve(D);
        for(int i = 0; i < D + 1; ++i) {
            if(ch->vertex(i) != t.infinite_vertex()) {
                vertices.push_back(ch->vertex(i));
                std::cout << ch->vertex(i)->data() << " ";
            }
        }
        std::cout << std::endl;
        for (int i = 0; i < D-1; i++) {
          for (int j = 0; j < D-1; j++) {

            if((i != j) && (vertices[i]->data() < vertices[j]->data())) {
                edges.insert(std::make_pair(vertices[i]->data(), vertices[j]->data()));
                }
          }
        }
    }
    for(auto pair : edges) {
        std::cout << "edge  between vertices " << pair.first << " and " << pair.second << std::endl;
    }

    std::cout << "Time taken to find edges: " << timer.time() << " seconds." << std::endl;

     return 0;
}