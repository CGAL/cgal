#include "Scene.h"
#include <QInputDialog>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Timer.h>

void Scene::benchmark_intersections(const double duration)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    // constructs tree
    std::cout << "Construct AABB tree...";
    CGAL::Timer timer;
    timer.start();
    Facet_tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    // generates random queries
    const int nb_queries = 1000000;
    std::cout << "Generates random queries...";
    std::vector<Ray> rays;
    std::vector<Line> lines;
    std::vector<Plane> planes;
    std::vector<Segment> segments;
    timer.start();
    srand(0);
    int i = 0;
    for(i=0; i<nb_queries; i++)
    {
        rays.push_back(random_ray(tree.bbox()));
        lines.push_back(random_line(tree.bbox()));
        planes.push_back(random_plane(tree.bbox()));
        segments.push_back(random_segment(tree.bbox()));
    }
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    // bench for all functions and query types
    bench_intersections(tree,duration,DO_INTERSECT,"do_intersect()",rays,lines,planes,segments,nb_queries);
    bench_intersections(tree,duration,ANY_INTERSECTED_PRIMITIVE,"any_intersected_primitive()",rays,lines,planes,segments,nb_queries);
    bench_intersections(tree,duration,ANY_INTERSECTION,"any_intersection()",rays,lines,planes,segments,nb_queries);
    bench_intersections(tree,duration,NB_INTERSECTIONS,"number_of_intersected_primitives()",rays,lines,planes,segments,nb_queries);
    bench_intersections(tree,duration,ALL_INTERSECTED_PRIMITIVES,"all_intersected_primitives()",rays,lines,planes,segments,nb_queries);
    bench_intersections(tree,duration,ALL_INTERSECTIONS,"all_intersections()",rays,lines,planes,segments,nb_queries);
}

void Scene::bench_intersections(Facet_tree& tree,
                                const double duration,
                                const int function,
                                const char *function_name,
                                const std::vector<Ray>& rays,
                                const std::vector<Line>& lines,
                                const std::vector<Plane>& planes,
                                const std::vector<Segment>& segments,
                                const int nb_queries)
{
    std::cout << "Benchmark " << function_name << std::endl;
    bench_intersection<Segment>(tree,function,duration,"segment",segments,nb_queries);
    bench_intersection<Ray>(tree,function,duration,"ray",rays,nb_queries);
    bench_intersection<Line>(tree,function,duration,"line",lines,nb_queries);
    bench_intersection<Plane>(tree,function,duration,"plane",planes,nb_queries);
}


void Scene::benchmark_distances(const double duration)
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    CGAL::Timer timer;
    timer.start();
    std::cout << "Construct AABB tree and internal KD tree...";
    Facet_tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    // benchmark
    bench_closest_point(tree,duration);
    bench_squared_distance(tree,duration);
    bench_closest_point_and_primitive(tree,duration);
}

std::size_t Scene::nb_digits(std::size_t value)
{
    std::size_t nb_digits = 0;
    while(value > 0)
    {
        nb_digits++;
        value /= 10;
    }
    return nb_digits;
}

// bench memory against number of facets in the tree
// the tree is reconstructed each timer in the mesh
// refinement loop
void Scene::bench_memory()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    std::cout << std::endl << "Benchmark memory" << std::endl;
    std::cout << "#Facets, Bytes, Mbytes, Bytes/primitive" << std::endl;
    while(m_pPolyhedron->size_of_facets() < 2000000)
    {
        // refines mesh at increasing speed
        Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
        std::size_t digits = nb_digits(m_pPolyhedron->size_of_facets());
        unsigned int nb_splits =
          static_cast<unsigned int>(0.2 * std::pow(10.0,(double)digits - 1.0));
        refiner.run_nb_splits(nb_splits);

        // constructs tree and measure memory before then after
        typedef CGAL::Memory_sizer::size_type size_type;
        size_type before = CGAL::Memory_sizer().virtual_size();
        Facet_tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
        tree.do_not_accelerate_distance_queries(); // 150 vs 61 bytes per primitive!

        size_type after = CGAL::Memory_sizer().virtual_size();
        size_type bytes = after - before; // in Bytes
        double mbytes = (double)bytes / (double)1048576; //  in MBytes
        double bpp = (double)bytes / (double)m_pPolyhedron->size_of_facets();
        std::cout << m_pPolyhedron->size_of_facets() << ", "
            << bytes << ", "
            << mbytes << ", "
            << bpp << std::endl;
    }
}

void Scene::bench_construction()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    std::cout << std::endl << "Benchmark construction" << std::endl;
    std::cout << "#Facets    alone (s)   with KD-tree (s)" << std::endl;

    while(m_pPolyhedron->size_of_facets() < 1000000)
    {
        // refines mesh at increasing speed
        Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
        std::size_t digits = nb_digits(m_pPolyhedron->size_of_facets());
        unsigned int nb_splits =
          static_cast<unsigned int>(0.2 * std::pow(10.0,(double)digits - 1.0));
        refiner.run_nb_splits(nb_splits);

        // constructs tree
        CGAL::Timer time1;
        time1.start();
        Facet_tree tree1(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second, *m_pPolyhedron);
        double duration_construction_alone = time1.time();

        CGAL::Timer time2;
        time2.start();
        Facet_tree tree2(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);
        double duration_construction_and_kdtree = time2.time();

        std::cout << m_pPolyhedron->size_of_facets() << "\t"
            << duration_construction_alone     << "\t"
            << duration_construction_and_kdtree << std::endl;
    }
}

void Scene::bench_intersections_vs_nbt()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    std::cout << std::endl << "Benchmark intersections against #triangles" << std::endl;
    std::cout << std::endl << "for random ray queries and all_intersections()" << std::endl;
    std::cout << "#Facets, #queries/s" << std::endl;

    // generates 10K random ray queries
    const int nb_queries = 10000;
    srand(0);
    std::vector<Ray> queries;
    for(int i=0;i<nb_queries;i++)
        queries.push_back(random_ray(m_bbox));

    while(m_pPolyhedron->size_of_facets() < 1000000)
    {
        // refines mesh at increasing speed
        Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
        std::size_t digits = nb_digits(m_pPolyhedron->size_of_facets());
        unsigned int nb_splits =
          static_cast<unsigned int>(0.2 * std::pow(10.0,(double)digits - 1.0));
        refiner.run_nb_splits(nb_splits);

        // constructs tree (out of timing)
        Facet_tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second,*m_pPolyhedron);

        // calls ray queries
        CGAL::Timer timer;
        timer.start();
        std::list<Facet_tree::Object_and_primitive_id> intersections;
        for(int i=0;i<nb_queries;i++)
            tree.all_intersections(queries[i],std::back_inserter(intersections));
        double duration = timer.time();
        int speed = (int)((double)nb_queries / (double)duration);

        std::cout << m_pPolyhedron->size_of_facets() << ", " << speed << std::endl;
    }
}

void Scene::bench_distances_vs_nbt()
{
    if(m_pPolyhedron == NULL)
    {
        std::cout << "Load polyhedron first." << std::endl;
        return;
    }

    std::cout << std::endl << "Benchmark distances against #triangles" << std::endl;
    std::cout << std::endl << "for random point queries and closest_point()" << std::endl;
    std::cout << "#Facets, #queries/s" << std::endl;

    // generates 10K random point queries
    const int nb_queries = 10000;
    std::vector<Point> queries;
    srand(0);
    for(int i=0;i<nb_queries;i++)
        queries.push_back(random_point(m_bbox));

    while(m_pPolyhedron->size_of_facets() < 1000000)
    {
        // refines mesh at increasing speed
        Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
        std::size_t digits = nb_digits(m_pPolyhedron->size_of_facets());
        unsigned int nb_splits =
          static_cast<unsigned int>(0.2 * std::pow(10.0,(double)digits - 1.0));
        refiner.run_nb_splits(nb_splits);

        // constructs tree (out of timing)
        Facet_tree tree(faces(*m_pPolyhedron).first, faces(*m_pPolyhedron).second, *m_pPolyhedron);

        // calls queries
        CGAL::Timer timer;
        timer.start();
        for(int i=0;i<nb_queries;i++)
            tree.closest_point(queries[i]);
        double duration = timer.time();
        int speed = (int)((double)nb_queries / (double)duration);

        std::cout << m_pPolyhedron->size_of_facets() << ", " << speed << std::endl;
    }
}

template <class Query>
void Scene::bench_intersection(Facet_tree& tree,
                               const int function,
                               const double duration,
                               const char *query_name,
                               const std::vector<Query>& queries,
                               const int nb_queries)
{
    CGAL::Timer timer;
    timer.start();
    unsigned int nb = 0;
    std::list<Facet_tree::Primitive_id> primitive_ids;
    std::list<Facet_tree::Object_and_primitive_id> intersections;
    while(timer.time() < duration)
    {
        const Query& query = queries[nb % nb_queries]; // loop over vector
        switch(function)
        {
        case DO_INTERSECT:
            tree.do_intersect(query);
            break;
        case ANY_INTERSECTION:
            tree.any_intersection(query);
            break;
        case NB_INTERSECTIONS:
            tree.number_of_intersected_primitives(query);
            break;
        case ALL_INTERSECTIONS:
            tree.all_intersections(query,std::back_inserter(intersections));
            break;
        case ANY_INTERSECTED_PRIMITIVE:
            tree.any_intersected_primitive(query);
            break;
        case ALL_INTERSECTED_PRIMITIVES:
            tree.all_intersected_primitives(query,std::back_inserter(primitive_ids));
        }
        nb++;
    }

    double speed = (double)nb / (double)timer.time();
    std::cout << speed << " queries/s with " << query_name << std::endl;
}



//////////////////////////////////////////////////////////////
// BENCH DISTANCES
//////////////////////////////////////////////////////////////

void Scene::bench_distance(Facet_tree& tree,
                           const int function,
                           const double duration)
{

    // generates 100K random point queries
    srand(0);
    unsigned int nb_queries = 100000;
    std::vector<Point> queries;
    for(unsigned int i=0;i<nb_queries;i++)
        queries.push_back(random_point(tree.bbox()));

    CGAL::Timer timer;
    timer.start();
    unsigned int nb = 0;
    while(timer.time() < duration)
    {
        const Point& query = queries[nb%nb_queries];
        switch(function)
        {
        case SQ_DISTANCE:
            tree.squared_distance(query);
            break;
        case CLOSEST_POINT:
            tree.closest_point(query);
            break;
        case CLOSEST_POINT_AND_PRIMITIVE_ID:
            tree.closest_point_and_primitive(query);
        }
        nb++;
    }

    double speed = (double)nb / (double)timer.time();
    std::cout << speed << " queries/s" << std::endl;
}

void Scene::bench_squared_distance(Facet_tree& tree,
                                   const double duration)
{
    std::cout << "Benchmark squared distance" << std::endl;
    bench_distance(tree,SQ_DISTANCE,duration);
}

void Scene::bench_closest_point(Facet_tree& tree,
                                const double duration)
{
    std::cout << "Benchmark closest point" << std::endl;
    bench_distance(tree,CLOSEST_POINT,duration);
}

void Scene::bench_closest_point_and_primitive(Facet_tree& tree,
                                              const double duration)
{
    std::cout << "Benchmark closest point and primitive" << std::endl;
    bench_distance(tree,CLOSEST_POINT_AND_PRIMITIVE_ID,duration);
}


