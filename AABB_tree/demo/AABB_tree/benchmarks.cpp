#include "Scene.h"
#include <QInputDialog>
#include <CGAL/Memory_sizer.h>

void Scene::benchmark_intersections(const int duration)
{
	// constructs tree
	std::cout << "Construct AABB tree...";
	QTime time;
	time.start();
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	// generates random queries
	const int nb_queries = 1000000;
	std::cout << "Generates random queries...";
	std::vector<Ray> rays;
	std::vector<Line> lines;
	std::vector<Plane> planes;
	std::vector<Segment> segments;
	time.start();
	srand(0);
	int i = 0;
	for(i=0; i<nb_queries; i++)
	{
		rays.push_back(random_ray(tree.bbox()));
		lines.push_back(random_line(tree.bbox()));
		planes.push_back(random_plane(tree.bbox()));
		segments.push_back(random_segment(tree.bbox()));
	}
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	// bench for all functions and query types
	bench_intersections(tree,duration,DO_INTERSECT,"do_intersect()",rays,lines,planes,segments,nb_queries);
	bench_intersections(tree,duration,ANY_INTERSECTED_PRIMITIVE,"any_intersected_primitive()",rays,lines,planes,segments,nb_queries);
	bench_intersections(tree,duration,ANY_INTERSECTION,"any_intersection()",rays,lines,planes,segments,nb_queries);
	bench_intersections(tree,duration,NB_INTERSECTIONS,"number_of_intersected_primitives()",rays,lines,planes,segments,nb_queries);
	bench_intersections(tree,duration,ALL_INTERSECTED_PRIMITIVES,"all_intersected_primitives()",rays,lines,planes,segments,nb_queries);
	bench_intersections(tree,duration,ALL_INTERSECTIONS,"all_intersections()",rays,lines,planes,segments,nb_queries);
}

void Scene::bench_intersections(Facet_tree& tree,
									              const int duration,
													      const int function,
													      char *function_name,
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


void Scene::benchmark_distances(const int duration)
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree and internal KD tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	tree.accelerate_distance_queries();
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	// benchmark
	bench_closest_point(tree,duration);
	bench_squared_distance(tree,duration);
	bench_closest_point_and_primitive(tree,duration);
}

unsigned int Scene::nb_digits(unsigned int value)
{
	unsigned int nb_digits = 0;
	while(value > 0)
	{
		nb_digits++;
		value /= 10;
	}
	return nb_digits;
}

// bench memory against number of facets in the tree
// the tree is reconstructed each time in the mesh 
// refinement loop
void Scene::bench_memory()
{
	std::cout << std::endl << "Benchmark memory" << std::endl;
	std::cout << "#Facets Bytes" << std::endl;
	while(m_pPolyhedron->size_of_facets() < 2000000)
	{
		// refines mesh at increasing speed
		Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
		unsigned int digits = nb_digits(m_pPolyhedron->size_of_facets());
		unsigned int nb_splits = (unsigned int)(0.2 * std::pow(10.0,(double)digits - 1.0));
		refiner.run_nb_splits(nb_splits);

		// constructs tree and measure memory before then after
		long before = CGAL::Memory_sizer().virtual_size();
		Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
		// tree.accelerate_distance_queries(); (100 vs 60 bytes per primitive!)

		long after = CGAL::Memory_sizer().virtual_size();
		long memory = after - before; // in Bytes
		std::cout << m_pPolyhedron->size_of_facets() << "\t" << memory << std::endl;
	}
}

void Scene::bench_construction()
{
	std::cout << std::endl << "Benchmark construction" << std::endl;
	std::cout << "#Facets    alone (ms)   with KD-tree (ms)" << std::endl;

	while(m_pPolyhedron->size_of_facets() < 1000000)
	{
		// refines mesh at increasing speed
		Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
		unsigned int digits = nb_digits(m_pPolyhedron->size_of_facets());
		unsigned int nb_splits = (unsigned int)(0.2 * std::pow(10.0,(double)digits - 1.0));
		refiner.run_nb_splits(nb_splits);

		// constructs tree
		QTime time1;
		time1.start();
		Facet_tree tree1(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
		double duration_construction_alone = time1.elapsed();

		QTime time2;
		time2.start();
		Facet_tree tree2(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
		tree2.accelerate_distance_queries();
		double duration_construction_and_kdtree = time2.elapsed();

		std::cout << m_pPolyhedron->size_of_facets() << "\t" 
			        << duration_construction_alone     << "\t" 
							<< duration_construction_and_kdtree << std::endl;
	}
}

void Scene::bench_intersections_vs_nbt()
{
	std::cout << std::endl << "Benchmark intersections against #triangles" << std::endl;
	std::cout << std::endl << "for random ray queries and all_intersections()" << std::endl;
	std::cout << "#Facets    #queries/s" << std::endl;

	while(m_pPolyhedron->size_of_facets() < 1000000)
	{
		// refines mesh at increasing speed
		Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
		unsigned int digits = nb_digits(m_pPolyhedron->size_of_facets());
		unsigned int nb_splits = (unsigned int)(0.2 * std::pow(10.0,(double)digits - 1.0));
		refiner.run_nb_splits(nb_splits);

		// constructs tree
		Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());

		// calls 10K random ray queries (neglects random generation of ray query)
		std::list<Object_and_primitive_id> intersections;
		QTime time;
		time.start();
		const int nb_queries = 10000;
		for(int i=0;i<nb_queries;i++)
		{
			Ray ray = random_ray(tree.bbox());
			tree.all_intersections(ray,std::back_inserter(intersections));
		}
		double duration = time.elapsed();
		int speed = (int)(1000 * nb_queries / duration);

		std::cout << m_pPolyhedron->size_of_facets() << "\t" << speed << std::endl;
	}
}

void Scene::bench_distances_vs_nbt()
{
	std::cout << std::endl << "Benchmark distances against #triangles" << std::endl;
	std::cout << std::endl << "for random point queries and closest_point()" << std::endl;
	std::cout << "#Facets    #queries/s" << std::endl;

	while(m_pPolyhedron->size_of_facets() < 1000000)
	{
		// refines mesh at increasing speed
		Refiner<Kernel,Polyhedron> refiner(m_pPolyhedron);
		unsigned int digits = nb_digits(m_pPolyhedron->size_of_facets());
		unsigned int nb_splits = (unsigned int)(0.2 * std::pow(10.0,(double)digits - 1.0));
		refiner.run_nb_splits(nb_splits);

		// constructs tree
		Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
		tree.accelerate_distance_queries();

		// calls 100K random point queries (neglects random generation of ray query)
		std::list<Object_and_primitive_id> intersections;
		QTime time;
		time.start();
		const int nb_queries = 10000;
		for(int i=0;i<nb_queries;i++)
		{
			Point query = random_point(tree.bbox());
			tree.closest_point(query);
		}
		double duration = time.elapsed();
		int speed = (int)(1000 * nb_queries / duration);

		std::cout << m_pPolyhedron->size_of_facets() << "\t" << speed << std::endl;
	}
}

template <class Query>
void Scene::bench_intersection(Facet_tree& tree,
									const int function,
									const int duration,
									const char *query_name,
									const std::vector<Query>& queries,
									const int nb_queries)
{
	QTime time;
	time.start();
	unsigned int nb = 0;
	std::list<Primitive_id> primitive_ids;
	std::list<Object_and_primitive_id> intersections;
	while(time.elapsed() < duration)
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

	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with " << query_name << std::endl;
}



//////////////////////////////////////////////////////////////
// BENCH DISTANCES
//////////////////////////////////////////////////////////////

void Scene::bench_distance(Facet_tree& tree,
						   const int function,
						   const int duration)
{
	QTime time;
	time.start();
	unsigned int nb = 0;
  srand(0);
	while(time.elapsed() < duration)
	{
		Point query = random_point(tree.bbox());
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

	// subtract random generation time to be more accurate
	double duration_random_and_tests = time.elapsed();
	time.start();
	for(unsigned int i=0;i<nb;i++)
		random_point(tree.bbox());
	double duration_random = time.elapsed();

	double elapsed = duration_random_and_tests - duration_random;
	double speed = 1000.0 * nb / elapsed;
	std::cout << speed << " queries/s" << std::endl;
}

void Scene::bench_squared_distance(Facet_tree& tree,
								   const int duration)
{
	std::cout << "Benchmark squared distance" << std::endl;
	bench_distance(tree,SQ_DISTANCE,duration);
}

void Scene::bench_closest_point(Facet_tree& tree,
								const int duration)
{
	std::cout << "Benchmark closest point" << std::endl;
	bench_distance(tree,CLOSEST_POINT,duration);
}

void Scene::bench_closest_point_and_primitive(Facet_tree& tree,
											  const int duration)
{
	std::cout << "Benchmark closest point and primitive" << std::endl;
	bench_distance(tree,CLOSEST_POINT_AND_PRIMITIVE_ID,duration);
}


