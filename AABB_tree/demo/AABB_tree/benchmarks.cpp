#include "Scene.h"
#include <QInputDialog>

void Scene::benchmark_intersections()
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	bench_do_intersect(tree);
	bench_nb_intersections(tree);
	bench_any_intersection(tree);
	bench_all_intersections(tree);
	bench_all_intersected_primitives(tree);
}

void Scene::benchmark_distances()
{
	QTime time;
	time.start();
	std::cout << "Construct AABB tree...";
	Facet_tree tree(m_pPolyhedron->facets_begin(),m_pPolyhedron->facets_end());
	tree.accelerate_distance_queries();
	std::cout << "done (" << time.elapsed() << " ms)" << std::endl;

	// TODO: check what the KD-tree is doing in there for large models.
	std::cout << "One single call to initialize KD-tree" << std::endl;
	tree.closest_point(CGAL::ORIGIN);

	// benchmark
	bench_closest_point(tree);
	bench_squared_distance(tree);
	bench_closest_point_and_primitive(tree);
}

//////////////////////////////////////////////////////////////
// BENCH INTERSECTIONS
//////////////////////////////////////////////////////////////

void Scene::bench_all_intersected_primitives(Facet_tree& tree)
{
	std::list<Primitive_id> primitive_ids;

	QTime time;
	time.start();
	std::cout << "Benchmark all_intersected_primitives" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Ray ray = random_ray();
		tree.all_intersected_primitives(ray,std::back_inserter(primitive_ids));
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;
	primitive_ids.clear();

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Line line = random_line();
		tree.all_intersected_primitives(line,std::back_inserter(primitive_ids));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;
	primitive_ids.clear();

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Segment segment = random_segment();
		tree.all_intersected_primitives(segment,std::back_inserter(primitive_ids));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Plane plane = random_plane();
		tree.all_intersected_primitives(plane,std::back_inserter(primitive_ids));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}


void Scene::bench_do_intersect(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark do_intersect" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Ray ray = random_ray();
		tree.do_intersect(ray);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Line line = random_line();
		tree.do_intersect(line);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Segment segment = random_segment();
		tree.do_intersect(segment);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Plane plane = random_plane();
		tree.do_intersect(plane);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}

void Scene::bench_nb_intersections(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark number_of_intersected_primitives" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Ray ray = random_ray();
		tree.number_of_intersected_primitives(ray);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Line line = random_line();
		tree.number_of_intersected_primitives(line);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Segment segment = random_segment();
		tree.number_of_intersected_primitives(segment);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Plane plane = random_plane();
		tree.number_of_intersected_primitives(plane);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}

void Scene::bench_any_intersection(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark any_intersection" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Ray ray = random_ray();
		tree.any_intersection(ray);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Line line = random_line();
		tree.any_intersection(line);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Segment segment = random_segment();
		tree.any_intersection(segment);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Plane plane = random_plane();
		tree.any_intersection(plane);
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}

void Scene::bench_all_intersections(Facet_tree& tree)
{
	std::list<Object_and_primitive_id> intersections;

	QTime time;
	time.start();
	std::cout << "Benchmark all_intersections" << std::endl;

	// with ray
	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Ray ray = random_ray();
		tree.all_intersections(ray,std::back_inserter(intersections));
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with ray" << std::endl;
	intersections.clear();

	// with line
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Line line = random_line();
		tree.all_intersections(line,std::back_inserter(intersections));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with line" << std::endl;
	intersections.clear();

	// with segment
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Segment segment = random_segment();
		tree.all_intersections(segment,std::back_inserter(intersections));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with segment" << std::endl;
	intersections.clear();

	// with plane
	nb = 0;
	time.start();
	while(time.elapsed() < 1000)
	{
		Plane plane = random_plane();
		tree.all_intersections(plane,std::back_inserter(intersections));
		nb++;
	}
	speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s with plane" << std::endl;
}


//////////////////////////////////////////////////////////////
// BENCH DISTANCES
//////////////////////////////////////////////////////////////

void Scene::bench_squared_distance(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark squared distance" << std::endl;

	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point query = random_point();
		tree.squared_distance(query);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s" << std::endl;
}


void Scene::bench_closest_point(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark closest point" << std::endl;

	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point query = random_point();
		tree.closest_point(query);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s" << std::endl;
}

void Scene::bench_closest_point_and_primitive(Facet_tree& tree)
{
	QTime time;
	time.start();
	std::cout << "Benchmark closest point and primitive" << std::endl;

	unsigned int nb = 0;
	while(time.elapsed() < 1000)
	{
		Point query = random_point();
		tree.closest_point_and_primitive(query);
		nb++;
	}
	double speed = 1000.0 * nb / time.elapsed();
	std::cout << speed << " queries/s" << std::endl;
}
