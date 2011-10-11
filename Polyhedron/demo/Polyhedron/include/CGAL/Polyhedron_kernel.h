// LP solver for kernel
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

// Taken from http://www.qhull.org/html/qhalf.htm

// If you do not know an interior point for the halfspaces, use linear programming
// to find one. Assume, n halfspaces defined by: aj*x1+bj*x2+cj*x3+dj>=0, j=1..n.
// Perform the following linear program:
//		max(x5) aj*x1+bj*x2+cj*x3+dj*x4-x5>=0, j=1..n

// Then, if [x1,x2,x3,x4,x5] is an optimal m_solution with x4,x5>0 we get:
//		aj*(x1/x4)+bj*(x2/x4)+cj*(x3/x4)+dj>=(x5/x4)>0, j=1..n
// and conclude that the point [x1/x4,x2/x4,x3/x4] is in the interior of all
// the halfspaces. Note that x5 is optimal, so this point is "way in" the
// interior (good for precision errors).

// After finding an interior point, the rest of the intersection algorithm is
// from Preparata & Shamos ['85, p. 316, "A simple case ..."]. Translate the
// halfspaces so that the interior point is the origin. Calculate the dual
// polytope. The dual polytope is the convex hull of the vertices dual to the
// original faces in regard to the unit sphere (i.e., halfspaces at distance
// d from the origin are dual to vertices at distance 1/d). Then calculate
// the resulting polytope, which is the dual of the dual polytope, and
// translate the origin back to the interior point [S. Spitz and S. Teller].

// NOTE: We changed this to max(x4) under constraints aj*x1+bj*x2+cj*x3+dj-x4>=0, j=1..n
//  i.e. aj*x1+bj*x2+cj*x3-x4 >= -dj, j=1..n
//  Then, if [x1,x2,x3,x4] is an optimal m_solution with x4 > 0 we pick
// the point [x1,x2,x3] as inside point.

template <class Kernel, class ET>
class Polyhedron_kernel
{
private:
	// basic types
	typedef typename Kernel::FT FT;
	typedef typename Kernel::Point_3 Point;
	typedef typename Kernel::Plane_3 Plane;
	typedef typename Kernel::Vector_3 Vector;
	typedef typename Kernel::Triangle_3 Triangle;

	// program and solution types
	typedef CGAL::Quadratic_program<double> LP;
	typedef typename CGAL::Quadratic_program_solution<ET> Solution;
	typedef typename Solution::Variable_value_iterator Variable_value_iterator;

	// linear program
	Solution m_solution;
	Point m_inside_point;

public:
	Polyhedron_kernel() {}
	~Polyhedron_kernel() {}

public:
	Point& inside_point() { return m_inside_point; }
	const Point& inside_point() const { return m_inside_point; }
	unsigned int number_of_iterations() { return m_solution.number_of_iterations(); }

public:

	template < class InputIterator >
	bool solve(InputIterator begin, InputIterator end)
	{
		// solve linear program with constraints Ax >= b
		LP lp(CGAL::LARGER,false);

		// column indices
		const int index_x1 = 0;
		const int index_x2 = 1;
		const int index_x3 = 2;
		const int index_x4 = 3;

		int j = 0;
		InputIterator it;
		for(it = begin; it != end; it++, j++)
		{
			const Triangle& triangle = *it;
			const Point& a = triangle[0];
			const Point& b = triangle[1];
			const Point& c = triangle[2];
			Plane plane(a,c,b); // note a c b to get inward normal

			Vector normal = unit_normal(plane);
			FT aj = normal.x();
			FT bj = normal.y();
			FT cj = normal.z();

			// constraint (line j): aj * x1 + bj * x2 + cj * x3 - x4 >= -dj
			lp.set_a(index_x1, j, aj);   // aj * x1
			lp.set_a(index_x2, j, bj);   // bj * x2
			lp.set_a(index_x3, j, cj);   // cj * x3
			lp.set_a(index_x4, j, -1.0); // -x4

			// right hand side (>= -dj)
			FT dj = distance_to_origin(plane) * CGAL::sign(plane.d());
			lp.set_b(j, -dj);  
		}

		// objective function -> max x4 (negative sign set because 
		// the lp solver always minimizes an objective function)
		lp.set_c(index_x4,-1.0);

		// add constraints over x4 (>= 0). somewhat redundant 
		// constraint as we later test for x4 > 0
		lp.set_l(index_x4,true,0.0); 

		// solve the linear program
		m_solution = CGAL::solve_linear_program(lp, ET());

		if(m_solution.is_infeasible())
			return false;

		if(!m_solution.is_optimal())
			return false;

		// get variables
		Variable_value_iterator X = m_solution.variable_values_begin();

		// Problem: X[0], X[1], ..., X[3] are proxy objects from
		// boost::facade_iterator. Those proxy objects are of type
		// operator_brackets_proxy<Variable_value_iterator>, which
		// is castable to a number type. In CGAL::to_double, the
		// automatic selection of a Real_embeddable_traits does not
		// work with such proxy objects (whose type is not
		// recognized). That is why we need to select the
		// Real_embeddable_traits manually, and create a functor
		// to_double manually. 	     -- Laurent Rineau, 2008/09/04
		typedef CGAL::Real_embeddable_traits<typename Variable_value_iterator::value_type> RE_traits;
		typename RE_traits::To_double to_double;

		// solution if x4 > 0
		double x4 = to_double(X[3]);
		if(x4 <= 0.0)
			return false;

		// define inside point as (x1;x2;x3)
		double x1 = to_double(X[0]); 
		double x2 = to_double(X[1]);
		double x3 = to_double(X[2]);
		m_inside_point = Point(x1,x2,x3);

		return true;
	}

	Vector unit_normal(const Plane& plane)
	{
		Vector n = plane.orthogonal_vector();
		return n / std::sqrt(n * n);
	}

	FT distance_to_origin(const Plane& plane)
	{
		return std::sqrt(CGAL::squared_distance(Point(CGAL::ORIGIN),plane));
	}

};

template <class Polyhedron, class OutputIterator>
void get_triangles(Polyhedron& polyhedron,
		   OutputIterator out)
{
	typedef typename Polyhedron::Facet_iterator Facet_iterator;
	typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
        typedef typename Polyhedron::Point_3 Point;
        typedef typename Polyhedron::Traits::Triangle_3 Triangle;

	for(Facet_iterator f = polyhedron.facets_begin();
		f != polyhedron.facets_end();
		f++)
	{
		Halfedge_handle he = f->halfedge();
		const Point& a = he->vertex()->point();
		const Point& b = he->next()->vertex()->point();
		const Point& c = he->next()->next()->vertex()->point();
		*out++ = Triangle(a,b,c);
	}
}
