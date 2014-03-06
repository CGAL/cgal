namespace CGAL {

/*!
\ingroup PkgGenerators

\brief computes a random convex polygon by writing its vertices (oriented
counterclockwise) in the list `l`, as the convex hull of \f$ n \$ random points in a disc centered in \f$0\f$ with radius \f$radius\$.
The generated polygon will have an average number of vertices \$f n^\frac{1}{3}(1+o(1))\f$.


\cgalHeading{Requirements}

- `Traits` is a model of the concept RandomPolygonTraits_2 
- `PointGenerator::value_type` is equivalent to 
`Traits::Point_2` and `OutputIterator::value_type`. 


The default traits class `Default_traits` is the kernel in which 
`Traits::Point_2` is defined. 

\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>` 
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>` 

\cgalHeading{Implementation}

The implementation is based on the method of eliminating self-intersections in 
a polygon by using so-called "2-opt" moves. Such a move eliminates an 
intersection between two edges by reversing the order of the vertices between 
the edges. No more than \f$ O(n^3)\f$ such moves are required to simplify a polygon 
defined on \f$ n\f$ points \cgalCite{ls-utstp-82}. 
Intersecting edges are detected using a simple sweep through the vertices 
and then one intersection is chosen at random to eliminate after each sweep. 
The worse-case running time is therefore \f$ O(n^4 \log n)\f$. 

\cgalHeading{Example}

The following program displays a random simple polygon with up to 100 
vertices, where the vertex coordinates are drawn uniformly from the 
unit square centered at the origin. 

\cgalExample{Generator/random_convex_polygon.cpp} 

*/

    template<class P,class GEN>
    void convex_random_polygon(size_t n,  typename Kernel_traits<P>::Kernel::FT radius, std::list<P> & l,GEN & gen, bool fast=true );
} /* namespace CGAL */
