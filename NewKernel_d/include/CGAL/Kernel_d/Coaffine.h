#ifndef CGAL_KD_COAFFINE_H
#define CGAL_KD_COAFFINE_H
#include <vector>
// DRAFT, DOESN'T WORK YET
namespace CGAL {
struct Orientation_in_subspace {
	std::vector<int> indices;
	// use indices.size() to know if the origin is added
	// an index of i means add the vector (0,...,1,...,0) with the 1 in i-th position
	// *OR* it means *don't* add this one...
};

template<class R_> struct Construct_orientation_in_subspace : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Construct_orientation_in_subspace)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::template Functor<Compute_cartesian_coordinate_tag>::type CCC;
	typedef typename R::LA LA;
	typedef typename LA::template Matrix<Dynamic_dimension,Dynamic_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;
	typedef typename R::template Functor<Point_dimension_tag>::type PD;
	typedef Orientation_in_subspace result_type;

	// This implementation is going to suck. Maybe we should push the
	// functionality into LA.
	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		PD pd (this->kernel());
		CCC ccc (this->kernel());
		int dim = pd(*f);
		Matrix coord (dim, dim);
		int col = 0;
		std::vector<int> proj;
		std::vector<int> rest; rest.reserve(dim);
		for(int i=0; i<dim; ++i) rest.push_back(i);
		for( ; f != e ; ++col, ++f ) {
			Point const&p=*f;
			for(int i=0; i<dim; ++i) coord(col, i) = ccc(p, i);
			int d = proj.size()+1;
			Matrix m (d, d);
			// Fill the matrix with what we already have
			for(int i=0; i<d; ++i)
			for(int j=0; j<d-1; ++j)
				m(i,j) = coord(i, proj[j]);
			// Try to complete with any other coordinate
			// TODO: iterate on rest by the end, or use a (forward_)list.
			for(std::vector<int>::iterator it=rest.begin();it!=rest.end();++it) {
				for(int i=0; i<d; ++i) m(i,d-1) = coord(i, *it);
				if(LA::sign_of_determinant(m)==0) continue;
				proj.push_back(*it);
				rest.erase(it);
			}
			CGAL_assertion(it!=rest.end());
		}
		Orientation_in_subspace o;
		// Check whether the origin is in the affine space (optimization)
		int d = proj.size();
		Matrix m (d+1, d+1);
		// Mu ? J'ai ecrit le code pour des vecteurs au lieu de points :-(
		for(int i=0; i<d; ++i)
			for(int j=0; j<d-1; ++j)
				m(i,j) = coord(i, proj[j]);
		return o;
	}
};

}
#endif
