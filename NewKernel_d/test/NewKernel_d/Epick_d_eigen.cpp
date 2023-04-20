#if  __cpp_aligned_new  >= 201606L
# define CGAL_NEWKERNEL_D_USE_EIGEN_VECTOR 1
# include "Epick_d.cpp"
#else
int main(){}
#endif
