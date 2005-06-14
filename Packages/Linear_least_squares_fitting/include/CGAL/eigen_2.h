#ifndef __EIGEN_2__
#define __EIGEN_2__

template <typename K>
void eigen_2(typename const K::FT *matrix,
             typename std::pair<K::Vector_2>& eigen_vectors,
             typename std::pair<K::FT>& eigen_values) 
{
}

// extract eigenvalues and eigenvectors from a 
// 2x2 symmetric matrix. Matrix numbering:
// a b
// b c 
template <typename K>
void eigen_symmetric_2(typename const K::FT *matrix,
                       typename std::pair<K::Vector_2>& eigen_vectors,
                       typename std::pair<K::FT>& eigen_values) 
{

}

#endif // __EIGEN_2__

