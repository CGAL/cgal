#ifndef __EIGEN_2__
#define __EIGEN_2__

// extract eigenvalues and eigenvectors from a 
// 2x2 symmetric matrix. Matrix numbering:
// a b
// b c 
// eigen values and vectors are sorted 
// in descendent order
template <typename K>
void eigen_symmetric_2(typename K::FT *matrix, // a b c
                       std::pair<typename K::Vector_2,
                                 typename K::Vector_2>& eigen_vectors,
                       std::pair<typename K::FT,
                                 typename K::FT>& eigen_values) 
{
  typedef typename K::FT FT;
  typedef typename K::Vector_2 Vector;
  // todo: replace
  FT a = matrix[0];
  FT b = matrix[1];
  FT c = matrix[2];
  FT p = c*c - 2*a*c + 4*b*b + a*a;
  if(p == 0) // isotropic case
    {
      eigen_values.first = (FT)1;
      eigen_values.second = (FT)1;
      // any vector is eigen vector
      eigen_vectors.first = Vector((FT)1,(FT)0);
      eigen_vectors.second = Vector((FT)0,(FT)1);
    }
  else
    {
      FT l1 = 0.5 * (-std::sqrt(p)+c+a);
      FT l2 = 0.5 * (std::sqrt(p)+c+a);
      if(l1 >= l2) // todo: check case where b = 0
	{
	  eigen_values.first = l1;
	  eigen_values.second = l2;
	  eigen_vectors.first = Vector((FT)1,-(std::sqrt(p)-c+a)/(b*2));
	  eigen_vectors.second = Vector((FT)1,(std::sqrt(p)+c-a)/(b*2));
	}
      else
	{
	  eigen_values.first = l2;
	  eigen_values.second = l1;
	  eigen_vectors.first = Vector((FT)1,(std::sqrt(p)+c-a)/(b*2));
	  eigen_vectors.second = Vector((FT)1,-(std::sqrt(p)-c+a)/(b*2));
	}
    }
}

#endif // __EIGEN_2__

