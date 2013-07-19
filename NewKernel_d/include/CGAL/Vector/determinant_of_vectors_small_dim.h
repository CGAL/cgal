#ifndef CGAL_VECTOR_DETVEC_SMALL_H
#define CGAL_VECTOR_DETVEC_SMALL_H
#include <CGAL/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/determinant_of_vectors.h>

#define CGAL_ALLOWED_INCLUSION 1

#define CGAL_CLASS Add_determinant_of_vectors_small_dim
#define CGAL_TAG Has_determinant_of_vectors_tag
#define CGAL_FUNC determinant_of_vectors
#define CGAL_SIGN_FUNC sign_of_determinant_of_vectors
#define CGAL_SHIFT 0

#include <CGAL/Vector/determinant_of_vectors_small_dim_internal.h>

#undef CGAL_CLASS
#undef CGAL_TAG
#undef CGAL_FUNC
#undef CGAL_SIGN_FUNC
#undef CGAL_SHIFT

#define CGAL_CLASS Add_determinant_of_vectors_omit_last_small_dim
#define CGAL_TAG Has_determinant_of_vectors_omit_last_tag
#define CGAL_FUNC determinant_of_vectors_omit_last
#define CGAL_SIGN_FUNC sign_of_determinant_of_vectors_omit_last
#define CGAL_SHIFT 1

#include <CGAL/Vector/determinant_of_vectors_small_dim_internal.h>

#undef CGAL_CLASS
#undef CGAL_TAG
#undef CGAL_FUNC
#undef CGAL_SIGN_FUNC
#undef CGAL_SHIFT

#undef CGAL_ALLOWED_INCLUSION

#endif
