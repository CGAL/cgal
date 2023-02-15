#include <CGAL/Linear_algebraHd.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/random_selection.h>
#include <CGAL/test_macros.h>
#include <CGAL/double.h>
#include <cstdlib>
#define VEC_DIM 10
#define MAT_DIM 10

#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>

#if defined( CGAL_USE_LEDA) || defined ( CGAL_USE_GMP )

typedef CGAL::Exact_integer RT;
typedef CGAL::Exact_rational FT;

#else

// The following are too slow :
// #include <CGAL/MP_Float.h>
// #include <CGAL/Quotient.h>
// typedef CGAL::MP_Float     RT;
// typedef CGAL::Quotient<RT> FT;
typedef double RT;
typedef double FT;

#endif

int main(int argc, char* argv[])
{
  CGAL_KD_SETDTHREAD(151);
  CGAL::IO::set_pretty_mode ( std::cerr );
  CGAL_TEST_START;
  {
    typedef RT NT;
    typedef CGAL::Linear_algebraHd<RT> LA;
    typedef LA::Matrix Matrix;
    typedef LA::Vector Vector;
    bool IOTEST = true;
    {
      int vec_dim;
      if (argc == 2) vec_dim = std::atoi(argv[1]);
      else           vec_dim = VEC_DIM;

      /* some construction and access ops */
      Vector v0(vec_dim), v1(vec_dim), v2(vec_dim);
      int F[] = { 1,2,3,4,5 };
      Vector v11(F,F+2), v12(F,F+3), v13(F,F+4), v14(F,F+5),
                 v15(v13.begin(),v13.end()),
                 v16(vec_dim,NT(1));
      CGAL_TEST(v13==v15){}
      CGAL_TEST(v0==v1){}
      for (int i = 0; i < vec_dim; i++) {
        v1[i] = i;
        v2[i] = vec_dim-i;
      }
      CGAL_TEST(v0!=v1){}
      Vector v3(v1);
      CGAL_TEST(v3==v1){}

      v1 += NT(2)*v2;
      v1 -= v2;
      v1 = -v1;
      CGAL_TEST((v1*v1 == NT(vec_dim*vec_dim*vec_dim))){}
      /* squared length of (v1 + v2) = v1*v1
         should be equal to dim*dim*dim */
      NT res1 = (v1*v1 - v2*v2);
      NT res2 = ((v1-v2)*(v1+v2));
      CGAL_TEST(res1==res2){}
      /* we test the third binomial formula */
      Vector v21(v13),v22(v13);
      v21 *= 13;
      CGAL_TEST(v21 == NT(13)*v22){}
      v21 /= 13;
      CGAL_TEST(v21 == v22){}
      if (IOTEST) CGAL_IO_TEST(v1,v2,CGAL::IO::ASCII);
    }

    {
      int mat_dim;
      if (argc == 2) mat_dim = std::atoi(argv[1]);
      else           mat_dim = MAT_DIM;

      /* some construction and access ops */
      Matrix::Identity ID;
      Matrix A(mat_dim,mat_dim), B(mat_dim),
             I(mat_dim,ID), One(mat_dim,mat_dim,NT(2));
      CGAL_TEST(A==B){}
      CGAL_TEST(One*I==One){}
      std::vector<Vector> F(mat_dim);

      int i,j;
      for (i = 0; i < mat_dim; i++) {
        Vector v(mat_dim);
        for (j = 0; j < mat_dim; j++) {
          A(i,j) = i;
          B(i,j) = j;
          v[j] = j;
        }
        F[i] = v;
      }
    #ifndef CGAL_SIMPLE_INTERFACE
      Matrix C(F), D(A), E(F.begin(),F.end());
    #else
      Matrix C(F), D(A), E(F);
    #endif
      CGAL_TEST(C==F && E==F){}
      /* A = (0,0,0,0...)
             (1,1,1,1...)
             (2,2,2,2...)
             ...
         B = (0,1,2,3...)
             (0,1,2,3...)
             (0,1,2,3...)
             ...
         C = A = D = E;
      */

      Matrix::iterator it;
      for (it = A.begin(), i = 0; it != A.end(); ++i,++it) {
        CGAL_TEST(A( i/mat_dim, i%mat_dim ) == *it){}
      }

      CGAL_TEST(A==C){}
      CGAL_TEST(A!=B){}
      CGAL_TEST(A==D){}
      CGAL_TEST(A.row_dimension()==mat_dim && A.column_dimension()==mat_dim){}
      CGAL_TEST(A.row(1)*B.column(1)==NT(mat_dim)){}

      /* some basic arithmetic testing */
      C = A; C += A;
      C -= NT(3)*A;
      C = -C;
      CGAL_TEST(C==A){}

      /* row sum test: */
      Vector ones(mat_dim), row_sum_vec(mat_dim);
      for (i = 0; i < mat_dim; i++) {
        ones[i] = 1;
        row_sum_vec[i] = i*mat_dim;
      }
      CGAL_TEST(A*ones==row_sum_vec){}

      C = A+B;
      /* matrix operations + ,* and |transpose|, |rank| */
      CGAL_TEST(C==LA::transpose(C)){}
      C = C-B;
      CGAL_TEST(C==A){}
      C = I+A;
      CGAL_TEST(LA::rank(C)==mat_dim){}

      /* matrix operations 2* , |determinant| */
      C = NT(2) * I;
      C = C * NT(2);
      Matrix L,U;
      Vector c;
      std::vector<int> q;
      NT det = LA::determinant(C, L, U, q, c); // det must be 2^{2*mat_dim}
      NT pot = 1; for (i=0; i<mat_dim; ++i) pot *= NT(4);
      CGAL_TEST(det == pot){}
      CGAL_TEST(LA::verify_determinant(C, det, L, U, q, c)){}
      CGAL_TEST(det == LA::determinant(C)){}
      CGAL_TEST(CGAL_NTS sign(det) == LA::sign_of_determinant(C)){}
      if (IOTEST) CGAL_IO_TEST(A,C,CGAL::IO::ASCII);
      /* a random linear solver task: */
      Vector b(mat_dim),x(mat_dim),e;
      NT denom;
      for (i = 0; i < mat_dim; i++) {
        for (j = 0; j < mat_dim; j++)
          C(i,j) = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
        b[i] = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
      }

      for (double f = 1.0; f > 0.0; f-=0.1) {
        Matrix E = C;
        for (i = 0; i< mat_dim; ++i)
          for (j = 0; j < mat_dim; ++j)
            if ( CGAL::get_default_random().get_double() > f ) E(i,j)=0;

        if (LA::linear_solver(E, b, x, denom, A, e)) {
          CGAL_TEST(E*x == b*denom){}
          LA::linear_solver(E,b,x,denom,e);
          CGAL_TEST(E*x == b*denom){}
          LA::linear_solver(E,b,x,denom);
          CGAL_TEST(E*x == b*denom){}
          CGAL_TEST(LA::is_solvable(E,b)){}
        }
        else {
          Vector null(mat_dim);
          CGAL_TEST(LA::transpose(E)*e==null && b*e!=NT(0)){}
          CGAL_TEST(!LA::is_solvable(E,b)){}
        }

        Matrix SV;
        if (LA::homogeneous_linear_solver(E,x)) {
          CGAL_TEST(E*x==Vector(mat_dim)){}
          int r = LA::homogeneous_linear_solver(E,SV);
          CGAL_TEST(r==LA::rank(SV)){}
        }

        // inverse
        if (LA::inverse(E,D,denom,c)) {
          CGAL_TEST(E*D==denom*Matrix(mat_dim,Matrix::Identity())){}
          CGAL_TEST(D==LA::inverse(E,denom)){}
        } else {
          CGAL_TEST(LA::transpose(E)*c==Vector(mat_dim)){}
        }
        LA::independent_columns(E,q);
      }

      if (mat_dim > 1) { // ueberbestimmt:
        Matrix N(mat_dim,mat_dim-1);
        for (i = 0; i < mat_dim; i++) {
          for (j = 0; j < mat_dim-1; j++)
            N(i,j) = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
          b[i] = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
        }
        if (LA::linear_solver(N,b,x,denom,e)) {
          CGAL_TEST(N*x == b*denom){}
        } else {
          Vector null(mat_dim-1);
          CGAL_TEST(LA::transpose(N)*e==null && b*e!=NT(0)){}
        }
      }

    }


  }
  {
    typedef FT NT;
    typedef CGAL::Linear_algebraCd<FT> LA;
    typedef LA::Matrix Matrix;
    typedef LA::Vector Vector;
    bool IOTEST = false;
    {
      int vec_dim;
      if (argc == 2) vec_dim = std::atoi(argv[1]);
      else           vec_dim = VEC_DIM;

      /* some construction and access ops */

      Vector v0(vec_dim), v1(vec_dim), v2(vec_dim);
      int F[] = { 1,2,3,4,5 };
      Vector v11(F,F+2), v12(F,F+3), v13(F,F+4), v14(F,F+5),
                 v15(v13.begin(),v13.end()),
                 v16(vec_dim,NT(1));
      CGAL_TEST(v13==v15){}
      CGAL_TEST(v0==v1){}
      for (int i = 0; i < vec_dim; i++) {
        v1[i] = i;
        v2[i] = vec_dim-i;
      }
      CGAL_TEST(v0!=v1){}
      Vector v3(v1);
      CGAL_TEST(v3==v1){}

      v1 += NT(2)*v2;
      v1 -= v2;
      v1 = -v1;
      CGAL_TEST((v1*v1 == NT(vec_dim*vec_dim*vec_dim))){}
      /* squared length of (v1 + v2) = v1*v1
         should be equal to dim*dim*dim */
      NT res1 = (v1*v1 - v2*v2);
      NT res2 = ((v1-v2)*(v1+v2));
      CGAL_TEST(res1==res2){}
      /* we test the third binomial formula */
      Vector v21(v13),v22(v13);
      v21 *= 13;
      CGAL_TEST(v21 == NT(13)*v22){}
      v21 /= 13;
      CGAL_TEST(v21 == v22){}

      if (IOTEST) CGAL_IO_TEST(v1,v2,CGAL::IO::ASCII);
    }

    {
      int mat_dim;
      if (argc == 2) mat_dim = std::atoi(argv[1]);
      else           mat_dim = MAT_DIM;

      /* some construction and access ops */
      Matrix::Identity ID;
      Matrix A(mat_dim,mat_dim), B(mat_dim),
             I(mat_dim,ID), One(mat_dim,mat_dim,NT(2));
      CGAL_TEST(A==B){}
      CGAL_TEST(One*I==One){}
      std::vector<Vector> F(mat_dim);

      int i,j;
      for (i = 0; i < mat_dim; i++) {
        Vector v(mat_dim);
        for (j = 0; j < mat_dim; j++) {
          A(i,j) = i;
          B(i,j) = j;
          v[j] = j;
        }
        F[i] = v;
      }
    #ifndef CGAL_SIMPLE_INTERFACE
      Matrix C(F), D(A), E(F.begin(),F.end());
    #else
      Matrix C(F), D(A), E(F);
    #endif
      CGAL_TEST(C==F && E==F){}
      /* A = (0,0,0,0...)
             (1,1,1,1...)
             (2,2,2,2...)
             ...
         B = (0,1,2,3...)
             (0,1,2,3...)
             (0,1,2,3...)
             ...
         C = A = D = E;
      */

      Matrix::iterator it;
      for (it = A.begin(), i = 0; it != A.end(); ++i,++it) {
        CGAL_TEST(A( i/mat_dim, i%mat_dim ) == *it){}
      }

      CGAL_TEST(A==C){}
      CGAL_TEST(A!=B){}
      CGAL_TEST(A==D){}
      CGAL_TEST(A.row_dimension()==mat_dim && A.column_dimension()==mat_dim){}
      CGAL_TEST(A.row(1)*B.column(1)==NT(mat_dim)){}

      /* some basic arithmetic testing */
      C = A; C += A;
      C -= NT(3)*A;
      C = -C;
      CGAL_TEST(C==A){}

      /* row sum test: */
      Vector ones(mat_dim), row_sum_vec(mat_dim);
      for (i = 0; i < mat_dim; i++) {
        ones[i] = 1;
        row_sum_vec[i] = i*mat_dim;
      }
      CGAL_TEST(A*ones==row_sum_vec){}

      C = A+B;
      /* matrix operations + ,* and |transpose|, |rank| */
      CGAL_TEST(C==LA::transpose(C)){}
      C = C-B;
      CGAL_TEST(C==A){}
      C = I+A;
      CGAL_TEST(LA::rank(C)==mat_dim){}

      /* matrix operations 2* , |determinant| */
      C = NT(2) * I;
      C = C * NT(2);
      Matrix L,U;
      Vector c;
      std::vector<int> q;
      NT det = LA::determinant(C, L, U, q, c); // det must be 2^{2*mat_dim}
      NT pot = 1; for (i=0; i<mat_dim; ++i) pot *= NT(4);
      CGAL_TEST(det == pot){}
      CGAL_TEST(LA::verify_determinant(C, det, L, U, q, c)){}
      CGAL_TEST(det == LA::determinant(C)){}
      CGAL_TEST(CGAL_NTS sign(det) == LA::sign_of_determinant(C)){}
      if (IOTEST) CGAL_IO_TEST(A,C,CGAL::IO::ASCII);
      // add binary test later: if (IOTEST) CGAL_IO_TEST(A,C,CGAL::IO::BINARY);
      /* a random linear solver task: */
      Vector b(mat_dim),x(mat_dim),e;
      NT denom;
      for (i = 0; i < mat_dim; i++) {
        for (j = 0; j < mat_dim; j++)
          C(i,j) = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
        b[i] = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
      }

      for (double f = 1.0; f > 0.0; f-=0.1) {
        Matrix E = C;
        for (i = 0; i< mat_dim; ++i)
          for (j = 0; j < mat_dim; ++j)
            if ( CGAL::get_default_random().get_double() > f ) E(i,j)=0;

        if (LA::linear_solver(E, b, x, denom, A, e)) {
          CGAL_TEST(E*x == b*denom){}
          LA::linear_solver(E,b,x,denom,e);
          CGAL_TEST(E*x == b*denom){}
          LA::linear_solver(E,b,x,denom);
          CGAL_TEST(E*x == b*denom){}
          CGAL_TEST(LA::is_solvable(E,b)){}
        }
        else {
          Vector null(mat_dim);
          CGAL_TEST(LA::transpose(E)*e==null && b*e!=NT(0)){}
          CGAL_TEST(!LA::is_solvable(E,b)){}
        }
        Matrix SV;
        if (LA::homogeneous_linear_solver(E,x)) {
          CGAL_TEST(E*x==Vector(mat_dim)){}
          int r = LA::homogeneous_linear_solver(E,SV);
          CGAL_TEST(r==LA::rank(SV)){}
        }
        // inverse
        if (LA::inverse(E,D,denom,c)) {
          CGAL_TEST(E*D==denom*Matrix(mat_dim,Matrix::Identity())){}
          CGAL_TEST(D==LA::inverse(E,denom)){}
        } else {
          CGAL_TEST(LA::transpose(E)*c==Vector(mat_dim)){}
        }
        LA::independent_columns(E,q);
      }

      if (mat_dim > 1) { // ueberbestimmt:
        Matrix N(mat_dim,mat_dim-1);
        for (i = 0; i < mat_dim; i++) {
          for (j = 0; j < mat_dim-1; j++)
            N(i,j) = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
          b[i] = CGAL::get_default_random().get_int(-mat_dim,mat_dim);
        }
        if (LA::linear_solver(N,b,x,denom,e)) {
          CGAL_TEST(N*x == b*denom){}
        } else {
          Vector null(mat_dim-1);
          CGAL_TEST(LA::transpose(N)*e==null && b*e!=NT(0)){}
        }
      }

    }


  }
  CGAL_TEST_END;
}


