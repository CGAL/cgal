#include <CGAL/basic.h>
#include <CGAL/Linear_algebra.h>
#include <CGAL/random_selection.h>
#include <CGAL/test_macros.h>
#define VEC_DIM 10
#define MAT_DIM 10

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

typedef CGAL::Linear_algebra<RT> LA;
typedef LA::Matrix Matrix;
typedef LA::Vector Vector;

int main(int argc, char* argv[])
{
  CGAL_TEST_START;
  { 
    int vec_dim;
    if (argc == 2) vec_dim = atoi(argv[1]);
    else           vec_dim = VEC_DIM;

    /* some construction and access ops */

    Vector v0(vec_dim), v1(vec_dim), v2(vec_dim);
    int F[] = { 1,2,3,4,5 };
    Vector v11(1,2), v12(1,2,3), v13(1,2,3,4), v14(F,F+5),
               v15(v13.begin(),v13.end()), 
               v16(vec_dim,Vector::Initialize(),1);
    CGAL_TEST(v13==v15);
    CGAL_TEST(v0==v1);
    for (int i = 0; i < vec_dim; i++) {
      v1[i] = i;
      v2[i] = vec_dim-i;
    }
    CGAL_TEST(v0!=v1);
    Vector v3(v1);
    CGAL_TEST(v3==v1);

    v1 += 2*v2;
    v1 -= v2;
    v1 = -v1;
    CGAL_TEST((v1*v1 == vec_dim*vec_dim*vec_dim));
    /* squared length of (v1 + v2) = v1*v1 
       should be equal to dim*dim*dim */
    RT res1 = (v1*v1 - v2*v2);
    RT res2 = ((v1-v2)*(v1+v2));
    CGAL_TEST(res1==res2);
    /* we test the third binomial formula */

    CGAL_IO_TEST(v1,v2);
  }

  { 
    int mat_dim;
    if (argc == 2) mat_dim = atoi(argv[1]);
    else           mat_dim = MAT_DIM;

    /* some construction and access ops */
    Matrix A(mat_dim,mat_dim), B(mat_dim),
           I(mat_dim,Matrix::Identity()),
           One(mat_dim,mat_dim,Matrix::Initialize(),2); 
    CGAL_TEST(A==B);
    CGAL_TEST(One*I==One);
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
    Matrix C(F), D(A), E(F.begin(),F.end());
    CGAL_TEST(C==F && E==F);
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

    CGAL_TEST(A==C);
    CGAL_TEST(A!=B);
    CGAL_TEST(A==D);
    CGAL_TEST(A.row_dimension()==mat_dim && A.column_dimension()==mat_dim);
    CGAL_TEST(A.row(1)*B.column(1)==mat_dim);
    CGAL_IO_TEST(A,C);

    /* some basic arithmetic testing */
    C += A; 
    C -= 3*A; 
    C = -C; 
    CGAL_TEST(C==A);

    /* row sum test: */
    Vector ones(mat_dim), row_sum_vec(mat_dim); 
    for (i = 0; i < mat_dim; i++) {
      ones[i] = 1; 
      row_sum_vec[i] = i*mat_dim;
    }
    CGAL_TEST(A*ones==row_sum_vec);

    C = A+B;
    /* matrix operations + ,* and |transpose|, |rank| */
    CGAL_TEST(C==LA::transpose(C));
    C = C-B;
    CGAL_TEST(C==A);

    C = I+A;
    CGAL_TEST(LA::rank(C)==mat_dim);
    
    /* matrix operations 2* , |determinant| */
    C = 2 * I;
    C = C * 2; 
    Matrix L,U; 
    Vector c; 
    std::vector<int> q; 
    RT det = LA::determinant(C, L, U, q, c); // det must be 2^{2*mat_dim}
    CGAL_TEST(log(det)==2*mat_dim);
    CGAL_TEST(LA::verify_determinant(C, det, L, U, q, c));
    CGAL_TEST(det == LA::determinant(C));
    CGAL_TEST(sign(det) == LA::sign_of_determinant(C));
    LA::independent_columns(C,q);

    /* a random linear solver task: */
    Vector b(mat_dim),x(mat_dim),e; 
    RT denom; 
    for (i = 0; i < mat_dim; i++) { 
      for (j = 0; j < mat_dim; j++) 
        E(i,j) = CGAL::default_random.get_int( - mat_dim,mat_dim); 
      b[i] = CGAL::default_random.get_int( - mat_dim,mat_dim); 
    }

    // linear solver
    if (LA::linear_solver(E, b, x, denom, A, e)) {
      CGAL_TEST(E*x == b*denom);
      LA::linear_solver(E,b,x,denom,e);
      CGAL_TEST(E*x == b*denom);
      LA::linear_solver(E,b,x,denom);
      CGAL_TEST(E*x == b*denom);
      CGAL_TEST(LA::is_solvable(E,b));
    }
    else {
      Vector null(mat_dim);
      CGAL_TEST(LA::transpose(E)*e==null && b*e!=0);
      CGAL_TEST(!LA::is_solvable(E,b));
    }

    Matrix SV;
    if (LA::homogeneous_linear_solver(E,x)) {
      CGAL_TEST(E*x==Vector(mat_dim));
      CGAL_TEST(LA::homogeneous_linear_solver(E,SV)==LA::rank(SV));
    }

    // inverse
    if (LA::inverse(E,D,denom,c)) {
      CGAL_TEST(E*D==denom*Matrix(mat_dim,Matrix::Identity()));
      CGAL_TEST(D==LA::inverse(E,denom));
    } else {
      CGAL_TEST(LA::transpose(E)*c==Vector(mat_dim));
    }
  }


  CGAL_TEST_END;
}


