#include <CGAL/determinant.h>
#include <CGAL/Random.h>

#include <boost/lexical_cast.hpp>

#include <iostream>

using CGAL::determinant;

template <typename RT>
RT det_77_alt(const RT& a00, const RT& a01, const RT& a02, const RT& a03, const RT& a04, const RT& a05, const RT& a06,
              const RT& a10, const RT& a11, const RT& a12, const RT& a13, const RT& a14, const RT& a15, const RT& a16,
              const RT& a20, const RT& a21, const RT& a22, const RT& a23, const RT& a24, const RT& a25, const RT& a26,
              const RT& a30, const RT& a31, const RT& a32, const RT& a33, const RT& a34, const RT& a35, const RT& a36,
              const RT& a40, const RT& a41, const RT& a42, const RT& a43, const RT& a44, const RT& a45, const RT& a46,
              const RT& a50, const RT& a51, const RT& a52, const RT& a53, const RT& a54, const RT& a55, const RT& a56,
              const RT& a60, const RT& a61, const RT& a62, const RT& a63, const RT& a64, const RT& a65, const RT& a66)
{
  const RT r1 = a06 * determinant(
                  a10, a11, a12, a13, a14, a15,
                                  a20, a21, a22, a23, a24, a25,
                                  a30, a31, a32, a33, a34, a35,
                                  a40, a41, a42, a43, a44, a45,
                                  a50, a51, a52, a53, a54, a55,
                                  a60, a61, a62, a63, a64, a65);

  const RT r2 = - a16 * determinant(a00, a01, a02, a03, a04, a05,

                                    a20, a21, a22, a23, a24, a25,
                                    a30, a31, a32, a33, a34, a35,
                                    a40, a41, a42, a43, a44, a45,
                                    a50, a51, a52, a53, a54, a55,
                                    a60, a61, a62, a63, a64, a65);

  const RT r3 = a26 * determinant(a00, a01, a02, a03, a04, a05,
                                  a10, a11, a12, a13, a14, a15,

                                  a30, a31, a32, a33, a34, a35,
                                  a40, a41, a42, a43, a44, a45,
                                  a50, a51, a52, a53, a54, a55,
                                  a60, a61, a62, a63, a64, a65);

  const RT r4 = - a36 * determinant(a00, a01, a02, a03, a04, a05,
                                    a10, a11, a12, a13, a14, a15,
                                    a20, a21, a22, a23, a24, a25,

                                    a40, a41, a42, a43, a44, a45,
                                    a50, a51, a52, a53, a54, a55,
                                    a60, a61, a62, a63, a64, a65);

  const RT r5 = a46 * determinant(a00, a01, a02, a03, a04, a05,
                                  a10, a11, a12, a13, a14, a15,
                                  a20, a21, a22, a23, a24, a25,
                                  a30, a31, a32, a33, a34, a35,

                                  a50, a51, a52, a53, a54, a55,
                                  a60, a61, a62, a63, a64, a65);

  const RT r6 = - a56 * determinant(a00, a01, a02, a03, a04, a05,
                                    a10, a11, a12, a13, a14, a15,
                                    a20, a21, a22, a23, a24, a25,
                                    a30, a31, a32, a33, a34, a35,
                                    a40, a41, a42, a43, a44, a45,

                                    a60, a61, a62, a63, a64, a65);

  const RT r7 = a66 * determinant(a00, a01, a02, a03, a04, a05,
                                  a10, a11, a12, a13, a14, a15,
                                  a20, a21, a22, a23, a24, a25,
                                  a30, a31, a32, a33, a34, a35,
                                  a40, a41, a42, a43, a44, a45,
                                  a50, a51, a52, a53, a54, a55

                                  );

  return r1 + r2 + r3 + r4 + r5 + r6 + r7;
}

int main(int, char**)
{
  assert(determinant<int>(4, 5, 1, 4, 6, 3, 1,
                          4, 3, 6, 4, 2, 7, 3,
                          6, 3, 3, 6, 2, 4, 5,
                          1, 4, 3, 5, 5, 6 ,1,
                          1, 3, 2, 7, 9, 6, 1,
                          7, 6, 5, 4, 6, 2, 2,
                          2, 3, 5, 7, 4, 3, 3) == 763);

  CGAL::Random rnd;
  std::cout << "Seed: " << rnd.get_seed() << std::endl;

  for(int k=0; k<100; ++k)
  {
    std::array<std::array<int, 7>, 7> mat;
    for(int i=0; i<7; ++i)
      for(int j=0; j<7; ++j)
        mat[i][j] = rnd.get_int(-18, 18);

    const int det_1 = determinant<int>(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[0][4], mat[0][5], mat[0][6],
                                       mat[1][0], mat[1][1], mat[1][2], mat[1][3], mat[1][4], mat[1][5], mat[1][6],
                                       mat[2][0], mat[2][1], mat[2][2], mat[2][3], mat[2][4], mat[2][5], mat[2][6],
                                       mat[3][0], mat[3][1], mat[3][2], mat[3][3], mat[3][4], mat[3][5], mat[3][6],
                                       mat[4][0], mat[4][1], mat[4][2], mat[4][3], mat[4][4], mat[4][5], mat[4][6],
                                       mat[0][0], mat[5][1], mat[5][2], mat[5][3], mat[5][4], mat[5][5], mat[5][6],
                                       mat[6][0], mat[6][1], mat[6][2], mat[6][3], mat[6][4], mat[6][5], mat[6][6]);

    const int det_2 = det_77_alt<int>(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[0][4], mat[0][5], mat[0][6],
                                      mat[1][0], mat[1][1], mat[1][2], mat[1][3], mat[1][4], mat[1][5], mat[1][6],
                                      mat[2][0], mat[2][1], mat[2][2], mat[2][3], mat[2][4], mat[2][5], mat[2][6],
                                      mat[3][0], mat[3][1], mat[3][2], mat[3][3], mat[3][4], mat[3][5], mat[3][6],
                                      mat[4][0], mat[4][1], mat[4][2], mat[4][3], mat[4][4], mat[4][5], mat[4][6],
                                      mat[0][0], mat[5][1], mat[5][2], mat[5][3], mat[5][4], mat[5][5], mat[5][6],
                                      mat[6][0], mat[6][1], mat[6][2], mat[6][3], mat[6][4], mat[6][5], mat[6][6]);

    std::cout << "dets: " << det_1 << " " << det_2 << std::endl;
    assert(det_1 == det_2);
  }

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
