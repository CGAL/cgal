inline
Sign
sign_of_determinant2x2_SAF(
    const Static_filter_error &a00,
    const Static_filter_error &a01,
    const Static_filter_error &a10,
    const Static_filter_error &a11,
    double & epsilon_0)
{
  typedef Static_filter_error FT;
 return static_cast<Sign>(static_cast<int>(CGAL::compare_SAF( a00*a11, a10*a01,
		epsilon_0))); }

inline
Sign
sign_of_determinant2x2_SAF(
    const Restricted_double &a00,
    const Restricted_double &a01,
    const Restricted_double &a10,
    const Restricted_double &a11,
    const double & epsilon_0)
{
  typedef Restricted_double FT;
 return static_cast<Sign>(static_cast<int>(CGAL::compare_SAF( a00*a11, a10*a01,
		epsilon_0))); }

inline
Sign
sign_of_determinant3x3_SAF(
    const Static_filter_error &a00,
    const Static_filter_error &a01,
    const Static_filter_error &a02,
    const Static_filter_error &a10,
    const Static_filter_error &a11,
    const Static_filter_error &a12,
    const Static_filter_error &a20,
    const Static_filter_error &a21,
    const Static_filter_error &a22,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::sign_SAF(det3x3_by_formula(a00, a01, a02,
                                      a10, a11, a12,
                                      a20, a21, a22),
		epsilon_0);
}

inline
Sign
sign_of_determinant3x3_SAF(
    const Restricted_double &a00,
    const Restricted_double &a01,
    const Restricted_double &a02,
    const Restricted_double &a10,
    const Restricted_double &a11,
    const Restricted_double &a12,
    const Restricted_double &a20,
    const Restricted_double &a21,
    const Restricted_double &a22,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::sign_SAF(det3x3_by_formula(a00, a01, a02,
                                      a10, a11, a12,
                                      a20, a21, a22),
		epsilon_0);
}

inline
Sign
sign_of_determinant4x4_SAF(
    const Static_filter_error &a00,
    const Static_filter_error &a01,
    const Static_filter_error &a02,
    const Static_filter_error &a03,
    const Static_filter_error &a10,
    const Static_filter_error &a11,
    const Static_filter_error &a12,
    const Static_filter_error &a13,
    const Static_filter_error &a20,
    const Static_filter_error &a21,
    const Static_filter_error &a22,
    const Static_filter_error &a23,
    const Static_filter_error &a30,
    const Static_filter_error &a31,
    const Static_filter_error &a32,
    const Static_filter_error &a33,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::sign_SAF(det4x4_by_formula(a00, a01, a02, a03,
                                      a10, a11, a12, a13,
                                      a20, a21, a22, a23,
                                      a30, a31, a32, a33),
		epsilon_0);
}

inline
Sign
sign_of_determinant4x4_SAF(
    const Restricted_double &a00,
    const Restricted_double &a01,
    const Restricted_double &a02,
    const Restricted_double &a03,
    const Restricted_double &a10,
    const Restricted_double &a11,
    const Restricted_double &a12,
    const Restricted_double &a13,
    const Restricted_double &a20,
    const Restricted_double &a21,
    const Restricted_double &a22,
    const Restricted_double &a23,
    const Restricted_double &a30,
    const Restricted_double &a31,
    const Restricted_double &a32,
    const Restricted_double &a33,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::sign_SAF(det4x4_by_formula(a00, a01, a02, a03,
                                      a10, a11, a12, a13,
                                      a20, a21, a22, a23,
                                      a30, a31, a32, a33),
		epsilon_0);
}

inline
Sign
sign_of_determinant5x5_SAF(
    const Static_filter_error &a00,
    const Static_filter_error &a01,
    const Static_filter_error &a02,
    const Static_filter_error &a03,
    const Static_filter_error &a04,
    const Static_filter_error &a10,
    const Static_filter_error &a11,
    const Static_filter_error &a12,
    const Static_filter_error &a13,
    const Static_filter_error &a14,
    const Static_filter_error &a20,
    const Static_filter_error &a21,
    const Static_filter_error &a22,
    const Static_filter_error &a23,
    const Static_filter_error &a24,
    const Static_filter_error &a30,
    const Static_filter_error &a31,
    const Static_filter_error &a32,
    const Static_filter_error &a33,
    const Static_filter_error &a34,
    const Static_filter_error &a40,
    const Static_filter_error &a41,
    const Static_filter_error &a42,
    const Static_filter_error &a43,
    const Static_filter_error &a44,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::sign_SAF(det5x5_by_formula(a00, a01, a02, a03, a04,
                                      a10, a11, a12, a13, a14,
                                      a20, a21, a22, a23, a24,
                                      a30, a31, a32, a33, a34,
                                      a40, a41, a42, a43, a44),
		epsilon_0);
}

inline
Sign
sign_of_determinant5x5_SAF(
    const Restricted_double &a00,
    const Restricted_double &a01,
    const Restricted_double &a02,
    const Restricted_double &a03,
    const Restricted_double &a04,
    const Restricted_double &a10,
    const Restricted_double &a11,
    const Restricted_double &a12,
    const Restricted_double &a13,
    const Restricted_double &a14,
    const Restricted_double &a20,
    const Restricted_double &a21,
    const Restricted_double &a22,
    const Restricted_double &a23,
    const Restricted_double &a24,
    const Restricted_double &a30,
    const Restricted_double &a31,
    const Restricted_double &a32,
    const Restricted_double &a33,
    const Restricted_double &a34,
    const Restricted_double &a40,
    const Restricted_double &a41,
    const Restricted_double &a42,
    const Restricted_double &a43,
    const Restricted_double &a44,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::sign_SAF(det5x5_by_formula(a00, a01, a02, a03, a04,
                                      a10, a11, a12, a13, a14,
                                      a20, a21, a22, a23, a24,
                                      a30, a31, a32, a33, a34,
                                      a40, a41, a42, a43, a44),
		epsilon_0);
}

inline
Sign
sign_of_determinant6x6_SAF(
    const Static_filter_error &a00,
    const Static_filter_error &a01,
    const Static_filter_error &a02,
    const Static_filter_error &a03,
    const Static_filter_error &a04,
    const Static_filter_error &a05,
    const Static_filter_error &a10,
    const Static_filter_error &a11,
    const Static_filter_error &a12,
    const Static_filter_error &a13,
    const Static_filter_error &a14,
    const Static_filter_error &a15,
    const Static_filter_error &a20,
    const Static_filter_error &a21,
    const Static_filter_error &a22,
    const Static_filter_error &a23,
    const Static_filter_error &a24,
    const Static_filter_error &a25,
    const Static_filter_error &a30,
    const Static_filter_error &a31,
    const Static_filter_error &a32,
    const Static_filter_error &a33,
    const Static_filter_error &a34,
    const Static_filter_error &a35,
    const Static_filter_error &a40,
    const Static_filter_error &a41,
    const Static_filter_error &a42,
    const Static_filter_error &a43,
    const Static_filter_error &a44,
    const Static_filter_error &a45,
    const Static_filter_error &a50,
    const Static_filter_error &a51,
    const Static_filter_error &a52,
    const Static_filter_error &a53,
    const Static_filter_error &a54,
    const Static_filter_error &a55,
    double & epsilon_0)
{
  typedef Static_filter_error FT;

  return CGAL::sign_SAF(det6x6_by_formula(a00, a01, a02, a03, a04, a05,
                                      a10, a11, a12, a13, a14, a15,
                                      a20, a21, a22, a23, a24, a25,
                                      a30, a31, a32, a33, a34, a35,
                                      a40, a41, a42, a43, a44, a45,
                                      a50, a51, a52, a53, a54, a55),
		epsilon_0);
}

inline
Sign
sign_of_determinant6x6_SAF(
    const Restricted_double &a00,
    const Restricted_double &a01,
    const Restricted_double &a02,
    const Restricted_double &a03,
    const Restricted_double &a04,
    const Restricted_double &a05,
    const Restricted_double &a10,
    const Restricted_double &a11,
    const Restricted_double &a12,
    const Restricted_double &a13,
    const Restricted_double &a14,
    const Restricted_double &a15,
    const Restricted_double &a20,
    const Restricted_double &a21,
    const Restricted_double &a22,
    const Restricted_double &a23,
    const Restricted_double &a24,
    const Restricted_double &a25,
    const Restricted_double &a30,
    const Restricted_double &a31,
    const Restricted_double &a32,
    const Restricted_double &a33,
    const Restricted_double &a34,
    const Restricted_double &a35,
    const Restricted_double &a40,
    const Restricted_double &a41,
    const Restricted_double &a42,
    const Restricted_double &a43,
    const Restricted_double &a44,
    const Restricted_double &a45,
    const Restricted_double &a50,
    const Restricted_double &a51,
    const Restricted_double &a52,
    const Restricted_double &a53,
    const Restricted_double &a54,
    const Restricted_double &a55,
    const double & epsilon_0)
{
  typedef Restricted_double FT;

  return CGAL::sign_SAF(det6x6_by_formula(a00, a01, a02, a03, a04, a05,
                                      a10, a11, a12, a13, a14, a15,
                                      a20, a21, a22, a23, a24, a25,
                                      a30, a31, a32, a33, a34, a35,
                                      a40, a41, a42, a43, a44, a45,
                                      a50, a51, a52, a53, a54, a55),
		epsilon_0);
}

