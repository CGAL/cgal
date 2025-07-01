double
determinant(
 double a00, double a01, double a02, double a03, double a04, double a05,
 double a10, double a11, double a12, double a13, double a14, double a15,
 double a20, double a21, double a22, double a23, double a24, double a25,
 double a30, double a31, double a32, double a33, double a34, double a35,
 double a40, double a41, double a42, double a43, double a44, double a45,
 double a50, double a51, double a52, double a53, double a54, double a55)
{
  double m01 = a00*a11 - a10*a01;
  double m02 = a00*a21 - a20*a01;
  double m03 = a00*a31 - a30*a01;
  double m04 = a00*a41 - a40*a01;
  double m05 = a00*a51 - a50*a01;
  double m12 = a10*a21 - a20*a11;
  double m13 = a10*a31 - a30*a11;
  double m14 = a10*a41 - a40*a11;
  double m15 = a10*a51 - a50*a11;
  double m23 = a20*a31 - a30*a21;
  double m24 = a20*a41 - a40*a21;
  double m25 = a20*a51 - a50*a21;
  double m34 = a30*a41 - a40*a31;
  double m35 = a30*a51 - a50*a31;
  double m45 = a40*a51 - a50*a41;

  double m012 = m01*a22 - m02*a12 + m12*a02;
  double m013 = m01*a32 - m03*a12 + m13*a02;
  double m014 = m01*a42 - m04*a12 + m14*a02;
  double m015 = m01*a52 - m05*a12 + m15*a02;
  double m023 = m02*a32 - m03*a22 + m23*a02;
  double m024 = m02*a42 - m04*a22 + m24*a02;
  double m025 = m02*a52 - m05*a22 + m25*a02;
  double m034 = m03*a42 - m04*a32 + m34*a02;
  double m035 = m03*a52 - m05*a32 + m35*a02;
  double m045 = m04*a52 - m05*a42 + m45*a02;
  double m123 = m12*a32 - m13*a22 + m23*a12;
  double m124 = m12*a42 - m14*a22 + m24*a12;
  double m125 = m12*a52 - m15*a22 + m25*a12;
  double m134 = m13*a42 - m14*a32 + m34*a12;
  double m135 = m13*a52 - m15*a32 + m35*a12;
  double m145 = m14*a52 - m15*a42 + m45*a12;
  double m234 = m23*a42 - m24*a32 + m34*a22;
  double m235 = m23*a52 - m25*a32 + m35*a22;
  double m245 = m24*a52 - m25*a42 + m45*a22;
  double m345 = m34*a52 - m35*a42 + m45*a32;

  double m0123 = m012*a33 - m013*a23 + m023*a13 - m123*a03;
  double m0124 = m012*a43 - m014*a23 + m024*a13 - m124*a03;
  double m0125 = m012*a53 - m015*a23 + m025*a13 - m125*a03;
  double m0134 = m013*a43 - m014*a33 + m034*a13 - m134*a03;
  double m0135 = m013*a53 - m015*a33 + m035*a13 - m135*a03;
  double m0145 = m014*a53 - m015*a43 + m045*a13 - m145*a03;
  double m0234 = m023*a43 - m024*a33 + m034*a23 - m234*a03;
  double m0235 = m023*a53 - m025*a33 + m035*a23 - m235*a03;
  double m0245 = m024*a53 - m025*a43 + m045*a23 - m245*a03;
  double m0345 = m034*a53 - m035*a43 + m045*a33 - m345*a03;
  double m1234 = m123*a43 - m124*a33 + m134*a23 - m234*a13;
  double m1235 = m123*a53 - m125*a33 + m135*a23 - m235*a13;
  double m1245 = m124*a53 - m125*a43 + m145*a23 - m245*a13;
  double m1345 = m134*a53 - m135*a43 + m145*a33 - m345*a13;
  double m2345 = m234*a53 - m235*a43 + m245*a33 - m345*a23;

  double m01234 = m0123*a44 - m0124*a34 + m0134*a24 - m0234*a14 + m1234*a04;
  double m01235 = m0123*a54 - m0125*a34 + m0135*a24 - m0235*a14 + m1235*a04;
  double m01245 = m0124*a54 - m0125*a44 + m0145*a24 - m0245*a14 + m1245*a04;
  double m01345 = m0134*a54 - m0135*a44 + m0145*a34 - m0345*a14 + m1345*a04;
  double m02345 = m0234*a54 - m0235*a44 + m0245*a34 - m0345*a24 + m2345*a04;
  double m12345 = m1234*a54 - m1235*a44 + m1245*a34 - m1345*a24 + m2345*a14;

  double m012345 = m01234*a55 - m01235*a45 + m01245*a35 - m01345*a25
                   + m02345*a15 - m12345*a05;
  return m012345;
}

int orientationC5(double p0, double p1, double p2, double p3, double p4, double p5,
                  double q0, double q1, double q2, double q3, double q4, double q5,
                  double r0, double r1, double r2, double r3, double r4, double r5,
                  double s0, double s1, double s2, double s3, double s4, double s5,
                  double t0, double t1, double t2, double t3, double t4, double t5,
                  double u0, double u1, double u2, double u3, double u4, double u5,
                  double v0, double v1, double v2, double v3, double v4, double v5)
group p0 q0 r0 t0 q0 u0 v0;
group p1 q1 r1 t1 q1 u1 v1;
group p2 q2 r2 t2 q2 u2 v2;
group p3 q3 r3 t3 q3 u3 v3;
group p4 q4 r4 t4 q4 u4 v4;
group p5 q5 r5 t5 q5 u5 v5;
{
    double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m05;
    m05 = (u0 - p0);
    double m06;
    m06 = (v0 - p0);

    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m15;
    m15 = (u1 - p1);
    double m16;
    m16 = (v1 - p1);

    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m25;
    m25 = (u2 - p2);
    double m26;
    m26 = (v2 - p2);

    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double m35;
    m35 = (u3 - p3);
    double m36;
    m36 = (v3 - p3);

    double m41;
    m41 = (q4 - p4);
    double m42;
    m42 = (r4 - p4);
    double m43;
    m43 = (s4 - p4);
    double m44;
    m44 = (t4 - p4);
    double m45;
    m45 = (u4 - p4);
    double m46;
    m46 = (v4 - p4);

    double m51;
    m51 = (q5 - p5);
    double m52;
    m52 = (r5 - p5);
    double m53;
    m53 = (s5 - p5);
    double m54;
    m54 = (t5 - p5);
    double m55;
    m55 = (u5 - p5);
    double m56;
    m56 = (v5 - p5);

    double det = determinant(m01, m02, m03, m04, m05, m06,
                             m11, m12, m13, m14, m15, m16,
                             m21, m22, m23, m24, m25, m26,
                             m31, m32, m33, m34, m35, m36,
                             m41, m42, m43, m44, m45, m46,
                             m51, m52, m53, m54, m55, m56);
    return sign(det);
}

//===========generated ==================

inline int orientationC5( double p0, double p1, double p2, double p3, double p4, double p5, double q0, double q1, double q2, double q3, double q4, double q5, double r0, double r1, double r2, double r3, double r4, double r5, double s0, double s1, double s2, double s3, double s4, double s5, double t0, double t1, double t2, double t3, double t4, double t5, double u0, double u1, double u2, double u3, double u4, double u5, double v0, double v1, double v2, double v3, double v4, double v5) {
    double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m05;
    m05 = (u0 - p0);
    double m06;
    m06 = (v0 - p0);
    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m15;
    m15 = (u1 - p1);
    double m16;
    m16 = (v1 - p1);
    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m25;
    m25 = (u2 - p2);
    double m26;
    m26 = (v2 - p2);
    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double m35;
    m35 = (u3 - p3);
    double m36;
    m36 = (v3 - p3);
    double m41;
    m41 = (q4 - p4);
    double m42;
    m42 = (r4 - p4);
    double m43;
    m43 = (s4 - p4);
    double m44;
    m44 = (t4 - p4);
    double m45;
    m45 = (u4 - p4);
    double m46;
    m46 = (v4 - p4);
    double m51;
    m51 = (q5 - p5);
    double m52;
    m52 = (r5 - p5);
    double m53;
    m53 = (s5 - p5);
    double m54;
    m54 = (t5 - p5);
    double m55;
    m55 = (u5 - p5);
    double m56;
    m56 = (v5 - p5);
    double det;
    det = determinant( m01, m02, m03, m04, m05, m06, m11, m12, m13, m14, m15, m16, m21, m22, m23, m24, m25, m26, m31, m32, m33, m34, m35, m36, m41, m42, m43, m44, m45, m46, m51, m52, m53, m54, m55, m56 );
    int int_tmp_result;
    double eps;
    double max1 = CGAL::abs(m01);
    double am = CGAL::abs(m02);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m11);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m12);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m21);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m22);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m31);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m32);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m41);
    if( (max1 < am) ) { max1 = am; }
    am = CGAL::abs(m42);
    if( (max1 < am) ) { max1 = am; }


    double max2 = CGAL::abs(m03);
    am = CGAL::abs(m11);

    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m12);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m13);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m21);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m22);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m23);

    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m31);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m32);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m33);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m41);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m42);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m43);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m51);
    if( (max2 < am) ) { max2 = am; }
    am = CGAL::abs(m52);
    if( (max2 < am) ) { max2 = am; }


    double max3 = CGAL::abs(m04);
    am = CGAL::abs(m14);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m24);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m34);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m44);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m54);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m05);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m15);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m25);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m35);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m45);
    if( (max3 < am) ) { max3 = am; }
    am = CGAL::abs(m55);
    if( (max3 < am) ) { max3 = am; }

    d
    double max5 = CGAL::abs(m06);
    am = CGAL::abs(m16);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m26);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m36);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m46);
    if( (max5 < am) ) { max5 = am; }
    am = CGAL::abs(m56);
    if( (max5 < am) ) { max5 = am; }


    double max6 = CGAL::abs(m13);
    am = CGAL::abs(m21);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m22);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m23);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m31);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m32);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m33);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m41);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m42);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m43);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m51);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m52);
    if( (max6 < am) ) { max6 = am; }
    am = CGAL::abs(m53);
    if( (max6 < am) ) { max6 = am; }


    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max6;
    upper_bound_1 = max6;
    if( (max5 < lower_bound_1) )
    {
        lower_bound_1 = max5;
    }
    else
    {
        if( (max5 > upper_bound_1) )
        {
            upper_bound_1 = max5;
        }
    }
    if( (max1 < lower_bound_1) )
    {
        lower_bound_1 = max1;
    }
    else
    {
        if( (max1 > upper_bound_1) )
        {
            upper_bound_1 = max1;
        }
    }
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    }
    else
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        }
    }
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    }
    else
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        }
    }
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    }
    else
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        }
    }
    if( (lower_bound_1 < 4.82472686053427214432e-50) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 1.87072209578355511223e+50) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.76403842114300859158e-12 * (((((max1 * max2) * max6) * max3) * max4) * max5));
        if( (det > eps) )
        {
            int_tmp_result = 1;
        }
        else
        {
            if( (det < -eps) )
            {
                int_tmp_result = -1;
            }
            else
            {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

