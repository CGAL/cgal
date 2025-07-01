double determinant(
 double a00,  double a01,  double a02,  double a03,
 double a10,  double a11,  double a12,  double a13,
 double a20,  double a21,  double a22,  double a23,
 double a30,  double a31,  double a32,  double a33)
{
  double m01 = a10*a01 - a00*a11;
  double m02 = a20*a01 - a00*a21;
  double m03 = a30*a01 - a00*a31;
  double m12 = a20*a11 - a10*a21;
  double m13 = a30*a11 - a10*a31;
  double m23 = a30*a21 - a20*a31;

  double m012 = m12*a02 - m02*a12 + m01*a22;
  double m013 = m13*a02 - m03*a12 + m01*a32;
  double m023 = m23*a02 - m03*a22 + m02*a32;
  double m123 = m23*a12 - m13*a22 + m12*a32;

  double m0123 = m123*a03 - m023*a13 + m013*a23 - m012*a33;
  return m0123;
}

int orientationC4(double p0, double p1, double p2, double p3,
                  double q0, double q1, double q2, double q3,
                  double r0, double r1, double r2, double r3,
                  double s0, double s1, double s2, double s3,
                  double t0, double t1, double t2, double t3)
group p0 q0 r0 t0 q0;
group p1 q1 r1 t1 q1;
group p2 q2 r2 t2 q2;
group p3 q3 r3 t3 q3;
{

    double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);

    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);

    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);

    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);

    double det = determinant(m01, m02, m03, m04,
                             m11, m12, m13, m14,
                             m21, m22, m23, m24,
                             m31, m32, m33, m34);
    return sign(det);


}

//===========generated ==================

inline int orientationC4( double p0, double p1, double p2, double p3, double q0, double q1, double q2, double q3, double r0, double r1, double r2, double r3, double s0, double s1, double s2, double s3, double t0, double t1, double t2, double t3) {
    double m01;
    m01 = (q0 - p0);
    double m02;
    m02 = (r0 - p0);
    double m03;
    m03 = (s0 - p0);
    double m04;
    m04 = (t0 - p0);
    double m11;
    m11 = (q1 - p1);
    double m12;
    m12 = (r1 - p1);
    double m13;
    m13 = (s1 - p1);
    double m14;
    m14 = (t1 - p1);
    double m21;
    m21 = (q2 - p2);
    double m22;
    m22 = (r2 - p2);
    double m23;
    m23 = (s2 - p2);
    double m24;
    m24 = (t2 - p2);
    double m31;
    m31 = (q3 - p3);
    double m32;
    m32 = (r3 - p3);
    double m33;
    m33 = (s3 - p3);
    double m34;
    m34 = (t3 - p3);
    double det;
    det = determinant( m01, m02, m03, m04, m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34 );
    int int_tmp_result;
    double eps;
    double max1 = fabs(m01);
    if( (max1 < fabs(m02)) )
    {
        max1 = fabs(m02);
    }
    if( (max1 < fabs(m03)) )
    {
        max1 = fabs(m03);
    }
    if( (max1 < fabs(m11)) )
    {
        max1 = fabs(m11);
    }
    if( (max1 < fabs(m12)) )
    {
        max1 = fabs(m12);
    }
    if( (max1 < fabs(m13)) )
    {
        max1 = fabs(m13);
    }
    if( (max1 < fabs(m23)) )
    {
        max1 = fabs(m23);
    }
    double max2 = fabs(m01);
    if( (max2 < fabs(m02)) )
    {
        max2 = fabs(m02);
    }
    if( (max2 < fabs(m11)) )
    {
        max2 = fabs(m11);
    }
    if( (max2 < fabs(m12)) )
    {
        max2 = fabs(m12);
    }
    if( (max2 < fabs(m21)) )
    {
        max2 = fabs(m21);
    }
    if( (max2 < fabs(m22)) )
    {
        max2 = fabs(m22);
    }
    if( (max2 < fabs(m23)) )
    {
        max2 = fabs(m23);
    }
    if( (max2 < fabs(m33)) )
    {
        max2 = fabs(m33);
    }
    double max3 = fabs(m04);
    if( (max3 < fabs(m14)) )
    {
        max3 = fabs(m14);
    }
    if( (max3 < fabs(m24)) )
    {
        max3 = fabs(m24);
    }
    if( (max3 < fabs(m34)) )
    {
        max3 = fabs(m34);
    }
    double max4 = fabs(m11);
    if( (max4 < fabs(m12)) )
    {
        max4 = fabs(m12);
    }
    if( (max4 < fabs(m21)) )
    {
        max4 = fabs(m21);
    }
    if( (max4 < fabs(m22)) )
    {
        max4 = fabs(m22);
    }
    if( (max4 < fabs(m31)) )
    {
        max4 = fabs(m31);
    }
    if( (max4 < fabs(m32)) )
    {
        max4 = fabs(m32);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
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
    if( (lower_bound_1 < 2.89273249588395272840e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    }
    else
    {
        if( (upper_bound_1 > 7.23700557733225900010e+75) )
        {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (3.17768858673611327578e-14 * (((max4 * max2) * max1) * max3));
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

