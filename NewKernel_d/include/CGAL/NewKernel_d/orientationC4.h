


inline int orientationC4(double p0, double p1, double p2, double p3,
                         double q0, double q1, double q2, double q3,
                         double r0, double r1, double r2, double r3,
                         double s0, double s1, double s2, double s3,
                         double t0, double t1, double t2, double t3)
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
    m01 = (q1 - p1);
    double m12;
    m02 = (r1 - p1);
    double m13;
    m03 = (s1 - p1);
    double m14;
    m04 = (t1 - p1);

    double m21;
    m01 = (q2 - p2);
    double m22;
    m02 = (r2 - p2);
    double m23;
    m03 = (s2 - p2);
    double m24;
    m04 = (t2 - p2);

    double m31;
    m01 = (q3 - p3);
    double m32;
    m02 = (r3 - p3);
    double m33;
    m03 = (s3 - p3);
    double m34;
    m04 = (t3 - p3);

    double det = determinant(m01, m02, m03, m04,
                             m11, m12, m13, m14,
                             m21, m22, m23, m24,
                             m31, m32, m33, m34);
    if (det > 0) {
        return 1; // positive orientation
    } else if (det < 0) {
        return -1; // negative orientation
    } else {
        return 0; // zero orientation
    }

}



