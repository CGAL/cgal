/**
 * @file   kernel/Vector3.h
 * @author Gernot Walzl
 * @date   2011-11-14
 */

#ifndef VECTOR3_H
#define VECTOR3_H

namespace kernel {

class Vector3 {
public:
    Vector3();
    Vector3(double x, double y, double z);
    Vector3(const Vector3& orig);
    virtual ~Vector3();

    double operator[](unsigned int i) const;
    double squared_length(void) const;
    double length(void) const;
    Vector3 normalize(void) const;
    double angle(const Vector3& v) const;
    Vector3 operator+(const Vector3& v) const;
    Vector3 operator-(const Vector3& v) const;
    Vector3 operator*(double s) const;
    Vector3 operator/(double s) const;

    /**
     * scalar product
     */
    double operator*(const Vector3& v) const;

    /**
     * cross product
     */
    Vector3 cross(const Vector3& v) const;

    bool operator==(const Vector3& v) const;

protected:
    double v_[3];
};

}

#endif /* VECTOR3_H */

