/**
 * @file   kernel/Vector2.h
 * @author Gernot Walzl
 * @date   2011-11-14
 */

#ifndef VECTOR2_H
#define VECTOR2_H

namespace kernel {

class Vector2 {
public:
    Vector2();
    Vector2(double x, double y);
    Vector2(const Vector2& orig);
    virtual ~Vector2();

    double operator[](unsigned int i) const;
    double squared_length(void) const;
    double length(void) const;
    Vector2 normalize(void) const;
    double angle(const Vector2& v) const;
    Vector2 operator+(const Vector2& v) const;
    Vector2 operator-(const Vector2& v) const;
    Vector2 operator*(double s) const;
    Vector2 operator/(double s) const;
    double operator*(const Vector2& v) const;
    bool operator==(const Vector2& v) const;

protected:
    double v_[2];
};

}

#endif /* VECTOR2_H */
