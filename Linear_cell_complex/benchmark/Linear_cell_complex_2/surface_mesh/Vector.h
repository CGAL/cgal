//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// SPDX-License-Identifier: GPL-2.0-only
//
//=============================================================================


#ifndef SURFACE_MESH_VECTOR_H
#define SURFACE_MESH_VECTOR_H


//== INCLUDES =================================================================

#include <iostream>
#include <assert.h>
#include <math.h>
#include <limits>
#include <cstring>

#ifdef WIN32
#undef min
#undef max
#endif


/// \addtogroup geometry geometry
/// @{



//== CLASS DEFINITION =========================================================


/** A vector class for an N-dimensional vector of scalar type T.
 Elements of a vector v can be accessed by v[0], v[1], ...
 For 3D vectors one can also use v.x, v.y and v.z.
 */
template <typename Scalar, int N>
class Vector
{
public:

    /// the scalar type of the vector
    typedef Scalar value_type;

    /// returns the dimension of the vector
    static int size() { return N; }


    /// default constructor, creates uninitialized values.
    Vector() {}

    /// construct from scalar s (fills all components with s)
    explicit Vector(const Scalar s)
    {
        for (int i=0; i<N; ++i)
            data_[i] = s;
    }

    /// construct from 2 scalars. only valid for 2D vectors.
    Vector(const Scalar x, const Scalar y)
    {
        assert(N==2);
        data_[0]=x;
        data_[1]=y;
    }

    /// construct from 3 scalars. only valid for 2D vectors.
    Vector(const Scalar x, const Scalar y, const Scalar z)
    {
        assert(N==3);
        data_[0]=x;
        data_[1]=y;
        data_[2]=z;
    }

    /// construct from 4 scalars. only valid for 2D vectors.
    Vector(const Scalar x, const Scalar y, const Scalar z, const Scalar w)
    {
        assert(N==4);
        data_[0]=x;
        data_[1]=y;
        data_[2]=z;
        data_[3]=w;
    }

    /// construct from vector of other scalar type
    template <typename OtherScalarType>
    explicit Vector(const Vector<OtherScalarType,N>& o)
    {
        for (int i=0; i<N; ++i)
            data_[i] = o.data_[i];
    }

    /// cast to vector of other scalar type
    template <typename OtherScalarType>
    operator Vector<OtherScalarType,N>()
    {
        Vector<OtherScalarType,N> v;
        for (int i=0; i<N; ++i)
            v[i] = data_[i];
        return v;
    }

    /// cast to Scalar array
    operator Scalar*()
    {
        return data_;
    }

    /// cast to const Scalar array
    operator const Scalar*() const
    {
        return data_;
    }

    /// cast to const Scalar array
    const Scalar* data() const
    {
        return data_;
    }

    /// get i'th element read-write
    Scalar& operator[](unsigned int i)
    {
        assert(i<N);
        return data_[i];
    }

    /// get i'th element read-only
    const Scalar operator[](unsigned int i) const
    {
        assert(i<N);
        return data_[i];
    }


    /// assign a scalar to all components
    Vector<Scalar,N>& operator=(const Scalar s)
    {
        for (int i=0; i<N; ++i)
            data_[i] = s;
        return *this;
    }

    /// assignment from a vector of different scalar type
    template<typename otherScalarType>
    Vector<Scalar,N>& operator=(const Vector<otherScalarType,N>& o)
    {
        for (int i=0; i<N; ++i)
            data_[i] = Scalar(o[i]);
        return *this;
    }

    /// component-wise comparison
    bool operator==(const Vector<Scalar,N>& other) const
    {
        for (int i=0; i<N; ++i)
            if (data_[i]!=other.data_[i])
                return false;
        return true;
    }

    /// component-wise comparison
    bool operator!=(const Vector<Scalar,N>& other) const
    {
        for (int i=0; i<N; ++i)
            if (data_[i]!=other.data_[i])
                return true;
        return false;
    }


    /// multiply vector by scalar s
    Vector<Scalar,N>& operator*=(const Scalar s)
    {
        for (int i=0; i<N; ++i)
            data_[i] *= s;
        return *this;
    }

    /// divide vector by scalar s
    Vector<Scalar,N>& operator/=(const Scalar s)
    {
        for (int i=0; i<N; ++i)
            data_[i] /= s;
        return *this;
    }

    /// subtract vector v
    Vector<Scalar,N>& operator-=(const Vector<Scalar,N>& v)
    {
        for (int i=0; i<N; ++i)
            data_[i] -= v.data_[i];
        return *this;
    }

    /// add vector v
    Vector<Scalar,N>& operator+=(const Vector<Scalar,N>& v)
    {
        for (int i=0; i<N; ++i)
            data_[i] += v.data_[i];
        return *this;
    }


    /// normalize vector, return normalized vector
    Vector<Scalar,N>& normalize()
    {
        Scalar n = norm(*this);
        if (n > (std::numeric_limits<Scalar>::min)())
            *this *= 1.0/n;
        return *this;
    }


    /// return vector with minimum of this and other in each component
    Vector<Scalar,N> minimize(const Vector<Scalar,N>& other)
    {
        for (int i = 0; i < N; ++i)
            if (other[i] < data_[i])
                data_[i] = other[i];
        return *this;
    }


    /// return vector with maximum of this and other in each component
    Vector<Scalar,N> maximize(const Vector<Scalar,N>& other)
    {
        for (int i = 0; i < N; ++i)
            if (other[i] > data_[i])
                data_[i] = other[i];
        return *this;
    }


public:
    /** The N values of type Scalar are the only data members
     of this class. This guarantees 100% compatibility with arrays of type
     Scalar and size N, allowing us to define the cast operators to and from
     arrays and array pointers */
    Scalar data_[N];
};



//== FUNCTIONS ================================================================


/// read the space-separated components of a vector from a stream
template <typename Scalar,int N>
inline std::istream& operator>>(std::istream& is, Vector<Scalar,N>& vec)
{
    for (int i=0; i<N; ++i)
        is >> vec[i];
    return is;
}


/// output a vector by printing its space-separated compontens
template <typename Scalar,int N>
inline std::ostream& operator<<(std::ostream& os, const Vector<Scalar,N>& vec)
{
    for (int i=0; i<N-1; ++i)
        os << vec[i] << " ";
    os << vec[N-1];
    return os;
}


/// negate vector
template <typename Scalar, int N>
inline Vector<Scalar,N> operator-(const Vector<Scalar,N>& v)
{
    Vector<Scalar,N> vv;
    for (int i=0; i<N; ++i)
        vv[i] = -v[i];
    return vv;
}


/// scalar * vector
template <typename Scalar, typename Scalar2, int N>
inline Vector<Scalar,N> operator*(const Scalar2 s, const Vector<Scalar,N>& v )
{
    return Vector<Scalar,N>(v) *= (Scalar)s;
}


/// vector * scalar
template <typename Scalar, typename Scalar2, int N>
inline Vector<Scalar,N> operator*(const Vector<Scalar,N>& v, const Scalar2 s)
{
    return Vector<Scalar,N>(v) *= (Scalar)s;
}


/// vector / scalar
template <typename Scalar, typename Scalar2, int N>
inline Vector<Scalar,N> operator/(const Vector<Scalar,N>& v, const Scalar2 s)
{
    return Vector<Scalar,N>(v) /= Scalar(s);
}


/// vector + vector
template <typename Scalar, int N>
inline Vector<Scalar,N> operator+(const Vector<Scalar,N>& v0, const Vector<Scalar,N>& v1)
{
    return Vector<Scalar,N>(v0) += v1;
}


/// vector - vector
template <typename Scalar, int N>
inline Vector<Scalar,N> operator-(const Vector<Scalar,N>& v0, const Vector<Scalar,N>& v1)
{
    return Vector<Scalar,N>(v0) -= v1;
}


/// compute the Euclidean norm of a vector
template <typename Scalar, int N>
inline Scalar norm(const Vector<Scalar,N>& v)
{
    Scalar s = v[0]*v[0];
    for (int i=1; i<N; ++i)
        s += v[i]*v[i];
    return (Scalar)sqrt(s);
}


/// compute the Euclidean norm of a vector
template <typename Scalar, int N>
inline Vector<Scalar,N> normalize(const Vector<Scalar,N>& v)
{
    return v/norm(v);
}


/// compute the squared Euclidean norm of a vector
template <typename Scalar, int N>
inline Scalar sqrnorm(const Vector<Scalar,N>& v)
{
    Scalar s = v[0]*v[0];
    for (int i=1; i<N; ++i)
        s += v[i]*v[i];
    return s;
}


/// compute the dot product of two vectors
template <typename Scalar, int N>
inline Scalar dot(const Vector<Scalar,N>& v0, const Vector<Scalar,N>& v1)
{
    Scalar p = v0[0]*v1[0];
    for (int i=1; i<N; ++i)
        p += v0[i]*v1[i];
    return p;
}


/// compute the Euclidean distance between two points
template <typename Scalar, int N>
inline Scalar distance(const Vector<Scalar,N>& v0, const Vector<Scalar,N>& v1)
{
    Scalar dist(0), d;
    for (int i=0; i<N; ++i)
    {
        d = v0[i] - v1[i];
        d *= d;
        dist += d;
    }
    return (Scalar)sqrt(dist);
}


/// compute the cross product of two vectors (only valid for 3D vectors)
template <typename Scalar>
inline Vector<Scalar,3> cross(const Vector<Scalar,3>& v0, const Vector<Scalar,3>& v1)
{
    return Vector<Scalar,3>(v0[1]*v1[2] - v0[2]*v1[1],
                            v0[2]*v1[0] - v0[0]*v1[2],
                            v0[0]*v1[1] - v0[1]*v1[0]);
}


//== TEMPLATE SPECIALIZATIONS FOR 3D ==========================================


#if 1

template <typename Scalar>
class Vector<Scalar,3>
{
public:

    typedef Scalar value_type;
    static int size() { return 3; }


    Vector() {}

    explicit Vector(const Scalar s) : x(s), y(s), z(s) {}

    Vector(const Scalar xx, const Scalar yy, const Scalar zz)
    : x(xx), y(yy), z(zz) {}

    template <typename OtherScalarType>
    explicit Vector(const Vector<OtherScalarType,3>& o)
    : x((Scalar)o.x), y((Scalar)o.y), z((Scalar)o.z) {}

    template <typename OtherScalarType>
    operator Vector<OtherScalarType,3>()
    {
        return Vector<OtherScalarType,3>(x,y,z);
    }

    operator Scalar*() { return &x;}
    operator const Scalar*() const { return &x; }
    const Scalar* data() const { return &x; }

    Scalar& operator[](unsigned int i)
    {
        assert(i<3);
        return (&x)[i];
    }

    const Scalar operator[](unsigned int i) const
    {
        assert(i<3);
        return (&x)[i];
    }

    Vector<Scalar,3>& operator=(const Scalar s)
    {
        x=y=z=s;
        return *this;
    }

    template <typename otherScalarType>
    Vector<Scalar,3>& operator=(const Vector<otherScalarType,3>& o)
    {
        x = (Scalar)o.x;
        y = (Scalar)o.y;
        z = (Scalar)o.z;
        return *this;
    }

    bool operator==(const Vector<Scalar,3>& o) const
    {
        if (x != o.x) return false;
        if (y != o.y) return false;
        if (z != o.z) return false;
        return true;
    }

    bool operator!=(const Vector<Scalar,3>& o) const
    {
        if (x != o.x) return true;
        if (y != o.y) return true;
        if (z != o.z) return true;
        return false;
    }

    Vector<Scalar,3>& operator*=(const Scalar s)
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    Vector<Scalar,3>& operator/=(const Scalar s)
    {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    Vector<Scalar,3>& operator-=(const Vector<Scalar,3>& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    Vector<Scalar,3>& operator+=(const Vector<Scalar,3>& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vector<Scalar,3>& normalize()
    {
        Scalar n = norm(*this);
        n = (n > (std::numeric_limits<Scalar>::min)()) ? 1.0/n : 0.0;
        x *= n;
        y *= n;
        z *= n;
        return *this;
    }

    Vector<Scalar,3> minimize(const Vector<Scalar,3>& o)
    {
        if (o.x < x) x = o.x;
        if (o.y < y) y = o.y;
        if (o.z < z) z = o.z;
        return *this;
    }

    Vector<Scalar,3> maximize(const Vector<Scalar,3>& o)
    {
        if (o.x > x) x = o.x;
        if (o.y > y) y = o.y;
        if (o.z > z) z = o.z;
        return *this;
    }


public:

    Scalar x,y,z;
};


template <typename Scalar>
inline Vector<Scalar,3> operator-(const Vector<Scalar,3>& v)
{
    return Vector<Scalar,3>(-v.x, -v.y, -v.z);
}

template <typename Scalar, typename Scalar2>
inline Vector<Scalar,3> operator*(const Scalar2 s, const Vector<Scalar,3>& v )
{
    return Vector<Scalar,3>(v.x*Scalar(s), v.y*Scalar(s), v.z*Scalar(s));
}

template <typename Scalar, typename Scalar2>
inline Vector<Scalar,3> operator*(const Vector<Scalar,3>& v, const Scalar2 s)
{
    return Vector<Scalar,3>(v.x*Scalar(s), v.y*Scalar(s), v.z*Scalar(s));
}

template <typename Scalar, typename Scalar2>
inline Vector<Scalar,3> operator/(const Vector<Scalar,3>& v, const Scalar2 s)
{
    return Vector<Scalar,3>(v.x/s, v.y/s, v.z/s);
}

template <typename Scalar>
inline Vector<Scalar,3> operator+(const Vector<Scalar,3>& v0, const Vector<Scalar,3>& v1)
{
    return Vector<Scalar,3>(v0.x+v1.x, v0.y+v1.y, v0.z+v1.z);
}

template <typename Scalar>
inline Vector<Scalar,3> operator-(const Vector<Scalar,3>& v0, const Vector<Scalar,3>& v1)
{
    return Vector<Scalar,3>(v0.x-v1.x, v0.y-v1.y, v0.z-v1.z);
}

template <typename Scalar>
inline Scalar norm(const Vector<Scalar,3>& v)
{
    Scalar s = v.x*v.x;
    s += v.y*v.y;
    s += v.z*v.z;
    return (Scalar)sqrt(s);
}

template <typename Scalar>
inline Scalar sqrnorm(const Vector<Scalar,3>& v)
{
    Scalar s = v.x*v.x;
    s += v.y*v.y;
    s += v.z*v.z;
    return s;
}

template <typename Scalar>
inline Vector<Scalar,3> normalize(const Vector<Scalar,3>& v)
{
    Scalar n = v.x*v.x;
    n += v.y*v.y;
    n += v.z*v.z;
    n = (Scalar)sqrt(n);
    return Vector<Scalar,3>(v.x/n, v.y/n, v.z/n);
}

template <typename Scalar>
inline Scalar dot(const Vector<Scalar,3>& v0, const Vector<Scalar,3>& v1)
{
    Scalar s = v0.x*v1.x;
    s += v0.y*v1.y;
    s += v0.z*v1.z;
    return s;
}

template <typename Scalar>
inline Scalar distance(const Vector<Scalar,3>& v0, const Vector<Scalar,3>& v1)
{
    Scalar dist(0), d;
    for (int i=0; i<3; ++i)
    {
        d = v0[i] - v1[i];
        d *= d;
        dist += d;
    }
    return (Scalar)sqrt(dist);
}

#endif



//== TYPEDEFS =================================================================


/** 2-byte signed vector */
typedef Vector<signed char,2> Vec2c;
/** 2-byte unsigned vector */
typedef Vector<unsigned char,2> Vec2uc;
/** 2-short signed vector */
typedef Vector<signed short int,2> Vec2s;
/** 2-short unsigned vector */
typedef Vector<unsigned short int,2> Vec2us;
/** 2-int signed vector */
typedef Vector<signed int,2> Vec2i;
/** 2-int unsigned vector */
typedef Vector<unsigned int,2> Vec2ui;
/** 2-float vector */
typedef Vector<float,2> Vec2f;
/** 2-double vector */
typedef Vector<double,2> Vec2d;

/** 3-byte signed vector */
typedef Vector<signed char,3> Vec3c;
/** 3-byte unsigned vector */
typedef Vector<unsigned char,3> Vec3uc;
/** 3-short signed vector */
typedef Vector<signed short int,3> Vec3s;
/** 3-short unsigned vector */
typedef Vector<unsigned short int,3> Vec3us;
/** 3-int signed vector */
typedef Vector<signed int,3> Vec3i;
/** 3-int unsigned vector */
typedef Vector<unsigned int,3> Vec3ui;
/** 3-float vector */
typedef Vector<float,3> Vec3f;
/** 3-double vector */
typedef Vector<double,3> Vec3d;

/** 4-byte signed vector */
typedef Vector<signed char,4> Vec4c;
/** 4-byte unsigned vector */
typedef Vector<unsigned char,4> Vec4uc;
/** 4-short signed vector */
typedef Vector<signed short int,4> Vec4s;
/** 4-short unsigned vector */
typedef Vector<unsigned short int,4> Vec4us;
/** 4-int signed vector */
typedef Vector<signed int,4> Vec4i;
/** 4-int unsigned vector */
typedef Vector<unsigned int,4> Vec4ui;
/** 4-float vector */
typedef Vector<float,4> Vec4f;
/** 4-double vector */
typedef Vector<double,4> Vec4d;


//=============================================================================
/// @}
//=============================================================================
#endif // VECTOR_H
//=============================================================================
