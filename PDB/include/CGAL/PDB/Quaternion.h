#ifndef CGAL_DSRPDB_QUATERNION_H_
#define CGAL_DSRPDB_QUATERNION_H_
#include <CGAL/PDB/basic.h>
#include <cmath>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Transform.h>
#include <iostream>

namespace CGAL { namespace PDB {
//! \cond
class Quaternion {
protected:
  typedef double NT;
  NT q_[4];

public:
  Quaternion(){
    q_[0]=0;
    q_[1]=0;
    q_[2]=0;
    q_[3]=0;
  };
  //! Construct it from the four doubles.
  Quaternion(NT a, NT b, NT c, NT d){
    q_[0]=a;
    q_[1]=b;
    q_[2]=c;
    q_[3]=d;
  };
    
    
  //! Construct from 3 doubles (the first quaternion value is implicit)
  Quaternion(NT b, NT c, NT d){
    q_[1]=b;
    q_[2]=c;
    q_[3]=d;
    q_[0]=std::sqrt(1-b*b-c*c-d*d);
  };
    
  //! Construct a quaternion for a rotation by angle around the axis
  Quaternion(const Vector &axis, double angle) {
    const double mypi= 3.14159265358979323846264338327950288;
    double f= angle/(2*mypi);
    double sf= std::sqrt(f);
    double n2= axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    double n= std::sqrt(n2);
    q_[0]= std::sqrt(1-f);
    q_[1]= sf*axis[0]/n;
    q_[2]= sf*axis[1]/n;
    q_[3]= sf*axis[2]/n;
    CGAL_assertion(q_[0]*q_[0]+ q_[1]*q_[1] + q_[2]*q_[2] + q_[3]*q_[3] < 1.1);
    CGAL_assertion(q_[0]*q_[0]+ q_[1]*q_[1] + q_[2]*q_[2] + q_[3]*q_[3] > 0.9);
  }

  template <class Arr>
  Quaternion(Arr va){
    for (int i=0; i<4; ++i){
      q_[i]=va[i];
    }
  }

  //Quaternion(const double m[3][3]);
  Quaternion(Transform t);

  NT& operator[](unsigned int i) {
    CGAL_assertion(i<4);
    return q_[i];
  }

  NT operator[](unsigned int i) const {
    CGAL_assertion(i<4);
    return q_[i];
  }

  Transform transform() const;
  Transform transform(Vector translation) const;
};

inline std::ostream &operator<<(std::ostream &out, const Quaternion &q) {
  out << "(" << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << ")";
  return out;
}

//! \endcond
}}
#endif
