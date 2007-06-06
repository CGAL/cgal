#ifndef DSRPDB_TRANSFORM_H
#define DSRPDB_TRANSFORM_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Quaternion.h>
#include <CGAL/PDB/Point.h>
#include <iostream>

CGAL_PDB_BEGIN_NAMESPACE
  //! This class provides a simple rigid transformation matrix
  /*!  The matrix has a rotational an a translational part. However,
    it is not enforced that the rotational part is actually just a
    rotation (this is slightly tricky due to possible numeric errors.
  */
  struct Transform {
    typedef double NT;
    Transform(){
      for (unsigned int i=0; i< 3; ++i){
	rot_[i][0]=0;
	rot_[i][1]=0;
	rot_[i][2]=0;
	rot_[i][i]=1;
	trans_[i]=0;
      }
    }

    //! Construct a transform from a rotation matrix and an offset vector.
    template <class TA, class TB>
    Transform(TA rot, TB trans){
      for (unsigned int i=0; i< 3; ++i){
	trans_[i]=trans[i];
	for (unsigned int j=0; j< 3; ++j){
	  rot_[i][j]=rot[i][j];
	}
      }
      
#ifndef NDEBUG
      NT m01 = rot[0][0]*rot[1][1] - rot[1][0]*rot[0][1];
      NT m02 = rot[0][0]*rot[2][1] - rot[2][0]*rot[0][1];
      NT m12 = rot[1][0]*rot[2][1] - rot[2][0]*rot[1][1];
      NT det = m01*rot[2][2] - m02*rot[1][2] + m12*rot[0][2];
      assert(det >0);
      assert(std::abs(1-det) < .25);
#endif
    }

    Transform(const Vector &trans, const Quaternion &q);


    //! Apply a transformation to a point
    template <class Point>
    Point operator()(const Point &pt) const {
      NT ix= pt.x();// + trans_[0];
      NT iy= pt.y();// + trans_[1];
      NT iz= pt.z();// + trans_[2];
      NT x= ix*rot_[0][0] + iy*rot_[0][1] + iz*rot_[0][2] + trans_[0];
      NT y= ix*rot_[1][0] + iy*rot_[1][1] + iz*rot_[1][2] + trans_[1];
      NT z= ix*rot_[2][0] + iy*rot_[2][1] + iz*rot_[2][2] + trans_[2];
      return Point(x,y,z);
    }



    //! Set the translation part of the transformation matrix.
    template <class Pt>
    void set_translation(Pt tr){
      trans_[0]=tr.x();
      trans_[1]=tr.y();
      trans_[2]=tr.z();
    }
    
    const NT* translation() const {
      return trans_;
    }

    Quaternion quaternion() const {
      return Quaternion(rot_);
    }

    void write(std::ostream &out) const;

    NT error(const Transform &o) const {
      NT n=0;
      for (unsigned int i=0; i< 3; ++i){
	n += std::abs(o.trans_[i]- trans_[i]);
	for (unsigned int j=0; j< 3; ++j){
	  n+= std::abs(o.rot_[i][j]- rot_[i][j]);
	}
      }
      return n;
    }

  private:
    NT rot_[3][3];
    NT trans_[3];
  };
  
  //! Write the transformation matrix to a stream.
  /*!
    The output format is 
    r r r t
    r r r t
    r r r t
    0 0 0 1
  */
  inline std::ostream &operator<<(std::ostream &out, const Transform &t) {
    t.write(out);
    return out;
  }
CGAL_PDB_END_NAMESPACE

#endif
