#ifndef HEXMESHING_PLANE_NORMAL_H
#define HEXMESHING_PLANE_NORMAL_H

namespace CGAL::Hexmeshing {
  /**
  * @brief Enumeration representing the normal direction of a plane in 3D space
  * 
  * This enum is used to specify the orientation of planes in the hexahedral mesh.
  * Each value represents a normal direction along one of the principal axes,
  * with an additional NONE value for unspecified or invalid cases.
  */
  enum PlaneNormal { 
    X,      ///< Normal direction along the X axis
    Y,      ///< Normal direction along the Y axis
    Z,      ///< Normal direction along the Z axis
    NONE = -1  ///< Represents no specific direction or invalid plane
  };
}


#endif