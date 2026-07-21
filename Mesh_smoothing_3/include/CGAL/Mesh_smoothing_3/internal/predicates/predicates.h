#pragma once
#include <Eigen/Eigen>
#include <array>


namespace Mesh_smoothing_3_internal {


namespace exact_predicates {

double orient3d(double pa[3], double pb[3], double pc[3], double pd[3]);

inline bool positive_tetrahedra(std::array<Eigen::Vector3d, 4> const &tetrahedra) {
    double pa[3] = {tetrahedra[0][0], tetrahedra[0][1], tetrahedra[0][2]};
    double pb[3] = {tetrahedra[1][0], tetrahedra[1][1], tetrahedra[1][2]};
    double pc[3] = {tetrahedra[2][0], tetrahedra[2][1], tetrahedra[2][2]};
    double pd[3] = {tetrahedra[3][0], tetrahedra[3][1], tetrahedra[3][2]};
    return orient3d(pa,pb,pc,pd) > 0; 
} 


#include "predicates_shewchuk.h"

}
}
