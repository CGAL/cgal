#include <CGAL/config.h> // needed for the LONGNAME flag

#if defined(CGAL_CFG_NO_LONGNAME_PROBLEM) || defined(_MSC_VER)
#define Quotient                        Qt
#define Homogeneous                     Hs
#define Cartesian                       Cn
#define Simple_cartesian                SC
#define Filtered_kernel                 FKl
#define Segment_2                       St
#define Point_2                         Pt2
#define Topological_map                 TM
#define Planar_map_2                    PMp
#define Arrangement_2                   Ar
#define I_Polyhedron_iterator           IPI
#define I_Polyhedron_const_iterator     IPCI
#define Arr_2_default_dcel              ADD
#define Arr_segment_exact_traits        AST
#define Arr_circles_real_traits         ACRT
#define Arr_polyline_traits             APT
#define Arr_base_node                   ABN
#define Arr_2_halfedge_base             AHB
#define Arr_2_vertex_base               AVB
#define Arr_2_face_base                 AFB
#define In_place_list_iterator          IPLI
#define In_place_list_const_iterator    IPLCI
#define allocator                       All
#define Planar_map_traits_wrap          PmTW
#define Td_X_trapezoid                  TXT
#define PL_X_curve_plus                 PXCP

// #define bidirectional_iterator_tag      BIT

#if defined(_MSC_VER)
// Has no effect, probably bug in MSVC
#pragma warning(disable:4503)
#endif

#endif
