#include <CGAL/config.h> // needed for the LONGNAME flag

#if defined(CGAL_CFG_LONGNAME_BUG)
#define Quotient                        Qt
#define Homogeneous                     Hs
#define Cartesian                       Cn
#define Simple_cartesian                SC
#define Cartesian_converter             CC
#define Filtered_kernel                 FKl
#define Topological_map                 TM
#define Pm_default_dcel                 PDD
#define Pm_halfedge_base                PHB
#define Pm_vertex_base                  PVB
#define Pm_face_base                    PFB
#define Planar_map_2                    PMp
// #define _Nonconst_traits                NTs
#define Td_traits                       TT
#define Pm_segment_traits_2             PST
#define Pm_traits_wrap_2                PMTW

#define Planar_map_with_intersections_2 PMWI
#define Arr_segment_traits_2            AST

#define Interval_converter              IC
#define NT_converter                    NC
#define Td_X_trapezoid                  TXT
#define Interval_nt                     IN
#define PL_X_curve_plus                 PXCP
#define Forward_circulator_tag          FCT

#define remove_in_face_interior         RIFI
#define I_Polyhedron_const_iterator     PCI
#define _Rb_tree_iterator               RTI

#define _Rb_tree                        RT
// #define Trapezoidal_decomposition_2     TD
// #define bidirectional_iterator_tag      BIT

#if defined(_MSC_VER)
// Has no effect, probably bug in MSVC
#pragma warning(disable:4503)
#endif

#endif
