
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\cgalRefines{Kernel::HasOn_3}
*/
class SphericalKernel::HasOn_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Returns `true` when `obj1` contains `obj2`, where Type1 and Type2 can be respectively:

<TABLE CELLSPACING=5 >
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
<TR>
<TD ALIGN=LEFT NOWRAP>
Type1
<TD ALIGN=LEFT NOWRAP>
Type2
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Sphere_3,
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Circular_arc_point_3
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Plane_3,
<TD ALIGN=LEFT NOWRAP>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Line_3,
<TD ALIGN=LEFT NOWRAP>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Circle_3,
<TD ALIGN=LEFT NOWRAP>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Line_arc_3, or
<TD ALIGN=LEFT NOWRAP>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Circular_arc_3
<TD ALIGN=LEFT NOWRAP>
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Plane_3,
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Line_arc_3
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Line_3
<TD ALIGN=LEFT NOWRAP>
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Sphere_3,
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Circular_arc_3
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Circle_3, or
<TD ALIGN=LEFT NOWRAP>
<TR>
<TD ALIGN=LEFT NOWRAP>
SphericalKernel::Plane_3
<TD ALIGN=LEFT NOWRAP>
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR>
</TABLE>


*/
bool
operator()(const Type1 &obj1, const Type2 &obj2);

/// @}

}; /* end SphericalKernel::HasOn_3 */

