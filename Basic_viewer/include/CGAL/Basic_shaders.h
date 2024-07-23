// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_BASIC_SHADERS_H
#define CGAL_BASIC_SHADERS_H

#include <CGAL/license/GraphicsView.h>

namespace CGAL
{

//------------------------------------------------------------------------------
const char VERTEX_SOURCE_COLOR[]=R"DELIM(
#version 150
in highp vec4 P;
in highp vec3 N;
in highp vec3 Color;

uniform highp mat4 u_Mvp;
uniform highp mat4 u_Mv;
uniform highp float u_PointSize;

out highp vec4 vs_fP; // view space position
out highp vec4 ls_fP; // local space position 
out highp vec3 fN;
out highp vec4 fColor;

void main(void)
{
  ls_fP = P;
  vs_fP = u_Mv * P;

  fN = mat3(u_Mv)* N;
  fColor = vec4(Color, 1.0);

  gl_Position = u_Mvp * P;
  gl_PointSize = u_PointSize;
}
)DELIM";

const char FRAGMENT_SOURCE_COLOR[]=R"DELIM(
#version 150
in highp vec4 vs_fP;
in highp vec4 ls_fP;
in highp vec3 fN;
in highp vec4 fColor;

uniform highp vec4  u_LightPos;
uniform highp vec4  u_LightDiff;
uniform highp vec4  u_LightSpec;
uniform highp vec4  u_LightAmb;
uniform highp float u_SpecPower;

uniform highp vec4  u_ClipPlane;
uniform highp vec4  u_PointPlane;
uniform highp float u_RenderingMode;
uniform highp float u_RenderingTransparency;

out highp vec4 out_color;

void main(void)
{
  highp vec3 L = u_LightPos.xyz - vs_fP.xyz;
  highp vec3 V = -vs_fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = vec4(max(dot(N,L), 0.0) * u_LightDiff.rgb * fColor.rgb, 1.0);
  highp vec4 ambient = vec4(u_LightAmb.rgb * fColor.rgb, 1.0);
  highp vec4 specular = pow(max(dot(R,V), 0.0), u_SpecPower) * u_LightSpec;

  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((ls_fP.xyz-u_PointPlane.xyz), u_ClipPlane.xyz));

  // rendering_mode == -1: draw all solid;
  // rendering_mode == 0: draw solid only;
  // rendering_mode == 1: draw transparent only;
  if (u_RenderingMode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  // draw corresponding part
  out_color = u_RenderingMode < 1 ? (diffuse + ambient) :
                      vec4(diffuse.rgb + ambient.rgb, u_RenderingTransparency);
}
)DELIM";

const char VERTEX_SOURCE_P_L[]=R"DELIM(
#version 150
in highp vec4 P;
in highp vec3 Color;

uniform highp mat4  u_Mvp;
uniform highp float u_PointSize;

out highp vec4 fColor;
out highp vec4 ls_fP; // local space 

void main(void)
{
  fColor = vec4(Color, 1.0);
  ls_fP = P;

  gl_PointSize = u_PointSize;
  gl_Position = u_Mvp * P;
}
)DELIM";

const char GEOMETRY_SOURCE_P_L[]=R"DELIM(
#version 150 
layout(points) in; 
layout(points, max_vertices = 1) out; 

in vec4 vColor[];
out vec4 fColor; 

void main(void)
{
  fColor = vColor[0];
  gl_Position = gl_in[0].gl_Position;
  EmitVertex();
  EndPrimitive();
}
)DELIM";

const char FRAGMENT_SOURCE_P_L[]=R"DELIM(
#version 150
in highp vec4 fColor;
in highp vec4 ls_fP;

uniform highp vec4  u_ClipPlane;
uniform highp vec4  u_PointPlane;
uniform highp float u_RenderingMode;

out highp vec4 out_color;

void main(void)
{
  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((ls_fP.xyz-u_PointPlane.xyz), u_ClipPlane.xyz));

  // rendering_mode == -1: draw both inside and outside;
  // rendering_mode == 0: draw inside only;
  // rendering_mode == 1: draw outside only;
  if (u_RenderingMode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  out_color = fColor;
}
)DELIM";

const char VERTEX_SOURCE_CLIPPING_PLANE[]=R"DELIM(
#version 150
in highp vec4 P;

uniform highp mat4 u_Vp;
uniform highp mat4 u_M;

void main(void)
{
  gl_Position = u_Vp * u_M * P;
}
)DELIM";

const char FRAGMENT_SOURCE_CLIPPING_PLANE[]=R"DELIM(
#version 150

out highp vec4 out_color;

void main(void)
{
  out_color = vec4(0.0, 0.0, 0.0, 1.0);
}
)DELIM";

const char VERTEX_SOURCE_LINE[]=R"DELIM(
#version 150
in highp vec3 P;
in highp vec3 Color;

uniform mat4 u_Mvp;

out VS_OUT {
  highp vec4 color;
} vs_out;

void main(void)
{
  vs_out.color = vec4(Color, 1.0);

  gl_Position = u_Mvp * vec4(P, 1.0);
}
)DELIM";

const char GEOMETRY_SOURCE_LINE[]=R"DELIM(
#version 150 
layout(lines) in; 
layout(line_strip, max_vertices = 2) out; 

in VS_OUT {
  highp vec4 color; 
} gs_in[];

out highp vec4 fColor; 

void main(void)
{
  fColor = gs_in[0].color;

  gl_Position = gl_in[0].gl_Position;
  EmitVertex();

  gl_Position = gl_in[1].gl_Position;
  EmitVertex();

  EndPrimitive();
}
)DELIM";

const char FRAGMENT_SOURCE_LINE[]=R"DELIM(
#version 150

in highp vec4 fColor;

out highp vec4 out_color;

void main(void)
{
  out_color = fColor;
}
)DELIM";

const char VERTEX_SOURCE_NORMAL[]=R"DELIM(
#version 150
in highp vec4 P;
in highp vec3 N;
 
uniform highp mat4  u_Mv;
uniform highp vec4  u_Color;
uniform highp float u_RenderingMode;

out VS_OUT {
  highp vec4  color;
  highp vec3  normal;
} vs_out;

void main(void)
{
  if (u_RenderingMode > 0.0)
  {
    vs_out.color = u_Color; 
  }
  else 
  {
    vs_out.color = vec4(abs(normalize(N)), 1); 
  }

  mat3 normalMatrix = mat3(transpose(inverse(u_Mv)));
  vs_out.normal = normalize(vec3(vec4(normalMatrix * N, 0.0)));
  
  gl_Position = u_Mv * P; 
}
)DELIM";

const char GEOMETRY_SOURCE_NORMAL[]=R"DELIM(
#version 150 
layout (triangles) in;
layout (line_strip, max_vertices = 6) out;

in VS_OUT {
  highp vec4  color;
  highp vec3  normal;
} gs_in[];

uniform highp mat4  u_Projection;
uniform highp float u_Factor;
uniform highp float u_SceneRadius;

out GS_OUT {
  highp vec4 color;
} gs_out;

void GenerateLine(int index)
{
  gs_out.color = gs_in[index].color; 
  
  gl_Position = u_Projection * gl_in[index].gl_Position;
  EmitVertex();

  gl_Position = u_Projection * (gl_in[index].gl_Position + vec4(gs_in[index].normal, 0.0) * u_SceneRadius * u_Factor);
  EmitVertex();

  EndPrimitive();
}

void main()
{
  GenerateLine(0); // first vertex normal
  GenerateLine(1); // second vertex normal
  GenerateLine(2); // third vertex normal
}
)DELIM";

const char FRAGMENT_SOURCE_NORMAL[]=R"DELIM(
#version 150

in GS_OUT {
  highp vec4 color;
} fs_in;

out highp vec4 out_color;

void main()
{
  out_color = fs_in.color;
}
)DELIM";

//------------------------------------------------------------------------------
//  compatibility shaders

const char VERTEX_SOURCE_COLOR_COMP[]=R"DELIM(
varying highp vec4 vertex;
varying highp vec3 normal;
varying highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
uniform highp float point_size;

varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 fColor;

void main(void)
{
  fP = mv_matrix * vertex;
  highp mat3 mv_matrix_3;
  mv_matrix_3[0] = mv_matrix[0].xyz;
  mv_matrix_3[1] = mv_matrix[1].xyz;
  mv_matrix_3[2] = mv_matrix[2].xyz;
  fN = mv_matrix_3* normal;
  fColor = vec4(color, 1.0);
  gl_PointSize = point_size;

  gl_Position = mvp_matrix * vertex;
}
)DELIM";

const char FRAGMENT_SOURCE_COLOR_COMP[]=R"DELIM(
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 fColor;

uniform highp vec4 u_LightPos;
uniform highp vec4 u_LightDiff;
uniform highp vec4 u_LightSpec;
uniform highp vec4 u_LightAmb;
uniform highp float u_SpecPower ;

void main(void)
{
  highp vec3 L = u_LightPos.xyz - fP.xyz;
  highp vec3 V = -fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = max(dot(N,L), 0.0) * u_LightDiff * fColor;
  highp vec4 specular = pow(max(dot(R,V), 0.0), u_SpecPower) * u_LightSpec;

  gl_FragColor = u_LightAmb*fColor + diffuse;
}
)DELIM";

const char VERTEX_SOURCE_P_L_COMP[]=R"DELIM(
varying highp vec4 vertex;
varying highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp float point_size;

varying highp vec4 fColor;

void main(void)
{
  gl_PointSize = point_size;
  fColor = vec4(color, 1.0);
  gl_Position = mvp_matrix * vertex;
}
)DELIM";

const char FRAGMENT_SOURCE_P_L_COMP[]=R"DELIM(
varying highp vec4 fColor;
void main(void)
{
  gl_FragColor = fColor;
}
)DELIM";

/* const char vertex_source_clipping_plane_comp[]=R"DELIM(
attribute highp vec4 vertex;

uniform highp mat4 vp_matrix;
uniform highp mat4 m_matrix;

void main(void)
{
  gl_Position = vp_matrix * m_matrix * vertex;
}
)DELIM";

const char fragment_source_clipping_plane_comp[]=R"DELIM(
out highp vec4 out_color;
void main(void)
{
  out_color = vec4(0.0, 0.0, 0.0, 1.0);
}
)DELIM";
*/

}

#endif // CGAL_BASIC_SHADERS_H
