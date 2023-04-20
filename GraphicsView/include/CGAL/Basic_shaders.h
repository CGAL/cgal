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
const char vertex_source_color[]=R"DELIM(
#version 150
in highp vec4 vertex;
in highp vec3 normal;
in highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
uniform highp float point_size;

out highp vec4 fP;
out highp vec3 fN;
out highp vec4 fColor;
out highp vec4 m_vertex;

void main(void)
{
  fP = mv_matrix * vertex;
  fN = mat3(mv_matrix)* normal;
  fColor = vec4(color, 1.0);
  gl_PointSize = point_size;

  m_vertex = vertex;

  gl_Position = mvp_matrix * vertex;
}
)DELIM";

const char fragment_source_color[]=R"DELIM(
#version 150
in highp vec4 fP;
in highp vec3 fN;
in highp vec4 fColor;
in highp vec4 m_vertex;

uniform highp vec4 light_pos;
uniform highp vec4 light_diff;
uniform highp vec4 light_spec;
uniform highp vec4 light_amb;
uniform highp float spec_power;

uniform highp vec4 clipPlane;
uniform highp vec4 pointPlane;
uniform highp float rendering_mode;
uniform highp float rendering_transparency;

out highp vec4 out_color;

void main(void)
{
  highp vec3 L = light_pos.xyz - fP.xyz;
  highp vec3 V = -fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = vec4(max(dot(N,L), 0.0) * light_diff.rgb * fColor.rgb, 1.0);
  highp vec4 ambient = vec4(light_amb.rgb * fColor.rgb, 1.0);
  highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;

  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((m_vertex.xyz-pointPlane.xyz), clipPlane.xyz));

  // rendering_mode == -1: draw all solid;
  // rendering_mode == 0: draw solid only;
  // rendering_mode == 1: draw transparent only;
  if (rendering_mode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  // draw corresponding part
  out_color = rendering_mode < 1 ? (diffuse + ambient) :
                      vec4(diffuse.rgb + ambient.rgb, rendering_transparency);
}
)DELIM";

const char vertex_source_p_l[]=R"DELIM(
#version 150
in highp vec4 vertex;
in highp vec3 color;

uniform highp mat4 mvp_matrix;
uniform highp float point_size;

out highp vec4 fColor;
out highp vec4 m_vertex;

void main(void)
{
  gl_PointSize = point_size;
  fColor = vec4(color, 1.0);
  m_vertex = vertex;
  gl_Position = mvp_matrix * vertex;
}
)DELIM";

const char fragment_source_p_l[]=R"DELIM(
#version 150
in highp vec4 fColor;
in highp vec4 m_vertex;

uniform highp vec4 clipPlane;
uniform highp vec4 pointPlane;
uniform highp float rendering_mode;

out highp vec4 out_color;

void main(void)
{
  // onPlane == 1: inside clipping plane, should be solid;
  // onPlane == -1: outside clipping plane, should be transparent;
  // onPlane == 0: on clipping plane, whatever;
  float onPlane = sign(dot((m_vertex.xyz-pointPlane.xyz), clipPlane.xyz));

  // rendering_mode == -1: draw both inside and outside;
  // rendering_mode == 0: draw inside only;
  // rendering_mode == 1: draw outside only;
  if (rendering_mode == (onPlane+1)/2) {
    // discard other than the corresponding half when rendering
    discard;
  }

  out_color = fColor;
}
)DELIM";

const char vertex_source_clipping_plane[]=R"DELIM(
#version 150
in highp vec4 vertex;

uniform highp mat4 vp_matrix;
uniform highp mat4 m_matrix;

void main(void)
{
  gl_Position = vp_matrix * m_matrix * vertex;
}
)DELIM";

const char fragment_source_clipping_plane[]=R"DELIM(
#version 150
out highp vec4 out_color;
void main(void)
{
  out_color = vec4(0.0, 0.0, 0.0, 1.0);
}
)DELIM";

//------------------------------------------------------------------------------
//  compatibility shaders

const char vertex_source_color_comp[]=R"DELIM(
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

const char fragment_source_color_comp[]=R"DELIM(
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 fColor;

uniform highp vec4 light_pos;
uniform highp vec4 light_diff;
uniform highp vec4 light_spec;
uniform highp vec4 light_amb;
uniform highp float spec_power ;

void main(void)
{
  highp vec3 L = light_pos.xyz - fP.xyz;
  highp vec3 V = -fP.xyz;

  highp vec3 N = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, N);
  highp vec4 diffuse = max(dot(N,L), 0.0) * light_diff * fColor;
  highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;

  gl_FragColor = light_amb*fColor + diffuse;
}
)DELIM";

const char vertex_source_p_l_comp[]=R"DELIM(
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

const char fragment_source_p_l_comp[]=R"DELIM(
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
