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
in highp   vec3 a_Pos;
in highp   vec3 a_Normal;
in mediump vec3 a_Color;

out highp   vec4 vs_fP; // view space position
out highp   vec4 ls_fP; // local space position 
out highp   vec3 fN;
out mediump vec4 fColor;

uniform highp   mat4  u_Mvp;
uniform highp   mat4  u_Mv;
uniform mediump float u_PointSize;
uniform mediump vec3  u_DefaultColor; 
uniform         bool  u_UseDefaultColor; 

void main(void)
{ 
  fColor = vec4(a_Color, 1.0); 
  if (u_UseDefaultColor)
  {
    fColor = vec4(u_DefaultColor, 1.0);
  }

  vec4 pos = vec4(a_Pos, 1.0); 

  ls_fP = pos;
  vs_fP = u_Mv * pos;

  fN = mat3(u_Mv)* a_Normal;

  gl_Position = u_Mvp * pos;
  gl_PointSize = u_PointSize;
}
)DELIM";

const char FRAGMENT_SOURCE_COLOR[]=R"DELIM(
#version 150
in highp   vec4 vs_fP;
in highp   vec4 ls_fP;
in highp   vec3 fN;
in mediump vec4 fColor;

out mediump vec4 out_color;

uniform highp   vec4  u_LightPos;
uniform mediump vec4  u_LightDiff;
uniform mediump vec4  u_LightSpec;
uniform mediump vec4  u_LightAmb;
uniform mediump float u_SpecPower;

uniform highp   vec4  u_ClipPlane;
uniform highp   vec4  u_PointPlane;
uniform mediump float u_RenderingMode;
uniform mediump float u_RenderingTransparency;

void main(void)
{
  highp vec3 L = u_LightPos.xyz - vs_fP.xyz;
  highp vec3 V = -vs_fP.xyz;

  highp vec3 a_Normal = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, a_Normal);
  highp vec4 diffuse = vec4(max(dot(a_Normal,L), 0.0) * u_LightDiff.rgb * fColor.rgb, 1.0);
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
in highp   vec3 a_Pos;
in mediump vec3 a_Color;

out mediump vec4 fColor; 
out highp   vec4 ls_fP; // local space 

uniform highp   mat4  u_Mvp;
uniform mediump float u_PointSize;
uniform         bool  u_IsOrthographic;
uniform mediump vec3  u_DefaultColor; 
uniform         bool  u_UseDefaultColor; 

bool EqualZero(float value)
{
  return abs(value) < 0.00001;
}

void main(void)
{
  fColor = vec4(a_Color, 1.0); 
  if (u_UseDefaultColor)
  {
    fColor = vec4(u_DefaultColor, 1.0);
  }

  vec4 pos = vec4(a_Pos, 1.0); 

  ls_fP  = pos;

  gl_Position  = u_Mvp * pos;

  float distance = gl_Position.w;
  if (u_IsOrthographic) 
  {
    distance = u_PointSize;
  } 

  float effectiveDistance = EqualZero(distance) ? 0.00001 : distance;
  gl_PointSize = u_PointSize / effectiveDistance * 5.0;
}
)DELIM";

const char VERTEX_SOURCE_SHAPE[]=R"DELIM(
#version 150
in highp   vec3 a_Pos;
in mediump vec3 a_Color;

out mediump vec4 gColor; 
out highp   vec4 ls_fP; 

uniform highp   mat4 u_Mvp;
uniform mediump vec3 u_DefaultColor;
uniform         bool u_UseDefaultColor;

void main(void)
{
  gColor = vec4(a_Color, 1.0);
  if (u_UseDefaultColor)
  {
    gColor = vec4(u_DefaultColor, 1.0);
  }

  gl_Position = vec4(a_Pos, 1.0);
}
)DELIM";

const char GEOMETRY_SOURCE_SPHERE[]=R"DELIM(
#version 150 
layout(points) in; 
layout(triangle_strip, max_vertices = 72) out; // max_vertices = (resolution+1) * 2 * latResolution 

#define PI 3.14159265358979323846

in mediump vec4 gColor[];

out mediump vec4 fColor;
out highp   vec4 ls_fP;

uniform highp   mat4 u_Mvp;
uniform mediump float u_Radius;

void drawSphere(in vec4 center, in float radius, in float resolution)
{
  float latResolution = resolution*0.5;
  float stepTheta = PI/latResolution;
  float stepPhi = 2*PI/resolution;
  for(int i=0; i<latResolution; ++i)
  {
    float theta1 = stepTheta*i;
    float theta2 = stepTheta*(i+1);
    for(int j=0; j<=resolution; ++j)
    {
      float phi = stepPhi*j;
      float x1 = center.x + radius * sin(theta1) * cos(phi);
      float y1 = center.y + radius * sin(theta1) * sin(phi);
      float z1 = center.z + radius * cos(theta1);
      ls_fP = vec4(x1, y1, z1, 1.0);
      gl_Position = u_Mvp * ls_fP;
      EmitVertex();

      float x2 = center.x + radius * sin(theta2) * cos(phi);
      float y2 = center.y + radius * sin(theta2) * sin(phi);
      float z2 = center.z + radius * cos(theta2); 
      ls_fP = vec4(x2, y2, z2, 1.0);
      gl_Position = u_Mvp * ls_fP;
      EmitVertex();
    }
    EndPrimitive();
  }
}

void main(void)
{
  fColor = gColor[0];

  int resolution = 8;
  vec4 center = gl_in[0].gl_Position; 

  drawSphere(center, u_Radius, resolution);
}
)DELIM";

const char GEOMETRY_SOURCE_CYLINDER[]=R"DELIM(
#version 150 
layout(lines) in; 
layout(triangle_strip, max_vertices = 22) out; 

#define PI 3.14159265358979323846

in mediump vec4 gColor[];

out mediump vec4 fColor;
out highp   vec4 ls_fP;

uniform highp   mat4 u_Mvp;
uniform mediump float u_Radius;

void drawCylinder(in vec3 u, in vec3 v, in vec4 bot, in vec4 top, in float radius, in float resolution)
{
  float step = 2*PI/resolution;
  for(int i=0; i<=resolution; ++i)
  {
    float theta = step*i;
    float cosf = radius*cos(theta);
    float sinf = radius*sin(theta);
    vec3 xAxis = cosf*u.xyz;
    vec3 yAxis = sinf*v.xyz;
    ls_fP = vec4(top.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    gl_Position = u_Mvp * ls_fP;
    EmitVertex();
    ls_fP = vec4(bot.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    gl_Position = u_Mvp * ls_fP;
    EmitVertex();
  }
  EndPrimitive(); 
}

void main(void)
{
  fColor = gColor[0];

  vec4 a = gl_in[0].gl_Position;
  vec4 b = gl_in[1].gl_Position;

  vec3 n = normalize(vec3(b.x-a.x, b.y-a.y, b.z-a.z)); // compute top normal 
  
  vec3 w = normalize(vec3(-n.z, n.x, n.y));

  // Axis vectors 
  vec3 u = normalize(cross(n, w));
  vec3 v = normalize(cross(n, u));

  int resolution = 10;

  drawCylinder(u, v, a, b, u_Radius, resolution);
}
)DELIM";

const char FRAGMENT_SOURCE_P_L[]=R"DELIM(
#version 150
in mediump vec4 fColor;
in highp   vec4 ls_fP;

out mediump vec4 out_color;

uniform highp   vec4  u_ClipPlane;
uniform highp   vec4  u_PointPlane;
uniform mediump float u_RenderingMode;

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
in highp vec3 a_Pos;

uniform highp mat4 u_Vp;
uniform highp mat4 u_M;

void main(void)
{
  gl_Position = u_Vp * u_M * vec4(a_Pos, 1.0);
}
)DELIM";

const char FRAGMENT_SOURCE_CLIPPING_PLANE[]=R"DELIM(
#version 150

out mediump vec4 out_color;

void main(void)
{
  out_color = vec4(0.0, 0.0, 0.0, 1.0);
}
)DELIM";

const char VERTEX_SOURCE_LINE[]=R"DELIM(
#version 150
in highp   vec3 a_Pos;
in mediump vec3 a_Color;

out VS_OUT {
  mediump vec4 color;
} vs_out; // vertex shader output

uniform mediump vec3 u_DefaultColor; 
uniform         bool u_UseDefaultColor; 

void main(void)
{
  vs_out.color = vec4(a_Color, 1.0);
  if (u_UseDefaultColor)
  {
    vs_out.color = vec4(u_DefaultColor, 1.0);
  }

  gl_Position = vec4(a_Pos, 1.0);
}
)DELIM";

const char GEOMETRY_SOURCE_ARROW[]=R"DELIM(
#version 150
layout(lines) in; 
layout(triangle_strip, max_vertices = 82) out; // max_vertices = resolution * 2 + 2 (cylinder) + resolution * 3 (disc) + resolution * 3 (cone)

#define PI 3.14159265358979323846

in VS_OUT {
  mediump vec4 color; 
} gs_in[]; // geometry shader input

out mediump vec4 fColor; 

uniform highp   mat4  u_Mvp;
uniform mediump float u_SceneRadius;

void drawTriangle(in vec4 v1, in vec4 v2, in vec4 v3)
{
  gl_Position = u_Mvp*v1;
  EmitVertex();
  gl_Position = u_Mvp*v2;
  EmitVertex();
  gl_Position = u_Mvp*v3;
  EmitVertex();
  EndPrimitive();
}

void drawTriangleFan(in vec3 u, in vec3 v, in vec4 center, in vec4 edge0, in float radius, in int resolution)
{
  float step = 2*PI/resolution;
  for(int i=0; i<resolution; ++i)
  {
    float theta = step*i;
    float cosf = radius*cos(theta);
    float sinf = radius*sin(theta);
    vec3 xAxis = cosf*u.xyz;
    vec3 yAxis = sinf*v.xyz;
    vec4 edge1 = vec4(edge0.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    theta = step*(i+1);
    cosf = radius*cos(theta);
    sinf = radius*sin(theta);
    xAxis = cosf*u.xyz;
    yAxis = sinf*v.xyz;
    vec4 edge2 = vec4(edge0.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    drawTriangle(center, edge1, edge2);
  } 
}

void drawDisc(in vec3 u, in vec3 v, in vec4 center, in float radius, in int resolution)
{
  drawTriangleFan(u, v, center, center, radius, resolution);
}

void drawCone(in vec3 u, in vec3 v, in vec3 n, in vec4 center, in float radius, in float height, in int resolution)
{
  drawTriangleFan(u, v, center, vec4(center.xyz-height*n.xyz, 1.0), radius, resolution);
}

void drawCylinder(in vec3 u, in vec3 v, in vec4 bot, in vec4 top, in float radius, in float resolution)
{
  float step = 2*PI/resolution;
  for(int i=0; i<=resolution; ++i)
  {
    float theta = step*i;
    float cosf = radius*cos(theta);
    float sinf = radius*sin(theta);
    vec3 xAxis = cosf*u.xyz;
    vec3 yAxis = sinf*v.xyz;
    gl_Position = u_Mvp * vec4(top.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    EmitVertex();
    gl_Position = u_Mvp * vec4(bot.xyz+xAxis.xyz+yAxis.xyz, 1.0);
    EmitVertex();
  }
  EndPrimitive(); 
}

void main(void)
{
  fColor = gs_in[0].color;

  vec4 a = gl_in[0].gl_Position;
  vec4 b = gl_in[1].gl_Position;

  vec3 n = normalize(vec3(b.x-a.x, b.y-a.y, b.z-a.z)); // compute top normal 
  vec3 w = normalize(vec3(-n.z, n.x, n.y));

  // Axis vectors 
  vec3 u = normalize(cross(n, w));
  vec3 v = normalize(cross(n, u));

  float radius = 0.013 * u_SceneRadius;
  float height = 0.035 * u_SceneRadius;
  int resolution = 10;

  vec4 c = vec4(b.xyz-height*n.xyz, 1.0);
  drawDisc(u, v, c, radius, resolution);
  drawCone(u, v, n, b, radius, height, resolution);
  drawCylinder(u, v, a, c, radius*0.5, resolution);
}
)DELIM";

const char GEOMETRY_SOURCE_LINE[]=R"DELIM(
#version 150
layout(lines) in; 
layout(line_strip, max_vertices = 2) out; 

in VS_OUT {
  mediump vec4 color;
} gs_in[]; // geometry shader input

out mediump vec4 fColor; 

uniform highp mat4 u_Mvp;

void main(void)
{
  fColor = gs_in[0].color;

  gl_Position = u_Mvp * gl_in[0].gl_Position;
  EmitVertex(); 
  gl_Position = u_Mvp * gl_in[1].gl_Position; 
  EmitVertex(); 

  EndPrimitive(); 
}
)DELIM";

const char FRAGMENT_SOURCE_LINE[]=R"DELIM(
#version 150

in mediump vec4 fColor;

out mediump vec4 out_color;

void main(void)
{
  out_color = fColor;
}
)DELIM";

const char VERTEX_SOURCE_NORMAL[]=R"DELIM(
#version 150
in highp vec3 a_Pos;
in highp vec3 a_Normal;

out VS_OUT {
  mediump vec4  color;
  highp   vec3  normal;
} vs_out; // vertex shader output
 
uniform highp   mat4 u_Mv;
uniform mediump vec3 u_DefaultColor;
uniform         bool u_UseDefaultColor;

void main(void)
{
  vs_out.color = vec4(abs(normalize(a_Normal)), 1.0); 
  if (u_UseDefaultColor)
  {
    vs_out.color = vec4(u_DefaultColor, 1.0); 
  }

  mat3 normalMatrix = mat3(transpose(inverse(u_Mv)));
  vs_out.normal = normalize(vec3(vec4(normalMatrix * a_Normal, 0.0)));
  
  gl_Position = vec4(a_Pos, 1.0); 
}
)DELIM";

const char GEOMETRY_SOURCE_NORMAL[]=R"DELIM(
#version 150 
layout (triangles) in;
layout (line_strip, max_vertices = 6) out;

in VS_OUT {
  mediump vec4  color;
  highp   vec3  normal;
} gs_in[]; // geometry shader input

out mediump vec4 fColor;
out highp   vec4 ls_fP;

uniform highp   mat4  u_Projection;
uniform mediump float u_Factor;
uniform mediump float u_SceneRadius;
uniform highp   mat4  u_Mv;
uniform bool          u_DisplayFaceNormal;

void GenerateLine(int index)
{
  fColor = gs_in[index].color; 
  
  ls_fP = gl_in[index].gl_Position;
  gl_Position = u_Projection * u_Mv * gl_in[index].gl_Position;
  EmitVertex();

  vec4 newPosition = u_Mv * gl_in[index].gl_Position + vec4(gs_in[index].normal, 0.0) * u_SceneRadius * u_Factor;
  ls_fP = inverse(u_Mv) * newPosition;
  gl_Position = u_Projection * newPosition;
  EmitVertex();

  EndPrimitive();
}

void DrawVerticesNormal()
{
  GenerateLine(0); // first vertex normal
  GenerateLine(1); // second vertex normal
  GenerateLine(2); // third vertex normal
}

void DrawFaceNormal()
{
  fColor = (gs_in[0].color + gs_in[1].color + gs_in[2].color) / 3; 
  vec4 center = (gl_in[0].gl_Position + gl_in[1].gl_Position + gl_in[2].gl_Position) / 3;
  ls_fP = center;
  gl_Position = u_Projection * u_Mv * center;
  EmitVertex();

  vec3 n = normalize((gs_in[0].normal.xyz + gs_in[1].normal.xyz + gs_in[2].normal.xyz) / 3);

  vec4 newPosition = u_Mv * center + vec4(n, 0.0) * u_SceneRadius * u_Factor;
  ls_fP = inverse(u_Mv) * newPosition;
  gl_Position = u_Projection * newPosition;
  EmitVertex();
}

void main()
{
  if (u_DisplayFaceNormal)
  {
    DrawFaceNormal();
  }
  else 
  {
    DrawVerticesNormal();
  }
}
)DELIM";

const char VERTEX_SOURCE_TRIANGLE[]=R"DELIM(
#version 150
in highp   vec3 a_Pos;

out VS_OUT {
  mediump vec4 color; 
  highp   vec4 ls_fP; 
} vs_out;

uniform highp mat4 u_Mvp;

void main(void)
{
  vec4 pos = vec4(a_Pos, 1.0);

  vs_out.color = vec4(0.85, 0.85, 0.85, 1.0);
  vs_out.ls_fP = pos;

  gl_Position = u_Mvp * pos;
}
)DELIM";

const char GEOMETRY_SOURCE_TRIANGLE[]=R"DELIM(
#version 150
layout (triangles) in;
layout (line_strip, max_vertices=4) out; 

in VS_OUT {
  mediump vec4 color; 
  highp   vec4 ls_fP; 
} gs_in[];

out mediump vec4 fColor; 
out highp   vec4 ls_fP; 

void main(void)
{
  fColor = gs_in[0].color; 
  ls_fP = gs_in[0].ls_fP;
  gl_Position = gl_in[0].gl_Position; 
  EmitVertex(); 

  ls_fP = gs_in[1].ls_fP;
  gl_Position = gl_in[1].gl_Position; 
  EmitVertex(); 

  ls_fP = gs_in[2].ls_fP;
  gl_Position = gl_in[2].gl_Position; 
  EmitVertex(); 

  ls_fP = gs_in[0].ls_fP;
  gl_Position = gl_in[0].gl_Position; 
  EmitVertex(); 
}
)DELIM";

const char VERTEX_SOURCE_LINE_WIDTH[]=R"DELIM(
#version 150

in highp   vec3 a_Pos; 
in mediump vec3 a_Color;
 
out VS_OUT {
  mediump float pointSize; 
  mediump vec4 color; 
  highp   vec4 ls_fP; 
} vs_out; 

uniform highp   mat4  u_Mvp;
uniform mediump float u_PointSize;
uniform         bool  u_IsOrthographic;
uniform mediump vec3  u_DefaultColor;
uniform         bool  u_UseDefaultColor;

bool EqualZero(float value)
{
  return abs(value) < 0.00001;
}

void main(void)
{
  vec4 pos = vec4(a_Pos, 1.0);

  vs_out.ls_fP = pos;
  vs_out.color = vec4(a_Color, 1.0); 
  if (u_UseDefaultColor)
  {
    vs_out.color = vec4(u_DefaultColor, 1.0);
  }

  gl_Position  = u_Mvp * pos;

  float distance = gl_Position.w;
  if (u_IsOrthographic) 
  {
    distance = u_PointSize;
  } 

  float effectiveDistance = EqualZero(distance) ? 0.00001 : distance;
  vs_out.pointSize = u_PointSize / effectiveDistance;
}
)DELIM";

const char GEOMETRY_SOURCE_LINE_WIDTH[]=R"DELIM(
#version 150
layout (lines) in;
layout (triangle_strip, max_vertices = 4) out;

in mediump vec4 g_Color[]; 

in VS_OUT {
  mediump float pointSize; 
  mediump vec4 color; 
  highp   vec4 ls_fP; 
} gs_in[]; 

out mediump vec4 fColor; 
out highp   vec4 ls_fP; 

uniform mediump float u_PointSize; 
uniform mediump vec2  u_Viewport;
uniform highp   mat4  u_Mvp;

vec2 ToScreenSpace(vec4 vertex)
{
  return vec2(vertex.xy / vertex.w) * u_Viewport; 
}

vec4 ToWorldSpace(vec4 vertex)
{
  return vec4((vertex.xy * vertex.w) / u_Viewport, vertex.zw); 
}

void main(void)
{
  vec2 p0 = ToScreenSpace(gl_in[0].gl_Position);
  vec2 p1 = ToScreenSpace(gl_in[1].gl_Position);
  vec2 v0 = normalize(p1 - p0);
  vec2 n0 = vec2(-v0.y, v0.x) * u_PointSize * 0.5;
  
  // line start
  gl_Position = ToWorldSpace(vec4(p0 - n0 * gs_in[0].pointSize, gl_in[0].gl_Position.zw)); 
  fColor = gs_in[0].color;
  ls_fP = inverse(u_Mvp) * gl_Position;
  EmitVertex();
  
  gl_Position = ToWorldSpace(vec4(p0 + n0 * gs_in[0].pointSize, gl_in[0].gl_Position.zw));
  fColor = gs_in[0].color;
  ls_fP = inverse(u_Mvp) * gl_Position;
  EmitVertex();
  
  // line end
  gl_Position = ToWorldSpace(vec4(p1 - n0 * gs_in[1].pointSize, gl_in[1].gl_Position.zw));
  fColor = gs_in[1].color;
  ls_fP = inverse(u_Mvp) * gl_Position;
  EmitVertex();
  
  gl_Position = ToWorldSpace(vec4(p1 + n0 * gs_in[1].pointSize, gl_in[1].gl_Position.zw));
  fColor = gs_in[1].color;
  ls_fP = inverse(u_Mvp) * gl_Position;
  EmitVertex();
}
)DELIM";

//------------------------------------------------------------------------------
//  compatibility shaders

const char VERTEX_SOURCE_COLOR_COMP[]=R"DELIM(
varying highp   vec3 a_Pos;
varying highp   vec3 a_Normal;
varying mediump vec3 a_Color;

varying highp   vec4 vs_fP; // view space position
varying highp   vec3 fN;
varying mediump vec4 fColor;

uniform highp   mat4  u_Mvp;
uniform highp   mat4  u_Mv;
uniform mediump float u_PointSize;
uniform mediump vec3  u_DefaultColor;
uniform         bool  u_UseDefaultColor;

void main(void)
{
  vec4 pos = vec4(a_Pos, 1.0);

  vs_fP = u_Mv * pos;
  highp mat3 mv_matrix_3;
  mv_matrix_3[0] = mv_matrix[0].xyz;
  mv_matrix_3[1] = mv_matrix[1].xyz;
  mv_matrix_3[2] = mv_matrix[2].xyz;
  fN = mv_matrix_3* a_Normal;

  fColor = vec4(a_Color, 1.0);
  if (u_UseDefaultColor)
  {
    fColor = vec4(u_DefaultColor, 1.0);
  }

  gl_PointSize = u_PointSize;

  gl_Position = u_Mvp * pos;
}
)DELIM";

const char FRAGMENT_SOURCE_COLOR_COMP[]=R"DELIM(
varying highp   vec4 vs_fP;
varying highp   vec3 fN;
varying mediump vec4 fColor;

uniform highp   vec4 u_LightPos;
uniform mediump vec4 u_LightDiff;
uniform mediump vec4 u_LightSpec;
uniform mediump vec4 u_LightAmb;
uniform mediump float u_SpecPower ;

void main(void)
{
  highp vec3 L = u_LightPos.xyz - vs_fP.xyz;
  highp vec3 V = -vs_fP.xyz;

  highp vec3 a_Normal = normalize(fN);
  L = normalize(L);
  V = normalize(V);

  highp vec3 R = reflect(-L, a_Normal);
  highp vec4 diffuse = max(dot(a_Normal,L), 0.0) * u_LightDiff * fColor;
  highp vec4 specular = pow(max(dot(R,V), 0.0), u_SpecPower) * u_LightSpec;

  gl_FragColor = u_LightAmb*fColor + diffuse;
}
)DELIM";

const char VERTEX_SOURCE_P_L_COMP[]=R"DELIM(
varying highp   vec3 a_Pos;
varying mediump vec3 a_Color;

varying mediump vec4 fColor;

uniform highp   mat4  u_Mvp;
uniform mediump float u_PointSize;
uniform mediump vec3  u_DefaultColor;
uniform         bool  u_UseDefaultColor;

void main(void)
{
  fColor = vec4(a_Color, 1.0);
  if (u_UseDefaultColor)
  {
    fColor = vec4(u_DefaultColor, 1.0);
  }

  gl_PointSize = u_PointSize;
  gl_Position = u_Mvp * vec4(a_Pos, 1.0);
}
)DELIM";

const char FRAGMENT_SOURCE_P_L_COMP[]=R"DELIM(
varying mediump vec4 fColor;

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
