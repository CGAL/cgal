// This file is part of GlSplat, a simple splatting C++ library
//
// Copyright (C) 2008-2009 Gael Guennebaud <g.gael@free.fr>
//
// GlSplat is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// GlSplat is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with GlSplat. If not, see <http://www.gnu.org/licenses/>.

//#extension GL_ARB_texture_rectangle : enable

#ifndef ES_DEPTH_INTERPOLATION
  #define ES_DEPTH_INTERPOLATION 0
#endif

#ifndef ES_OUTPUT_DEPTH
  #define ES_OUTPUT_DEPTH 0
#endif

// avoid an annoying bug with the nvidia driver 87XX serie.
#define epsilon 0.000001

uniform vec4 viewport;

#ifndef ES_DEFERRED_SHADING

uniform sampler2DRect ColorWeight;
#if (ES_OUTPUT_DEPTH==1)
uniform sampler2DRect Depth;
#endif

void Finalization(void)
{
  vec4 color = texture2DRect(ColorWeight, gl_FragCoord.st - viewport.xy + epsilon);
  #if (ES_OUTPUT_DEPTH==1)
  gl_FragDepth = texture2DRect(Depth, gl_FragCoord.st + epsilon).x;
  #endif
  if (color.w<0.001)
    discard;
  gl_FragColor = color/color.w;
  gl_FragColor.a = 1.;
}

#else

uniform vec2 unproj;

uniform sampler2DRect ColorWeight;
uniform sampler2DRect NormalWeight;

#if ( (ES_DEPTH_INTERPOLATION==0) || (ES_OUTPUT_DEPTH==1))
uniform sampler2DRect Depth;
#endif

void Finalization(void)
{
  vec4 color = texture2DRect(ColorWeight, gl_FragCoord.st - viewport.xy + epsilon);

  if (color.w<0.001)
    discard;


  if(color.w>0.001)
    color.xyz /= color.w;

  vec3 viewVec = normalize(gl_TexCoord[0].xyz);
  vec4 normaldepth = texture2DRect(NormalWeight, gl_FragCoord.st + epsilon);

  normaldepth.xyz = normaldepth.xyz/normaldepth.w;

  #if (ES_OUTPUT_DEPTH==1)
  gl_FragDepth = texture2DRect(Depth, gl_FragCoord.st + epsilon).x;
  #endif

  #if ES_DEPTH_INTERPOLATION==2
    float depth = -normaldepth.z;
  #elif ES_DEPTH_INTERPOLATION==1
    float depth = unproj.y/(2.0*normaldepth.z+unproj.x-1.0);
  #else
    float depth = texture2DRect(Depth, gl_FragCoord.st + epsilon).x;
    depth = unproj.y/(2.0*depth+unproj.x-1.0);
  #endif

  vec3 normal = normaldepth.xyz;
  #if ES_DEPTH_INTERPOLATION!=0
    normal.z = sqrt(1. - dot(vec3(normal.xy,0),vec3(normal.xy,0)));
  #endif
  normal = normalize(normal);
  vec3 eyePos = gl_TexCoord[0].xyz * depth;

  gl_FragColor = meshlabLighting(color, eyePos, normal);
  gl_FragColor.a = 1.0;
}

#endif



