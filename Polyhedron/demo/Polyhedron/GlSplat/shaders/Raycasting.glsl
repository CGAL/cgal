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

#pragma optimize(on)

#ifndef ES_EWA_HINT
    #define ES_EWA_HINT 0
#endif

#ifndef ES_DEPTH_INTERPOLATION
    #define ES_DEPTH_INTERPOLATION 0
#endif

//--------------------------------------------------------------------------------
// shared variables
//--------------------------------------------------------------------------------

// custom vertex attributes
//attribute float radius;

#ifdef CLIPPED_SPLAT
attribute vec3 secondNormal;
varying vec4 clipLine;
#endif

// standard uniforms
uniform float expeRadiusScale;
uniform float expePreComputeRadius;
uniform float expeDepthOffset;

// varying
varying vec4 covmat;
varying vec3 fragNormal;

varying vec3 fragNoverCdotN;
varying vec3 fragCenter;
varying float scaleSquaredDistance;

#ifdef ES_ATI_WORKAROUND
varying vec4 fragCenterAndRadius;
#endif

#ifdef ES_DEPTH_CORRECTION
varying float depthOffset;
#endif

uniform vec2 halfVp;
uniform float oneOverEwaRadius;


#ifdef ES_BACKFACE_SHADING
  #undef ES_EARLY_BACK_FACE_CULLING
  //#define ES_EWA_HINT 2
#endif

//--------------------------------------------------------------------------------
// Visibility Splatting
//    Vertex Shader
//--------------------------------------------------------------------------------

#ifdef GLSPLAT__VisibilityVP__
varying vec2 scaledFragCenter2d;
void VisibilityVP(void)
{
    vec3 normal = normalize(gl_NormalMatrix * gl_Normal);
    // Point in eye space
    vec4 ePos = gl_ModelViewMatrix * gl_Vertex;

    float dotpn = dot(normal.xyz,ePos.xyz);

    vec4 oPos;

    #ifdef ES_EARLY_BACK_FACE_CULLING
    // back_face culling
    oPos = vec4(0,0,1,0);
    if(dotpn<0.)
    {
    #endif

    float radius = gl_MultiTexCoord2.x * expeRadiusScale;

    vec4 pointSize;
    pointSize.x = radius * expePreComputeRadius / ePos.z;
    gl_PointSize = max(1.0, pointSize.x);

    scaleSquaredDistance = 1.0 / (radius * radius);
    //fragNormal = normal;
    fragCenter = ePos.xyz;
    fragNoverCdotN = normal/dot(ePos.xyz,normal);

    #ifndef ES_DEPTH_CORRECTION
    ePos.xyz += normalize(ePos.xyz) * expeDepthOffset * radius;
    #else
    //ePos.xyz += normalize(ePos.xyz) * expeDepthOffset * radius;
    depthOffset = expeDepthOffset * radius;
    #endif

    oPos = gl_ProjectionMatrix * ePos;

    #if (ES_EWA_HINT>0)
    scaledFragCenter2d = 0.5*((oPos.xy/oPos.w)+1.0)*halfVp*oneOverEwaRadius;
    #endif

    #ifdef ES_ATI_WORKAROUND
    fragCenterAndRadius.xyz = (oPos.xyz/oPos.w) + 1.0;
    fragCenterAndRadius.xy = fragCenterAndRadius.xy*halfVp;
    fragCenterAndRadius.z = fragCenterAndRadius.z*0.5;
    fragCenterAndRadius.w = pointSize.x;
    #endif

    #ifndef ES_EARLY_BACK_FACE_CULLING
    oPos.w = oPos.w * (dotpn<0.0 ? 1.0 : 0.0);
    #else
    }
    #endif

    gl_Position = oPos;
}

#endif

//--------------------------------------------------------------------------------
// Visibility Splatting
//    Fragment Shader
//--------------------------------------------------------------------------------

#ifdef GLSPLAT__VisibilityFP__
varying vec2 scaledFragCenter2d;
uniform vec3 rayCastParameter1;
uniform vec3 rayCastParameter2;
uniform vec2 depthParameterCast;

void VisibilityFP(void)
{
    #ifdef ES_ATI_WORKAROUND
    vec3 fragCoord;
    fragCoord.xy = fragCenterAndRadius.xy + (gl_TexCoord[0].st-0.5) * fragCenterAndRadius.w;
    fragCoord.z = fragCenterAndRadius.z;
    #else
    vec3 fragCoord = gl_FragCoord.xyz;
    #endif
    // compute q in object space
    vec3 qOne = rayCastParameter1 * fragCoord + rayCastParameter2;
    float oneOverDepth = dot(qOne,-fragNoverCdotN);
    float depth = (1.0/oneOverDepth);
    vec3 diff = fragCenter + qOne * depth;
    float r2 = dot(diff,diff);

    #if (ES_EWA_HINT>0)
    vec2 d2 = oneOverEwaRadius*gl_FragCoord.xy - scaledFragCenter2d;
    float r2d = dot(d2,d2);
    gl_FragColor = vec4(min(r2d,r2*scaleSquaredDistance));
    #else
    gl_FragColor = vec4(r2*scaleSquaredDistance);
    #endif

    #ifdef ES_DEPTH_CORRECTION
    oneOverDepth = 1.0/(-depth+depthOffset);
    gl_FragDepth = depthParameterCast.x * oneOverDepth + depthParameterCast.y;
    #endif
}

#endif

#ifdef GLSPLAT__AttributeVP__

varying vec2 scaledFragCenter2d;

void AttributeVP(void)
{
    // transform normal
    vec3 normal = normalize(gl_NormalMatrix * gl_Normal);
    // Point in eye space
    vec4 ePos = gl_ModelViewMatrix * gl_Vertex;

    float dotpn = dot(normal.xyz,ePos.xyz);

    vec4 oPos;

    #ifdef ES_EARLY_BACK_FACE_CULLING
    // back_face culling
    oPos = vec4(0,0,1,0);
    if(dotpn<0.)
    {
    #endif

    #ifdef ES_BACKFACE_SHADING
    if(dotpn>0.)
    {
        dotpn = -dotpn;
        normal = -normal;
    }
    #endif

    float radius = gl_MultiTexCoord2.x * expeRadiusScale * 1.1;

    vec4 pointSize;
    pointSize.x = radius * expePreComputeRadius / ePos.z;

    #if (ES_EWA_HINT>0)
    gl_PointSize = max(2.0, pointSize.x);
    #else
    gl_PointSize = max(1.0, pointSize.x);
    #endif

    scaleSquaredDistance = 1. / (radius * radius);
    fragNormal = normal;
    fragCenter = ePos.xyz;
    fragNoverCdotN = normal/dot(ePos.xyz,normal);

    // Output color
    #ifdef ES_DEFERRED_SHADING
        fragNormal.xyz = normal.xyz;
        gl_FrontColor = gl_Color;
    #else
        // Output color
        #ifdef ES_LIGHTING
        gl_FrontColor = expeLighting(gl_Color, ePos.xyz, normal.xyz, 1.);
        #else
        gl_FrontColor = meshlabLighting(gl_Color, ePos.xyz, normal.xyz);
        #endif
    #endif

    oPos = gl_ModelViewProjectionMatrix * gl_Vertex;

    #ifdef ES_ATI_WORKAROUND
    fragCenterAndRadius.xyz = (oPos.xyz/oPos.w) + 1.0;
    fragCenterAndRadius.xy = fragCenterAndRadius.xy*halfVp;
    fragCenterAndRadius.z = fragCenterAndRadius.z*0.5;
    fragCenterAndRadius.w = pointSize.x;
    #endif

    #if (ES_EWA_HINT>0)
    scaledFragCenter2d = ((oPos.xy/oPos.w)+1.0)*halfVp*oneOverEwaRadius;
    #endif

    #ifndef ES_EARLY_BACK_FACE_CULLING
    oPos.w = oPos.w * (dotpn<0. ? 1.0 : 0.0);
    #else
    }
    #endif

    gl_Position = oPos;
}

#endif

//--------------------------------------------------------------------------------
// EWA Splatting
//    Fragment Shader
//--------------------------------------------------------------------------------

#ifdef GLSPLAT__AttributeFP__
// this sampler is only used by this fragment shader

varying vec2 scaledFragCenter2d;
uniform vec3 rayCastParameter1;
uniform vec3 rayCastParameter2;
uniform vec2 depthParameterCast;

// uniform sampler1D Kernel1dMap;

void AttributeFP(void)
{
    #ifdef ES_ATI_WORKAROUND
    vec3 fragCoord;
    fragCoord.xy = fragCenterAndRadius.xy + (gl_TexCoord[0].st-0.5) * fragCenterAndRadius.w;
    fragCoord.z = fragCenterAndRadius.z;
    #else
    vec3 fragCoord = gl_FragCoord.xyz;
    #endif

#if 1
    vec3 qOne = rayCastParameter1 * fragCoord + rayCastParameter2;
    float oneOverDepth = dot(qOne,fragNoverCdotN);
    float depth = (1.0/oneOverDepth);
    vec3 diff = fragCenter - qOne * depth;
    float r2 = dot(diff,diff);

    #if (ES_EWA_HINT>0)
    vec2 d2 = oneOverEwaRadius*gl_FragCoord.xy - scaledFragCenter2d;
    float r2d = dot(d2,d2);
//     float weight = texture1D(Kernel1dMap, min(r2d,r2*scaleSquaredDistance)).a;
    float weight = min(r2d,r2*scaleSquaredDistance);
    weight = clamp(1.-weight,0,1);
    weight = weight*weight;
    #else
    //float weight = texture1D(Kernel1dMap, r2*scaleSquaredDistance).a;
    float weight = clamp(1.-r2*scaleSquaredDistance,0.0,1.0);
    weight = weight*weight;
    #endif
    weight *= 0.1; // limits overflow

    #ifdef ES_DEPTH_CORRECTION
    gl_FragDepth = depthParameterCast.x * oneOverDepth + depthParameterCast.y;
    #endif

    #ifdef ES_DEFERRED_SHADING
    gl_FragData[0].rgb = gl_Color.rgb;
    gl_FragData[1].xyz = fragNormal.xyz;
    gl_FragData[1].w = weight;
    gl_FragData[0].w = weight;

    #if ES_DEPTH_INTERPOLATION==2 // linear space
        gl_FragData[1].z = -depth;
    #elif ES_DEPTH_INTERPOLATION==1 // window space
        #ifdef ES_DEPTH_CORRECTION
        gl_FragData[1].z = gl_FragDepth;
        #else
        gl_FragData[1].z = fragCoord.z;
        #endif
    #endif

    #else
    gl_FragColor.rgb = gl_Color.rgb;
    gl_FragColor.w = weight;
    #endif
#endif
}

#endif
