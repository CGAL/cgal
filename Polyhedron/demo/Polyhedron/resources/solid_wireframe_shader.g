#version 150

layout (lines) in;
layout (triangle_strip, max_vertices = 6) out;

in VS_OUT
{
  vec4 out_color;
} gs_in[2];

uniform mat4 viewport_matrix;
uniform mat4 viewport_matrix_inv;
uniform float width;

out GS_OUT
{
   vec4 color;
} gs_out;

void main(void)
{

  vec4 p0 = gl_in[0].gl_Position;
  vec4 p1 = gl_in[1].gl_Position;
  float w0 = p0.w;
  float w1 = p1.w;

  vec3 P0 = vec3(viewport_matrix * vec4(p0.x/w0, p0.y/w0, p0.z/w0, p0.w/w0));
  vec3 P1 = vec3(viewport_matrix * vec4(p1.x/w1, p1.y/w1, p1.z/w1, p1.w/w1));

  vec3 n = cross(vec3(0,0,1), P1-P0);
  vec3 N;
  if(n == vec3(0,0,0) )
    N = vec3(0,0,0);
  else
    N = normalize(n);

  vec3 A1 = P0+0.5*width*N;
  vec3 A2 = P0-0.5*width*N;
  vec3 B1 = P1+0.5*width*N;
  vec3 B2 = P1-0.5*width*N;
  vec4 a1 = viewport_matrix_inv*vec4(A1.x * w0, A1.y*w0, A1.z*w0, w0);
  vec4 a2 = viewport_matrix_inv*vec4(A2.x * w0, A2.y*w0, A2.z*w0, w0);
  vec4 b1 = viewport_matrix_inv*vec4(B1.x * w1, B1.y*w1, B1.z*w1, w1);
  vec4 b2 = viewport_matrix_inv*vec4(B2.x * w1, B2.y*w1, B2.z*w1, w1);


  gl_Position = a1;
  gs_out.color = gs_in[0].out_color;
  EmitVertex();
  
  gl_Position = a2;
  gs_out.color = gs_in[0].out_color;
  EmitVertex();

  gl_Position = b1;
  gs_out.color = gs_in[1].out_color;
  EmitVertex();
  EndPrimitive();
  
gl_Position = b2;
gs_out.color = gs_in[1].out_color;
EmitVertex();

gl_Position = b1;
gs_out.color = gs_in[1].out_color;
EmitVertex();

gl_Position = a2;
gs_out.color = gs_in[0].out_color;
EmitVertex();
EndPrimitive();
}
