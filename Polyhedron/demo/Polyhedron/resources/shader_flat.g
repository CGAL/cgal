#version 430 core

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in VS_OUT
{
  vec4 fP;
  vec3 normal;
  vec4 out_color;
  float dist[6];
} gs_in[3];

out GS_OUT
{
  vec4 fP;
  flat vec3 normal;
  flat vec4 color;
  float dist[6];
} gs_out;

uniform mat4 mv_matrix;

void main(void)
{
  gl_Position = gl_in[0].gl_Position;
  gs_out.fP = gs_in[0].fP;

for(int i=0; i< 6; ++i)
  gs_out.dist[i] = gs_in[0].dist[i];
  EmitVertex();

  gl_Position = gl_in[1].gl_Position;
  gs_out.fP = gs_in[1].fP;
  for(int i=0; i< 6; ++i)
    gs_out.dist[i] = gs_in[1].dist[i];
EmitVertex();

  gl_Position = gl_in[2].gl_Position;
  gs_out.fP = gs_in[2].fP;

  // We're only writing the output color for the last
  // vertex here because they're flat attributes,
  // and the last vertex is the provoking vertex by default
  vec3 normal = vec3(0.0,0.0,0.0);
  for(int i=0; i<3; ++i)
  {
    normal+=gs_in[i].normal/3;
  }
  gs_out.normal = mat3(mv_matrix)* normal;
  gs_out.color = gs_in[2].out_color;

  for(int i=0; i< 6; ++i)
    gs_out.dist[i] = gs_in[2].dist[i];

  EmitVertex();

  EndPrimitive();
}
