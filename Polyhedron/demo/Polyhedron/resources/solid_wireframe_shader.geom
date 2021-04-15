#version 150

layout (lines) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_OUT
{
  vec4 out_color;
  float dist[6];
} gs_in[2];

uniform vec2 viewport;
uniform float width;
uniform float near;
uniform float far;

out GS_OUT
{
   vec4 color;
   float dist[6];
} gs_out;


void main(void)
{
  ;
  //linearized arbitrary Z offset to keep the widelines in front of the edges as much as possible 
  //Problem if coef is not 0, edges might be visible through faces.
  float coef = 0.0f;
  float z_offset = 0;//(2.0 * near) / (far + near - coef * (far - near)) - (2.0 * near) / (far + near);
  vec3 ndc0 = gl_in[0].gl_Position.xyz / gl_in[0].gl_Position.w;
  vec3 ndc1 = gl_in[1].gl_Position.xyz / gl_in[1].gl_Position.w;
  
  vec2 polyline = normalize(ndc1.xy - ndc0.xy);
  vec2 ortho_dir = vec2(-polyline.y, polyline.x);
  vec2 offset = (vec2(width) / viewport) * ortho_dir;
  
  vec4 cpos = gl_in[0].gl_Position;
  gl_Position = vec4(cpos.xy + offset*cpos.w, cpos.z-z_offset, cpos.w);
  gs_out.color = gs_in[0].out_color;
  gs_out.dist = gs_in[0].dist;
  EmitVertex();
  
  cpos = gl_in[0].gl_Position;
  gl_Position = vec4(cpos.xy - offset*cpos.w, cpos.z-z_offset, cpos.w);
  gs_out.color = gs_in[0].out_color;
  gs_out.dist = gs_in[0].dist;
  EmitVertex();
  
  cpos = gl_in[1].gl_Position;
  gl_Position = vec4(cpos.xy + offset*cpos.w, cpos.z-z_offset, cpos.w);
  gs_out.color = gs_in[1].out_color;
  gs_out.dist = gs_in[1].dist;
  EmitVertex();
  
  cpos = gl_in[1].gl_Position;
  gl_Position = vec4(cpos.xy - offset*cpos.w, cpos.z-z_offset, cpos.w);
  gs_out.color = gs_in[1].out_color;
  gs_out.dist = gs_in[1].dist;
  EmitVertex();
  EndPrimitive();

}
