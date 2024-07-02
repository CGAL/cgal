
#version 330

in vec3 vpos[];
out vec4 vCol;

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

void main() {
  const vec3 lightDir = normalize(vec3(1,.5,.5));

  // compute the normal for the current triangle
  vec3 triNormal = normalize(cross(vpos[1]-vpos[0], vpos[2]-vpos[0]));
  //float c = clamp(dot(lightDir,triNormal), 0, 1);
  float c = abs(dot(lightDir,triNormal));
  vCol = vec4(.2, .2,0,1) + vec4(c,c,0,1);

  gl_Position = gl_in[0].gl_Position; EmitVertex();
  gl_Position = gl_in[1].gl_Position; EmitVertex();
  gl_Position = gl_in[2].gl_Position; EmitVertex();
  EndPrimitive();
}
