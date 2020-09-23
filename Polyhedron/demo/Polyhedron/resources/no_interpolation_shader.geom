#version 150

layout (lines_adjacency) in;
layout (triangle_strip, max_vertices = 6) out;

in VS_OUT
{
  vec4 fP;
  vec3 fN;
  vec4 out_color;
} gs_in[4];

out GS_OUT
{
  vec4 fP;
  vec3 fN;
  flat vec4 color[4];
  vec2 uv;
  flat vec4 prob[4];
} gs_out;

void main(void)
{
  gl_Position = gl_in[0].gl_Position;
  gs_out.fN = gs_in[0].fN;
  gs_out.fP = gs_in[0].fP;
  gs_out.uv = vec2(0.0, 0.0);

  EmitVertex();

  gl_Position = gl_in[1].gl_Position;
  gs_out.fN = gs_in[1].fN;
  gs_out.fP = gs_in[1].fP;
  gs_out.uv = vec2(0.0, 1.0);
  EmitVertex();

  gl_Position = gl_in[3].gl_Position;
  gs_out.fN = gs_in[3].fN;
  gs_out.fP = gs_in[3].fP;
  gs_out.uv = vec2(1.0, 0.0);

  // We're only writing the output color for the last
  // vertex here because they're flat attributes,
  // and the last vertex is the provoking vertex by default
  gs_out.color[0] = gs_in[0].out_color;
  gs_out.color[1] = gs_in[1].out_color;
  gs_out.color[2] = gs_in[2].out_color;
  gs_out.color[3] = gs_in[3].out_color;
  for(int i=0; i<4; ++i)
   for(int j=0;j<4;++j)
     gs_out.prob[i][j]=0;
  for(int i=0; i<4; ++i)
  {
    for(int j=0; j<4; ++j)
    {
      if(gs_out.color[i] == gs_out.color[j])
      {
        gs_out.prob[i][j] = 1.0;
      }
    }
  }
  EmitVertex();

  EndPrimitive();

  gl_Position = gl_in[1].gl_Position;
  gs_out.fN = gs_in[1].fN;
  gs_out.fP = gs_in[1].fP;
  gs_out.uv = vec2(0.0, 1.0);
  EmitVertex();

  gl_Position = gl_in[2].gl_Position;
  gs_out.fN = gs_in[2].fN;
  gs_out.fP = gs_in[2].fP;
  gs_out.uv = vec2(1.0, 1.0);
  EmitVertex();

  gl_Position = gl_in[3].gl_Position;
  gs_out.fN = gs_in[3].fN;
  gs_out.fP = gs_in[3].fP;
  gs_out.uv = vec2(1.0, 0.0);
  // Again, only write the output color for the last vertex
  gs_out.color[0] = gs_in[0].out_color;
  gs_out.color[1] = gs_in[1].out_color;
  gs_out.color[2] = gs_in[2].out_color;
  gs_out.color[3] = gs_in[3].out_color;
  for(int i=0; i<4; ++i)
    for(int j=0;j<4;++j)
      gs_out.prob[i][j]=0;
  for(int i=0; i<4; ++i)
  {
    for(int j=0; j<4; ++j)
    {
      if(gs_out.color[i] == gs_out.color[j])
      {
        gs_out.prob[i][j] = 1.0;
      }
    }
  }
  EmitVertex();

  EndPrimitive();
}
