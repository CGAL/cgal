#version 150

in VS_OUT
{
  vec4 fP;
  vec4 out_color;
  float dist[6];
  vec4 vertex;
  vec3 normal;
  flat float stat_value;
}fs_in;

uniform bool is_two_side;
uniform vec4 light_pos;
uniform vec4 light_diff;
uniform vec4 light_spec;
uniform vec4 light_amb;
uniform float spec_power ;
uniform bool is_clipbox_on;
uniform float min_threshold;
uniform float max_threshold;

out vec4 out_color;

void main(void)
{
  if(fs_in.stat_value < min_threshold
     || fs_in.stat_value > max_threshold)
    discard;
  if(is_clipbox_on)
    if(fs_in.dist[0]>0.0 ||
       fs_in.dist[1]>0.0 ||
       fs_in.dist[2]>0.0 ||
       fs_in.dist[3]>0.0 ||
       fs_in.dist[4]>0.0 ||
       fs_in.dist[5]>0.0)
      discard;
  vec3 L = light_pos.xyz - fs_in.fP.xyz;
  vec3 V = -fs_in.fP.xyz;

  vec3 N;
  if(fs_in.normal == vec3(0.0,0.0,0.0))
       N = vec3(0.0,0.0,0.0);
   else
       N = fs_in.normal;
  L = normalize(L);
  V = normalize(V);

  vec3 R = reflect(-L, N);
  vec4 diffuse;
  if(!is_two_side)
      diffuse = max(dot(N,L),0) * light_diff*fs_in.out_color;
  else
      diffuse = max(abs(dot(N,L)),0) * light_diff*fs_in.out_color;
  vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;
  out_color = vec4((fs_in.out_color*light_amb).xyz + diffuse.xyz + specular.xyz,1.0);
}
