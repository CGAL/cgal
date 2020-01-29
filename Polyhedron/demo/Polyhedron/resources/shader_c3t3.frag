#version 150
in vec4 color;
in vec4 fP; 
in vec3 fN; 
uniform vec4 light_pos;  
uniform vec4 light_diff; 
uniform vec4 light_spec; 
uniform vec4 light_amb;  
uniform float spec_power ;
uniform int is_two_side; 
uniform bool is_selected;
uniform float near;
uniform float far;
uniform float width;
uniform float height;
uniform bool comparing;
uniform bool writing;
uniform sampler2D sampler;
uniform float alpha;
uniform bool is_surface;
out  vec4 out_color;

float depth(float z)
{
  return (2 * near) / (far + near - z * (far - near));
}

void main(void) {
  float d = depth(gl_FragCoord.z);
  float test = texture(sampler, vec2(gl_FragCoord.x/width, gl_FragCoord.y/height)).r;
  if(comparing && d <= test)
  discard;
  if(writing)
    out_color = vec4(d,d,d,1.0);
  else
  {
    if(color.w<0 || is_surface)
    {

      vec4 my_color = vec4(color.xyz, 1.);
      vec3 L = light_pos.xyz - fP.xyz; 
      vec3 V = -fP.xyz; 
      vec3 N; 
      if(fN == vec3(0.0,0.0,0.0))
      {
        out_color = my_color;
        return;
      }
      else 
          N = normalize(fN); 
      L = normalize(L); 
      V = normalize(V); 
      vec3 R = reflect(-L, N); 
      vec4 diffuse; 
      if(is_two_side == 1) 
          diffuse = abs(dot(N,L)) * light_diff * my_color; 
      else 
          diffuse = max(dot(N,L), 0.0) * light_diff * my_color; 
      vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; 
      vec4 ret_color = vec4((my_color*light_amb).xyz + diffuse.xyz + specular.xyz,1); 
      if(is_selected)
          out_color = vec4(ret_color.r+70.0/255.0, ret_color.g+70.0/255.0, ret_color.b+70.0/255.0, alpha);
      else
          out_color = vec4(ret_color.xyz, alpha);
    }
    else
      discard;
    }
}


