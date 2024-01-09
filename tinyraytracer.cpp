#include "geometry.h" // NOTICE: the class of vectors lives inside the geometry.h file
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector> // NOTICE: Use of cpp vector

// Create the material struct
struct Material {
  Material(const Vec2f &a, const Vec3f &color, const float &spec)
      : albedo(a), diffuse_color(color), specular_exponent(spec) {}
  Material() : albedo(1, 0), diffuse_color(), specular_exponent() {}
  Vec2f albedo;
  Vec3f diffuse_color;
  float specular_exponent;
};

// Create the light struct
//
struct Light {
  Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
  Vec3f position;
  float intensity;
};

// Create a Sphere struct which implements functions to understand ray
// intersections We need only four numbers to describe a sphere in memeory a
// three-dimensional vector for the centre and a scale for the radius
//
//
struct Sphere {
  Vec3f center;
  float radius;
  Material material;

  Sphere(const Vec3f &c, const float &r, const Material &m)
      : center(c), radius(r), material(m) {}

  bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) {
    // This function is based of the Ray-Sphere Intersection algorithm
    //
    Vec3f L = center - orig;
    float tca = L * dir;
    float d2 = L * L - tca * tca;
    if (d2 > radius * radius)
      return false;
    float thc = sqrtf(radius * radius - d2);
    t0 = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0)
      t0 = t1;
    if (t0 < 0)
      return false;
    return true;
  }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) { return I - N * 2.f * (I * N); }

bool scene_interesect(const Vec3f &orig, const Vec3f &dir,
                      std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N,
                      Material &material) {
  float spheres_dist = std::numeric_limits<float>::max();
  for (size_t i = 0; i < spheres.size(); i++) {
    float dist_i;
    if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
      spheres_dist = dist_i;
      hit = orig + dir * dist_i;
      N = (hit - spheres[i].center).normalize();
      material = spheres[i].material;
    }
  }
  return spheres_dist < 1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir,
               std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
  Vec3f point, N;
  Material material;

  if (!scene_interesect(orig, dir, spheres, point, N, material)) {
    return Vec3f(0.2, 0.7, 0.8); // Background color
  }

  float diffuse_light_intensity = 0, specular_light_intensity = 0;
  for (size_t i = 0; i < lights.size(); i++) {
    Vec3f light_dir = (lights[i].position - point).normalize();
    //SHADOWS:
    float light_distance = (lights[i].position - point).norm(); 
    Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; //Checking if the point lies in the shadow of the lights[i]
    Vec3f shadow_pt, shadow_N; 
    Material tmp_material; 
    if (scene_interesect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmp_material) && (shadow_pt - shadow_orig).norm() < light_distance) continue; //In other words, we don't render the diffusion and specular light if this is the case
    //diffuse light 
    diffuse_light_intensity +=
        lights[i].intensity * std::max(0.f, light_dir * N);
    //Speculatar light
    specular_light_intensity +=
        pow(std::max(0.f, -reflect(-light_dir, N) * dir),
            material.specular_exponent) *
        lights[i].intensity;
  }

  return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity*material.albedo[1];
}

void render(std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
  const int width = 1024;
  const int height = 768;
  const float fov = M_PI / 2;
  std::vector<Vec3f> framebuffer(width * height);

  for (size_t j = 0; j < height; j++) {
    for (size_t i = 0; i < width; i++) {
      float x = (2 * (i + 0.5) / (float)width - 1) * tan(fov / 2.) * width /
                (float)height;
      (float)height;
      float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
      Vec3f dir = Vec3f(x, y, -1).normalize();
      framebuffer[i + j * width] =
          cast_ray(Vec3f(0, 0, 0), dir, spheres, lights);
    }
  }

  std::ofstream ofs; // save the framebuffer to file
  ofs.open("./out.ppm");
  ofs << "P6\n" << width << " " << height << "\n255\n";
  for (size_t i = 0; i < height * width; ++i) {
    for (size_t j = 0; j < 3; j++) {
      Vec3f &c = framebuffer[i];
      float max = std::max(c[0], std::max(c[1], c[2]));
      if (max > 1)
        c = c * (1. / max);
      ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
    }
  }
  ofs.close();
}

int main() {
  Material ivory(Vec2f(0.6, 0.3), Vec3f(0.4, 0.4, 0.3), 50.);
  Material red_rubber(Vec2f(0.9, 0.1), Vec3f(0.3, 0.1, 0.1), 10.);

  // Create our example spheres
  std::vector<Sphere> spheres;
  spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
  spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber));
  spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red_rubber));
  spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, ivory));

  // Create the light
  std::vector<Light> lights;
  lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));
  lights.push_back(Light(Vec3f(30, 50, -25), 1.8));
  lights.push_back(Light(Vec3f(30, 20, 30), 1.7));

  // Rendering:
  render(spheres, lights);

  std::cout << "Rendering finished.\n";

  return 0;
}
