#include <iostream>
#include "random.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_video.h>
#include <thread>
#include <mutex>
#include "vector.h"
#include "world.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"

/// Semaphore
struct Semaphore {
private:
  size_t value = 0;
  std::mutex mut;

public:
  explicit Semaphore(size_t val): value(val) {}
  
  void post(size_t val = 1) {
    std::unique_lock<std::mutex> lk(mut);
    value += val;
  }
  
  size_t get_value() {
    std::unique_lock<std::mutex> lk(mut);
    return value;
  }
  
  bool try_wait(size_t* val) {
    std::unique_lock<std::mutex> lk(mut);
    if (value == 0) {
      *val = value;
      return false;
    } else {
      value--;
      *val = value;
      return true;
    }
  }
};

Vec3<> color(const Ray& r, const World& world, const int depth) {
  Hit hit;
  if (world.hit(r, 0.001, MAXFLOAT, hit)) {
    Ray scattered;
    Vec3<> attenuation;
    Vec3<> emission = hit.mat->emitted(hit.u, hit.v, hit.p);
    if (depth < 50 && hit.mat->scatter(r, hit, attenuation, scattered)) {
      // TODO: Shadow rays
      return emission + attenuation*color(scattered, world, depth + 1);
    } else {
      return emission;
    }
  } else {
    return Vec3<>{0.0};
    Vec3<> dir = r.direction().normalized();
    double t = 0.5*(dir.y + 1.0);
    return (1.0 - t) * Vec3<>{1.0, 1.0, 1.0} + t * Vec3<>{0.5, 0.7, 1.0};
  }
}

struct Region {
  const size_t nx, ny;
  size_t x0, x1;
  size_t y0, y1;
};

struct ThreadArgs {
  const World& world;
  const Camera& cam;
  /// Screen buffer
  uint32_t* pixels;
  /// Number of samples per pixels
  int ns;
  /// All work items (a.k.a Regions of the screen)
  const std::vector<Region>& regions;
  /// Number of work items remaining
  Semaphore& sem;
  /// Number of threads done
  Semaphore& sem_done;
};

void thread_work(ThreadArgs args) {
  size_t row = 0;
  while (args.sem.try_wait(&row)) {
    const Region& reg = args.regions[row];
    for (size_t j = reg.y0; j <= reg.y1; j++) {
      for (size_t i = reg.x0; i <= reg.x1; i++) {
        Vec3<> t_color{};
        for (int s = 0; s < args.ns; s++) {
          double u = double(i + drand48()) / double(reg.nx);
          double v = double(j + drand48()) / double(reg.ny);
          Ray r = args.cam.get_ray(u, v);
          t_color += color(r, args.world, 0);
        }
        if (t_color == NAN) { std::cerr << "NAN" << std::endl; }
        if (t_color == INFINITY) { std::cerr << "INF" << std::endl; }
        t_color /= double(args.ns);
        t_color = {std::sqrt(t_color.x), std::sqrt(t_color.y), std::sqrt(t_color.z)}; // Gamma-2 correction
        auto ir = uint32_t(t_color.x * 255);
        auto ig = uint32_t(t_color.y * 255);
        auto ib = uint32_t(t_color.z * 255);
        auto ia = uint32_t(1);
        uint32_t pixel = 0;
        pixel += (ia << (8 * 3));
        pixel += (ir << (8 * 2));
        pixel += (ig << (8 * 1));
        pixel += (ib << (8 * 0));
        args.pixels[((reg.ny - j) * reg.nx) + i] = pixel;
      }
    }
  }
  args.sem_done.post();
}

int main() {
  SDL_Init(SDL_INIT_EVERYTHING);

#ifdef __LINUX__
  init_random_generator();
#endif

  const size_t nx = 720;
  const size_t ny = 400;
  const size_t ns = 10; // Number of samples per pixel
  
  SDL_Window* window = SDL_CreateWindow("RayTracer", 0, 0, nx, ny, 0);
  SDL_Surface* scr = SDL_GetWindowSurface(window);
  uint32_t* pixels = (uint32_t*) scr->pixels;
  
  Vec3<> lookfrom = {-2, 2, 1};
  Vec3<> lookat = {0, 0, -1};
  double dist_to_focus = (lookfrom - lookat).length();
  Camera cam{lookfrom, lookat, 60, double(nx) / double(ny), 0.1, dist_to_focus};
  
  World world;
  world.default_scene();
  world.bake_world();
  
  const size_t num_threads = std::thread::hardware_concurrency() == 0 ? 4 : std::thread::hardware_concurrency();
  std::cout << "Starting " << num_threads << " number of threads." << std::endl;
  std::vector<std::thread> threads{};
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<Region> regions{};
  const size_t num_rows = ny / num_threads;
  std::cout << " - Number of rows per thread: " << num_rows << std::endl;
  for (size_t i = 0; i < size_t(ny / num_rows); i++) {
    regions.emplace_back(Region{nx, ny, 0, nx, i * num_rows, (i + 1) * num_rows});
  }
  Semaphore sem{0};
  Semaphore sem_done{0};
  sem.post(regions.size());
  for (size_t i = 0; i < num_threads; i++) {
    ThreadArgs args{world, cam, pixels, ns, regions, sem, sem_done};
    threads.emplace_back(std::thread{thread_work, args});
  }
  
  bool quit = false;
  bool finished = false;
  while (!quit) {
    SDL_Event event{};
    while (SDL_PollEvent(&event)) {
      if (event.key.keysym.sym == SDLK_ESCAPE || event.type == SDL_QUIT) {
          quit = true;
          break;
      }
      switch (event.type) {
        case SDL_KEYDOWN:
          switch (event.key.keysym.sym) {
            case SDLK_a:
              // TODO: Camera movemments
              break;
          }
          break;
      }
    }
    SDL_UpdateWindowSurface(window);
    
    if (sem_done.get_value() == num_threads && !finished) {
      finished = true;
      auto end = std::chrono::high_resolution_clock::now();
      auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
      std::cout << " - " << diff << " ms, " << (ns*nx*ny)/(diff/1000.0)/1'000'000 << " Mrays/s" << std::endl;
    }
  }
  // FIXME: No way to async terminate a thread (except via pthread APIs).
  for (auto& t : threads) { t.join(); }
  return EXIT_SUCCESS;
}