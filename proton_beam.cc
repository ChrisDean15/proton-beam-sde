#include "material.cc"
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <vector>

#ifndef PB
#define PB

struct proton_path {

  proton_path(const double e0, const std::vector<double> x0,
              const std::vector<double> w0, const int n)
      : energy(n, e0), s(n, 0), x(n, x0), omega(n, w0), u(3, 0), z(3, 0),
        w(3, 0) {}

  void reset(const double e0, const std::vector<double> x0,
             const std::vector<double> w0) {
    x[0] = x0;
    omega[0] = w0;
    energy[0] = e0;
    s[0] = 0;
    return;
  }

  double log_a(const int k, const int m, const double theta) const {
    double ret = log(theta + 2 * k - 1) + gsl_sf_lnpoch(theta + m, k - 1) -
                 gsl_sf_lnfact(m) - gsl_sf_lnfact(k - m);
    return ret;
  }

  double b(const int k, const int m, const double t, const double theta) const {
    double ret = 1;
    if (k > 0) {
      ret = exp(log_a(k, m, theta) - k * (k + theta - 1) * t / 2);
    }
    return ret;
  }

  int c(const int m, const double t, const double theta) const {
    int i = 0;
    double b_curr = b(m, m, t, theta);
    double b_next = b(m + 1, m, t, theta);
    while (b_next >= b_curr) {
      i++;
      b_curr = b_next;
      b_next = b(i + m + 1, m, t, theta);
    }
    return i;
  }

  int number_of_blocks(const double t, gsl_rng *gen) const {
    int m = 0;
    double theta = 1;
    if (t < 0.07) {
      double mu = 2 / t;
      double sigma = sqrt(2 / (3 * t));
      m = round(mu + sigma * gsl_ran_gaussian_ziggurat(gen, 1));
    } else {
      std::vector<int> k(1, 0);
      bool proceed = true;
      double u = gsl_rng_uniform_pos(gen);
      double smin = 0, smax = 0, increment = 0;
      while (proceed) {
        k[m] = ceil(c(m, t, theta) / 2.0);
        for (int i = 0; i < k[m]; i++) {
          increment = b(m + 2 * i, m, t, theta) - b(m + 2 * i + 1, m, t, theta);
          smin += increment;
          smax += increment;
        }
        increment = b(m + 2 * k[m], m, t, theta);
        smin += increment - b(m + 2 * k[m] + 1, m, t, theta);
        smax += increment;
        while (smin < u && u < smax) {
          for (int i = 0; i <= m; i++) {
            k[i]++;
            increment = b(i + 2 * k[i], i, t, theta);
            smax = smin + increment;
            smin += increment - b(i + 2 * k[i] + 1, i, t, theta);
          }
        }
        if (smin > u) {
          proceed = false;
        } else {
          k.push_back(0);
          m++;
        }
      }
    }
    return m;
  }

  double wright_fisher(const double r, gsl_rng *gen) const {
    double y;
    if (r > 1e-9) {
      int m = number_of_blocks(r, gen);
      y = gsl_ran_beta(gen, 1, 1 + m);
    } else {
      y = r / 2;
      y = fabs(gsl_ran_gaussian_ziggurat(gen, sqrt(r * y * (1 - y))));
    }
    return y;
  }

  double dist(const std::vector<double> &y, const std::vector<double> &z) {
    double ret = 0;
    for (unsigned int i = 0; i < y.size(); i++) {
      ret += (z[i] - y[i]) * (z[i] - y[i]);
    }
    return sqrt(ret);
  }
  void spherical_bm(const double time_increment, int &ix, gsl_rng *gen,
                    const Material &mat) {
    z[0] = sin(omega[ix - 1][0]) * cos(omega[ix - 1][1]);
    z[1] = sin(omega[ix - 1][0]) * sin(omega[ix - 1][1]);
    z[2] = cos(omega[ix - 1][0]);
    double y = wright_fisher(
        pow(mat.multiple_scattering_sd(energy[ix - 1], time_increment), 2),
        gen);

    double theta = 2 * M_PI * gsl_rng_uniform(gen);
    // Set up defaults for when z is near (0, 0, 1)
    u[0] = 1;
    u[1] = 1;
    u[2] = 1;
    if (z[0] > 0) {
      u[0] = -1;
    }
    if (z[1] > 0) {
      u[1] = -1;
    }
    double denom = sqrt(z[0] * z[0] + z[1] * z[1] + (z[2] - 1) * (z[2] - 1));
    if (denom > 1e-10) {
      u[0] = -z[0] / denom;
      u[1] = -z[1] / denom;
      u[2] = (1 - z[2]) / denom;
    }
    z[0] = 2 * sqrt(y * (1 - y)) * cos(theta);
    z[1] = 2 * sqrt(y * (1 - y)) * sin(theta);
    z[2] = 1 - 2 * y;
    w[0] = (1 - 2 * u[0] * u[0]) * z[0] - 2 * u[0] * u[1] * z[1] -
           2 * u[0] * u[2] * z[2];
    w[1] = (1 - 2 * u[1] * u[1]) * z[1] - 2 * u[0] * u[1] * z[0] -
           2 * u[1] * u[2] * z[2];
    w[2] = (1 - 2 * u[2] * u[2]) * z[2] - 2 * u[0] * u[2] * z[0] -
           2 * u[1] * u[2] * z[1];
    if (ix == int(omega.size())) {
      omega.resize(2 * omega.size(), omega.back());
      energy.resize(2 * energy.size(), energy.back());
      x.resize(2 * x.size(), x.back());
      s.resize(2 * s.size(), s.back());
    }
    omega[ix][0] = acos(w[2]);
    omega[ix][1] = atan2(w[1], w[0]);
    // if(omega[ix][1]<0) {
    //  omega[ix][1]+=2*M_PI;
    // }
    double v0 = omega[ix - 1][0];
    double v1 = omega[ix][0];
    double w0 = omega[ix - 1][1];
    double w1 = omega[ix][1];
    double dt = time_increment;

    // Check for division by zero in x and y position updates
    double denom_xy = (v0 - v1 + w0 - w1) * (v0 - v1 - w0 + w1);
    if (fabs(denom_xy) > 1e-9) {
      x[ix][0] = x[ix - 1][0] -
                 dt *
                     ((v0 - v1) * (cos(v0) * cos(w0) - cos(v1) * cos(w1)) +
                      (w0 - w1) * (sin(v0) * sin(w0) - sin(v1) * sin(w1))) /
                     denom_xy;
      x[ix][1] = x[ix - 1][1] +
                 dt *
                     ((w0 - w1) * (cos(w0) * sin(v0) - cos(w1) * sin(v1)) -
                      (v0 - v1) * (cos(v0) * sin(w0) - cos(v1) * sin(w1))) /
                     denom_xy;
    } else {
      // Linear approximation when denominator is too small
      x[ix][0] = x[ix - 1][0] + dt * sin(v0) * cos(w0);
      x[ix][1] = x[ix - 1][1] + dt * sin(v0) * sin(w0);
    }

    // Z position update
    if (fabs(v0 - v1) > 1e-9) {
      x[ix][2] = x[ix - 1][2] + (sin(v0) - sin(v1)) * dt / (v0 - v1);
    } else {
      x[ix][2] = x[ix - 1][2] - dt * (cos(v0) + cos(v1)) / 2;
    }
    energy[ix] = energy[ix - 1] -
                 fmax(mat.bethe_bloch(energy[ix - 1]) * dt +
                          sqrt(dt) * mat.energy_straggling_sd(energy[ix - 1]) *
                              gsl_ran_gaussian_ziggurat(gen, 1),
                      0);
    energy[ix] = fmax(energy[ix], 0);
    s[ix] = energy[ix - 1] - energy[ix];
    ix++;
    return;
  }

  int simulate(const double dt, const double absorption_energy,
               const double air_gap, gsl_rng *gen,
               std::vector<Material> &materials) {
    double nonelastic_jump_rate, rutherford_elastic_jump_rate, alpha, u;
    int ix = 1;
    int material_index = 0;
    int crossed = 0;
    while (energy[ix - 1] > absorption_energy) {
      if (crossed == 0 && x[ix - 1][0] >= air_gap) {
        material_index++;
        crossed = 1;
      }
      spherical_bm(dt, ix, gen, materials[material_index]);
      if (energy[ix - 1] > absorption_energy) {
        nonelastic_jump_rate =
            materials[material_index].nonelastic_rate(energy[ix - 1]);
        // 2 * materials[material_index].multiple_scattering_sd(energy[ix - 1],
        //        dt);

        rutherford_elastic_jump_rate =
            materials[material_index].rutherford_and_elastic_rate(
                energy[ix - 1]);
        alpha = rutherford_elastic_jump_rate + nonelastic_jump_rate;
        if (gsl_rng_uniform(gen) < 1 - exp(-alpha * dt)) {
          u = gsl_rng_uniform(gen);
          if (u < rutherford_elastic_jump_rate / alpha) {
            materials[material_index].rutherford_elastic_scatter(
                omega[ix - 1], energy[ix - 1], gen);
          } else {
            materials[material_index].nonelastic_scatter(omega[ix - 1],
                                                         energy[ix - 1], gen);
            s[ix - 1] = (energy[ix - 2] - energy[ix - 1]);
          }
        }
      }
    }
    return ix;
  }

  void print() {
    for (unsigned int i = 0; i < x.size(); i++) {
      std::cout << x[i][0] << " " << x[i][1] << " " << x[i][2] << " "
                << omega[i][0] << " " << omega[i][1] << " " << energy[i] << " "
                << s[i] << std::endl;
    }
    return;
  }

  std::vector<double> energy, s;
  std::vector<std::vector<double>> x, omega;
  // Dummy vectors for spherical BM
  std::vector<double> u, z, w;
};

#endif
