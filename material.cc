#include "cross_sections.cc"
#include <cstdlib>
#include <vector>

#ifndef MAT
#define MAT

struct Atom {

  Atom(const int a0, const int z0, const std::string el_r,
       const std::string ne_r, const std::string ne_y, const std::string el_a,
       const std::string ne_a, const std::string ne_e)
      : a(a0), z(z0), el_rate(el_r), ne_rate(ne_r), ne_yield(ne_y),
        el_angle_cdf(el_a), ne_angle_cdf(ne_a), ne_energy_cdf(ne_e) {}

  // Constructor for zero non-elastic rate for hydrogen
  Atom(const int a0, const int z0, const std::string el_r,
       const std::string el_a)
      : a(a0), z(z0), el_rate(el_r), ne_rate(), ne_yield(), el_angle_cdf(el_a),
        ne_angle_cdf(), ne_energy_cdf() {}

  Atom(const Atom &other)
      : a(other.a), z(other.z), el_rate(other.el_rate), ne_rate(other.ne_rate),
        ne_yield(other.ne_yield), el_angle_cdf(other.el_angle_cdf),
        ne_angle_cdf(other.ne_angle_cdf), ne_energy_cdf(other.ne_energy_cdf) {}

  double cm_to_lab_frame(const double ang, const double e0,
                         const double e_delta) const {
    double mp = 938.346; // mass of proton * c^2, MeV
    double mn = mp * a;  // mass of colliding nucleus * c^2, MeV
    double u = sqrt(e0 * (e0 + mp)) / (e0 + mp + mn);
    double g = 1 / sqrt(1 - u * u);
    double p = sqrt((e0 - e_delta) * (e0 - e_delta + 2 * mp));
    double e = e0 - e_delta + mp;
    double v_ratio = u * (e - u * p) / (p - u * e);
    return atan(sin(ang) / (g * (cos(ang) + v_ratio)));
  }

  const int a, z;
  CS_1d el_rate, ne_rate, ne_yield;
  CS_2d el_angle_cdf, ne_angle_cdf;
  CS_3d ne_energy_cdf;
};

struct Material {

  Material(std::vector<Atom> &atoms, const std::vector<int> &id,
           const std::vector<double> &x0, const double d0, const double I0)
      : density(d0), I(I0 / 1e6), x(x0), at() {
    for (unsigned int i = 0; i < id.size(); i++) {
      at.push_back(atoms[id[i]]);
    }
  }

  Material() : density(), I(), x(), at() {}

  void read_material(const std::string filename, std::vector<Atom> &atoms) {
    std::ifstream file;
    file.open(filename + ".txt");
    std::string line, token;
    getline(file, line);
    density = atof(line.c_str());
    getline(file, line);
    I = atof(line.c_str()) / 1e6;
    while (getline(file, line)) {
      std::stringstream iss;
      iss << line;
      getline(iss, token, ' ');
      at.push_back(atoms[atoi(token.c_str())]);
      getline(iss, token, ' ');
      x.push_back(atof(token.c_str()));
    }
    file.close();
    return;
  }

  double bethe_bloch(const double e) const {
    double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double num = 0;
    double denom = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      num += x[i] * 0.3072 * at[i].z * density *
             (log(2 * mecsq * betasq / I * (1 - betasq)) - betasq) /
             betasq; // MeV / cm
      denom += at[i].a * x[i];
    }
    return num / denom;
  }

  double multiple_scattering_sd(const double e, const double dt) const {
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double pv = (2 * mpcsq + e) * e / (mpcsq + e);
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double c = 29979245800; // speed of light
    double vel = sqrt(betasq) * c;
    double p = pv / sqrt(betasq); // momentum in MeV / c.
    // effective chi_c_sq is just the sum of individual elements
    double chi_c_sq = 0;
    std::vector<double> chi_a_sq_vec(at.size());
    for (unsigned int i = 0; i < at.size(); i++) {
      chi_c_sq += x[i] * at[i].z * (at[i].z + 1.0) / at[i].a;
      chi_a_sq_vec[i] = 2.007e-5 * pow(at[i].z, 2 / 3) *
                        (1 + 3.34 * pow(at[i].z / (137 * vel), 2)) / (p * p);
    }
    chi_c_sq *= 0.157 * dt * density / (pv * pv);
    // effective chi_a_sq is a weighted average on the log-scale
    double chi_a_sq = 0;
    double denom = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      chi_a_sq +=
          x[i] * at[i].z * (at[i].z + 1) * log(chi_a_sq_vec[i]) / at[i].a;
      denom += x[i] * at[i].z * (at[i].z + 1.0) / at[i].a;
    }
    chi_a_sq = exp(chi_a_sq / denom);
    double omega = chi_c_sq / chi_a_sq;
    double F = 0.98;
    double v = omega / (2 * (1 - F));
    double ret = sqrt(chi_c_sq * ((1 + v) * log(1 + v) / v - 1)) / (1 + F * F);
    ret *= dt / sqrt(ret * ret + dt * dt);
    return ret;
  }

  double energy_straggling_sd(const double e) const {
    double alpha = 1 / 137.0;
    double log_hbar = -21 * log(10) + log(4.136) - log(2 * M_PI); // MeV * s
    double log_c = log(29979245800);                              // cm / s
    double log_avogadro = log(6) + 23 * log(10);
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    double a = 0;
    double z = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a; // average molar mass
      z += x[i] * at[i].z; // electrons per average molecule
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    double ret =
        4 * M_PI * z * (1 - betasq / 2) / sqrt(1 - betasq) *
        exp(2 * (log(alpha) + log_hbar + log_c) + log_molecule_density);
    return sqrt(ret);
  }

  double nonelastic_rate(const double e) const {
    double log_avogadro = log(6) + 23 * log(10);
    double log_barns_to_cmsq = -24 * log(10);
    double a = 0;
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a; // average molar mass
      ret += at[i].a * x[i] * at[i].ne_rate.evaluate(e);
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density) / a;
    return ret; // rate per cm
  }

  double elastic_rate(const double e) const {
    double log_avogadro = log(6) + 23 * log(10);
    double log_barns_to_cmsq = -24 * log(10);
    double a = 0;
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a; // average molar mass
      ret += at[i].a * x[i] * at[i].el_rate.evaluate(e);
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density) / a;
    return ret; // rate per cm
  }

  double rutherford_rate(const double e, const double lb) const {
    double log_ahbarc = log(197.3 / 137) - 13 * log(10);
    double log_avogadro = log(6) + 23 * log(10);
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double pv = (2 * mpcsq + e) * e / (mpcsq + e);
    double a = 0;
    double az = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      a += x[i] * at[i].a; // average molar mass
      az += x[i] * at[i].a * at[i].z * at[i].z;
    }
    double log_molecule_density =
        log(density) + log_avogadro - log(a); // molecules / cm^3
    double sig =
        exp(2 * (log_ahbarc + log(cos(lb / 2)) - log(pv) - log(sin(lb / 2))) +
            log_molecule_density) *
        M_PI;
    double ret = sig * az / a;
    return ret;
  }

  void nonelastic_scatter(std::vector<double> &ang, double &e, gsl_rng *gen) {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double rate = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      rate += at[i].a * x[i] * at[i].ne_rate.evaluate(e);
    }
    double u = gsl_rng_uniform(gen);
    double ind = 0;
    double tmp = at[ind].a * x[ind] * at[ind].ne_rate.evaluate(e) / rate;
    while (tmp < u) {
      ind++;
      tmp += at[ind].a * x[ind] * at[ind].ne_rate.evaluate(e) / rate;
    }
    double alpha = 0;
    double e_old = e;
    if (gsl_rng_uniform(gen) > at[ind].ne_yield.evaluate(e)) {
      // proton absorbed & track ends
      e = 0;
    } else {
      alpha = at[ind].ne_angle_cdf.sample(e, gen);
      e = at[ind].ne_energy_cdf.sample(e, alpha, gen);
      alpha = at[ind].cm_to_lab_frame(alpha, e_old, e_old - e);
      if (0 <= beta && beta < M_PI / 2) {
        ang[0] -= atan(sin(beta) * tan(alpha));
        ang[1] -= atan(cos(beta) * tan(alpha));
      } else if (M_PI / 2 <= beta && beta < M_PI) {
        ang[0] -= atan(sin(M_PI - beta) * tan(alpha));
        ang[1] += atan(cos(M_PI - beta) * tan(alpha));
      } else if (M_PI <= beta && beta < 3 * M_PI / 2) {
        ang[0] += atan(sin(3 * M_PI / 2 - beta) * tan(alpha));
        ang[1] += atan(cos(3 * M_PI / 2 - beta) * tan(alpha));
      } else {
        ang[0] += atan(sin(2 * M_PI - beta) * tan(alpha));
        ang[1] -= atan(cos(2 * M_PI - beta) * tan(alpha));
      }
    }
    return;
  }

  void elastic_scatter(std::vector<double> &ang, const double e, gsl_rng *gen) {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double rate = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      rate += at[i].a * x[i] * at[i].el_rate.evaluate(e);
    }
    double u = gsl_rng_uniform(gen);
    double ind = 0;
    double tmp = at[ind].a * x[ind] * at[ind].el_rate.evaluate(e) / rate;
    while (tmp < u) {
      ind++;
      tmp += at[ind].a * x[ind] * at[ind].el_rate.evaluate(e) / rate;
    }
    double alpha = at[ind].el_angle_cdf.sample(e, gen);
    alpha = at[ind].cm_to_lab_frame(alpha, e, 0);
    if (0 <= beta && beta < M_PI / 2) {
      ang[0] -= atan(sin(beta) * tan(alpha));
      ang[1] -= atan(cos(beta) * tan(alpha));
    } else if (M_PI / 2 <= beta && beta < M_PI) {
      ang[0] -= atan(sin(M_PI - beta) * tan(alpha));
      ang[1] += atan(cos(M_PI - beta) * tan(alpha));
    } else if (M_PI <= beta && beta < 3 * M_PI / 2) {
      ang[0] += atan(sin(3 * M_PI / 2 - beta) * tan(alpha));
      ang[1] += atan(cos(3 * M_PI / 2 - beta) * tan(alpha));
    } else {
      ang[0] += atan(sin(2 * M_PI - beta) * tan(alpha));
      ang[1] -= atan(cos(2 * M_PI - beta) * tan(alpha));
    }
    return;
  }

  void rutherford_scatter(std::vector<double> &ang, const double lb,
                          gsl_rng *gen) {
    double beta = 2 * M_PI * gsl_rng_uniform(gen);
    double u = gsl_rng_uniform(gen);
    double alpha = acos((cos(lb) - u * pow(cos(lb / 2), 2)) /
                        (1 - u * pow(cos(lb / 2), 2)));
    if (0 <= beta && beta < M_PI / 2) {
      ang[0] -= atan(sin(beta) * tan(alpha));
      ang[1] -= atan(cos(beta) * tan(alpha));
    } else if (M_PI / 2 <= beta && beta < M_PI) {
      ang[0] -= atan(sin(M_PI - beta) * tan(alpha));
      ang[1] += atan(cos(M_PI - beta) * tan(alpha));
    } else if (M_PI <= beta && beta < 3 * M_PI / 2) {
      ang[0] += atan(sin(3 * M_PI / 2 - beta) * tan(alpha));
      ang[1] += atan(cos(3 * M_PI / 2 - beta) * tan(alpha));
    } else {
      ang[0] += atan(sin(2 * M_PI - beta) * tan(alpha));
      ang[1] -= atan(cos(2 * M_PI - beta) * tan(alpha));
    }
    return;
  }

  void print() const {
    std::cout << "density = " << density << std::endl;
    std::cout << "I = " << I << std::endl;
    for (unsigned int i = 0; i < at.size(); i++) {
      std::cout << x[i] << " " << at[i].z << " " << at[i].a << std::endl;
    }
    return;
  }

  double density; // density, g / cm^3
  double I;       // mean excitation energy, MeV
  std::vector<double> x;
  std::vector<Atom> at;
};

#endif
