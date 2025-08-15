#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef CS
#define CS
// necessary coef for proton-proton is ia = oa = 1 iz = oz = 0,
// this class is incident/outgoing particle dependent generally,
// preparing setup for general case
// defining here instead of materials as
// appears in ENDF cross section struct as input to member function
struct AtomConsts {
  const int a, z;
  const double Mi, Mo, ANU, BNU, iNU, oNU, Si, So;
  double AngleScatterCoef(const int AC, const int NC, const int ZC,
                          const int AX, const int NX, const int ZX,
                          const double sep) const {
    double temp = 15.68 * (AC - AX) -
                  (28.07 * ((pow(NC - ZC, 2) / static_cast<double>(AC)) -
                            (pow(NX - ZX, 2) / static_cast<double>(AX)))) -
                  (18.56 * (pow(AC, 2.0 / 3.0) - pow(AX, 2.0 / 3.0))) +
                  (33.22 * ((pow(NC - ZC, 2) / pow(AC, 4.0 / 3.0)) -
                            (pow(NX - ZX, 2) / pow(AX, 4.0 / 3.0)))) -
                  (0.717 * ((pow(ZC, 2) / pow(AC, 1.0 / 3.0)) -
                            (pow(ZX, 2) / pow(AX, 1.0 / 3.0)))) +
                  (1.211 * ((pow(ZC, 2) / static_cast<double>(AC)) -
                            (pow(ZX, 2) / static_cast<double>(AX)))) -
                  sep;
    return temp;
  }
  AtomConsts(const int a0, const int z0, const int ia, const int iz,
             const int oa, const int oz, const int mi, const int mo,
             const double sepi, const double sepo)
      : a(a0), z(z0), Mi(mi), Mo(mo),
        ANU(((a0 - z0) * 1.0087 + z0 * 1.0073) / 1.0087),
        BNU(((a0 + ia - oa - z0 - iz + oz) * 1.0087 + (z0 + iz - oz) * 1.0073) /
            1.0087),
        iNU(((ia - iz) * 1.0087 + iz * 1.0073) / 1.0087),
        oNU(((oa - oz) * 1.0087 + oz * 1.0073) / 1.0087),
        Si(AngleScatterCoef(a0 + ia, a0 + ia - z0 - iz, z0 + iz, a0, a0 - z0,
                            z0, sepi)),
        So(AngleScatterCoef(a0 + oa, a0 + oa - z0 - oz, z0 + oz, a0, a0 - z0,
                            z0, sepo)) {}

  AtomConsts(const AtomConsts &other)
      : a(other.a), z(other.z), Mi(other.Mi), Mo(other.Mo), ANU(other.ANU),
        BNU(other.BNU), iNU(other.iNU), oNU(other.oNU), Si(other.Si),
        So(other.So) {}
  AtomConsts() : a(), z(), Mi(), Mo(), ANU(), BNU(), iNU(), oNU(), Si(), So() {}
};
struct CS_1d {

  CS_1d(const std::string filename) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    std::stringstream iss;

    getline(file, line);
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }

    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      rate.push_back(atof(token.c_str()));
    }
  }

  CS_1d(const CS_1d &other) : energy(other.energy), rate(other.rate) {}

  CS_1d() : energy(), rate() {}

  void print() const {
    for (unsigned int i = 0; i < energy.size(); i++) {
      std::cout << energy[i] << " ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < rate.size(); i++) {
      std::cout << rate[i] << " ";
    }
    std::cout << std::endl;
    return;
  }

  double evaluate(const double e) const {
    int r;
    double ret = 0;
    if (energy.size() > 0) {
      if (e <= energy[0]) {
        ret = rate[0];
      } else if (e >= energy.back()) {
        ret = rate.back();
      } else {
        r = 1;
        while (e > energy[r]) {
          r++;
        }
        ret = ((energy[r] - e) * rate[r - 1] + (e - energy[r - 1]) * rate[r]) /
              (energy[r] - energy[r - 1]);
      }
    }
    return ret;
  }

  std::vector<double> energy, rate;
};

struct CS_3d {

  CS_3d(const std::string filename) : energy(), exit_energy(), cdf(), rvalue() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_energy.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss5;
      iss5 << line;
      while (getline(iss5, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      rvalue.push_back(tmp_vec);
    }
    file.close();
  }
  CS_3d(const CS_3d &other)
      : energy(other.energy), exit_energy(other.exit_energy), cdf(other.cdf),
        rvalue(other.rvalue) {}

  CS_3d() : energy(), exit_energy(), cdf(), rvalue() {}

  int find_energy_index(const double e) const {
    int ret = 0;
    while (ret < int(energy.size()) && e > energy[ret]) {
      ret++;
    }
    return ret;
  }
  std::vector<double> sample(const double e, gsl_rng *gen,
                             AtomConsts &coef) const {
    int enindex = find_energy_index(e);
    int densindex = 0;
    double densrng = gsl_rng_uniform(gen);
    if (enindex == 0) {
      while (cdf[enindex][densindex] < densrng) {
        densindex++;
      }
    } else {
      double densrng2 = gsl_rng_uniform(gen);
      double difcheck =
          (e - energy[enindex - 1]) / (energy[enindex] - energy[enindex - 1]);
      if (densrng2 < difcheck) {
        while (cdf[enindex][densindex] < densrng) {
          densindex++;
        }
      } else {
        enindex--;
        while (cdf[enindex][densindex] < densrng) {
          densindex++;
        }
      }
    }
    double outenergycm;
    double outrvalue;
    if (densindex == 0) {
      outenergycm = exit_energy[enindex][densindex];
      outrvalue = rvalue[enindex][densindex];
    } else {
      double difference =
          (densrng - cdf[enindex][densindex - 1]) /
          (cdf[enindex][densindex] - cdf[enindex][densindex - 1]);
      outenergycm = difference * exit_energy[enindex][densindex] +
                    (1 - difference) * exit_energy[enindex][densindex];
      outrvalue = difference * rvalue[enindex][densindex] +
                  (1 - difference) * rvalue[enindex][densindex];
    }

    double ea = e * (coef.ANU / (coef.ANU + coef.iNU));
    double eb = outenergycm * ((coef.BNU + coef.oNU) / (coef.BNU));
    double Ea = ea + coef.Si;
    double Eb = eb + coef.So;
    double X1 = fmin(Ea, 130) * Eb / Ea;
    double X3 = fmin(Ea, 41) * Eb / Ea;
    double aval = (0.04 * X1) + (1.8 * 1e-6 * pow(X1, 3)) +
                  (6.7 * 1e-7 * coef.Mi * coef.Mo * pow(X3, 4));
    double cdfc2 = outrvalue * cosh(aval) - sinh(aval);
    double cdfc1 = 2 * sinh(aval);
    double u2 = gsl_rng_uniform(gen);
    double z = cdfc1 * u2 + cdfc2;
    double z2 = (z + pow(pow(z, 2) - pow(outrvalue, 2) + 1, 1.0 / 2.0)) /
                (outrvalue + 1);
    double outanglecm = log(z2) / aval;
    outanglecm = fmax(outanglecm, -1);
    outanglecm = fmin(outanglecm, 1);
    double outenergylab =
        outenergycm + (e * coef.iNU * coef.oNU / pow(coef.ANU + coef.iNU, 2)) +
        ((2 * pow(coef.iNU * coef.oNU * outenergycm * e, 1.0 / 2.0) *
          outanglecm) /
         (coef.ANU + coef.iNU));
    double outanglelab =
        (pow(outenergycm / outenergylab, 1.0 / 2.0) * outanglecm) +
        (pow(coef.iNU * coef.oNU * e / outenergylab, 1.0 / 2.0) /
         (coef.ANU + coef.iNU));
    if (fabs(outanglelab) > 1) {
      outanglelab = fmin(1, outanglelab);
      outanglelab = fmax(-1, outanglelab);
    }
    std::vector<double> returnarray = {outenergylab, outanglelab};
    return returnarray;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_energy, cdf, rvalue;
};

struct CS_2d {

  CS_2d(const std::string filename) : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_angle.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
  }
  CS_2d(const CS_2d &other)
      : energy(other.energy), exit_angle(other.exit_angle), cdf(other.cdf) {}

  CS_2d() : energy(), exit_angle(), cdf() {}

  int find_energy_index(const double e) const {
    int ret = 0;
    while (ret < int(energy.size()) && e > energy[ret]) {
      ret++;
    }
    return ret;
  }
  double sample(const double e, gsl_rng *gen) const {
    int enindex = find_energy_index(e);
    int densindex = 0;
    double densrng = gsl_rng_uniform(gen);
    if (enindex == 0) {
      while (cdf[enindex][densindex] < densrng) {
        densindex++;
      }
    } else {
      double densrng2 = gsl_rng_uniform(gen);
      double difcheck =
          (e - energy[enindex - 1]) / (energy[enindex] - energy[enindex - 1]);
      if (densrng2 < difcheck) {
        while (cdf[enindex][densindex] < densrng) {
          densindex++;
        }
      } else {
        enindex--;
        while (cdf[enindex][densindex] < densrng) {
          densindex++;
        }
      }
    }
    double outanglecm;
    if (densindex == 0) {
      outanglecm = exit_angle[enindex][densindex];
    } else {
      double difference =
          (densrng - cdf[enindex][densindex - 1]) /
          (cdf[enindex][densindex] - cdf[enindex][densindex - 1]);
      outanglecm = difference * exit_angle[enindex][densindex] +
                   (1 - difference) * exit_angle[enindex][densindex - 1];
    }

    return outanglecm;
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_angle, cdf;
};

#endif
