#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef GRID
#define GRID

struct Grid {

  Grid(const int n, const double dx_)
      : dx(dx_), x_pp(n), x_pm(n), x_mp(n), x_mm(n) {
    std::vector<double> tmp(1, 0);
    std::vector<std::vector<double>> tmp2(1, tmp);
    for (int i = 0; i < n; i++) {
      x_pp[i] = tmp2;
      x_pm[i] = tmp2;
      x_mp[i] = tmp2;
      x_mm[i] = tmp2;
    }
  }

  void print(const std::string filename) {
    std::ofstream outfile(filename);
    for (unsigned int ix = 0; ix < x_pp.size(); ix++) {
      for (unsigned int j = 0; j < x_pp[ix].size(); j++) {
        for (unsigned int k = 0; k < x_pp[ix][j].size(); k++) {
          if (x_pp[ix][j][k] > 0) {
            outfile << dx * ix << " " << dx * j << " " << dx * k << " "
                    << x_pp[ix][j][k] << std::endl;
          }
        }
      }
      for (unsigned int j = 0; j < x_pm[ix].size(); j++) {
        for (unsigned int k = 0; k < x_pm[ix][j].size(); k++) {
          if (x_pm[ix][j][k] > 0) {
            outfile << dx * ix << " " << dx * j << " " << -dx * (k + 1) << " "
                    << x_pm[ix][j][k] << std::endl;
          }
        }
      }
      for (unsigned int j = 0; j < x_mp[ix].size(); j++) {
        for (unsigned int k = 0; k < x_mp[ix][j].size(); k++) {
          if (x_mp[ix][j][k] > 0) {
            outfile << dx * ix << " " << -dx * (j + 1) << " " << dx * k << " "
                    << x_mp[ix][j][k] << std::endl;
          }
        }
      }
      for (unsigned int j = 0; j < x_mm[ix].size(); j++) {
        for (unsigned int k = 0; k < x_mm[ix][j].size(); k++) {
          if (x_mm[ix][j][k] > 0) {
            outfile << dx * ix << " " << -dx * (j + 1) << " " << -dx * (k + 1)
                    << " " << x_mm[ix][j][k] << std::endl;
          }
        }
      }
    }
    outfile.close();
    return;
  }

  void add(const std::vector<std::vector<double>> &y,
           const std::vector<double> &s, const int len) {
    int ix;
    unsigned int iy, iz;
    std::vector<double> tmp(1, 0);
    std::vector<std::vector<double>> tmp2(1, tmp);
    for (int i = 0; i < len; i++) {
      ix = floor(y[i][0] / dx);
      iy = floor(fabs(y[i][1]) / dx);
      iz = floor(fabs(y[i][2]) / dx);
      if (ix >= 0) {
        if (ix >= int(x_pp.size())) {
          x_pp.resize(ix + 1, tmp2);
          x_pm.resize(ix + 1, tmp2);
          x_mp.resize(ix + 1, tmp2);
          x_mm.resize(ix + 1, tmp2);
        }
        if (y[i][1] >= 0 && y[i][2] >= 0) {
          if (iy >= x_pp[ix].size()) {
            x_pp[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_pp[ix][iy].size()) {
            x_pp[ix][iy].resize(iz + 1, 0);
          }
          x_pp[ix][iy][iz] += s[i];
        } else if (y[i][1] >= 0 && y[i][2] < 0) {
          if (iy >= x_pm[ix].size()) {
            x_pm[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_pm[ix][iy].size()) {
            x_pm[ix][iy].resize(iz + 1, 0);
          }
          x_pm[ix][iy][iz] += s[i];
        } else if (y[i][1] < 0 && y[i][2] >= 0) {
          if (iy >= x_mp[ix].size()) {
            x_mp[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_mp[ix][iy].size()) {
            x_mp[ix][iy].resize(iz + 1, 0);
          }
          x_mp[ix][iy][iz] += s[i];
        } else {
          if (iy >= x_mm[ix].size()) {
            x_mm[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_mm[ix][iy].size()) {
            x_mm[ix][iy].resize(iz + 1, 0);
          }
          x_mm[ix][iy][iz] += s[i];
        }
      }
    }
    return;
  }
  double interp_t_calc(double &x, double v) {
    if (v>0){   
      if (((fabs(ceil(x/dx)*dx - x))/v)<1e-8){
        x+=1e-7;
        return (ceil(x/dx)*dx+dx - x)/v;
      } else {
        return (ceil(x/dx)*dx - x)/v;
      }
    } else if (v<0){
      if (((fabs(floor(x/dx)*dx - x))/v)<1e-8){
        x+= -1e-7;
        return (floor(x/dx)*dx - dx - x)/v;
      } else {
        return (floor(x/dx)*dx - x)/v;
      }
    } else {
      return 1;
    }
  }
  void add_interp(const std::vector<std::vector<double>> &y,
           const std::vector<double> &s, const int len) {
    int ix;
    unsigned int iy, iz;
    double interp_x, interp_y, interp_z, t_interp_x, t_interp_y, t_interp_z, t_min,
        tmp_start_x, tmp_start_y, tmp_start_z, tmp_end_x, tmp_end_y, tmp_end_z,tmp_t_left;
    for (int i = 1; i < len; i++) {
      tmp_t_left = 1.0;
      tmp_start_x = y[i-1][0];
      tmp_start_y = y[i-1][1];
      tmp_start_z = y[i-1][2];
      tmp_end_x = y[i][0];
      tmp_end_y = y[i][1];
      tmp_end_z = y[i][2];
      interp_x = tmp_end_x - tmp_start_x;
      interp_y = tmp_end_y - tmp_start_y;
      interp_z = tmp_end_z - tmp_start_z;
      while (tmp_t_left > 0) {
      t_interp_x = interp_t_calc(tmp_start_x,interp_x);
      t_interp_y = interp_t_calc(tmp_start_y,interp_y);
      t_interp_z = interp_t_calc(tmp_start_z,interp_z);
      t_min = fmin(fmin(t_interp_x,t_interp_y),t_interp_z);
      ix = floor(tmp_start_x / dx);
      iy = floor(fabs(tmp_start_y) / dx);
      iz = floor(fabs(tmp_start_z) / dx);
      add_data(ix,iy,iz,s[i]*fmin(t_min,tmp_t_left),tmp_start_y,tmp_start_z);
      tmp_t_left -= t_min;
      tmp_start_x += interp_x * t_min;
      tmp_start_y += interp_y * t_min;
      tmp_start_z += interp_z * t_min;
      }
    }
    return;
  }
  void add_data(int ix, unsigned int iy, unsigned int iz, double value, double y, double z) {
    std::vector<double> tmp(1, 0);
    std::vector<std::vector<double>> tmp2(1, tmp);
    if (ix >= 0) {
        if (ix >= int(x_pp.size())) {
          x_pp.resize(ix + 1, tmp2);
          x_pm.resize(ix + 1, tmp2);
          x_mp.resize(ix + 1, tmp2);
          x_mm.resize(ix + 1, tmp2);
        }
        if (y >= 0 && z >= 0) {
          if (iy >= x_pp[ix].size()) {
            x_pp[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_pp[ix][iy].size()) {
            x_pp[ix][iy].resize(iz + 1, 0);
          }
          x_pp[ix][iy][iz] += value;
        } else if (y >= 0 && z < 0) {
          if (iy >= x_pm[ix].size()) {
            x_pm[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_pm[ix][iy].size()) {
            x_pm[ix][iy].resize(iz + 1, 0);
          }
          x_pm[ix][iy][iz] += value;
        } else if (y < 0 && z >= 0) {
          if (iy >= x_mp[ix].size()) {
            x_mp[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_mp[ix][iy].size()) {
            x_mp[ix][iy].resize(iz + 1, 0);
          }
          x_mp[ix][iy][iz] += value;
        } else {
          if (iy >= x_mm[ix].size()) {
            x_mm[ix].resize(iy + 1, tmp);
          }
          if (iz >= x_mm[ix][iy].size()) {
            x_mm[ix][iy].resize(iz + 1, 0);
          }
          x_mm[ix][iy][iz] += value;
        }
      }
  }

  double dx;
  std::vector<std::vector<std::vector<double>>> x_pp, x_pm, x_mp, x_mm;
};

#endif
