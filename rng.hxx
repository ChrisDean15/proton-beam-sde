#pragma once
#include <memory>
#include <random>

class RNG
{
public:
  RNG() = delete;

  static std::mt19937 &get_generator()
  {
    static std::mt19937 gen;
    return gen;
  }
  static void SetSeed(long unsigned int seed){
    get_generator().seed(seed);
  }

  void operator=(const RNG &) = delete;
  RNG(RNG &other) = delete;
};