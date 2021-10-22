#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP
#include "ParticleType.h"

class ResonanceType : public ParticleType {
public:
  ResonanceType(std::string name, double mass, int charge, double width);
  double GetWidth() const override;
  void Print() const override;

private:
  double const fWidth;
};

#endif