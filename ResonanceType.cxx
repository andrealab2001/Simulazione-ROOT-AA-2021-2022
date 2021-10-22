#include "ResonanceType.h"

ResonanceType::ResonanceType(std::string name, double mass, int charge,
                             double width)
    : ParticleType(name, mass, charge), fWidth{width} {}

double ResonanceType::GetWidth() const { return fWidth; }

void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << "Resonance width: " << fWidth << '\n';
}