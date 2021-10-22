#include "ParticleType.h"

ParticleType::ParticleType(std::string name, double mass, int charge)
    : fName{name}, fMass{mass}, fCharge{charge} {}

std::string ParticleType::GetName() const { return fName; }

double ParticleType::GetMass() const { return fMass; }

int ParticleType::GetCharge() const { return fCharge; }

double ParticleType::GetWidth() const { return 0; }

void ParticleType::Print() const {
  std::cout << "Particle's parameters:" << '\n';
  std::cout << "Name: " << fName << '\n';
  std::cout << "Mass: " << fMass << '\n';
  std::cout << "Charge: " << fCharge << "\n\n";
}
