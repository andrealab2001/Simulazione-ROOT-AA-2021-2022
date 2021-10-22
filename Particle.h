#ifndef PARTICLE_H
#define PARTICLE_H
#include "ParticleType.h"
#include "ResonanceType.h"
#include <cmath>
#include <cstdlib>

class Particle {
public:
  Particle(std::string name, double pX = 0., double pY = 0., double pZ = 0.);
  Particle();
  static int GetFailureValue();
  static void AddParticleType(std::string name, double mass, int charge,
                              double width = 0.);
  static void PrintParticleType();
  int GetIndex() const;
  double GetPx() const;
  double GetPy() const;
  double GetPz() const;
  double GetNorm2P() const;
  double GetMass() const;
  double GetCharge() const;
  double GetEnergy() const;
  double InvMass(Particle &p) const;
  void SetIndex(int indexIn);
  void SetIndex(std::string name);
  void SetP(double pX, double pY, double pZ);
  void Print() const;
  int Decay2body(Particle &dau1, Particle &dau2) const;

private:
  static const int fMaxNumParticleType = 10;
  static int fNParticleType;
  static int fFindFailure;
  static ParticleType *fParticleType[fMaxNumParticleType];
  int fIndex;
  double fPx, fPy, fPz;
  static int FindParticle(std::string particleName);
  void Boost(double bx, double by, double bz);
};

#endif