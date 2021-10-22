#include "Particle.h"

int Particle::fNParticleType{};

int Particle::fFindFailure{fMaxNumParticleType};

ParticleType *Particle::fParticleType[fMaxNumParticleType];

int Particle::FindParticle(std::string particleName) {
  for (int i = 0; i < fNParticleType; ++i) {
    if (fParticleType[i]->GetName() == particleName)
      return i;
  }
  std::cout << "No match found for " << particleName << '\n';
  std::cout << "Returned failure index " << fFindFailure << "\n\n";
  return fFindFailure;
}

Particle::Particle(std::string name, double pX, double pY, double Pz)
    : fPx{pX}, fPy{pY}, fPz{Pz} {
  fIndex = FindParticle(name);
}

Particle::Particle() : fIndex{fFindFailure}, fPx{0.}, fPy{0.}, fPz{0.} {}

int Particle::GetFailureValue() { return fFindFailure; }

void Particle::AddParticleType(std::string name, double mass, int charge,
                               double width) {
  int index = FindParticle(name);
  if (index != fFindFailure || fNParticleType == fMaxNumParticleType)
    return;
  ++fNParticleType;
  if (width == 0.) {
    fParticleType[fNParticleType - 1] = new ParticleType{name, mass, charge};
  } else {
    fParticleType[fNParticleType - 1] =
        new ResonanceType{name, mass, charge, width};
  }
}

void Particle::PrintParticleType() {
  for (int i = 0; i < fNParticleType; ++i) {
    fParticleType[i]->Print();
  }
}

int Particle::GetIndex() const { return fIndex; }

double Particle::GetPx() const { return fPx; }

double Particle::GetPy() const { return fPy; }

double Particle::GetPz() const { return fPz; }

double Particle::GetNorm2P() const { return fPx * fPx + fPy * fPy + fPz * fPz; }

double Particle::GetMass() const {
  return fIndex != fFindFailure ? fParticleType[fIndex]->GetMass() : 0.;
}

double Particle::GetCharge() const {
  return fIndex != fFindFailure ? fParticleType[fIndex]->GetCharge() : 100;
}

double Particle::GetEnergy() const {
  double energy = sqrt(GetMass() * GetMass() + GetNorm2P());
  return energy;
}

double Particle::InvMass(Particle &p) const {
  double pX = fPx + p.GetPx(), pY = fPy + p.GetPy(), pZ = fPz + p.GetPz();
  double invMass =
      sqrt((GetEnergy() + p.GetEnergy()) * (GetEnergy() + p.GetEnergy()) -
           (pX * pX + pY * pY + pZ * pZ));
  return invMass;
}

void Particle::SetIndex(int indexIn) {
  if (indexIn >= fNParticleType && indexIn != fFindFailure)
    return;
  else
    fIndex = indexIn;
}

void Particle::SetIndex(std::string name) {
  int index = FindParticle(name);
  if (index != fFindFailure)
    fIndex = index;
}

void Particle::SetP(double pX, double pY, double pZ) {
  fPx = pX;
  fPy = pY;
  fPz = pZ;
}

void Particle::Print() const {
  std::cout << "Index: " << fIndex << '\n';
  std::cout << "Momentum components: "
            << "( " << fPx << ',' << fPy << ',' << fPz << " )"
            << "\n\n";
}

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetMass() == 0.0) {
    std::cout << "Decayment cannot be preformed if mass is zero" << '\n';
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIndex < fFindFailure) { // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleType[fIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cout
        << "Decayment cannot be preformed because mass is too low in this "
           "channel"
        << '\n';
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
            pout * cos(theta));
  dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
            -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz) {

  double energy = GetEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}