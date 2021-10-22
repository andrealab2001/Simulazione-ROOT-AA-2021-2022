#include "Particle.h"
#include "ParticleType.h"
#include "ResonanceType.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"

R__LOAD_LIBRARY(ParticleType_cxx.so)
R__LOAD_LIBRARY(ResonanceType_cxx.so)
R__LOAD_LIBRARY(Particle_cxx.so)

void setStyle() {
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(1110);
}

int main() {

  setStyle();

  Particle::AddParticleType("Pi+", 0.13957, 1);
  Particle::AddParticleType("Pi-", 0.13957, -1);
  Particle::AddParticleType("K+", 0.49367, 1);
  Particle::AddParticleType("K-", 0.49367, -1);
  Particle::AddParticleType("P+", 0.93827, 1);
  Particle::AddParticleType("P-", 0.93827, -1);
  Particle::AddParticleType("K*", 0.89166, 0, 0.050);

  TH1F *hParticleType = new TH1F("hParticleType", "Particle Types", 7, 0, 7);
  TH1F *hAzimuthalAngle = new TH1F(
      "hAzimuthalAngle", "Azimuthal Angle distribution", 1000, 0., 2 * M_PI);
  TH1F *hPolarAngle =
      new TH1F("hPolarAngle", "Polar Angle distribution", 1000, 0., M_PI);
  TH1F *hMomentum =
      new TH1F("hMomentum", "Momentum distribution (module)", 1000, 0., 4.);
  TH1F *hMomTrasv =
      new TH1F("hMomTrasv", "Trasverse momentum (module)", 1000, 0., 1.4);
  TH1F *hEnergy = new TH1F("hEnergy", "Particles energy", 1000, 0., 4.);
  TH1F *hInvMassTot =
      new TH1F("hInvMass", "Invariant Mass distribution", 1000, 0., 7.);
  TH1F *hInvMassDisc =
      new TH1F("hInvMassDisc",
               "Invariant Mass distribution (discordant charge)", 1000, 0., 7.);
  TH1F *hInvMassConc =
      new TH1F("hInvMassConc",
               "Invariant Mass distribution (concordant charge)", 1000, 0., 7.);
  TH1F *hInvMass_PiNeg_KPos = new TH1F(
      "hInvMass_PiNeg_KPos",
      "Invariant Mass distribution for #pi^{-} and K^{+}", 1000, 0., 7.);
  TH1F *hInvMass_PiPos_KNeg = new TH1F(
      "hInvMass_PiPos_KNeg",
      "Invariant Mass distribution for #pi^{+} and K^{-}", 1000, 0., 7.);
  TH1F *hInvMass_PiNeg_KNeg = new TH1F(
      "hInvMass_PiNeg_KNeg",
      "Invariant Mass distribution for #pi^{-} and K^{-}", 1000, 0., 7.);
  TH1F *hInvMass_PiPos_KPos = new TH1F(
      "hInvMass_PiPos_KPos",
      "Invariant Mass distribution for #pi^{+} and K^{+}", 1000, 0., 7.);
  TH1F *hInvMassDecay = new TH1F(
      "hInvMassDecay", "Invariant Mass distribution for K* decay products",
      1000, 0., 2.);
  int const x = 20;
  Particle particles[100 + x];

  for (int i = 0; i < 1e5; ++i) {
    int k = 0;
    for (int j = 0; j < 100; ++j) {
      double phi = gRandom->Uniform(0., 2 * M_PI);
      double theta = gRandom->Uniform(0., M_PI);
      hAzimuthalAngle->Fill(phi);
      hPolarAngle->Fill(theta);
      double p = gRandom->Exp(1.);
      hMomentum->Fill(p);
      double pX = p * sin(theta) * cos(phi);
      double pY = p * sin(theta) * sin(phi);
      double pZ = p * cos(theta);
      hMomTrasv->Fill(sqrt(pX * pX + pY * pY));
      particles[j].SetP(pX, pY, pZ);
      double prob = gRandom->Uniform(0., 1.);
      double type = gRandom->Uniform(0., 1.);

      if (prob < 0.8) {
        if (type < 0.5)
          particles[j].SetIndex("Pi+");
        else
          particles[j].SetIndex("Pi-");
      } else if (prob < 0.9) {
        if (type < 0.5)
          particles[j].SetIndex("K+");
        else
          particles[j].SetIndex("K-");
      } else if (prob < 0.99) {
        if (type < 0.5)
          particles[j].SetIndex("P+");
        else
          particles[j].SetIndex("P-");
      } else {
        particles[j].SetIndex("K*");
      }

      hParticleType->Fill(particles[j].GetIndex());
      hEnergy->Fill(particles[j].GetEnergy());

      if (particles[j].GetIndex() == 6) {
        double prob2 = gRandom->Uniform(0, 1);
        for (; k < x; ++k) {
          if (particles[100 + k].GetIndex() == Particle::GetFailureValue())
            break;
        }

        if (k == 20)
          break;

        if (prob2 < 0.5) {
          particles[100 + k].SetIndex("Pi+");
          particles[100 + k + 1].SetIndex("K-");
        } else {
          particles[100 + k].SetIndex("Pi-");
          particles[100 + k + 1].SetIndex("K+");
        }
        particles[j].Decay2body(particles[100 + k], particles[100 + k + 1]);
        hInvMassDecay->Fill(particles[100 + k].InvMass(particles[100 + k + 1]));
      }
    }
    // ciclo su tutte le particelle

    for (int r = 0; r < 100 + k + 1; ++r) {
      for (int s = r + 1; s < 100 + k + 1; ++s) {

        hInvMassTot->Fill(particles[r].InvMass(particles[s]));

        if (particles[r].GetCharge() * particles[s].GetCharge() < 0)
          hInvMassDisc->Fill(particles[r].InvMass(particles[s]));

        if (particles[r].GetCharge() * particles[s].GetCharge() > 0)
          hInvMassConc->Fill(particles[r].InvMass(particles[s]));

        if ((particles[r].GetIndex() == 0 && particles[s].GetIndex() == 3) ||
            (particles[r].GetIndex() == 3 && particles[s].GetIndex() == 0))
          hInvMass_PiPos_KNeg->Fill(particles[r].InvMass(particles[s]));

        if ((particles[r].GetIndex() == 1 && particles[s].GetIndex() == 2) ||
            (particles[r].GetIndex() == 2 && particles[s].GetIndex() == 1))
          hInvMass_PiNeg_KPos->Fill(particles[r].InvMass(particles[s]));

        if ((particles[r].GetIndex() == 1 && particles[s].GetIndex() == 3) ||
            (particles[r].GetIndex() == 3 && particles[s].GetIndex() == 1))
          hInvMass_PiNeg_KNeg->Fill(particles[r].InvMass(particles[s]));

        if ((particles[r].GetIndex() == 0 && particles[s].GetIndex() == 2) ||
            (particles[r].GetIndex() == 2 && particles[s].GetIndex() == 0))
          hInvMass_PiPos_KPos->Fill(particles[r].InvMass(particles[s]));
      }
    }

    // devo ripulire la "coda" dell'array, altrimenti all'iterazione successiva
    // il controllo troverà occupato

    for (int r = 0; r <= k + 1; ++r) {
      particles[100 + r].SetIndex(Particle::GetFailureValue());
    }
  }

  // cosmetica e scrittura su file

  TH1F *histos[14] = {
      hParticleType,       hAzimuthalAngle,     hPolarAngle,
      hMomentum,           hMomTrasv,           hEnergy,
      hInvMassTot,         hInvMassConc,        hInvMassDisc,
      hInvMass_PiNeg_KNeg, hInvMass_PiNeg_KPos, hInvMass_PiPos_KNeg,
      hInvMass_PiPos_KPos, hInvMassDecay};

  const TString titlesX[7] = {"Particle Types",
                              "Azimuthal Angle (rad)",
                              "Polar Angle (rad)",
                              "Momentum module (GeV/c)",
                              "Trasverse Momentum module (GeV/c)",
                              "Energy (GeV)",
                              "Invariant Mass (GeV/c^{2})"};

  const TString titleY = "Occurrencies";

  const TString types[7] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}",
                            "P^{+}",   "P^{-}",   "K*"};

  TFile *file = new TFile("generation.root", "RECREATE");

  for (int i = 1; i < 8; ++i) {
    hParticleType->GetXaxis()->SetBinLabel(i, types[i - 1]);
  }

  for (int i = 0; i < 14; ++i) {
    histos[i]->GetYaxis()->SetTitleOffset(1.3);
    histos[i]->GetYaxis()->SetTitle(titleY);
    if (i < 6)
      histos[i]->GetXaxis()->SetTitle(titlesX[i]);
    else
      histos[i]->GetXaxis()->SetTitle(titlesX[6]);
    histos[i]->SetLineColor(9);
    histos[i]->Write();
  }

  file->Close();
}
