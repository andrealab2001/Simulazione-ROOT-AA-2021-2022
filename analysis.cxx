void setStyle2() {
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptStat(112210);
}

void analysis() {

  gStyle->SetOptFit(111);

  TFile *gen = new TFile("generation.root");

  TH1F *hTypes = (TH1F *)gen->Get("hParticleType");

  cout << "Proportions of generated particles for each bin:" << endl;
  cout << "================================================" << endl;

  double total = 0.;
  for (int i = 1; i <= hTypes->GetNbinsX(); ++i) {
    double proportion = (hTypes->GetBinContent(i) / hTypes->GetEntries()) * 100;
    cout << "Bin " << i << ": " << proportion << " % " << endl;
    total += proportion;
  }
  cout << "================================================" << endl;
  cout << "Total: " << total << " % " << endl;

  TH1F *hPolAng = (TH1F *)gen->Get("hPolarAngle");
  TH1F *hAzimAng = (TH1F *)gen->Get("hAzimuthalAngle");
  TH1F *hMom = (TH1F *)gen->Get("hMomentum");

  TCanvas *angles1 = new TCanvas("Pol. angle", "Pol. angle fit");
  hPolAng->Fit("pol0");
  hPolAng->GetFunction("pol0")->SetLineColor(kRed);
  hPolAng->Draw();
  TLegend *legPol = new TLegend(.1, .7, .3, .9, "Legend");
  legPol->AddEntry(hPolAng, "Punti sperimentali");
  legPol->AddEntry(hPolAng->GetFunction("pol0"), "Funzione di fit");
  legPol->Draw("SAME");

  TCanvas *angles2 = new TCanvas("Az. angle", "Az. angle fit");
  hAzimAng->Fit("pol0");
  hAzimAng->GetFunction("pol0")->SetLineColor(kRed);
  hAzimAng->Draw();
  TLegend *legAz = new TLegend(.1, .7, .3, .9, "Legend");
  legAz->AddEntry(hAzimAng, "Punti sperimentali");
  legAz->AddEntry(hAzimAng->GetFunction("pol0"), "Funzione di fit");
  legAz->Draw("SAME");

  TCanvas *momentum = new TCanvas("Momentum", "Mom. intensity fit");
  hMom->Fit("expo");
  hMom->GetFunction("expo")->SetLineColor(kRed);
  hMom->Draw();
  TLegend *legMom = new TLegend(.1, .7, .3, .9, "Legend");
  legMom->AddEntry(hMom, "Punti sperimentali");
  legMom->AddEntry(hMom->GetFunction("expo"), "Funzione di fit");
  legMom->Draw("SAME");

  TH1F *hPiKConc = (TH1F *)gen->Get("hInvMass_Pi_KConc");
  TH1F *hPiKDisc = (TH1F *)gen->Get("hInvMass_Pi_KDisc");
  TH1F *hSubtract =
      new TH1F("hSubtract", "Discordant minus Concordant charge #pi and K",
               1000, 0., 7.);
  hSubtract->Add(hPiKDisc, hPiKConc, 1., -1.);
  hSubtract->GetXaxis()->SetRange(30, 300);
  hSubtract->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
  hSubtract->GetYaxis()->SetTitleOffset(1.2);
  hSubtract->GetYaxis()->SetTitle("Occurencies");

  TCanvas *subtraction = new TCanvas(
      "subtraction", "Discordant minus Concordant charge #pi and K");
  hSubtract->Fit("gaus", "", "", 0.6, 1.4);
  hSubtract->GetFunction("gaus")->SetLineColor(kRed);
  hSubtract->Draw();
  TLegend *legSub = new TLegend(.1, .7, .3, .9, "Legend");
  legSub->AddEntry(hSubtract, "Punti sperimentali");
  legSub->AddEntry(hSubtract->GetFunction("gaus"), "Funzione di fit");
  legSub->Draw("SAME");

  TH1F *hSubtract2 = new TH1F(
      "hSubtract2", "Discordant minus Concordant charge (total)", 1000, 0., 7.);
  TH1F *hConc = (TH1F *)gen->Get("hInvMassConc");
  TH1F *hDisc = (TH1F *)gen->Get("hInvMassDisc");
  hSubtract2->Add(hDisc, hConc, 1., -1.);
  hSubtract2->GetXaxis()->SetRange(30, 300);
  hSubtract2->GetYaxis()->SetTitleOffset(1.2);
  hSubtract2->GetYaxis()->SetTitle("Occurrencies");
  hSubtract2->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
  TCanvas *subtraction2 =
      new TCanvas("subtraction2", "Discordant minus Concordant charge (total)");

  hSubtract2->Fit("gaus", "", "", 0.6, 1.4);
  hSubtract2->GetFunction("gaus")->SetLineColor(kRed);
  hSubtract2->Draw();
  TLegend *legSub2 = new TLegend(.1, .7, .3, .9, "Legend");
  legSub2->AddEntry(hSubtract2, "Punti sperimentali");
  legSub2->AddEntry(hSubtract->GetFunction("gaus"), "Funzione di fit");
  legSub2->Draw("SAME");

  subtraction->Print("canvas/subtractionPartial.pdf");
  subtraction->Print("canvas/subtractionPartial.C");
  subtraction->Print("canvas/subtractionPartial.root");

  subtraction2->Print("canvas/subtractionTotal.pdf");
  subtraction2->Print("canvas/subtractionTotal.C");
  subtraction2->Print("canvas/subtractionTotal.root");
}