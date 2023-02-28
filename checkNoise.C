#if !defined(__CLING__) || defined(__ROOTCLING__)
// ROOT header
#include <TROOT.h>
#include <TChain.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TLine.h>
// O2 header
#include "DataFormatsTRD/TriggerRecord.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/Constants.h"

#include <string>
#include <vector>
#endif

using namespace o2::trd;
using namespace o2::trd::constants;

struct ChannelInfo
{
  ChannelInfo() = default;
  ChannelInfo(const ChannelInfo&) = default;
  ChannelInfo& operator=(const ChannelInfo& rhs) = default;

  uint32_t det{0};
  uint32_t sec{0};
  uint32_t stack{0};
  uint32_t layer{0};
  uint32_t row{0};
  int col{0};
	int rob{0};
	int mcmcol{0};
  int channel{0};
  uint32_t channelGlb{0}; // global channel number (0..167)
  uint32_t colGlb{0};
  float adcMean{0}; // mean ADC value for this pad
	int noisydigits{0};
  uint64_t adcSum{0}; // sum of ADC_i values
  uint64_t adcSumSquared{0}; // sum of ADC_i^2
  float adcSquaredMean{0};
  float variance{0};
  uint32_t nEntries{0}; // number of ADC values stored
  float adcRMS{0}; // the sum of (ADC_i - ADC_mean)^2
	//float StdDeviation{0}; // rms= noise***
  int init{0}; // flag whether this channel sent some ADC values (1-yes, 0-no)
  ClassDefNV(ChannelInfo, 1);
};

void drawTrdGrid(TLine* line)
{
  line->DrawLine(15.5, 0, 15.5, 2592);
  line->DrawLine(31.5, 0, 31.5, 2592);
  line->DrawLine(43.5, 0, 43.5, 2592);
  line->DrawLine(59.5, 0, 59.5, 2592);
  for (int iSec = 1; iSec < 18; ++iSec) {
    float yPos = iSec * 144 - 0.5;
    line->DrawLine(0, yPos, 76, yPos);
  }
}

std::vector<TH2F*> createTrdPadHistsPerLayer()
{
  std::vector<TH2F*> hLayers;
  for (int iLayer = 0; iLayer < 6; ++iLayer) {
    hLayers.push_back(new TH2F(Form("layer%i", iLayer), Form("Noisy Channels count per pad in layer %i;stack;sector", iLayer), 76, -0.5, 75.5, 2592, -0.5, 2591.5));
    auto xax = hLayers.back()->GetXaxis();
    xax->SetBinLabel(8, "0");
    xax->SetBinLabel(24, "1");
    xax->SetBinLabel(38, "2");
    xax->SetBinLabel(52, "3");
    xax->SetBinLabel(68, "4");
    xax->SetTicks("-");
    xax->SetTickSize(0.01);
    xax->SetLabelSize(0.045);
    xax->SetLabelOffset(0.01);
    xax->SetTitleOffset(-1);
    auto yax = hLayers.back()->GetYaxis();
    for (int iSec = 0; iSec < 18; ++iSec) {
      auto lbl = std::to_string(iSec);
      yax->SetBinLabel(iSec * 144 + 72, lbl.c_str());
    }
    yax->SetTicks("-");
    yax->SetTickSize(0.01);
    yax->SetLabelSize(0.045);
    yax->SetLabelOffset(0.01);
    yax->SetTitleOffset(1.4);
    hLayers.back()->SetStats(0);
  }
  return hLayers;
}

// std::vector<TH2F*> createNoisyDetectorHistsPerLayer()
// {
//   std::vector<TH2F*> hnoisyDet;
//   for (int iLayer = 0; iLayer < 6; ++iLayer) {
//     hnoisyDet.push_back(new TH2F(Form("layer%i", iLayer), Form("Noisy Digit count per pad in layer %i;info.stack;info.sector", iLayer), 76, -0.5, 75.5, 290304, -0.5, 290303.5));
//     auto xax = hnoisyDet.back()->GetXaxis();
//     xax->SetBinLabel(8, "0");
//     xax->SetBinLabel(24, "1");
//     xax->SetBinLabel(38, "2");
//     xax->SetBinLabel(52, "3");
//     xax->SetBinLabel(68, "4");
//     xax->SetTicks("-");
//     xax->SetTickSize(0.01);
//     xax->SetLabelSize(0.045);
//     xax->SetLabelOffset(0.01);
//     xax->SetTitleOffset(-1);
//     auto yax = hnoisyDet.back()->GetYaxis();
//     for (int iSec = 0; iSec < 18; ++iSec) {
//       auto lbl = std::to_string(iSec);
//       yax->SetBinLabel(iSec * 290304 + 72, lbl.c_str());
//     }
//     yax->SetTicks("-");
//     yax->SetTickSize(0.01);
//     yax->SetLabelSize(0.045);
//     yax->SetLabelOffset(0.01);
//     yax->SetTitleOffset(1.4);
//     hnoisyDet.back()->SetStats(0);
//   }
//   return hnoisyDet;
// }

void checkNoise()
{
  TChain chain("o2sim");
  chain.AddFile("/home/archita/Downloads/trddigits (1).root");
  std::vector<Digit> digits, *digitInPtr{&digits};
  chain.SetBranchAddress("TRDDigit", &digitInPtr);
  std::cout << "Total number of entries in the tree: " << chain.GetEntries() << std::endl;

  std::vector<TriggerRecord> *trigRecs = nullptr;
  chain.SetBranchAddress("TriggerRecord", &trigRecs);
  std::cout << "Total number of entries in the tree: " << chain.GetEntries() << std::endl;

  auto fOut = new TFile("digitsOutput.root", "recreate");
//  auto fOut = new TFile("checkNoise.root", "recreate");
//  auto hDet = new TH1F("det", "Detector number for digit;detector;counts", 540, -0.5, 539.5);
//  auto hnoisyrms = new TH2F("hnoisyrms",";Padrow; Channels",16*4+12, 0, 16*4+12, 168, 0, 168);
  TH2 *h2padHadAdc = new TH2F("h2padHadAdc", ";Chamber;Pad", 540, 0, 539, 144 * 16, 0, 144 * 16 - 1);

  auto hLayers = createTrdPadHistsPerLayer();
//  auto hnoisyDet = createNoisyDetectorHistsPerLayer();

  //auto hLayer0 = new TH2F("layer0", "Digit count per pad in layer 0;rowGlb;colGlb", 76, -0.5, 75.5, 2592, -0.5, 2591.5);
  int nChannelsPerRow = NMCMROBINCOL * 2 * NADCMCM; // 4*2*21 = 168
  int nChannelsPerChamberC1 = NROWC1 * nChannelsPerRow; // 16*168
  int nChannelsPerChamberC0 = NROWC0 * nChannelsPerRow; //  12*168
  int nChannelsTotal = NSECTOR * NLAYER * (NSTACK - 1) * nChannelsPerChamberC1; // for stacks 0, 1, 3 and 4
  nChannelsTotal += NSECTOR * NLAYER * nChannelsPerChamberC0; // for stack 2
  int nChannelsPerSector = nChannelsTotal / NSECTOR;
  int nChannelsPerLayer = nChannelsPerSector / NLAYER;

  // create a vector with one ChannelInfo object per TRD readout channel
  printf("Creating channel information vector of size %i\n", nChannelsTotal);
  std::vector<ChannelInfo> channelInfos(nChannelsTotal);
  ChannelInfo channelInfo;


  auto tree = new TTree("padValues", "Mean and RMS ADC information per pad");
  tree->Branch("padInfo", &channelInfo);

  printf("Looping over all provided digits and filling channel information\n");
  for (int iEntry = 0; iEntry < chain.GetEntries(); ++iEntry) {
    chain.GetEntry(iEntry); // for each TimeFrame there is one tree entry

//    for (const auto& digit : digits) {
  for(const auto &trigRec : *trigRecs) {
    int nDigits = trigRec.getNumberOfDigits();   		//Checks the no. of digits for every trigger

    for(int iDigit = 0; iDigit <nDigits; ++iDigit) {
	  Digit digit = (digits)[iDigit + trigRec.getFirstDigit()];

      //   if (skipSharedDigits && digit.isSharedDigit()) {
      //   continue;
      // }
//      int det = digit.getdetector();
//      int layer = det % 6;
//      hdet->fill(det);
//      int stack = (det % 30) / 6;
//      int rowglb = stack < 3 ? digit.getpadrow() + stack * 16 : digit.getpadrow() + 44 + (stack - 3) * 16; // pad row within whole sector
//      int sec = det / 30;
//      int colglb = digit.getpadcol() + sec * 144; // pad column number from 0 to nsectors * 144
//      hlayers[layer]->fill(rowglb, colglb);

	int iadcSum = digit.getADCsum();		//ADC sum
	int padrow = digit.getPadRow();			// 0-15
	int rob  = digit.getROB();					// read out board within chamber [0-7] [0-5] depending on C0 or C1
	int mcm = digit.getMCM();					// 0-7
	int padcol		= digit.getPadCol();
	int channel = digit.getChannel();  		// 0-20
  int mcmcol = (rob % 2) ? mcm % 4 + NMCMROBINCOL : mcm % 4;
  int channelGlb = mcmcol * NADCMCM + 20 - channel;
	int detector = digit.getDetector();		// 0-539
	int sector = HelperMethods::getSector(detector);
	int layer = HelperMethods::getLayer(detector);
	int stack = HelperMethods::getStack(detector);
	ArrayADC adcArray = digit.getADC();		//ADC values
	int stackOffset = (stack < 3) ? stack * nChannelsPerChamberC1 : (stack - 1) * nChannelsPerChamberC1 + nChannelsPerChamberC0;
  int colGlb = digit.getPadCol() + sector * 144; // pad column number from 0 to NSectors * 144
  int rowGlb = (stack < 3) ? digit.getPadRow() + stack * 16 : digit.getPadRow() + 44 + (stack - 3) * 16; // pad row within whole sector

  int index = sector * nChannelsPerSector + layer * nChannelsPerLayer + stackOffset+ nChannelsPerRow * padrow + channelGlb;

  //hLayers[layer]->Fill(rowGlb, colGlb, channelGlb);
	auto& info = channelInfos[index];

        if (!info.init) {
        // the first time we see data from this channel we fill the detector information etc.
        info.det = detector;
        info.sec = sector;
        info.stack = stack;
        info.layer = layer;
        info.row = rowGlb;
	      info.rob = rob;
	      info.mcmcol = mcmcol;
        info.col = colGlb;
        info.channel = channel;
        info.channelGlb = channelGlb;
	      info.init = 1;
//        hDet->Fill(info.det);
	      }

	      for (int i = 0; i < TIMEBINS; ++i) {
      	  auto adc = digit.getADC()[i];
      	  info.adcSum += adc;
	        info.adcSumSquared += (adc*adc);
	        float meanCurrent = info.adcMean;
      	  info.adcMean += (adc - meanCurrent) / (info.nEntries + 1);
          info.variance += (info.nEntries * (info.nEntries + 1) * (info.adcMean - meanCurrent) * (info.adcMean - meanCurrent));

          info.nEntries += 1;
          info.adcRMS = TMath::Sqrt(info.variance/(info.nEntries-1));


        }
        //hLayers[info.layer]->Fill(info.row, info.channelGlb);
        //if(info.init>0&& info.nEntries>1&&info.adcRMS >=8) {
         hLayers[layer]->Fill(rowGlb, colGlb, channelGlb);
        //}
    }
  }
}

printf("Done reading the input from the digits\n");
int countFilledChannels = 0;
for (auto& info : channelInfos) {
	if (info.init) {
		++countFilledChannels;
	}
//    printf("%i_%i_%i: Row %i, Pad  %i, meanADC(%f), variance(%f), nEntries(%i)\n", info.sec, info.stack, info.layer, info.row, info.col, info.adcMean, info.variance, info.nEntries);

    channelInfo = info;
  	tree->Fill();

    // if(info.det==205 && info.adcRMS>=6 && info.adcRMS<=7.5) {
    //   //cout << "Sector: " << info.sec << endl;
    // hnoisyDet[info.layer]->Fill(info.row, info.channelGlb,info.adcRMS);
    // }
  //if(info.adcRMS >=8){
    // //hDet->Fill(info.det);
    // hLayers[info.layer]->Fill(info.row, info.channel);
  //}
}

//   auto c1 = new TCanvas("c1", "c1");
//   hDet->Draw();
// //  c1->SaveAs("digits-det.png");
//   c1->SaveAs("adcRMS-det.png");

  auto c2 = new TCanvas("c2", "c2", 1400, 1000);
  auto line = new TLine();
  line->SetLineStyle(kDashed);
  c2->Divide(3, 2);

  for (int iLayer = 0; iLayer < 6; ++iLayer) {
    auto pad = c2->cd(iLayer + 1);
    pad->SetRightMargin(0.15);
    hLayers[iLayer]->Draw("colz");
    hLayers[iLayer]->SetStats(kTRUE);
    drawTrdGrid(line);
    pad->SetLogz();
  }
  c2->Update();
//  c2->SaveAs("digits-pad.png");
  c2->SaveAs("529829_adcRMS-pad.png");

//   auto c3 = new TCanvas("c3", "c3", 800, 600);
//   auto line = new TLine();
//   line->SetLineStyle(kDashed);
//   c3->Divide(2,3);
//
// //for (int iSector =0; iSector <18; iSector++) {
//   for (int iLayer = 0; iLayer < 6; ++iLayer) {
//     auto pad = c3->cd(iLayer + 1);
//     pad->SetRightMargin(0.15);
//     hnoisyDet[iLayer]->Draw("colz");
//     hnoisyDet[iLayer]->SetStats(kTRUE);
//     drawTrdGrid(line);
//     pad->SetLogz();
//   }
// //}
//   c3->Update();
//   c3->SaveAs("NoiseMap_det.png");

  tree->Write();
  delete tree;
  //fOut->Close();


/*
  hDet->Write();
  delete hDet;
  for (auto& h : hLayers) {
    h->Write();
    delete h;
  }
  fOut->Close();
*/
}
