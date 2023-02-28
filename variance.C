#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TROOT.h>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include <TAxis.h>
#include <TLine.h>
#include <TStyle.h>
#include <TChain.h>
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include <fairlogger/Logger.h>
#include "TRDBase/Geometry.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/HelperMethods.h"
#include "DataFormatsTRD/Constants.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#endif

using namespace std;
using namespace o2::trd;
using namespace o2::trd::constants;
using std::ofstream;

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
  uint32_t rowGlb{0};
  int col{0};
	int rob{0};
	int mcmcol{0};
  int channel{0};
  uint32_t channelGlb{0}; // global channel number (0..167)
  uint32_t colGlb{0};
  //std::vector<int> adcvalues(1378944);    //adc values
  float adcMean{0}; // mean ADC value for this pad
	int noisydigits{0};
  int timebins{0};
  //uint64_t adc{0};
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

std::vector<TH2F*> createNoisyChannelHistsPerLayer()
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


void variance()         //(bool skipSharedDigits = true)
{

// prepare to read the TRD digits
TChain chain("o2sim");
chain.AddFile("/home/archita/Downloads/trddigits (1).root");
std::vector<Digit> digits, *digitInPtr{&digits};
chain.SetBranchAddress("TRDDigit", &digitInPtr);
std::vector<TriggerRecord> *trigRecs = nullptr;
chain.SetBranchAddress("TriggerRecord", &trigRecs);
std::cout << "Total number of entries in the tree: " << chain.GetEntries() << std::endl;

//o2::trd::Geometry *mGeo;
int nRows = 16*4+ 12;
// TH1 *plots0[18][6];				//for meanADC histogram
// TH2 *plots1[18][6];		  	//for mean ADC histogram wrt channels
// TH2 *plots2[18][6]; 	 		//for mean ADC histogram wrt padrows
//TH2 *plots3[18][6]; 			//for SD (noise map)
// TH1 *plots4[18][6];	 			//for RMS (or SD)ADC histogram

// auto hnoisyDet = new TH2F("hnoisyDet", ";detector;adcRMS", 540, -0.5, 539.5, 1700, -1, 180);


int NEvents =0;
int maxPadRowInSector = 4 * constants::NROWC1 + constants::NROWC0;  // 16*4 + 12*1 = 76 padrows
int NChannels;

int nChannelsPerRow = NMCMROBINCOL * 2 * NADCMCM; // 4*2*21 = 168
int nChannelsPerChamberC1 = NROWC1 * nChannelsPerRow; // 16*168
int nChannelsPerChamberC0 = NROWC0 * nChannelsPerRow; //  12*168
int nChannelsTotal = NSECTOR * NLAYER * (NSTACK - 1) * nChannelsPerChamberC1; // for stacks 0, 1, 3 and 4
nChannelsTotal += NSECTOR * NLAYER * nChannelsPerChamberC0; // for stack 2
int nChannelsPerSector = nChannelsTotal / NSECTOR;
int nChannelsPerLayer = nChannelsPerSector / NLAYER;

//std::vector<int> adcvalues(1378944);

// create a vector with one ChannelInfo object per TRD readout channel
printf("Creating channel information vector of size %i\n", nChannelsTotal);
std::vector<ChannelInfo> channelInfos(nChannelsTotal);
ChannelInfo channelInfo;
auto fOut = new TFile("trdNoise.root", "recreate");
auto hDet190 = new TH2F("det_190", "Detector 190: Noise Map;Pad Row;Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);
auto hDet204 = new TH2F("det_204", "Detector 204: Noise Map;Pad Row;Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);
auto hDet52 = new TH2F("det_52", "Detector 52: Noise Map;Pad Row;Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);
auto hDet189 = new TH2F("det_189", "Detector 189: Noise Map;Pad Row;Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);
auto hbaseline = new TH1F("adcMean", "Baseline plot for adcRMS<8; adcMean", 1500,5,20);
auto hpulseheight = new TH1F("pulse-height", "pulse-height plot for Det 190", 30, -0.5, 29.5);
auto hLayers = createNoisyChannelHistsPerLayer();

auto tree = new TTree("padValues", "Mean and RMS ADC information per pad");
tree->Branch("padInfo", &channelInfo);

printf("Looping over all provided digits and filling channel information\n");
for (int iEntry = 0; iEntry < chain.GetEntries(); ++iEntry)   //tf
{
chain.GetEntry(iEntry);

	for(const auto &trigRec : *trigRecs)
	{
		int nDigits = trigRec.getNumberOfDigits();   		//Checks the no. of digits for every trigger

		for(int iDigit = 0; iDigit <nDigits; ++iDigit)
		{
			Digit digit = (digits)[iDigit + trigRec.getFirstDigit()];
      // if (skipSharedDigits && digit.isSharedDigit()) {
      //   continue;
      // }
			int iadcSum		= digit.getADCsum();		//ADC sum
			int padrow		= digit.getPadRow();			// 0-15
			int rob 			= digit.getROB();					// read out board within chamber [0-7] [0-5] depending on C0 or C1
			int mcm 			= digit.getMCM();					// 0-7
			//int padcol		= digit.getPadCol();
			int channel		= digit.getChannel();  		// 0-20

      int mcmcol 		= (rob % 2) ? mcm % 4 + NMCMROBINCOL : mcm % 4;
			int channelGlb = mcmcol * NADCMCM + 20 - channel;
			int detector	= digit.getDetector();		// 0-539
			int sector		= detector / 30;
			int layer			= detector % 6;
			int stack			= (detector % 30) / 6;

			ArrayADC adcArray = digit.getADC();		//ADC values

			int stackOffset = (stack < 3) ? stack * nChannelsPerChamberC1 : (stack - 1) * nChannelsPerChamberC1 + nChannelsPerChamberC0;

			int index = sector * nChannelsPerSector + layer * nChannelsPerLayer + stackOffset+ nChannelsPerRow * padrow + channelGlb;
      int rowGlb = (stack < 3) ? digit.getPadRow() + stack * 16 : digit.getPadRow() + 44 + (stack - 3) * 16; // pad row within whole sector
      int colGlb = digit.getPadCol() + sector * 144; // pad column number from 0 to NSectors * 144

  		auto& info = channelInfos[index];

      if (!info.init) {
        // the first time we see data from this channel we fill the detector information etc.
        info.det = detector;
        info.sec = sector;
        info.stack = stack;
        info.layer = layer;
        info.row = padrow;
        info.rowGlb = rowGlb;
				info.rob = rob;
				info.mcmcol = mcmcol;
        info.col = digit.getPadCol();
        info.colGlb = colGlb;
        info.channel = channel;
        info.channelGlb = channelGlb;
        //info.adcvalues[TIMEBINS]= adcArray;
			  info.init = 1;
			}

      for (int i = 0; i < TIMEBINS; ++i)
      {
       auto adc = digit.getADC()[i];
       info.adcSum += adc;
		   info.adcSumSquared += (adc*adc);
		// mean, variance and standard deviation are calculated recursively

      float meanCurrent = info.adcMean;
      info.adcMean += (adc - meanCurrent) / (info.nEntries + 1);
      info.variance += (info.nEntries * (info.nEntries + 1) * (info.adcMean - meanCurrent) * (info.adcMean - meanCurrent));
      info.nEntries += 1;
      info.adcRMS = TMath::Sqrt(info.variance/(info.nEntries-1));


    }
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
//   if(info.adcRMS >=8){
//   hLayers[info.layer]->Fill(info.rowGlb, info.colGlb, info.adcRMS);
// } else continue;
 //  if (info.det==204){
 //    hDet204->Fill(info.rowGlb,info.channelGlb, info.adcRMS);
 // } else continue;

  // if (info.adcRMS<8){
  //   hbaseline->Fill(info.adcMean);
  // } else continue;

  // if (info.det==190)
  //
  // hDet190->Fill(info.det,info.rowGlb, info.adcRMS);
  //
  // if (info.det==204.6)
  // hDet204_6->Fill(info.det,info.rowGlb, info.adcRMS);


//  hnoisyDet->Fill(info.det, info.adcRMS);

  // if(info.det==205 && info.adcRMS>=6 && info.adcRMS<=7.5) {
  //   //cout << "Sector: " << info.sec << endl;
  //   hnoisyDet[info.layer]->Fill(info.row, info.channelGlb,info.adcRMS);
  // }

/*  if (info.det==190) {   //&& (info.row==11) && (info.col==125) && (info.channel==2)){

for(int i=0; i<TIMEBINS; i++){
      //for(const auto &adc: adcArray)
      //tree->Scan("det:row:col", "init>0 && adcRMS>6 && det==190");
      //cout << info.adc << endl;
      hpulseheight->Fill(i, info.adcvalues[i]); */ //This piece isn't working here but in pulseheightplot.C
    }
//  }
//}
//}
printf("Got information for %i out of %i channels\n", countFilledChannels, nChannelsTotal);

// auto c1 = new TCanvas("Det190_pulseheightplot", "Det190_pulseheightplot");
// hpulseheight->Draw("colz");
// auto c1 = new TCanvas("c1", "c1");
// hDet204->Draw("colz");

//auto c2 = new TCanvas("c2", "c2");
//hbaseline->Draw("colz");
// auto c2 = new TCanvas("c2", "c2");
// hDet52->Draw("colz");
//
// auto c3 = new TCanvas("c3", "c3");
 //hDet52->Draw("colz");
//hDet->GetZaxis()->SetRangeUser(0, 10);
//c1->SaveAs("529829_digits-det.png");
// hnoisyDet->Draw();
// hnoisyDet->SetStats(kTRUE);
// c1->SaveAs("hnoisyDet_adcRMS.png");
/*
auto c3 = new TCanvas("c3", "c3", 1400, 1000);
auto line = new TLine();
line->SetLineStyle(kDashed);
c3->Divide(3,2);
for (int iLayer = 0; iLayer < 6; ++iLayer) {
  auto pad = c3->cd(iLayer + 1);
  pad->SetRightMargin(0.15);
  hLayers[iLayer]->Draw("colz");
  //hLayers[iLayer]->SetStats(kTRUE);
  drawTrdGrid(line);
  pad->SetLogz();
}
//}
c3->Update();
c3->SaveAs("529969_NoisyPads-location_adcRMS>8.png");
*/
// c1->Update();
// hpulseheight->Write();
tree->Write();
// delete tree;
fOut->Close();


}

// TCanvas *c1 = new TCanvas("c1", "TRDPadNoise", 800, 600);
// hnoisyADC->Draw("colz");
// hnoisyADC->SetStats(kTRUE);
// //std::string plotTitle = Form("Sector S%d: L%d",iSector, iLayer);
// //hnoisyADC->SetTitle(plotTitle.data());
// hnoisyADC->GetXaxis()->CenterTitle();
// hnoisyADC->GetYaxis()->CenterTitle();
// hnoisyADC->GetXaxis()->SetLabelSize(0.05);
// hnoisyADC->GetYaxis()->SetLabelSize(0.05);
// hnoisyADC->GetYaxis()->SetTitleSize(0.05);
// hnoisyADC->GetXaxis()->SetTitleSize(0.05);
// hnoisyADC->GetZaxis()->SetRangeUser(0, 3);
// hnoisyADC->GetZaxis()->SetLabelSize(0.05);
