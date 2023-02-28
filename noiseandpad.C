#include <vector>
#include <TROOT.h>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include <TLine.h>
#include <TStyle.h>
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include <fairlogger/Logger.h>
#include "TRDBase/Geometry.h"
#include "DataFormatsTRD/Tracklet64.h"
//#include "DataFormatsTRD/HelperMethods.h"

using namespace o2::trd;
using namespace o2::trd::constants;
// void drawTrdGrid(TLine* line)
// {
//   line->DrawLine(15.5, 0, 15.5, 3024);
//   line->DrawLine(31.5, 0, 31.5, 3024);
//   line->DrawLine(43.5, 0, 43.5, 3024);
//   line->DrawLine(59.5, 0, 59.5, 3024);
//   for (int iSec = 1; iSec < 18; ++iSec) {
//     float yPos = iSec * 168 - 0.5;
//     line->DrawLine(0, yPos, 72, yPos);
//   }
// }
//
// std::vector<TH2F*> createTrdPadHistsPerLayer()
// {
//   std::vector<TH2F*> hLayers;
//   for (int iLayer = 0; iLayer < 6; ++iLayer) {
//     hLayers.push_back(new TH2F(Form("layer%i", iLayer), Form("Noisy Digit count per pad in layer %i;stack;sector", iLayer), 76, -0.5, 75.5, 3024, -0.5, 3023.5));
//     auto xax = hLayers.back()->GetXaxis();
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
//     auto yax = hLayers.back()->GetYaxis();
//     for (int iSec = 0; iSec < 18; ++iSec) {
//       auto lbl = std::to_string(iSec);
//       yax->SetBinLabel(iSec * 168 + 72, lbl.c_str());
//     }
//     yax->SetTicks("-");
//     yax->SetTickSize(0.01);
//     yax->SetLabelSize(0.045);
//     yax->SetLabelOffset(0.01);
//     yax->SetTitleOffset(1.4);
//     hLayers.back()->SetStats(0);
//   }
//   return hLayers;
// }

void noiseandpad(string digitFile = "/home/archita/Downloads/trddigits (2).root")
{
	TFile *fin = TFile::Open(digitFile.data());	//Load digits
	TTree *digitTree = (TTree *)fin->Get("o2sim");

	std::vector<Digit> *digits = nullptr;
	digitTree->SetBranchAddress("TRDDigit", &digits);

  std::vector<TriggerRecord> *trigRecs = nullptr;
  digitTree->SetBranchAddress("TriggerRecord", &trigRecs);

  // get total number of time frames in the data
  int nev = digitTree->GetEntries();
  LOG(info) << nev << " TF entries found";


  o2::trd::Geometry *mGeo;
  int nRows = 16 * 4 + 12;
  TH1 *plots0[18][6];		//for meanADC histogram
	TH2 *plots3[18][6]; 		//for SD (noise map)
	TH2 *plots4[18][6];		//for RMS (or SD)ADC histogram
	int NEvents =0;

	std::vector<int> adcSum(1378944);
	std::vector<int> adcSqSum(1378944);
	std::vector<int> nEntries(1378944);
	int maxPadRowInSector = 4 * constants::NROWC1 + constants::NROWC0;  // 16*4 + 12*1 = 76 padrows

	int nChannelsPerRow = NMCMROBINCOL * 2 * NADCMCM; // 4*2*21 = 168
	int nChannelsPerChamberC1 = NROWC1 * nChannelsPerRow; // 16*168
	int nChannelsPerChamberC0 = NROWC0 * nChannelsPerRow; //  12*168
	int nChannelsTotal = NSECTOR * NLAYER * (NSTACK - 1) * nChannelsPerChamberC1; // for stacks 0, 1, 3 and 4
	nChannelsTotal += NSECTOR * NLAYER * nChannelsPerChamberC0; // for stack 2
	int nChannelsPerSector = nChannelsTotal / NSECTOR;
	int nChannelsPerLayer = nChannelsPerSector / NLAYER;

	TH2 *hnoisyDet = new TH2F("hnoisyDet", ";Padrow; Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);


	for(int iSector =0; iSector <18; iSector++)
	{
		for (int iLayer = 0; iLayer < 6; iLayer++)
		{
			// std::string h0 = Form("meanADC_histogram_%d_%d", iSector, iLayer);
			// plots0[iSector][iLayer] = new TH1F(h0.data(), ";meanADC", 1500 , 5, 20);

			std::string h3 = Form("hadcsumsq_%d_%d", iSector, iLayer);
			plots3[iSector][iLayer] = new TH2F(h3.data(), ";Padrow; Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);

			// std::string h4 = Form("rmsADC_histogram_%d_%d", iSector, iLayer);
			// plots4[iSector][iLayer] = new TH2F(h4.data(), ";Padrow; Channels", 16*4+12, 0, 16*4+12, 168, 0, 168);


		}
	}

	for(int iev = 0; iev< nev; ++iev) //loop over each TF
	{
		int ndigits = digitTree->GetEntry(iev);

		for(const auto &trigRec: *trigRecs)
		{
			int nDigits = trigRec.getNumberOfDigits();   		//Checks the no. of digits for every trigger
//			cout << nDigits << " No. of digits per trigger " << endl;
//		  Printf("Line %d",__LINE__);
//			std::vector<int> nevents;
//			nevents.push_back(nev);
//			NEvents = NEvents + nevents.size();
//			std::cout << " No. of events" << NEvents << endl;

			for(int iDigit = 0; iDigit <nDigits; ++iDigit)
			{
				Digit digit = (*digits)[iDigit + trigRec.getFirstDigit()];

				int iadcSum		= digit.getADCsum();		//ADC sum
				int padrow		= digit.getPadRow();			// 0-15
				int rob 			= digit.getROB();					// read out board within chamber [0-7] [0-5] depending on C0 or C1
				int mcm 			= digit.getMCM();					// 0-7
				//int padcol		= digit.getPadCol();
				int channel		= digit.getChannel();  		// 0-20

	      int mcmcol 		= (digit.getROB() % 2) ? digit.getMCM() % 4 + NMCMROBINCOL : digit.getMCM() % 4;
				int channelGlb = mcmcol * NADCMCM + 20 - channel;
				int detector	= digit.getDetector();		// 0-539

				//if(digit.getDetector() == 190)
			//	{
				int sector		= HelperMethods::getSector(detector);
				int layer			= HelperMethods::getLayer(detector);
				int stack			= HelperMethods::getStack(detector);
				ArrayADC adcArray = digit.getADC();		//ADC values
				int stackOffset = (stack < 3) ? stack * nChannelsPerChamberC1 : (stack - 1) * nChannelsPerChamberC1 + nChannelsPerChamberC0;
				int index = sector * nChannelsPerSector + layer * nChannelsPerLayer + stackOffset+ nChannelsPerRow * padrow + channelGlb;
        int colGlb = digit.getPadCol() + sector * 144; // pad column number from 0 to NSectors * 144
        int rowGlb = (stack < 3) ? digit.getPadRow() + stack * 16 : digit.getPadRow() + 44 + (stack - 3) * 16; // pad row within whole sector


        //   cout << digit << endl;
					//cout << detector << endl;
				for (int i = 0; i < TIMEBINS; ++i)
				{
				for(const auto &adc: adcArray)
				{
					//cout<< adc << endl;
					//if(digit.getDetector() == 190) {
					adcSum[index] += adc;
					//cout << adcSum[index] << " ADC Sum " << endl;

					adcSqSum[index] += (adc*adc);
					nEntries[index]++;
				}
			 }
		 //} else continue;
		}
	}
}

	TCanvas *c1 = new TCanvas("c1", "TRDNoiseMap", 800, 600);
	c1->Divide(2,3);

	for(int iSector =0; iSector< 1; iSector++)
	{
		// for(int iStack = 0; iStack<5; iStack++)
		// { int stackOffset = (iStack < 3) ? iStack * nChannelsPerChamberC1 : (iStack - 1) * nChannelsPerChamberC1 + nChannelsPerChamberC0;
		for(int iLayer =0; iLayer<6; iLayer++)
		{
			for(int iPadrow=0; iPadrow<maxPadRowInSector; iPadrow++)
			{
				for(int iChannel=0; iChannel<nChannelsPerRow; iChannel++)
				{

					int index = iSector * nChannelsPerSector + iLayer * nChannelsPerLayer + nChannelsPerRow * iPadrow + iChannel;


								float adcsum 	= adcSum[index];
								//cout << adcsum << " ADC Sum " << endl;

								float adcsqsum =  adcSqSum[index];
								//cout << adcsqsum << " ADC Sq Sum " << endl;

								float n   	= nEntries[index];
								//cout << n << " No. of Entries " << endl;

								if (n <1)
								continue;

								float adcmean	= adcsum/n;				//Dividing by "n" takes care of the normalisation factor.
//								cout << adcmean << " ADC mean values " << endl;	//mean = baseline***

								float adcsqmean	= adcsqsum/n;
								//cout << adcsqmean << " ADC sq mean values " << endl;
//
								float var 	= (adcsqmean-(adcmean* adcmean));
								//cout << var << " ADC Variance values " << endl;

								float SD = TMath::Sqrt(var);				//rms= noise***
//								cout << SD << "Standard Deviation" << endl;
                //if(SD>=8){
                //  cout << SD << "Standard Deviation" << endl;
                // cout << digit << endl;
                 plots3[iSector][iLayer]->Fill(iPadrow, iChannel, SD);
              //}
						}
					}
				}
			}

	for(int iSector =0; iSector<1; iSector++)
	 {
		for (int iLayer = 0; iLayer < 6; iLayer++)
	    	{
	        	c1->cd(iLayer + 1);
						//c1->cd();
						//c1->cd();
						c1->SetLogy();
						plots3[iSector][iLayer]->Draw("colz");
	       	  plots3[iSector][iLayer]->SetStats(kTRUE);
		        plots3[iSector][iLayer]->SetTitle(Form(" Layers %i", iLayer));
	       	 	plots3[iSector][iLayer]->GetXaxis()->CenterTitle();
	        	plots3[iSector][iLayer]->GetYaxis()->CenterTitle();
	        	plots3[iSector][iLayer]->GetXaxis()->SetLabelSize(0.05);
	        	plots3[iSector][iLayer]->GetYaxis()->SetLabelSize(0.05);
	        	plots3[iSector][iLayer]->GetYaxis()->SetTitleSize(0.05);
	        	plots3[iSector][iLayer]->GetXaxis()->SetTitleSize(0.05);
						//plots3[iSector][iLayer]->GetZaxis()->SetRangeUser(0, 3);
	        	plots3[iSector][iLayer]->GetZaxis()->SetLabelSize(0.05);
	        }
	}




}
