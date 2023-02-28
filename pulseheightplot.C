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
#include <TProfile.h>
#include <TProfile2D.h>
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

void pulseheightplot()
{
TChain chain("o2sim");
chain.AddFile("/home/archita/Downloads/trddigits (1).root");    //Noise run: 529969
std::vector<Digit> digits, *digitInPtr{&digits};
chain.SetBranchAddress("TRDDigit", &digitInPtr);
std::vector<TriggerRecord> *trigRecs = nullptr;
chain.SetBranchAddress("TriggerRecord", &trigRecs);
std::cout << "Total number of entries in the tree: " << chain.GetEntries() << std::endl;

std::vector<int> adcvalues(1378944);

auto hpulseheight = new TH2F("pulse-height", "pulse-height plot", 30, 0, 30, 1024, 0, 1024);
//auto hpulseheight = new TH1F("pulse-height", "pulse-height plot for Det 190", 30, -0.5, 29.5);
//auto hpulseheight = new TH2F("pulse-height", "pulse-height plot for Det 190; Timebin; Channel", 30, -0.5, 29.5, 168, 0, 168);

auto *hprof_tbadc = new TProfile("hprof_tbadc","Profile of ADC counts versus Timebins ", 30, 0, 30, 0, 1024);
//auto *hprof_tbadc = new TProfile("hprof_tbadc","Profile of ADC counts versus Timebins ", 30, -0.5, 29.5,0, 1024);

int n_ph = 0;

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

			int iadcSum		= digit.getADCsum();		//ADC sum
			int padrow		= digit.getPadRow();			// 0-15
			int padcol		= digit.getPadCol();
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

		//GOOD CHANNELS
      //if((detector==95) && (padrow==0) && (padcol==2)) {


		//BAD CHANNELS
			//if((detector==190) && (padrow==11) && (padcol==125) && (channel==2))
      //if((detector==204) && (padrow==1) && (padcol==99) && (channel==10))
			//if((detector==246) && (padrow==0) && (padcol==81) && (channel==10))
			//if((detector==264) && (padrow==3) && (padcol==125) && (channel==20))
			//if((detector==447) && (padrow==0) && (padcol==15) && (channel==4))
			//if((detector==447) && (padrow==0) && (padcol==16) && (channel==3))
			if((detector==135) && (padrow==6) && (padcol==65) && (channel==8))
			{
        //n_ph++;
        for(int tb = 0; tb<TIMEBINS; tb++)
        {
          adcvalues[tb]= adcArray[tb];
          hpulseheight->Fill(tb, adcvalues[tb]);
					//hpulseheight->Fill(tb, channelGlb, adcvalues[tb]);
					hprof_tbadc->Fill(tb,adcvalues[tb]);
        }

			}
    }
  }

}
//auto c1 = new TCanvas("Det190_pulseheightplot", "Det190_pulseheightplot");
auto c1 = new TCanvas("pulseheightplot", "pulseheightplot");

//hpulseheight->Draw("colz");
hpulseheight->ProfileX();
hpulseheight->Draw("colz");
hprof_tbadc->Draw("same");
c1->Update();

}
