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
			//if((detector==447) && (padrow==0) && (padcol==17) && (channel==2))
			//if((detector==447) && (padrow==0) && (padcol==18) && (channel==1))
			//if((detector==447) && (padrow==0) && (padcol==19) && (channel==0))
			//if((detector==135) && (padrow==6) && (padcol==65) && (channel==8))
			//if((detector==498) && (padrow==14) && (padcol==33) && (channel==4))
			//if((detector==498) && (padrow==14) && (padcol==34) && (channel==3))
			//if((detector==498) && (padrow==14) && (padcol==35) && (channel==2))
			//if((detector==498) && (padrow==14) && (padcol==36) && (channel==1))
			//if((detector==498) && (padrow==14) && (padcol==37) && (channel==0))
			//if((detector==27) && (padrow==11) && (padcol==71) && (channel==20))
			//if((detector==11) && (padrow==1) && (padcol==109) && (channel==0))
			//if((detector==23) && (padrow==13) && (padcol==73) && (channel==0))
			//if((detector==29) && (padrow==6) && (padcol==55) && (channel==0))
			//if((detector==54) && (padrow==11) && (padcol==37) && (channel==0))
			//if((detector==233) && (padrow==2) && (padcol==73) && (channel==0))
			//if((detector==365) && (padrow==14) && (padcol==73) && (channel==0))
			//if((detector==497) && (padrow==11) && (padcol==35) && (channel==20))
			//if((detector==21) && (padrow==11) && (padcol==70) && (channel==3))
			//if((detector==47) && (padrow==4) && (padcol==72) && (channel==1))
			//if((detector==47) && (padrow==5) && (padcol==72) && (channel==1))
			//if((detector==47) && (padrow==5) && (padcol==73) && (channel==0))
			//if((detector==135) && (padrow==6) && (padcol==57) && (channel==16))
			//if((detector==204) && (padrow==1) && (padcol==99) && (channel==10))
			//if((detector==222) && (padrow==3) && (padcol==61) && (channel==12))
			//if((detector==218) && (padrow==11) && (padcol==91) && (channel==0))
			//if((detector==213) && (padrow==11) && (padcol==91) && (channel==0))
			//if((detector==233) && (padrow==2) && (padcol==72) && (channel==1))
			//if((detector==233) && (padrow==5) && (padcol==71) && (channel==20))
			//if((detector==352) && (padrow==13) && (padcol==71) && (channel==20))
			//if((detector==59) && (padrow==12) && (padcol==36) && (channel==1))
			//if((detector==233) && (padrow==5) && (padcol==73) && (channel==0))
			//if((detector==233) && (padrow==8) && (padcol==71) && (channel==20))
			//if((detector==242) && (padrow==14) && (padcol==71) && (channel==20))
			//if((detector==441) && (padrow==11) && (padcol==18) && (channel==1))
			//if((detector==489) && (padrow==15) && (padcol==72) && (channel==1))
			//if((detector==20) && (padrow==14) && (padcol==72) && (channel==1))
			//if((detector==58) && (padrow==9) && (padcol==141) && (channel==4))
			//if((detector==135) && (padrow==6) && (padcol==55) && (channel==18))
			//if((detector==192) && (padrow==3) && (padcol==0) && (channel==19))
			//if((detector==200) && (padrow==4) && (padcol==0) && (channel==19))
			//if((detector==195) && (padrow==3) && (padcol==0) && (channel==19))
			//if((detector==195) && (padrow==8) && (padcol==0) && (channel==19))
			//if((detector==233) && (padrow==2) && (padcol==71) && (channel==20))
			//if((detector==233) && (padrow==5) && (padcol==72) && (channel==1))
			//if((detector==246) && (padrow==12) && (padcol==53) && (channel==20))
			//if((detector==260) && (padrow==11) && (padcol==35) && (channel==20))
			//if((detector==249) && (padrow==11) && (padcol==73) && (channel==0))
			if((detector==512) && (padrow==7) && (padcol==73) && (channel==0))
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
