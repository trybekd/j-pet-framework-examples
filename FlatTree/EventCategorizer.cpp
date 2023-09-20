/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.cpp
 */

#include "EventCategorizer.h"
#include "EventCategorizerTools.h"
#include "SourcePos.h"
#include <JPetMCHit/JPetMCHit.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include <TSystem.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using namespace jpet_options_tools;
using namespace std;

const double EventCategorizer::kUnknownEventType = 66666666;

EventCategorizer::EventCategorizer(const char* name): JPetUserTask(name){}

EventCategorizer::~EventCategorizer() {
  std::cout << "in destructor" << std::endl;
  
  if(fOutFile){
    //fOutFile->Delete("FlatTree;1");
  }

  if (pTree22) {
    
    //delete pTree22;
    //pTree22 = nullptr;
  }
}

bool EventCategorizer::init()
{
  //al below should be put into a function
  //try to read input file
  opts = fParams.getOptions();
  fInName = jpet_options_tools::getInputFile(opts);
  
  // srand((unsigned) time(0));
  // int randomNumber = rand();
  // auto uniqueName = std::to_string(randomNumber);


  std::string delimiter = "/";

  size_t pos = 0;
  std::string token;
  while ((pos = fInName.find(delimiter)) != std::string::npos) {
    token = fInName.substr(0, pos);
    //    std::cout << token << std::endl;
    fInName.erase(0, pos + delimiter.length());
  }

  std::string outname = "flatTree_"+fInName;
  if(gSystem->AccessPathName(outname.c_str())){
    std::cout << "writting output file " << outname << std::endl;
    fOutFile = std::make_unique<TFile>(outname.c_str(),"CREATE");

  } else {
    std::cout << "file " << outname << " exists: ABORTING... " << std::endl;
    ERROR("File exist.");
    exit(-1);
  }

  pTree22 = new TTree("FlatTree","flat tree");
  pTree22->Branch("timeWindowNumber",&fTimeWindowNumber,"timeWindowNumber/I");
  pTree22->Branch("eventNumber",&fEventNumber,"eventNumber/I");
  pTree22->Branch("numberOfHits",&fNumberOfHits,"numberOfHits/I");
  pTree22->Branch("x",&fPosX);
  pTree22->Branch("y",&fPosY);
  pTree22->Branch("z",&fPosZ);
  pTree22->Branch("energy",&fEnergy);
  pTree22->Branch("time",&fTime);
  pTree22->Branch("hitType",&fHitType);
  pTree22->Branch("vtxIndex",&fVtxIndex);
  pTree22->Branch("recoOrthoVtxX",&fRecoOrthoVtxPosX, "recoOrthoVtxX/F");
  pTree22->Branch("recoOrthoVtxY",&fRecoOrthoVtxPosY, "recoOrthoVtxY/F");
  pTree22->Branch("recoOrthoVtxZ",&fRecoOrthoVtxPosZ, "recoOrthoVtxZ/F");
  pTree22->Branch("isAcc",&fIsAcc, "isAcc/O");
  pTree22->Branch("isSecondary",&fIsSecondary, "isSecondary/O");
  pTree22->Branch("isScattered",&fIsScattered, "isScattered/O");
  pTree22->Branch("isOPs",&fIsOPs, "isOPs/O");
  pTree22->Branch("isPickOff",&fIsPickOff, "isPickOff/O");
  pTree22->Branch("containsPrompt",&fContainsPrompt, "containsPrompt/O");
  pTree22->Branch("timeDiff",&fTimeDiff, "timeDiff/F");
  
  INFO("Event categorization started.");
  // Parameter for back to back categorization
  if (isOptionSet(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey)){
    fB2BSlotThetaDiff = getOptionAsFloat(fParams.getOptions(), kBack2BackSlotThetaDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kBack2BackSlotThetaDiffParamKey.c_str(), fB2BSlotThetaDiff
    ));
  }
  // Parameter for scattering determination
  if (isOptionSet(fParams.getOptions(), kScatterTOFTimeDiffParamKey)) {
    fScatterTOFTimeDiff = getOptionAsFloat(fParams.getOptions(), kScatterTOFTimeDiffParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kScatterTOFTimeDiffParamKey.c_str(), fScatterTOFTimeDiff
    ));
  }
  // Parameters for deexcitation TOT cut
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMinParamKey)) {
    fDeexTOTCutMin = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMinParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMinParamKey.c_str(), fDeexTOTCutMin
    ));
  }
  if (isOptionSet(fParams.getOptions(), kDeexTOTCutMaxParamKey)) {
    fDeexTOTCutMax = getOptionAsFloat(fParams.getOptions(), kDeexTOTCutMaxParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kDeexTOTCutMaxParamKey.c_str(), fDeexTOTCutMax
    ));
  }
  if (isOptionSet(fParams.getOptions(), kMaxTimeDiffParamKey)) {
    fMaxTimeDiff = getOptionAsFloat(fParams.getOptions(), kMaxTimeDiffParamKey);
  } else {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kMaxTimeDiffParamKey.c_str(), fMaxTimeDiff));
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }
  if (isOptionSet(fParams.getOptions(), kTOTCalculationType)) {
    fTOTCalculationType = getOptionAsString(fParams.getOptions(), kTOTCalculationType);
  } else {
    WARNING("No TOT calculation option given by the user. Using standard sum.");
  }


  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  // Initialise hisotgrams
  if(fSaveControlHistos) initialiseHistograms();
  return true;
}

bool EventCategorizer::exec()
{

    fEventNumber = 0;
    fNumberOfHits = 0;
  
    if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    auto timeWindowMC = dynamic_cast<const JPetTimeWindowMC* const>(fEvent); 
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
        
      fEventNumber = i;

      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));

      // Check types of current event
      bool is2Gamma = EventCategorizerTools::checkFor2Gamma(
        event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff, fMaxTimeDiff
      );
      bool is3Gamma = EventCategorizerTools::checkFor3Gamma(
        event, getStatistics(), fSaveControlHistos
      );
      // bool isPrompt = EventCategorizerTools::checkForPrompt(
      //   event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType
      // );
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType
      );

      fEnergy.clear();
      fTime.clear();
      fPosX.clear();
      fPosY.clear();
      fPosZ.clear();
      fHitType.clear();
      fVtxIndex.clear();
      fRecoOrthoVtxPosX = 0;
      fRecoOrthoVtxPosY = 0;
      fRecoOrthoVtxPosZ = 0;

      fNumberOfHits = event.getHits().size();
      auto vhits = event.getHits();
      
      int count = 0;
      auto firstHit = vhits[0];
      for(const auto& hit : event.getHits()){

  count++;
	fPosX.push_back(hit.getPosX());
        fPosY.push_back(hit.getPosY());
        fPosZ.push_back(hit.getPosZ()); 
        fEnergy.push_back(hit.getEnergy()); 
        fTime.push_back(hit.getTime()); 
        if (timeWindowMC) {
          auto mcHit = timeWindowMC->getMCHit<JPetMCHit>(hit.getMCindex());
          fHitType.push_back(mcHit.getGenGammaMultiplicity());
          fVtxIndex.push_back(mcHit.getMCVtxIndex());

          const TVector3& pos = SourcePos::idSourcePosMap[mcHit.getMCVtxIndex()];

          fRecoOrthoVtxPosX = pos.x();
          fRecoOrthoVtxPosY = pos.y();
          fRecoOrthoVtxPosZ = pos.z();
          
        } else {
          fHitType.push_back(kUnknownEventType);
          fVtxIndex.push_back(0);

          fRecoOrthoVtxPosX = 0;
          fRecoOrthoVtxPosY = 0;
          fRecoOrthoVtxPosZ = 0;
        }
      }
      bool isPrompt = kFALSE;
      if (false) {
        auto hits = event.getHits();
        std::vector<float> energyDeex;
        // I need to identify the prompt out of 4 hits and remove from hits to set the gammas for the GPS IP calculation
        //  alternatively I can create a new container (vector) as a copy of the hits so this is not altered
        // I need to know also if isOPs or isPickOff
        if (hits.size() > 4)
          continue;
        else
        {
          //	  std::cout << hits.size() << std::endl;
          for (int i = 0; i < hits.size(); i++)
          {
            getStatistics().fillHistogram("Deex_Ene_mult4", hits.at(i).getEnergy());
            if (hits.at(i).getEnergy() >= 650 && hits.at(i).getEnergy() <= 1200)
            {
              getStatistics().fillHistogram("Deex_Ene_energyWindow", hits.at(i).getEnergy());
              energyDeex.push_back(hits.at(i).getEnergy());
            }
          }
          getStatistics().fillHistogram("mult_prompt", energyDeex.size());
          auto result = std::max_element(energyDeex.begin(), energyDeex.end());

          if (result != energyDeex.end())
          {
            getStatistics().fillHistogram("Deex_Ene", *result);
            isPrompt = kTRUE;
          }
          energyDeex.clear();
	}
      }

      JPetEvent newEvent = event;
      
      //reset event type
      newEvent.setEventType(JPetEventType::kUnknown);
      
      if(isPrompt){
	newEvent.addEventType(JPetEventType::kPrompt);
      }else{
	// if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
	// if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
	
	// if(isScattered) newEvent.addEventType(JPetEventType::kScattered);
      }

      if ( std::equal(fVtxIndex.begin() + 1, fVtxIndex.end(), fVtxIndex.begin()) )
      	fIsAcc = kFALSE;
      else{
      	fIsAcc = kTRUE;
      }

      for(auto i : fHitType){
	auto number = i%10;
	auto n = 1;
	if(i%100!=0 && number!=0 && i>10)
	  fIsScattered=kTRUE;
	if(number==1)
	  fContainsPrompt=kTRUE;
	else if(number==2)
	  fIsPickOff=kTRUE;
	else if(number==3)
	  fIsOPs=kTRUE;
	else if(number==0){
	  fIsSecondary=kTRUE;
	  while(number==0&&i>0){
	    n = n*10;	    
	    number = (i/n)%10;
	    if(i%n*100 == 0 && number!=0 && i/n>10)
	      fIsScattered=kTRUE;
	    if(number==1)
	      fContainsPrompt=kTRUE;
	    else if(number==2)
	      fIsPickOff=kTRUE;
	    else if(number==3)
	      fIsOPs=kTRUE;
	  }
	}
      }//end auto i:fHitType
    
      pTree22->Fill();

      fIsSecondary = kFALSE;
      fIsScattered = kFALSE;
      fContainsPrompt = kFALSE;
      fIsPickOff = kFALSE;
      fIsOPs = kFALSE;

      if(fSaveControlHistos){
        for(auto hit : event.getHits()){

          getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
        }
      }
      events.push_back(newEvent);
    }
    saveEvents(events);
  } else { 
    fTimeWindowNumber++;
    return false; 
  }
  fTimeWindowNumber++;
  return true;
}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  if(!pTree22) std::cout << "wtf??" << std::endl;
  fOutFile->cd();
  pTree22->Write();
  //pTree22->Print();
  fOutFile->Write("",TObject::kWriteDelete);
  
  fOutFile->Close();
  //delete fOutFile;
  return true;
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void EventCategorizer::initialiseHistograms(){
  
  // General histograms
  getStatistics().createHistogramWithAxes(
    new TH2D("All_XYpos", "Hit position XY", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Hit X position [cm]", "Hit Y position [cm]"
  );
  
  // Histograms for 2Gamma category
  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Zpos", "Z-axis position of 2 gamma hits", 201, -50.25, 50.25),
    "Z axis position [cm]", "Number of Hits"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_DLOR", "Delta LOR distance", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_ThetaDiff", "Angle difference of 2 gamma hits ", 181, -0.5, 180.5),
    "Hits theta diff [deg]", "Counts"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_TimeDiff", "Time difference of 2 gamma hits", 200, -10100.0, 99900.0),
    "Time Difference [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("2Gamma_Dist", "B2B hits distance", 150, -0.5, 149.5),
    "Distance [cm]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Annih_TOF", "Annihilation pairs Time of Flight", 201, -3015.0, 3015.0),
    "Time of Flight [ps]", "Number of Annihilation Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_XY", "XY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "X position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_ZX", "ZX position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "X position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("AnnihPoint_ZY", "ZY position of annihilation point", 240, -60.25, 59.75, 240, -60.25, 59.75),
    "Z position [cm]", "Y position [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH1D("Annih_DLOR", "Delta LOR distance of annihilation photons", 100, -0.25, 49.25),
    "Delta LOR [cm]", "Counts"
  );

  // Histograms for 3Gamama category
  getStatistics().createHistogramWithAxes(
    new TH2D("3Gamma_Angles", "Relative angles - transformed", 250, -0.5, 249.5, 20, -0.5, 199.5),
    "Relative angle 1-2", "Relative angle 2-3"
  );

  // Histograms for scattering category
  getStatistics().createHistogramWithAxes(
    new TH1D("ScatterTOF_TimeDiff", "Difference of Scatter TOF and hits time difference",
    3.0*fScatterTOFTimeDiff, -0.5, 3.0*fScatterTOFTimeDiff-0.5),
    "Scat_TOF - time diff [ps]", "Number of Hit Pairs"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_PrimaryTOT", "Angle of scattering vs. TOT of primary hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of primary hit [ps]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("ScatterAngle_ScatterTOT", "Angle of scattering vs. TOT of scattered hits",
    200, -0.5, 199.5, 200, -100.0, 39900.0),
    "Scattering Angle", "TOT of scattered hit [ps]"
  );

  // Histograms for deexcitation
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_TOT_cut", "Energy all hits", 200, 0.0, 2000.0),
    "Energy [keV]", "Number of Hits"
  );
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_Ene_mult4", "Energy 4 hits multiplicity", 200, 0.0, 2000.0),
    "Energy [keV]", "Number of Hits"
  );
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_Ene_energyWindow", "Energy 4 hits multiplicity between [650, 1200] [kEv]", 200, 0.0, 2000.0),
    "Energy [keV]", "Number of Hits"
  );
  getStatistics().createHistogramWithAxes(
    new TH1D("Deex_Ene", "Energy Prompt 4 hits multiplicity between [650, 1200] [kEv]", 200, 0.0, 2000.0),
    "Energy [keV]", "Number of Hits"
  );
  getStatistics().createHistogramWithAxes(
    new TH1D("mult_prompt", "Multiplicity prompt candidates in Energy Window [650, 1200] [kEv]", 10, -.5, 9.5),
    "Energy [keV]", "Number of Hits"
  );
  getStatistics().createHistogramWithAxes(
    new TH1D("htimeDiff", "Absolute time difference between hits", 250, -.5, 300000),
    "Abs(t_{i} - t_{j})", "Number of Hits"
  );
}

