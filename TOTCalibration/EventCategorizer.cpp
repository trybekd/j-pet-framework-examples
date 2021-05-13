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

#include <JPetAnalysisTools/JPetAnalysisTools.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "HitFinderTools.h"
#include "EventCategorizerTools.h"
#include "EventCategorizer.h"
#include <iostream>

using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char* name): JPetUserTask(name) {}

EventCategorizer::~EventCategorizer() {}

bool EventCategorizer::init()
{
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
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    vector<JPetEvent> events;
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));

      // Check types of current event
      bool is2Gamma = EventCategorizerTools::checkFor2Gamma(
        event, getStatistics(), fSaveControlHistos, fB2BSlotThetaDiff, fMaxTimeDiff
      );
      bool is3Gamma = EventCategorizerTools::checkFor3Gamma(
        event, getStatistics(), fSaveControlHistos
      );
      bool isPrompt = EventCategorizerTools::checkForPrompt(
        event, getStatistics(), fSaveControlHistos, fDeexTOTCutMin, fDeexTOTCutMax, fTOTCalculationType
      );
      bool isScattered = EventCategorizerTools::checkForScatter(
        event, getStatistics(), fSaveControlHistos, fScatterTOFTimeDiff, fTOTCalculationType
								);

      JPetEvent newEvent = event;
      if(is2Gamma) newEvent.addEventType(JPetEventType::k2Gamma);
      if(is3Gamma) newEvent.addEventType(JPetEventType::k3Gamma);
      if(isPrompt) newEvent.addEventType(JPetEventType::kPrompt);
      if(isScattered) newEvent.addEventType(JPetEventType::kScattered);

      auto mult = event.getHits().size();
      for (auto hit : event.getHits()) {
	if(fSaveControlHistos){
	  auto tot = HitFinderTools::calculateTOT(hit, HitFinderTools::getTOTCalculationType(fTOTCalculationType));
	  getStatistics().fillHistogram("TOT_all_hits", tot);
	  //EPR TOT per scint                                                                                                                                                                             
	  getStatistics().fillHistogram("tot_per_scin",
					tot, (float)(hit.getScintillator().getID()));                                                                                                                               
	  getStatistics().fillHistogram("tot_per_scin_zpos",
					tot, (float)(hit.getScintillator().getID()), hit.getPosZ());
	  if(mult==1){
	    getStatistics().fillHistogram("TOT_all_hits_mult1", tot);
	    getStatistics().fillHistogram("tot_per_scin_mult1",
					  tot, (float)(hit.getScintillator().getID()));
	  }
	  if(mult==2){
	    getStatistics().fillHistogram("TOT_all_hits_mult2", tot);
	    getStatistics().fillHistogram("tot_per_scin_mult2",
					  tot, (float)(hit.getScintillator().getID()));
	  }
	  if(mult==3){
	    getStatistics().fillHistogram("TOT_all_hits_mult3", tot);
	    getStatistics().fillHistogram("tot_per_scin_mult3",
					  tot, (float)(hit.getScintillator().getID()));
	  }
	  if(mult==4){
	    getStatistics().fillHistogram("TOT_all_hits_mult4", tot);
	    getStatistics().fillHistogram("tot_per_scin_mult4",
					  tot, (float)(hit.getScintillator().getID()));
	  }
	  if(mult>4){
	    getStatistics().fillHistogram("TOT_all_hits_multgt4", tot);
	    getStatistics().fillHistogram("tot_per_scin_multgt4",
					  tot, (float)(hit.getScintillator().getID()));
	  }
      
	  getStatistics().fillHistogram("All_XYpos", hit.getPosX(), hit.getPosY());
        }//end savehistos
      }//end hits iter
      events.push_back(newEvent);
    }
    saveEvents(events);
  } else { return false; }
  return true;
}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
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
    new TH1D("Deex_TOT_cut", "TOT of all hits with deex cut (30,50) ns", 200, 24950.0, 54950.0),
    "TOT [ps]", "Number of Hits"
  );

  //calibration plots (already in HitFinder)
  
  getStatistics().createHistogramWithAxes(
    new TH2D("tot_per_scin", "Hit TOT per Scintillator ID",
    250, -255., 99750.0, 192, 0.5, 192.5),
    "TOT hit", "ID of Scintillator"
  );

  getStatistics().createHistogramWithAxes(
    new TH3D("tot_per_scin_zpos", "Hit TOT per Scintillator ID and Z position",
	     250, -255., 99750.0, 192, 0.5, 192.5, 200, -49.75, 50.25),
    "TOT hit", "ID of Scintillator","Z [cm]"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("tot_per_scin_mult1", "Hit TOT per Scintillator ID multiplicity 1",
    250, -255., 99750.0, 192, 0.5, 192.5),
    "TOT hit with mult==1", "ID of Scintillator"
  );
  getStatistics().createHistogramWithAxes(
    new TH2D("tot_per_scin_mult2", "Hit TOT per Scintillator ID multiplicity 1",
    250, -255., 99750.0, 192, 0.5, 192.5),
    "TOT hit with mult==1", "ID of Scintillator"
  );
  getStatistics().createHistogramWithAxes(
    new TH2D("tot_per_scin_mult3", "Hit TOT per Scintillator ID multiplicity 1",
    250, -255., 99750.0, 192, 0.5, 192.5),
    "TOT hit with mult==1", "ID of Scintillator"
  );
  getStatistics().createHistogramWithAxes(
    new TH2D("tot_per_scin_mult4", "Hit TOT per Scintillator ID multiplicity 1",
    250, -255., 99750.0, 192, 0.5, 192.5),
    "TOT hit with mult==1", "ID of Scintillator"
  );
  getStatistics().createHistogramWithAxes(
    new TH2D("tot_per_scin_multgt4", "Hit TOT per Scintillator ID multiplicity 1",
    250, -255., 99750.0, 192, 0.5, 192.5),
    "TOT hit with mult==1", "ID of Scintillator"
  );

  getStatistics().createHistogramWithAxes(
    new TH2D("hit_pos_per_scin", "Hit Position per Scintillator ID",
    200, -49.75, 50.25, 192, 0.5, 192.5),
    "Hit z position [cm]", "ID of Scintillator"
  );

  // TOT calculating for all hits and reco flags
  getStatistics().createHistogramWithAxes(
    new TH1D("TOT_all_hits", "TOT of all hits", 200, -250.0, 99750.0),
    "Time over Threshold [ps]", "Number of Hits"
  );
    getStatistics().createHistogramWithAxes(
					    new TH1D("TOT_all_hits_mult1", "TOT of all hits", 200, -250.0, 99750.0),
					    "Time over Threshold [ps]", "Number of Hits"
					    );
    getStatistics().createHistogramWithAxes(
					    new TH1D("TOT_all_hits_mult2", "TOT of all hits", 200, -250.0, 99750.0),
					    "Time over Threshold [ps]", "Number of Hits"
					    );
    getStatistics().createHistogramWithAxes(
					    new TH1D("TOT_all_hits_mult3", "TOT of all hits", 200, -250.0, 99750.0),
					    "Time over Threshold [ps]", "Number of Hits"
					    );
    getStatistics().createHistogramWithAxes(
					    new TH1D("TOT_all_hits_mult4", "TOT of all hits", 200, -250.0, 99750.0),
					    "Time over Threshold [ps]", "Number of Hits"
					    );
    getStatistics().createHistogramWithAxes(
					    new TH1D("TOT_all_hits_multgt4", "TOT of all hits", 200, -250.0, 99750.0),
					    "Time over Threshold [ps]", "Number of Hits"
					    );

}
