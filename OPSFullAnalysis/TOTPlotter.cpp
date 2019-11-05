/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
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
 *  @file TOTPlotter.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "TOTPlotter.h"
#include "OPSAnalysisTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetMCHit/JPetMCHit.h>

using namespace jpet_options_tools;
using namespace ops_analysis_tools;
using namespace std;

TOTPlotter::TOTPlotter(const char* name): JPetUserTask(name) {}

JPetHit TOTPlotter::calculateTOT(const JPetHit& hit){

  double tot = 0.;
  for(auto& phys_sig: {hit.getSignalA(), hit.getSignalB()}){

    const JPetRawSignal& raw_sig = phys_sig.getRecoSignal().getRawSignal();
    for(auto& p: raw_sig.getTimesVsThresholdValue(JPetSigCh::Leading)){
      // fill multiplicity of particular threshold values at leading edge
      getStatistics().getHisto1D("thr_mult")->Fill(p.second.first);
    }

    // calculate single-signal TOT
    auto leading_points = raw_sig.getTimesVsThresholdNumber(JPetSigCh::Leading);
    auto trailing_points = raw_sig.getTimesVsThresholdNumber(JPetSigCh::Trailing);
      
    for(int i=1;i<5;i++){
      auto lead_search = leading_points.find(i);
      auto trail_search = trailing_points.find(i);
      if (lead_search != leading_points.end()
          && trail_search != trailing_points.end()){
        tot += (trail_search->second - lead_search->second) / 1000.; // in ns
      }
    }
    
  }
  
  JPetHit new_hit = hit;
  new_hit.setEnergy(tot);
  return new_hit;
}

bool TOTPlotter::init()
{

  INFO("TOT plotting started");

  fOutputEvents = new JPetTimeWindow("JPetHit");

  /************************************************************************/
  /* Histograms for data                                                  */
  /************************************************************************/
  for(int i=1; i <= getParamBank().getBarrelSlotsSize(); ++i){

    getStatistics().createHistogram(
				    new TH1F(Form("tot_strip_good_%d", getParamBank().getBarrelSlot(i).getID()),
                                             "TOT for good hits; TOT [ns]; counts", 1000., 0., 100.)
				    );
    getStatistics().createHistogram(
				    new TH1F(Form("tot_strip_corrupted_%d", getParamBank().getBarrelSlot(i).getID()),
                                             "TOT for corrupted hits; TOT [ns]; counts", 1000., 0., 100.)
				    );

  }
  
  getStatistics().createHistogram(new TH1F("thr_mult", "threshold value multiplicity in raw signals",
                                           100, 0., 500.)
                                  );

  /**************************************************************************/
  /* Histograms for MC                                                      */
  /**************************************************************************/
  getStatistics().createHistogram(
                                  new TH1F("mc_edep_2g", 
                                           "MC deposited energy - 2g anh. photons; E [keV]; counts",
                                           1300., 0., 1300.)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("mc_edep_3g", 
                                           "MC deposited energy - 3g anh. photons; E [keV]; counts",
                                           1300., 0., 1300.)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("mc_edep_prompt", 
                                           "MC deposited energy - prompt photons; E [keV]; counts",
                                           1300., 0., 1300.)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("mc_edep_scat", 
                                           "MC deposited energy - secondary scatterings; E [keV]; counts",
                                           1300., 0., 1300.)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("mc_edep_unknown", 
                                           "MC deposited energy - secondary scatterings; E [keV]; counts",
                                           1300., 0., 1300.)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("n_prompt_in_tw", "No. prompt photon candidates in time window",
                                           100, -0.5, 99.5)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("n_anh_in_tw", "No. prompt photon candidates in time window",
                                           100, -0.5, 99.5)
                                  );
  
  
  // load limits of TOT/Edep cuts
  for(auto key : fTOTcutKeys){
    if (isOptionSet(fParams.getOptions(), key)){
      fTOTcuts.push_back(getOptionAsFloat(fParams.getOptions(), key));
    }else{
      ERROR(Form("TOT cut value (%s) not provided by the user!", key.c_str()));
    }
  }
  if(fTOTcuts.size() == fTOTcutKeys.size()){
    INFO(Form("Loaded TOT cut values: (%lf, %lf) and (%lf, %lf).",
              fTOTcuts[0], fTOTcuts[1], fTOTcuts[2], fTOTcuts[3]));
  }
  
  for(auto key : fEdepMCcutKeys){
    if (isOptionSet(fParams.getOptions(), key)){
      fMCEdepCuts.push_back(getOptionAsFloat(fParams.getOptions(), key));
    }else{
      ERROR(Form("MC Edep cut value (%s) not provided by the user!", key.c_str()));
    }
  }
  if(fMCEdepCuts.size() == fEdepMCcutKeys.size()){
    INFO(Form("Loaded Edep cut values for MC: (%lf, %lf) and (%lf, %lf).",
              fMCEdepCuts[0], fMCEdepCuts[1], fMCEdepCuts[2], fMCEdepCuts[3]));  
  }
    
  return true;
}

bool TOTPlotter::exec()
{
  const JPetTimeWindowMC* time_window_mc = nullptr;
  if (time_window_mc = dynamic_cast<const JPetTimeWindowMC*>(fEvent)) {
    fIsMC = true;    
  }

  auto time_window = dynamic_cast<const JPetTimeWindow*>(fEvent);

  uint n = time_window->getNumberOfEvents();

  uint n_anh_hits = 0;
  uint n_prompt_hits = 0;
  
  for(uint i=0;i<n;++i){

    const JPetHit & hit =  dynamic_cast<const JPetHit&>(time_window->operator[](i));
    
    JPetHit new_hit;
    if(!fIsMC){
      new_hit = calculateTOT(hit);
    }else{
      new_hit = hit;
      new_hit.setRecoFlag(JPetHit::Good);
    }

    // count annihilation and prompt photon candidates
    HitCandidateType type = identifyHitType(new_hit, fTOTcuts, fIsMC, fMCEdepCuts);
    if(type==HitCandidateType::Annihilation){
      n_anh_hits++;
    }else if(type==HitCandidateType::Prompt){
      n_prompt_hits++;
    }
    
    if(fIsMC){ // fill MC histograms of deposited energy
      
      auto& mc_hit = time_window_mc->getMCHit<JPetMCHit>(new_hit.getMCindex());
      int hit_type = mc_hit.getGenGammaMultiplicity();
      std::string histo_name = "mc_edep_unknown";
      if(hit_type == 1){
        histo_name = "mc_edep_prompt";
      }else if(hit_type == 2){
        histo_name = "mc_edep_2g";
      }else if(hit_type == 3){
        histo_name = "mc_edep_3g";
      }else if(hit_type >= 100){
        histo_name = "mc_edep_scat";
      }

      getStatistics().getHisto1D(histo_name.c_str())->Fill(new_hit.getEnergy());
      
    }else{ // fill TOT histograms for data

      double tot = new_hit.getEnergy();

      if(new_hit.getRecoFlag() == JPetHit::Good){
        getStatistics().getHisto1D(Form("tot_strip_good_%d", new_hit.getBarrelSlot().getID()))->Fill(tot);
      }else{
        getStatistics().getHisto1D(Form("tot_strip_corrupted_%d", new_hit.getBarrelSlot().getID()))->Fill(tot);
      }
    }
    
    fOutputEvents->add<JPetHit>(new_hit);
  }

  getStatistics().getHisto1D("n_prompt_in_tw")->Fill(n_prompt_hits);
  getStatistics().getHisto1D("n_anh_in_tw")->Fill(n_anh_hits);

  return true;
}

bool TOTPlotter::terminate()
{
  INFO("Done plotting TOT-s");
  return true;
}

