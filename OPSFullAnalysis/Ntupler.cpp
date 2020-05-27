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
 *  @file Ntupler.cpp
 */

#include <iostream>
#include "Ntupler.h"
#include "OPSAnalysisTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;
using namespace ops_analysis_tools;
using namespace std;

Ntupler::Ntupler(const char* name): JPetUserTask(name) {}

bool Ntupler::init()
{
  INFO("Started reduction of data to ntuples.");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  if (isOptionSet(fParams.getOptions(), fThresholdValuesKey)){
    fThresholdValues = getOptionAsVectorOfInts(fParams.getOptions(), fThresholdValuesKey);

    INFO(Form("Loaded the following threshold values: (%d, %d, %d, %d)",
              fThresholdValues[0],
              fThresholdValues[1],
              fThresholdValues[2],
              fThresholdValues[3]
              ));

  }else{
    WARNING("TOT values were not provided by the user!");
  }

  if(isOptionSet(fParams.getOptions(), "file_std::vector<std::string>")){
    fOutFileName = getOptionAsVectorOfStrings(fParams.getOptions(), "file_std::vector<std::string>").front();    
  }
  
  // initialize output tree
  fOutFileName.append(".root");
  fOutFile = new TFile(fOutFileName.c_str(), "RECREATE");
  fOutTree = new TTree("T", "o-Ps event candidates");
  
  fOutTree->Branch("nhits", &fNumberOfHits, "nhits/b");
  fOutTree->Branch("times", &fHitTimes);
  fOutTree->Branch("pos", &fHitPos);
  fOutTree->Branch("tots_flat", &fHitTOTsFlat);
  fOutTree->Branch("tots_proportional", &fHitTOTsProportional);  
  fOutTree->Branch("scins", &fHitScinIDs);

  return true;
}

bool Ntupler::exec()
{    
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    int n_events = timeWindow->getNumberOfEvents();
    
    for (int entry = 0; entry < n_events; ++entry) {
      const JPetEvent &event =
          dynamic_cast<const JPetEvent &>(timeWindow->operator[](entry));

      const auto &hits = event.getHits();
      fNumberOfHits = event.getHits().size();

      for (auto &hit : hits) {

        double tot_proportional = calculateTOTproportional(hit);        
        double tot_flat = hit.getEnergy();

        // reject hits with badly reconstructed TOT
        if (tot_proportional <= 0. ||
            tot_flat <= 0.) {
          continue;
        }

        fHitPos.push_back(hit.getPos());
        fHitTimes.push_back(hit.getTime() / 1000.);
        fHitScinIDs.push_back(hit.getBarrelSlot().getID());
        fHitTOTsFlat.push_back(tot_flat);        
        fHitTOTsProportional.push_back(tot_proportional);
      }

      fOutTree->Fill();
      resetRow();
    }
    
  } else {
    return false;
  }
  return true;
}

bool Ntupler::terminate()
{
  fOutTree->Write();
  fOutFile->Close();

  // delete fOutTree;
  // delete fOutFile;
  
  INFO("Finished reduction of data to ntuples.");  
  return true;
}

double Ntupler::calculateTOTproportional(const JPetHit& hit) const {

  auto getThresholdWeight = [&](int i) -> double {
    if (i <= 1) {
      return 1.;
    }
    return 1.*(fThresholdValues[i-1] - fThresholdValues[i-2]) / fThresholdValues[0];
  };
  
  double tot = 0.;
  
  for(auto& phys_sig: {hit.getSignalA(), hit.getSignalB()}){
    
    const JPetRawSignal& raw_sig = phys_sig.getRecoSignal().getRawSignal();
    
    auto leading_points = raw_sig.getTimesVsThresholdNumber(JPetSigCh::Leading);
    auto trailing_points = raw_sig.getTimesVsThresholdNumber(JPetSigCh::Trailing);

    for (int i = 1; i <= 4; i++) {
      auto lead_search = leading_points.find(i);
      auto trail_search = trailing_points.find(i);
      if (lead_search != leading_points.end()
          && trail_search != trailing_points.end()){
        tot += getThresholdWeight(i) *
          (trail_search->second - lead_search->second) / 1000.; // in ns
      }
    }
  }

  
  return tot;
}

void Ntupler::resetRow() {

  fNumberOfHits = 0;
  fHitScinIDs.clear();
  fHitPos.clear();
  fHitTimes.clear();
  fHitTOTsFlat.clear();
  fHitTOTsProportional.clear();
  
}
