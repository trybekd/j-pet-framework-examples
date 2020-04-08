/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
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
 *  @file OPSCleaner.cpp
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <functional>
#include "OPSCleaner.h"
#include "OPSAnalysisTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetMCHit/JPetMCHit.h>

using namespace jpet_options_tools;
using namespace ops_analysis_tools;

using namespace std;



OPSCleaner::OPSCleaner(const char* name): JPetUserTask(name) {}

void OPSCleaner::bookHisto(TH1* h){

  std::string base_name = std::string(h->GetName());
  for(auto& tt: fHistoSuffices){
    std::string name = base_name + std::string("_init") + tt.second;
    getStatistics().createHistogram(h->Clone(name.c_str()));
    name = base_name + std::string("_sel") + tt.second;
    getStatistics().createHistogram(h->Clone(name.c_str()));
  }
  
}

TH1F* OPSCleaner::getHisto1D(std::string name, MCEventType type, bool selected = false){
  
  std::string full_name = name
    + (selected ? std::string("_sel") : std::string("_init") )
    + fHistoSuffices.at(type);
  return getStatistics().getHisto1D(full_name.c_str());
}

TH2F* OPSCleaner::getHisto2D(std::string name, MCEventType type, bool selected = false){
  
  std::string full_name = name
    + (selected ? std::string("_sel") : std::string("_init") )
    + fHistoSuffices.at(type);
  return getStatistics().getHisto2D(full_name.c_str());
}

bool OPSCleaner::init()
{
  INFO("Cleaning of o-Ps->3g candidate sample started.");

  fOutputEvents = new JPetTimeWindow("JPetOpsEvent");

  if (isOptionSet(fParams.getOptions(), fAngleSumCutKey)){
    fAngleSumCut = getOptionAsFloat(fParams.getOptions(), fAngleSumCutKey);
  }else{
    WARNING("Angles sum cut value not provided by the user!");
  }

  // initialize counters
  for(int i=0;i<10;++i){
    fEventCouters[i] = 0;
  }

  /************************************************************************/
  /* New histograms, to be filled after every step of the cleaning        */
  /************************************************************************/

  bookHisto(
            new TH2F("E2_E1",
                     "E2 vs E1;E1 [keV];E2 [keV]",
                     600, -0.5, 599.5,
                     600, -0.5, 599.5)
            );
  
  bookHisto(
            new TH2F("E3_E1",
                     "E3 vs E1;E1 [keV];E3 [keV]",
                     600, -0.5, 599.5,
                     600, -0.5, 599.5)
            );
  
  bookHisto(
            new TH2F("E3_E2",
                     "E3 vs E2;E2 [keV];E3 [keV]",
                     600, -0.5, 599.5,
                     600, -0.5, 599.5)
            );
  
  bookHisto(
            new TH2F("th2_th1",
                     "th2 vs th1;th1 [deg];th2 [deg]",
                     200, -0.5, 189.5,
                     200, -0.5, 189.5)
            );

  bookHisto(
            new TH2F("th3_th1",
                     "th3 vs th1;th1 [deg];th3 [deg]",
                     200, -0.5, 189.5,
                     200, -0.5, 189.5)
            );
  
  bookHisto(
            new TH2F("th3_th2",
                     "th3 vs th2;th2 [deg];th3 [deg]",
                     200, -0.5, 189.5,
                     200, -0.5, 189.5)
            );
  
  bookHisto(
            new TH2F("3_hit_angles",
                     "3 Hit 3D angles difference;"
                     "Smallest angle + Second smallest angle [deg];"
                     "Second smallest angle - Smallest angle [deg]",
                     360, -0.5, 359.5,
                     360, -0.5, 359.5)
            );
  
  bookHisto(
            new TH2F("anh_XY",
                     "transverse position of the o-Ps->3g decay point;"
                     "Y [cm]; X [cm]",
                     100, -50.5, 49.5,
                     100, -50.5, 49.5
                     )
            );
  
  bookHisto(
            new TH2F("anh_XZ",
                     "position of the o-Ps->3g decay point in XZ;"
                     "Z [cm]; X [cm]",
                     100, -50., 50.,
                     100, -50., 50.
                     )
            );  

  bookHisto(new TH1F("anh_R",
                     "Transverse radius of the o-Ps->3g decay point; R [cm]",
                     100, 0., 50.
                     )
            );

  bookHisto(new TH1F("anh_R_jacobian",
                     "Transverse radius of the o-Ps->3g decay point - Jacobian weighted; R [cm]",
                     100, 0., 50.
                     )
            );
  
  bookHisto(
            new TH1F("lifetime",
                     "Time between deexcitation and annihilation;"
                     "#Delta t [ns]",
                     2000, -260., 260.)
            );  

  bookHisto(
            new TH1F("Sk1",
                     "|S|#bullet |k_{1}|;|S|#bullet |k_{1}|",
                     1000, -1.1, 1.1)
            );

  bookHisto(
            new TH1F("Sk2",
                     "|S|#bullet |k_{2}|;|S|#bullet |k_{2}|",
                     1000, -1.1, 1.1)
            );

  bookHisto(
            new TH1F("Sk3",
                     "|S|#bullet |k_{3}|;|S|#bullet |k_{3}|",
                     1000, -1.1, 1.1)
            );

  bookHisto(
            new TH1F("Sk1xk2",
                     "|S|#bullet |k_{1}#times k_{2}|;|S|#bullet|k_{1}#times k_{2}|",
                     1000, -1.1, 1.1)
            );

  bookHisto(
            new TH1F("Sk1.Sk1xk2",
                     "(|S|#bullet |k_{1}|)(|S|#bullet |k_{1}#times k_{2}|);"
                     "(|S|#bullet |k_{1}|)(|S|#bullet |k_{1}#times k_{2}|)",
                     1000, -1.1, 1.1)
            );  

  bookHisto(
            new TH2F("dvts",
                     "d-ct;|d-ct|_{1} [cm]; |d-ct|_{2} [cm]",
                     400, 0., 100.,
                     400, 0., 100.)
            );  

  // single histos
  getStatistics().createHistogram(new TH1F("pair_mc_mult", "MCGenMult of hit pair", 120, -0.5, 119.5));
  getStatistics().createHistogram(new TH1F("triple_mc_mult", "MCGenMult of hit triple", 120, -0.5, 119.5));
  getStatistics().createHistogram(new TH1F("all_mults_different", "Were all hit MC mults  different?", 2, -0.5, 1.5));

  getStatistics().createHistogram(new TH2F("Mults_pair_single", "MCGenMult of pair vs single hit;single;pair",
                                           120, -0.5, 119.5,
                                           120, -0.5, 119.5));
  
    return true;
}

bool OPSCleaner::exec()
{
  const JPetTimeWindowMC* time_window_mc = nullptr;
  if (time_window_mc = dynamic_cast<const JPetTimeWindowMC*>(fEvent)) {
    fIsMC = true;    
  }
  
  if (auto time_window = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n_events = time_window->getNumberOfEvents();
    for(int i=0;i<n_events;++i){

      const JPetOpsEvent & event = time_window->getEvent<JPetOpsEvent>(i);

      EventInfo evt_info;
      evt_info.event_odd = (i%2==1);
      evt_info.kinematics = calcKinematics(event);
      evt_info.angles_xy = calcSmallestAnglesXY(event);
      auto& angles = std::get<0>(evt_info.kinematics); 
      evt_info.angles_3d = std::make_pair(angles[0], angles[1]);
      // shuffle the angles for the 2D relative plots
      std::random_shuffle(angles.begin(), angles.end());

      auto& energies = std::get<1>(evt_info.kinematics); 
      
      // calculate the operators
      evt_info.operators = calcOperators(event, evt_info.kinematics);

      // Calculate the LOR distances
      evt_info.distances = calcLORdistances(event);

      // Perform scatter tests
      evt_info.scatter_tests = calcScatterTests(event);
      
      // in case of MC, identify event type
      if(fIsMC){
        evt_info.type = OTHER;

        std::array<JPetMCHit, 3> mc_hits;
        for(int k=0; k<3; ++k){
          mc_hits[k] = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(k).getMCindex());        
        }
      
        if( mc_hits[0].getMCVtxIndex() == mc_hits[1].getMCVtxIndex() &&
            mc_hits[0].getMCVtxIndex() == mc_hits[2].getMCVtxIndex() ){
          // hits from same simulated event

          bool triple = false;
          bool pair = false;
          for(int j=0; j<3; ++j){
            for(int k=j+1; k<3; ++k){
              if(mc_hits[j].getGenGammaMultiplicity() == mc_hits[k].getGenGammaMultiplicity()){
                // there was a pair with the same multiplicity
                pair = true;
                int pair_mult = mc_hits[k].getGenGammaMultiplicity();
                int single_mult = mc_hits[ 3 - j - k].getGenGammaMultiplicity();
                getStatistics().getHisto1D("pair_mc_mult")->Fill( pair_mult );

                if( single_mult == pair_mult ){
                  triple = true;
                  break;
                }else{ // the third hit had a different multiplicity
                  getStatistics().getHisto2D("Mults_pair_single")->Fill(single_mult, pair_mult);

                  // cases B2B_SCAT, B2B_PROMPT
                  if(pair_mult == 2){ // there was a back-to-back event
                    if( single_mult = 1 ){
                      evt_info.type = B2B_PROMPT;
                    }
                    if( single_mult >= 100 ){
                      evt_info.type = B2B_SCAT;
                    }
                  }
                
                }
              
              }
            }
          }

          // true if there was not even a single pair
          getStatistics().getHisto1D("all_mults_different")->Fill(!pair);
        
          if(triple){ // three hits with the same MC multiplicity
            getStatistics().getHisto1D("triple_mc_mult")->Fill(mc_hits[0].getGenGammaMultiplicity());

            if(mc_hits[0].getGenGammaMultiplicity() == 3){
              evt_info.type = SIGNAL; // signal (3-photon) event!
            }

            if(mc_hits[0].getGenGammaMultiplicity() == 0){
              evt_info.type = POSSIBLE_SIGNAL; // signal (3-photon) event!
            }

          
          }
        
        }else{
          // random coincidence
          evt_info.type = RANDOM;
          int index = mc_hits[0].getMCVtxIndex();
          std::cout << "R = ("
                    << mc_hits[0].getMCVtxIndex() - index << " "
                    << mc_hits[0].getGenGammaMultiplicity() << ") ("
                    << mc_hits[1].getMCVtxIndex() - index << " "
                    << mc_hits[1].getGenGammaMultiplicity() << ") ( "
                    << mc_hits[2].getMCVtxIndex() - index << " "
                    << mc_hits[2].getGenGammaMultiplicity() << ")" << std::endl;
        }

        if(evt_info.type == OTHER){
          std::cout << "O = " << mc_hits[0].getGenGammaMultiplicity() << " "
                    << mc_hits[1].getGenGammaMultiplicity() << " "
                    << mc_hits[2].getGenGammaMultiplicity() << std::endl;
        }
      
        /********************************************************************/
        /* Study of particular cuts starts here                             */
        /********************************************************************/
        MCEventType evt_type = evt_info.type;
        fillHistos(event, evt_info, false);
        // also fill total signal/background histos
        if( evt_info.type != SIGNAL && evt_info.type != POSSIBLE_SIGNAL ){
          evt_info.type = ALL_BCG;
          fillHistos(event, evt_info, false);
        }
        // and fill histos for any kinds of events
        evt_info.type = ANY;
        fillHistos(event, evt_info, false);
        // restore the original type
        evt_info.type = evt_type;
      
        if( evt_info.scatter_tests.at(0) > 15.0 ){
        
          double r = event.getAnnihilationPoint().Perp();
          if( r > 4.0 && r < 20.0 ){
            fillHistos(event, evt_info, true);
            // also fill total signal/background histos
            if( evt_info.type != SIGNAL && evt_info.type != POSSIBLE_SIGNAL ){
              evt_info.type = ALL_BCG;
              fillHistos(event, evt_info, true);
            }
            // and fill histos for any kinds of events
            evt_info.type = ANY;
            fillHistos(event, evt_info, true);
            // restore the original type
            evt_info.type = evt_type;

          }
        }
      }
    } // end loop over events

  } else {
    return false;
  }
  return true;
}

bool OPSCleaner::terminate()
{
  INFO("Cleaning of o-Ps->3g candidate sample finished.");
  
  return true;
}

void OPSCleaner::fillHistos(const JPetOpsEvent& event, EventInfo evt_info, bool selected){

  auto& angles = std::get<0>(evt_info.kinematics); 
  auto& energies = std::get<1>(evt_info.kinematics); 

  // fill angle histograms
  getHisto2D("3_hit_angles", evt_info.type, selected)->Fill(evt_info.angles_3d.second + evt_info.angles_3d.first,
                                          evt_info.angles_3d.second - evt_info.angles_3d.first);

  if(evt_info.event_odd){
    getHisto2D("th2_th1", evt_info.type, selected)->Fill(angles[0], angles[1]);
    getHisto2D("th3_th2", evt_info.type, selected)->Fill(angles[1], angles[2]);
    getHisto2D("th3_th1", evt_info.type, selected)->Fill(angles[0], angles[2]);
  }else{
    getHisto2D("th2_th1", evt_info.type, selected)->Fill(angles[1], angles[0]);
    getHisto2D("th3_th2", evt_info.type, selected)->Fill(angles[2], angles[1]);
    getHisto2D("th3_th1", evt_info.type, selected)->Fill(angles[2], angles[0]);
  }
    
  // fill energy histograms
  if(evt_info.event_odd){
    getHisto2D("E2_E1", evt_info.type, selected)->Fill(energies[0], energies[1]);
    getHisto2D("E3_E2", evt_info.type, selected)->Fill(energies[1], energies[2]);
    getHisto2D("E3_E1", evt_info.type, selected)->Fill(energies[0], energies[2]);
  }else{
    getHisto2D("E2_E1", evt_info.type, selected)->Fill(energies[1], energies[0]);
    getHisto2D("E3_E2", evt_info.type, selected)->Fill(energies[2], energies[1]);
    getHisto2D("E3_E1", evt_info.type, selected)->Fill(energies[2], energies[0]);
  }
  
  // lifetime
  if(event.hasPrompt()){
    getHisto1D("lifetime", evt_info.type, selected)->Fill(event.getLifeTime() / 1000.);
  }

  // annihilation point location
  getHisto2D("anh_XZ", evt_info.type, selected)->Fill(event.getAnnihilationPoint().Z(), event.getAnnihilationPoint().X());

  if( fabs(event.getAnnihilationPoint().Z()) > 4.0 ){
    getHisto2D("anh_XY", evt_info.type, selected)->Fill(event.getAnnihilationPoint().Y(), event.getAnnihilationPoint().X());

    double r = event.getAnnihilationPoint().Perp();
    getHisto1D("anh_R", evt_info.type, selected)->Fill(r);
    getHisto1D("anh_R_jacobian", evt_info.type, selected)->Fill(r, 1./r);
  }
      
  // operators
  getHisto1D("Sk1", evt_info.type, selected)->Fill(evt_info.operators.at(0));
  getHisto1D("Sk1xk2", evt_info.type, selected)->Fill(evt_info.operators.at(1));
  getHisto1D("Sk1.Sk1xk2", evt_info.type, selected)->Fill(evt_info.operators.at(2));

  // scatter tests
  getHisto2D("dvts", evt_info.type, selected)->Fill(evt_info.scatter_tests.at(0), evt_info.scatter_tests.at(1));
  
}




