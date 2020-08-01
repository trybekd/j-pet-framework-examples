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
 *  @file OPSReconstructor.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "OPSReconstructor.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetGeomMapping/JPetGeomMapping.h>
#include <JPetMCHit/JPetMCHit.h>
#include "JPetOpsEvent.h"
#include <fstream>
using namespace jpet_options_tools;

using namespace std;

OPSReconstructor::OPSReconstructor(const char* name): JPetUserTask(name) {}

bool OPSReconstructor::init()
{

  INFO("Reconstruction of o-Ps->3g decays started.");

  fOutputEvents = new JPetTimeWindow("JPetOpsEvent");

  fReconstructor = Reconstructor::GetInstance();
  double scin_length = getParamBank().getScintillators().begin()->second->getScinSize(JPetScin::kLength);
  scin_length *= 0.1; // in cm
  fReconstructor->setBarrelLength( scin_length );
  fReconstructor->setChamberRadius( 10.0 );

  getStatistics().createHistogram(new TH1F("anh_hits_in_evt", "No. annihilation hits in event", 6, -0.5, 5.5));
  getStatistics().createHistogram(new TH1F("evt_rej_z", "Was event rejected by Z requirement?", 2, -0.5, 1.5));
  getStatistics().createHistogram(new TH1F("reco_error", "Was reconstruction error?", 2, -0.5, 1.5));

  
  // create histograms for annihilation position
  getStatistics().createHistogram(new TH2F("decay point XY",
					   "transverse position of the o-Ps->3g decay point;"
					   "X [cm]; Y [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );

  getStatistics().createHistogram(new TH2F("decay point XZ",
					   "position of the o-Ps->3g decay point in XZ;"
					   "Z [cm]; X [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );

  // create histograms for annihilation position - MC signal events only
  getStatistics().createHistogram(new TH2F("decay point XY - signal",
					   "transverse position of the o-Ps->3g decay point;"
					   "X [cm]; Y [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );

  getStatistics().createHistogram(new TH2F("decay point XZ - signal",
					   "position of the o-Ps->3g decay point in XZ;"
					   "Z [cm]; X [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );
  
  getStatistics().createHistogram(new TH1F("generated multiplicity", "Generated event multiplicity",
                                           4, -0.5, 3.5
                                           )
                                  );
  
  
  getStatistics().createHistogram(new TH1F("z_res", "Hit Z resolution", 100, -5.0, 5.0));
  getStatistics().createHistogram(new TH1F("t_res", "Hit T resolution", 100, -1.0, 1.0));
  
  // optionally read a text file with additional per-strip time calibration
  if (isOptionSet(fParams.getOptions(), fTimeRecalibFileKey)){
    fTimeRecalibFile = getOptionAsString(fParams.getOptions(), fTimeRecalibFileKey);
    fShouldRecalibrateTimes = true;
  }else{
    WARNING("No time recalibration data provided by the user. Hit times will not be recalibrated.");
  }

  if(fShouldRecalibrateTimes){

    INFO("Time recalibration was requested by the user. Loading time recalibration data.");

    ifstream in_file(fTimeRecalibFile, ios::in);
    double calib_constant;
    while(in_file >> calib_constant){
      fTimeRecalibConstants.push_back(calib_constant);
    }
    in_file.close();
    
  }
  
  return true;
}

bool OPSReconstructor::exec()
{
  const JPetTimeWindowMC* time_window_mc = nullptr;
  if (time_window_mc = dynamic_cast<const JPetTimeWindowMC*>(fEvent)) {
    fIsMC = true;    
  }
  
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n = timeWindow->getNumberOfEvents();

    for(uint i=0;i<n;++i){
      const JPetEvent & event = timeWindow->getEvent<JPetEvent>(i);

      // if not 3gamma event, do not attempt to reconstruct
      if( !event.isTypeOf(JPetEventType::k3Gamma) ){
        JPetOpsEvent ops_event(event);
	fOutputEvents->add<JPetOpsEvent>(ops_event);        
        continue;
      }
      
      // remove deexcitation photon candidate
      auto hits = event.getHits();
      hits.erase(
		 remove_if(hits.begin(), hits.end(),
			   [](const JPetHit& hit){
			     if( hit.getQualityOfEnergy() > 1.2 || hit.getQualityOfEnergy() < 0.8 ){
			       return true;
			     }
			     return false;
			   }),
		 hits.end());


      // temporary - skip event if more than 3 annihilation hit candiadtes
      getStatistics().getHisto1D("anh_hits_in_evt")->Fill(hits.size());
      if( hits.size() > 3){
        continue;
      }

      /********************************************************************/
      /* MC analysis                                                      */
      /********************************************************************/
      if(fIsMC){

        // check z and t resolution
        for(auto & hit : hits){
          auto& mc_hit = time_window_mc->getMCHit<JPetMCHit>(hit.getMCindex());
          getStatistics().getHisto1D("z_res")->Fill(hit.getPos().Z() - mc_hit.getPos().Z());
          getStatistics().getHisto1D("t_res")->Fill((hit.getTime() - mc_hit.getTime())/1000.);
        }

      }
      /**************************************************************************/
      /* End of MC analysis                                                     */
      /**************************************************************************/
      
      //
      bool z_ok = true;

      std::vector<JPetHit*> annih_hits;
      
      for(auto & hit : hits){
	if( fabs(hit.getPos().Z()) > 23.0 ){
	  z_ok = false;
	}
      }

      // skip events with bad Z of hits
      getStatistics().getHisto1D("evt_rej_z")->Fill(!z_ok);
      if(!z_ok){
        continue;
      }

      // if time recalibration was requested by the user,
      // apply the calibration to hit times
      if(fShouldRecalibrateTimes){
        for(auto & hit : hits){
          int bs_id = hit.getBarrelSlot().getID();
          hit.setTime(hit.getTime() - fTimeRecalibConstants.at(bs_id-1)*1000.);
        }
      }

      // test with no smearing at all
      /*
      if(fIsMC){

        // use the exact hit positions from MCHits
        for(auto & hit : hits){
          auto& mc_hit = time_window_mc->getMCHit<JPetMCHit>(hit.getMCindex());
          hit.setPosX( mc_hit.getPosX() );
          hit.setPosY( mc_hit.getPosY() );
        }

      }
      */      

      // setup reconstruction
      fReconstructor->setGamma(0, hits.at(0));
      fReconstructor->setGamma(1, hits.at(1));
      fReconstructor->setGamma(2, hits.at(2));

      int error = 0;
      TVector3 sol[2];
      double t[2];
      error = fReconstructor->getSolution(sol[0], t[0], 0);
      error = fReconstructor->getSolution(sol[1], t[1], 1);
      
      getStatistics().getHisto1D("reco_error")->Fill(error!=0);
      
      if(error == 0){
	getStatistics().getHisto2D("decay point XY")->Fill(sol[1].X(), sol[1].Y());
	getStatistics().getHisto2D("decay point XZ")->Fill(sol[1].Z(), sol[1].X());

        // test only for events where all hits came from the same annihilation
        if(fIsMC){
          auto& mc_hit0 = time_window_mc->getMCHit<JPetMCHit>(hits.at(0).getMCindex());
          auto& mc_hit1 = time_window_mc->getMCHit<JPetMCHit>(hits.at(1).getMCindex());
          auto& mc_hit2 = time_window_mc->getMCHit<JPetMCHit>(hits.at(2).getMCindex());
        
          if( mc_hit0.getMCVtxIndex() == mc_hit1.getMCVtxIndex() &&
              mc_hit0.getMCVtxIndex() == mc_hit2.getMCVtxIndex() ){

            getStatistics().getHisto1D("generated multiplicity")->Fill(mc_hit0.getGenGammaMultiplicity());
          
            getStatistics().getHisto2D("decay point XY - signal")->Fill(sol[1].X(), sol[1].Y());
            getStatistics().getHisto2D("decay point XZ - signal")->Fill(sol[1].Z(), sol[1].X());
          }
        }
        
	JPetOpsEvent ops_event(event);
	ops_event.setAnnihilationPoint(sol[1]);
	ops_event.setAnnihilationTime(t[1]*1000.); // time resulting from reconstructor is in [ns]
	
	fOutputEvents->add<JPetOpsEvent>(ops_event);
	
      }
    } // end loop over events in a time window
      
  } else {
    return false;
  }
  return true;
}

bool OPSReconstructor::terminate()
{
  INFO("Reconstruction of o-Ps->3g decays done.");
  return true;
}

