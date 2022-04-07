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
 *  @file HitFinder.cpp
 */

using namespace std;

#include <JPetAnalysisTools/JPetAnalysisTools.h>
#include <JPetGeomMapping/JPetGeomMapping.h>
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>

#include "HitFinder.h"
#include "HitFinderTools.h"
#include "ToTEnergyConverterFactory.h"
#include "UniversalFileLoader.h"

#include <map>
#include <string>
#include <vector>

using namespace tot_energy_converter;
using namespace jpet_options_tools;
using namespace boost::property_tree;

HitFinder::HitFinder(const char *name) : JPetUserTask(name) {}

HitFinder::~HitFinder() {}

bool HitFinder::init() {
  INFO("Hit finding Started");
  fOutputEvents = new JPetTimeWindow("JPetHit");

  // Reading values from the user options if available
  // Getting bool for using bad signals
  if (isOptionSet(fParams.getOptions(), kUseCorruptedSignalsParamKey)) {
    fUseCorruptedSignals =
        getOptionAsBool(fParams.getOptions(), kUseCorruptedSignalsParamKey);
    if (fUseCorruptedSignals) {
      WARNING("Hit Finder is using Corrupted Signals, as set by the user");
    } else {
      WARNING("Hit Finder is NOT using Corrupted Signals, as set by the user");
    }
  } else {
    WARNING("Hit Finder is not using Corrupted Signals (default option)");
  }
  // Allowed time difference between signals on A and B sides
  if (isOptionSet(fParams.getOptions(), kABTimeDiffParamKey)) {
    fABTimeDiff = getOptionAsFloat(fParams.getOptions(), kABTimeDiffParamKey);
  }
  // Getting velocities file from user options
  auto velocitiesFile = std::string("dummyCalibration.txt");
  if (isOptionSet(fParams.getOptions(), kVelocityFileParamKey)) {
    velocitiesFile =
        getOptionAsString(fParams.getOptions(), kVelocityFileParamKey);
  } else {
    WARNING(
        "No path to the file with velocities was provided in user options.");
  }
  // Getting number of Reference Detector Scintillator ID
  if (isOptionSet(fParams.getOptions(), kRefDetScinIDParamKey)) {
    fRefDetScinID = getOptionAsInt(fParams.getOptions(), kRefDetScinIDParamKey);
  } else {
    WARNING(Form("No value of the %s parameter provided by the user, "
                 "indicating that Reference Detector was not used.",
                 kRefDetScinIDParamKey.c_str()));
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos =
        getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }

  // Use of velocities file
  JPetGeomMapping mapper(getParamBank());
  auto tombMap = mapper.getTOMBMapping();
  fVelocities =
      UniversalFileLoader::loadConfigurationParameters(velocitiesFile, tombMap);
  if (fVelocities.empty()) {
    ERROR("Velocities map seems to be empty");
  }

  // Loading parameters for conversion to ToT to energy
  if (isOptionSet(fParams.getOptions(), kConvertToTParamKey)) {
    fConvertToT = getOptionAsBool(fParams.getOptions(), kConvertToTParamKey);
    if (fConvertToT) {
      INFO("Hit finder performs conversion of ToT to deposited energy with "
           "provided params.");
      fToTConverterFactory.loadConverterOptions(fParams.getOptions());
    } else {
      INFO("Hit finder will not convert ToT to deposited energy since no user "
           "parameters are provided.");
    }
  }

  // Loading parameters for TOT synchronizationi
  if (isOptionSet(fParams.getOptions(), kUseToTSyncParamKey)) {
    fSyncToT = getOptionAsBool(fParams.getOptions(), kUseToTSyncParamKey);
    if (fSyncToT) {
      if (isOptionSet(fParams.getOptions(), kTOTConstantsFileParamKey)) {
        INFO("Hit finder will perform ToT synchronization.");
	std::string kSync = getOptionAsString(fParams.getOptions(), kTOTConstantsFileParamKey);
	read_json(kSync, fConstantsTree);
      } else {
        ERROR("No file for TOT synchronization provided. No synchroniztion "
              "applied.");
      }
    } else {
      WARNING("Hit finder won't performe the ToT synchronization.");
    }
  }

  // if (isOptionSet(fParams.getOptions(), kTOTConstantsFileParamKey)) {
  //   read_json(getOptionAsString(fParams.getOptions(),
  //   kTOTConstantsFileParamKey), fTOTConstantsTree);
  // } else {
  //   INFO("No file for TOT synchronization provided. No synchroniztion
  //   applied");
  // }

  if (isOptionSet(fParams.getOptions(), kTOTCalculationType)) {
    fTOTCalculationType =
        getOptionAsString(fParams.getOptions(), kTOTCalculationType);
  } else {
    WARNING("No TOT calculation option given by the user. Using standard sum.");
  }

  // Control histograms
  if (fSaveControlHistos) {
    initialiseHistograms();
  }
  return true;
}

bool HitFinder::exec() {
  if (auto &timeWindow = dynamic_cast<const JPetTimeWindow *const>(fEvent)) {
    auto signalsBySlot =
        HitFinderTools::getSignalsBySlot(timeWindow, fUseCorruptedSignals);
    auto totConverter = fToTConverterFactory.getEnergyConverter();
    auto allHits = HitFinderTools::matchAllSignals(
        signalsBySlot, fVelocities, fABTimeDiff, fRefDetScinID, fConvertToT,
        totConverter, getStatistics(), fSaveControlHistos);
    if (fSaveControlHistos) {
      getStatistics().fillHistogram("hits_per_time_slot", allHits.size());
    }
    saveHits(allHits);
  } else
    return false;
  return true;
}

bool HitFinder::terminate() {
  INFO("Hit finding ended");
  return true;
}

void HitFinder::saveHits(const std::vector<JPetHit> &hits) {
  auto sortedHits = JPetAnalysisTools::getHitsOrderedByTime(hits);
  getStatistics().fillHistogram("mult_hits", sortedHits.size());
  auto mult = sortedHits.size();

  for (const auto &hit : sortedHits) {
    if (fSaveControlHistos) {

      // for the sides
      auto type = HitFinderTools::getTOTCalculationType(fTOTCalculationType);
      // std::map<int, double> thrToTOT_sideA =
      // hit.getSignalA().getRecoSignal().getRawSignal().getTOTsVsThresholdValue();
      // std::map<int, double> thrToTOT_sideB =
      // hit.getSignalB().getRecoSignal().getRawSignal().getTOTsVsThresholdValue();
      // // for (auto it = thrToTOT_sideA.begin(); it != thrToTOT_sideA.end();
      // it++) {
      // // 	std::cout << "Thresholds " << it->first << " " << it->second <<
      // std::endl;
      // // }
      // auto totsideA = HitFinderTools::calculateTOTside(hit, thrToTOT_sideA,
      // type); auto totsideB = HitFinderTools::calculateTOTside(hit,
      // thrToTOT_sideB, type);
      // getStatistics().fillHistogram("totSideA_per_scin",
      // 			  totsideA,
      // (float)(hit.getScintillator().getID()));
      // getStatistics().fillHistogram("totSideB_per_scin",
      // 			  totsideB,
      // (float)(hit.getScintillator().getID()));
      // calculate tot for Zposition studies and correction for different
      // thresholds
      float tot;
      float tot1;
      float tot2;
      float tot3;
      // synchronization ==> needs to be moved somewhere else for official
      if (fSyncToT) {
	
	tot = HitFinderTools::syncTOT(hit, tot, fConstantsTree);
	tot1 = HitFinderTools::syncTOT(hit, tot1, fConstantsTree);
        tot2 = HitFinderTools::syncTOT(hit, tot2, fConstantsTree);
        tot3 = HitFinderTools::syncTOT(hit, tot3, fConstantsTree);
      }
      else{
	tot = HitFinderTools::calculateTOT(hit, getStatistics(), type);
	tot1 = HitFinderTools::calculateTOTPlot(hit, 30, type);	       
	tot2 = HitFinderTools::calculateTOTPlot(hit, 80, type);
	tot3 = HitFinderTools::calculateTOTPlot(hit, 190, type);
      }
      // EPR TOT per scint
      getStatistics().fillHistogram("TOT_all_hits", tot);
      getStatistics().fillHistogram("tot_per_scin", tot,
                                    (float)(hit.getScintillator().getID()));

      getStatistics().fillHistogram("tot_per_scin_zpos", tot,
                                    (float)(hit.getScintillator().getID()),
                                    hit.getPosZ());

      getStatistics().fillHistogram("tot_vs_zpos", tot, hit.getPosZ());
      getStatistics().fillHistogram("tot_vs_zpos1", tot1, hit.getPosZ());
      getStatistics().fillHistogram("tot_vs_zpos2", tot2, hit.getPosZ());
      getStatistics().fillHistogram("tot_vs_zpos3", tot3, hit.getPosZ());

      for (int scinID = 1; scinID <= 192; scinID++) {
        if (TMath::Abs(hit.getPosZ()) <= 23 &&
            hit.getScintillator().getID() == scinID)
          getStatistics().fillHistogram(Form("tot_zPos_scinID%i", scinID), tot,
                                        hit.getPosZ());
      }
      if (mult == 1) {
        getStatistics().fillHistogram("TOT_all_hits_mult1", tot);
        getStatistics().fillHistogram("tot_per_scin_mult1", tot,
                                      (float)(hit.getScintillator().getID()));
      }
      if (mult == 2) {
        getStatistics().fillHistogram("TOT_all_hits_mult2", tot);
        getStatistics().fillHistogram("tot_per_scin_mult2", tot,
                                      (float)(hit.getScintillator().getID()));
      }
      if (mult == 3) {
        getStatistics().fillHistogram("TOT_all_hits_mult3", tot);
        getStatistics().fillHistogram("tot_per_scin_mult3", tot,
                                      (float)(hit.getScintillator().getID()));
      }
      if (mult == 4) {
        getStatistics().fillHistogram("TOT_all_hits_mult4", tot);
        getStatistics().fillHistogram("tot_per_scin_mult4", tot,
                                      (float)(hit.getScintillator().getID()));
      }
      if (mult > 4) {
        getStatistics().fillHistogram("TOT_all_hits_multgt4", tot);
        getStatistics().fillHistogram("tot_per_scin_multgt4", tot,
                                      (float)(hit.getScintillator().getID()));
      }
      // end EPR
      if (hit.getRecoFlag() == JPetHit::Good) {
        getStatistics().fillHistogram("TOT_good_hits", tot);
      } else if (hit.getRecoFlag() == JPetHit::Corrupted) {
        getStatistics().fillHistogram("TOT_corr_hits", tot);
      }
    }
    fOutputEvents->add<JPetHit>(hit);
  }
}

void HitFinder::initialiseHistograms() {

  getStatistics().createHistogramWithAxes(
      new TH1D("good_vs_bad_hits", "Number of good and corrupted Hits created",
               3, 0.5, 3.5),
      "Quality", "Number of Hits");

  std::vector<std::pair<unsigned, std::string>> binLabels;
  binLabels.push_back(std::make_pair(1, "GOOD"));
  binLabels.push_back(std::make_pair(2, "CORRUPTED"));
  binLabels.push_back(std::make_pair(3, "UNKNOWN"));
  getStatistics().setHistogramBinLabel(
      "good_vs_bad_hits", getStatistics().AxisLabel::kXaxis, binLabels);

  getStatistics().createHistogramWithAxes(
      new TH1D("hits_per_time_slot", "Number of Hits in Time Window", 100, -0.5,
               99.5),
      "Hits in Time Slot", "Number of Time Slots");

  getStatistics().createHistogramWithAxes(
      new TH2D("time_diff_per_scin",
               "Signals Time Difference per Scintillator ID",
               4 * fABTimeDiff / 10, -2 * fABTimeDiff, 2 * fABTimeDiff, 192,
               0.5, 192.5),
      "A-B time difference", "ID of Scintillator");
  // TOT needs to have the "correct" boundaries for all mehods
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_per_scin", "Hit TOT per Scintillator ID", 2 * 250, -255.,
               199500.0, 192, 0.5, 192.5),
      "TOT hit", "ID of Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH2D("totSideA_per_scin", "Hit TOT side A per Scintillator ID",
               2 * 250, -255., 199500.0, 192, 0.5, 192.5),
      "TOT hit", "ID of Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH2D("totSideB_per_scin", "Hit TOT side B per Scintillator ID",
               2 * 250, -255., 199500.0, 192, 0.5, 192.5),
      "TOT hit", "ID of Scintillator");
  // be careful with the scint numers
  for (int scinID = 1; scinID <= 192; scinID++) {
    getStatistics().createHistogramWithAxes(
        new TH2D(Form("tot_zPos_scinID%i", scinID),
                 "Hit TOT vs Z position per single scint.", 2 * 250, -255.,
                 199500.0, 47, -23.5, 23.5),
        "TOT hit", "Z [cm]");
  }
  getStatistics().createHistogramWithAxes(
      new TH3D("tot_per_scin_zpos",
               "Hit TOT per Scintillator ID and Z position", 250, -255.,
               99750.0, 192, 0.5, 192.5, 47, -23.5, 23.5),
      "TOT hit", "ID of Scintillator", "Z [cm]");
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_vs_zpos", "Hit TOT along Z position", 2 * 250, -255.,
               199500.0, 47, -23.5, 23.5),
      "TOT hit", "ID of Scintillator", "Z [cm]");
  // for thresholds
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_vs_zpos1", "Hit TOT along Z position threshold 1, 2, 3, 4",
               2 * 250, -255., 199500.0, 47, -23.5, 23.5),
      "TOT hit", "ID of Scintillator", "Z [cm]");
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_vs_zpos2", "Hit TOT along Z position threshold 2, 3, 4",
               175, -255., 69825.0, 47, -23.5, 23.5),
      "TOT hit", "ID of Scintillator", "Z [cm]");

  getStatistics().createHistogramWithAxes(
      new TH2D("tot_vs_zpos3", "Hit TOT along Z position threshold 3, 4", 125,
               -255., 49875.0, 47, -23.5, 23.5),
      "TOT hit", "ID of Scintillator", "Z [cm]");
  // end thresholds z position tot

  getStatistics().createHistogramWithAxes(
      new TH2D("tot_per_scin_mult1",
               "Hit TOT per Scintillator ID multiplicity 1", 2 * 250, -255.,
               199500.0, 192, 0.5, 192.5),
      "TOT hit with mult==1", "ID of Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_per_scin_mult2",
               "Hit TOT per Scintillator ID multiplicity 1", 2 * 250, -255.,
               199500.0, 192, 0.5, 192.5),
      "TOT hit with mult==1", "ID of Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_per_scin_mult3",
               "Hit TOT per Scintillator ID multiplicity 1", 2 * 250, -255.,
               199500.0, 192, 0.5, 192.5),
      "TOT hit with mult==1", "ID of Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_per_scin_mult4",
               "Hit TOT per Scintillator ID multiplicity 1", 2 * 250, -255.,
               199500.0, 192, 0.5, 192.5),
      "TOT hit with mult==1", "ID of Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH2D("tot_per_scin_multgt4",
               "Hit TOT per Scintillator ID multiplicity 1", 2 * 250, -255.,
               199500.0, 192, 0.5, 192.5),
      "TOT hit with mult==1", "ID of Scintillator");
  // remember active material covers only 23 cm each side from center
  getStatistics().createHistogramWithAxes(
      new TH2D("hit_pos_per_scin", "Hit Position per Scintillator ID", 200,
               -49.75, 50.25, 192, 0.5, 192.5),
      "Hit z position [cm]", "ID of Scintillator");
  // multiplicity ordered hits
  getStatistics().createHistogramWithAxes(
      new TH1D("mult_hits", "multiplicity hits", 50, -.5, 49.5), "multiplicity",
      "Number of Hits");

  // TOT calculating for all hits and reco flags
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_all_hits", "TOT of all hits", 2 * 250, -255.0, 199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_all_hits_mult1", "TOT of all hits", 2 * 250, -255.0,
               199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_all_hits_mult2", "TOT of all hits", 2 * 250, -255.0,
               199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_all_hits_mult3", "TOT of all hits", 2 * 250, -255.0,
               199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_all_hits_mult4", "TOT of all hits", 2 * 250, -255.0,
               199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_all_hits_multgt4", "TOT of all hits", 2 * 250, -255.0,
               199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_good_hits", "TOT of hits with GOOD flag", 2 * 250, -255.0,
               199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("TOT_corr_hits", "TOT of hits with CORRUPTED flag", 2 * 250,
               -255.0, 199500.0),
      "Time over Threshold [ps]", "Number of Hits");
  getStatistics().createHistogramWithAxes(
      new TH1D("remain_signals_per_scin",
               "Number of Unused Signals in Scintillator", 192, 0.5, 192.5),
      "ID of Scintillator", "Number of Unused Signals in Scintillator");
  getStatistics().createHistogramWithAxes(
      new TH1D("remain_signals_tdiff",
               "Time Diff of an unused signal and the consecutive one", 200,
               fABTimeDiff - 125.0, 49875.0 + fABTimeDiff),
      "Time difference [ps]", "Number of Signals");

  if (fConvertToT) {
    auto converterRange = fToTConverterFactory.getEnergyConverter().getRange();
    auto totConverter = fToTConverterFactory.getEnergyConverter();

    auto minToT = converterRange.first;
    auto maxToT = converterRange.second;
    auto minEDep = totConverter(converterRange.first);
    auto maxEDep = totConverter(converterRange.second);

    getStatistics().createHistogramWithAxes(
        new TH1D("conv_tot_range",
                 "TOT of hits in range of conversion function", 200, minToT,
                 maxToT),
        "Time over Threshold [ps]", "Number of Hits");
    getStatistics().createHistogramWithAxes(
        new TH1D(
            "conv_dep_energy",
            "Deposited energy of hits, converted from ToT with provied formula",
            200, minEDep, maxEDep),
        "Deposited energy [keV]", "Number of Hits");
    getStatistics().createHistogramWithAxes(
        new TH2D("conv_dep_energy_vs_tot",
                 "Deposited energy of hits, converted from ToT with provied "
                 "formula vs. input ToT",
                 200, minEDep, maxEDep, 200, minToT, maxToT),
        "Deposited energy [keV]", "ToT of Hit [ps]");
  }

  // getStatistics().createHistogramWithAxes(
  //     new TH1D("TOTdiff", "TOT(i) - TOT(i-1)", 250, -5000.0, 5000.0),
  //     "TOTdiff [ps]", "Number of Hits");
  // getStatistics().createHistogramWithAxes(
  //     new TH1D(
  //         "reversed_tot_thr2",
  //         "Number hit with TOTi < TOTi-1 == 1 else TOTi > TOTi-1 with thi = 80",
  //         2, -0.5, 1.5),
  //     "True/False", "Number of Hits");
  // getStatistics().createHistogramWithAxes(
  //     new TH1D("reversed_tot_thr3",
  //              "Number hit with TOTi < TOTi-1 == 1 else TOTi > TOTi-1 with thi "
  //              "= 190",
  //              2, -0.5, 1.5),
  //     "True/False", "Number of Hits");
  // getStatistics().createHistogramWithAxes(
  //     new TH1D("reversed_tot_thr4",
  //              "Number hit with TOTi < TOTi-1 == 1 else TOTi > TOTi-1 with thi "
  //              "= 300",
  //              2, -0.5, 1.5),
  //     "True/False", "Number of Hits");
  // getStatistics().createHistogramWithAxes(
  //     new TH2D("thr_vs_tot", "thrNr vs TOT of the hits", 250, -5000.0,
  //              30000.0, 6, -0.5, 5.5),
  //     "ToT of Hit [ps]", "THR [#number]");
}
