/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file HitFinderTools.cpp
 */

#include "HitFinderTools.h"
#include "UniversalFileLoader.h"
#include <TMath.h>
#include <cmath>
#include <map>
#include <vector>

using namespace tot_energy_converter;
using namespace std;

/**
 * Helper method for sotring signals in vector
 */
void HitFinderTools::sortByTime(vector<JPetPhysSignal> &sigVec) {
  sort(sigVec.begin(), sigVec.end(),
       [](const JPetPhysSignal &sig1, const JPetPhysSignal &sig2) {
         return sig1.getTime() < sig2.getTime();
       });
}

/**
 * Method distributing Signals according to Scintillator they belong to
 */
map<int, vector<JPetPhysSignal>>
HitFinderTools::getSignalsBySlot(const JPetTimeWindow *timeWindow,
                                 bool useCorrupts) {
  map<int, vector<JPetPhysSignal>> signalSlotMap;
  if (!timeWindow) {
    WARNING("Pointer of Time Window object is not set, returning empty map");
    return signalSlotMap;
  }
  const unsigned int nSignals = timeWindow->getNumberOfEvents();
  for (unsigned int i = 0; i < nSignals; i++) {
    auto physSig =
        dynamic_cast<const JPetPhysSignal &>(timeWindow->operator[](i));
    if (!useCorrupts && physSig.getRecoFlag() == JPetBaseSignal::Corrupted) {
      continue;
    }
    int slotID = physSig.getBarrelSlot().getID();
    auto search = signalSlotMap.find(slotID);
    if (search == signalSlotMap.end()) {
      vector<JPetPhysSignal> tmp;
      tmp.push_back(physSig);
      signalSlotMap.insert(pair<int, vector<JPetPhysSignal>>(slotID, tmp));
    } else {
      search->second.push_back(physSig);
    }
  }
  return signalSlotMap;
}

/**
 * Loop over all Scins invoking matching procedure
 */
vector<JPetHit> HitFinderTools::matchAllSignals(
    map<int, vector<JPetPhysSignal>> &allSignals,
    const map<unsigned int, vector<double>> &velocitiesMap, double timeDiffAB,
    int refDetScinId, bool convertToT, const ToTEnergyConverter &totConverter,
    JPetStatistics &stats, bool saveHistos) {
  vector<JPetHit> allHits;
  for (auto &slotSigals : allSignals) {
    // Loop for Reference Detector ID
    if (slotSigals.first == refDetScinId) {
      for (auto refSignal : slotSigals.second) {
        auto refHit = createDummyRefDetHit(refSignal);
        allHits.push_back(refHit);
      }
      continue;
    }
    // Loop for other slots than reference one
    auto slotHits = matchSignals(slotSigals.second, velocitiesMap, timeDiffAB,
                                 convertToT, totConverter, stats, saveHistos);
    allHits.insert(allHits.end(), slotHits.begin(), slotHits.end());
  }
  return allHits;
}

/**
 * Method matching signals on the same Scintillator
 */
vector<JPetHit> HitFinderTools::matchSignals(
    vector<JPetPhysSignal> &slotSignals,
    const map<unsigned int, vector<double>> &velocitiesMap, double timeDiffAB,
    bool convertToT, const ToTEnergyConverter &totConverter,
    JPetStatistics &stats, bool saveHistos) {
  vector<JPetHit> slotHits;
  vector<JPetPhysSignal> remainSignals;
  sortByTime(slotSignals);
  while (slotSignals.size() > 0) {
    auto physSig = slotSignals.at(0);
    if (slotSignals.size() == 1) {
      remainSignals.push_back(physSig);
      break;
    }
    for (unsigned int j = 1; j < slotSignals.size(); j++) {
      if (slotSignals.at(j).getTime() - physSig.getTime() < timeDiffAB) {
        if (physSig.getPM().getSide() != slotSignals.at(j).getPM().getSide()) {
          auto hit = createHit(physSig, slotSignals.at(j), velocitiesMap,
                               convertToT, totConverter, stats, saveHistos);
          slotHits.push_back(hit);
          slotSignals.erase(slotSignals.begin() + j);
          slotSignals.erase(slotSignals.begin() + 0);
          break;
        } else {
          if (j == slotSignals.size() - 1) {
            remainSignals.push_back(physSig);
            slotSignals.erase(slotSignals.begin() + 0);
            break;
          } else {
            continue;
          }
        }
      } else {
        if (saveHistos &&
            physSig.getPM().getSide() != slotSignals.at(j).getPM().getSide()) {
          stats.fillHistogram("remain_signals_tdiff",
                              slotSignals.at(j).getTime() - physSig.getTime());
        }
        remainSignals.push_back(physSig);
        slotSignals.erase(slotSignals.begin() + 0);
        break;
      }
    }
  }
  if (remainSignals.size() > 0 && saveHistos) {
    stats.fillHistogram("remain_signals_per_scin",
                        (float)(remainSignals.at(0).getPM().getScin().getID()),
                        remainSignals.size());
  }
  return slotHits;
}

/**
 * Method for Hit creation - setting all fields, that make sense here
 */
JPetHit HitFinderTools::createHit(
    const JPetPhysSignal &signal1, const JPetPhysSignal &signal2,
    const map<unsigned int, vector<double>> &velocitiesMap, bool convertToT,
    const ToTEnergyConverter &totConverter, JPetStatistics &stats,
    bool saveHistos) {
  JPetPhysSignal signalA;
  JPetPhysSignal signalB;
  if (signal1.getPM().getSide() == JPetPM::SideA) {
    signalA = signal1;
    signalB = signal2;
  } else {
    signalA = signal2;
    signalB = signal1;
  }
  auto radius = signalA.getPM().getBarrelSlot().getLayer().getRadius();
  auto theta = TMath::DegToRad() * signalA.getPM().getBarrelSlot().getTheta();
  auto velocity = UniversalFileLoader::getConfigurationParameter(
      velocitiesMap, getProperChannel(signalA));
  checkTheta(theta);

  JPetHit hit;
  hit.setSignalA(signalA);
  hit.setSignalB(signalB);
  hit.setTime((signalA.getTime() + signalB.getTime()) / 2.0);
  hit.setQualityOfTime(-1.0);
  hit.setTimeDiff(signalB.getTime() - signalA.getTime());
  hit.setQualityOfTimeDiff(-1.0);
  hit.setScintillator(signalA.getPM().getScin());
  hit.setBarrelSlot(signalA.getPM().getBarrelSlot());
  hit.setPosX(radius * cos(theta));
  hit.setPosY(radius * sin(theta));
  hit.setPosZ(velocity * hit.getTimeDiff() / 2000.0);

  if (convertToT) {
    auto tot = calculateTOT(hit, stats);
    /// Checking if provided conversion function accepts calculated value of ToT
    if (tot > totConverter.getRange().first &&
        tot < totConverter.getRange().second) {
      auto energy = totConverter(tot);
      if (!isnan(energy)) {
        hit.setEnergy(energy);
        stats.fillHistogram("conv_tot_range", tot);
        stats.fillHistogram("conv_dep_energy", energy);
        stats.fillHistogram("conv_dep_energy_vs_tot", energy, tot);
      } else {
        hit.setEnergy(-1.0);
      }
    } else {
      hit.setEnergy(-1.0);
    }
  } else {
    hit.setEnergy(-1.0);
  }
  hit.setQualityOfEnergy(-1.0);

  if (signalA.getRecoFlag() == JPetBaseSignal::Good &&
      signalB.getRecoFlag() == JPetBaseSignal::Good) {
    hit.setRecoFlag(JPetHit::Good);
    if (saveHistos) {
      stats.fillHistogram("good_vs_bad_hits", 1);
      stats.fillHistogram("time_diff_per_scin", hit.getTimeDiff(),
                          (float)(hit.getScintillator().getID()));
      stats.fillHistogram("hit_pos_per_scin", hit.getPosZ(),
                          (float)(hit.getScintillator().getID()));
    }
  } else if (signalA.getRecoFlag() == JPetBaseSignal::Corrupted ||
             signalB.getRecoFlag() == JPetBaseSignal::Corrupted) {
    hit.setRecoFlag(JPetHit::Corrupted);
    if (saveHistos) {
      stats.fillHistogram("good_vs_bad_hits", 2);
    }
  } else {
    hit.setRecoFlag(JPetHit::Unknown);
    if (saveHistos) {
      stats.fillHistogram("good_vs_bad_hits", 3);
    }
  }
  return hit;
}

/**
 * Method for Hit creation in case of reference detector.
 * Setting only some necessary fields.
 */
JPetHit HitFinderTools::createDummyRefDetHit(const JPetPhysSignal &signalB) {
  JPetHit hit;
  hit.setSignalB(signalB);
  hit.setTime(signalB.getTime());
  hit.setQualityOfTime(-1.0);
  hit.setTimeDiff(0.0);
  hit.setQualityOfTimeDiff(-1.0);
  hit.setEnergy(-1.0);
  hit.setQualityOfEnergy(-1.0);
  hit.setScintillator(signalB.getPM().getScin());
  hit.setBarrelSlot(signalB.getPM().getBarrelSlot());
  return hit;
}

/**
 * Helper method for getting TOMB channel
 */
int HitFinderTools::getProperChannel(const JPetPhysSignal &signal) {
  auto someSigCh = signal.getRecoSignal().getRawSignal().getPoints(
      JPetSigCh::Leading, JPetRawSignal::ByThrValue)[0];
  return someSigCh.getTOMBChannel().getChannel();
}

/**
 * Helper method for checking if theta is in radians
 */
void HitFinderTools::checkTheta(const double &theta) {
  if (theta > 2 * TMath::Pi()) {
    WARNING("Probably wrong values of Barrel Slot theta - conversion to "
            "radians failed. Check please.");
  }
}

/**
 * Calculation of the total TOT of the hit - Time over Threshold:
 * the weighted sum of the TOTs on all of the thresholds (1-4) and on the both
 * sides (A,B)
 */

HitFinderTools::TOTCalculationType
HitFinderTools::getTOTCalculationType(const std::string &type) {
  if (type == "rectangular") {
    return TOTCalculationType::kThresholdRectangular;
  } else if (type == "trapeze") {
    return TOTCalculationType::kThresholdTrapeze;
  } else if (type == "standard") {
    return TOTCalculationType::kSimplified;
  } else if (type == "elena") {
    return TOTCalculationType::kThresholdElena;
  } else if (type == "triangle") {
    return TOTCalculationType::kThresholdTriangle;
  } else {
    WARNING("Undefinied type for the calculation of the TOTs. Probably missing "
            "option in userParams, typo, or mistake.");
    return TOTCalculationType::kSimplified;
  }
}

double HitFinderTools::calculateTOT(const JPetHit &hit, JPetStatistics &stats,
                                    TOTCalculationType type) {
  double tot = 0.0;

  std::map<int, double> thrToTOT_sideA =
      hit.getSignalA().getRecoSignal().getRawSignal().getTOTsVsThresholdValue();
  std::map<int, double> thrToTOT_sideB =
      hit.getSignalB().getRecoSignal().getRawSignal().getTOTsVsThresholdValue();

  // auto someSigCh = signal.getRecoSignal().getRawSignal()
  // .getPoints(JPetSigCh::Leading, JPetRawSignal::ByThrValue)[0]; //this should
  // return the same at the end

  std::map<int, double> leadingPointsA =
      hit.getSignalA().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(
          JPetSigCh::Leading);
  std::map<int, double> leadingPointsB =
      hit.getSignalB().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(
          JPetSigCh::Leading);
  std::map<int, double> trailingPointsA =
      hit.getSignalA().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(
          JPetSigCh::Trailing);
  std::map<int, double> trailingPointsB =
      hit.getSignalB().getRecoSignal().getRawSignal().getTimesVsThresholdNumber(
          JPetSigCh::Trailing);
  // std::vector<JPetSigCh> ff =
  // hit.getSignalB().getRecoSignal().getRawSignal().getPoints(JPetSigCh::Leading,
  // JPetRawSignal::ByThrValue); std::vector<JPetSigCh> ll =
  // hit.getSignalB().getRecoSignal().getRawSignal().getPoints(JPetSigCh::Trailing,
  // JPetRawSignal::ByThrValue); std::cout << " times " << std::endl; for(auto
  // it: ff){

  //   // for (auto ii = leadingPointsB.begin(); ii != leadingPointsB.end();
  //   ii++){
  //   // 	std::cout << "joer " <<ii->second << std::endl;
  //   // }
  //   std::cout << "leading " << it.getThresholdNumber()<< " " <<
  //   it.getValue()<<std::endl;

  // }
  // for(auto it2: ll){
  //   std::cout << "trailing " << it2.getThresholdNumber()<< " " <<
  //   it2.getValue()<<std::endl;
  // std::cout << "TOT " << it2.getValue() - it.getValue() << std::endl;
  //  }
  // // for (auto it = std::next(ff.begin(), 1); it != ff.end(); it++) {
  // //   std::cout << "Trailing " << it->first << " " << it->second <<
  // std::endl;
  // // }
  // // for (auto it = std::next(leadingPointsA.begin(), 1); it !=
  // leadingPointsA.end(); it++) {
  // //   std::cout << "Leading " << it->first << " " << it->second <<
  // std::endl;
  // // }

  tot += calculateTOTside(hit, thrToTOT_sideA, leadingPointsA, stats, type);
  tot += calculateTOTside(hit, thrToTOT_sideB, leadingPointsB, stats, type);
  return tot;
}

double HitFinderTools::calculateTOTPlot(const JPetHit &hit, int thr,
                                        TOTCalculationType type) {
  double tot = 0.0;

  std::map<int, double> thrToTOT_sideA =
      hit.getSignalA().getRecoSignal().getRawSignal().getTOTsVsThresholdValue();
  std::map<int, double> thrToTOT_sideB =
      hit.getSignalB().getRecoSignal().getRawSignal().getTOTsVsThresholdValue();

  tot += calculateTOTsidePlot(hit, thrToTOT_sideA, thr, type);
  tot += calculateTOTsidePlot(hit, thrToTOT_sideB, thr, type);
  return tot;
}

double HitFinderTools::calculateTOTside(
    const JPetHit &hit, const std::map<int, double> &thrToTOT_side,
    std::map<int, double> trailingPoints, JPetStatistics &stats,
    TOTCalculationType type) {
  double tot = 0., weight = 1.;
  int thrNr = 0;
  auto areaTip = 0;

  if (!thrToTOT_side.empty()) {
    double firstThr = thrToTOT_side.begin()->first;
    tot += weight * thrToTOT_side.begin()->second;

    if (thrToTOT_side.begin()->first == 30)
      stats.fillHistogram("thr_vs_tot", thrToTOT_side.begin()->second, 1);
    else if (thrToTOT_side.begin()->first == 80)
      stats.fillHistogram("thr_vs_tot", thrToTOT_side.begin()->second, 2);
    else if (thrToTOT_side.begin()->first == 190)
      stats.fillHistogram("thr_vs_tot", thrToTOT_side.begin()->second, 3);
    else
      stats.fillHistogram("thr_vs_tot", thrToTOT_side.begin()->second, 4);

    if (thrToTOT_side.size() > 1) {
      for (auto it = std::next(thrToTOT_side.begin(), 1);
           it != thrToTOT_side.end(); it++) {

        if (it->first == 80) {
          if (it->second < std::prev(it, 1)->second)
            stats.fillHistogram("reversed_tot_thr2", 1);
          else
            stats.fillHistogram("reversed_tot_thr2", 0);
        } else if (it->first == 190) {
          if (it->second < std::prev(it, 1)->second)
            stats.fillHistogram("reversed_tot_thr3", 1);
          else
            stats.fillHistogram("reversed_tot_thr3", 0);
        } else {
          if (it->second < std::prev(it, 1)->second)
            stats.fillHistogram("reversed_tot_thr4", 1);
          else
            stats.fillHistogram("reversed_tot_thr4", 0);
        }

        thrNr++;
        stats.fillHistogram("TOTdiff", it->second - std::prev(it, 1)->second);
        if (it->first == 30)
          stats.fillHistogram("thr_vs_tot", it->second, 1);
        else if (it->first == 80)
          stats.fillHistogram("thr_vs_tot", it->second, 2);
        else if (it->first == 190)
          stats.fillHistogram("thr_vs_tot", it->second, 3);
        else
          stats.fillHistogram("thr_vs_tot", it->second, 4);

        switch (type) {
        case TOTCalculationType::kSimplified:
          weight = 1.;
          break;
        case TOTCalculationType::kThresholdRectangular:
          weight = (it->first - std::prev(it, 1)->first) / firstThr;
          break;
        case TOTCalculationType::kThresholdTrapeze:
          weight = (it->first - std::prev(it, 1)->first) / firstThr;
          tot += weight * fabs(it->second - std::prev(it, 1)->second) / 2;
          break;
          // this actually does the trapeze Area = (B+b)*h / 2   B=base b=top
        case TOTCalculationType::kThresholdElena:
          weight = (it->first - std::prev(it, 1)->first) / (2 * firstThr);
          tot += weight * std::prev(it, 1)->second;
          break;
        case TOTCalculationType::kThresholdTriangle:
          if ((it != thrToTOT_side.end()) &&
              (next(it) == thrToTOT_side.end() & thrNr > 2 &&
               it->first == 300)) {
            // points at the last element
            // std::cout << "going to calculate vcp " << thrNr << " " <<
            // it->first << " " << it->second << std::endl;
            areaTip = calculateAreaTriangle(hit, thrToTOT_side, trailingPoints,
                                            thrNr, it->first);
            tot += areaTip /
                   firstThr; // remember to normalize to the first threshold
          }
          weight = (it->first - std::prev(it, 1)->first) / (2 * firstThr);
          tot += weight * std::prev(it, 1)->second;
          break;
        }
        tot += weight * it->second;
      }
    }
  } else
    return 0;
  return tot;
}

double
HitFinderTools::calculateTOTsidePlot(const JPetHit &hit,
                                     const std::map<int, double> &thrToTOT_side,
                                     int thr, TOTCalculationType type) {
  double tot = -666., weight = 1.;

  int cc = 0.;
  if (!thrToTOT_side.empty()) {
    for (auto it1 = std::next(thrToTOT_side.begin(), 0);
         it1 != thrToTOT_side.end(); it1++) {
      double firstThr = it1->first;
      cc++;
      if (firstThr != thr) {
        continue;
      }
      tot += weight * it1->second;
      if (thrToTOT_side.size() > 1) {
        for (auto it = std::next(thrToTOT_side.begin(), cc);
             it != thrToTOT_side.end(); it++) {
          switch (type) {
          case TOTCalculationType::kSimplified:
            weight = 1.;
            break;
          case TOTCalculationType::kThresholdRectangular:
            weight = (it->first - std::prev(it, 1)->first) / firstThr;
            break;
          case TOTCalculationType::kThresholdTrapeze:
            weight = (it->first - std::prev(it, 1)->first) / firstThr;
            tot += weight * fabs(it->second - std::prev(it, 1)->second) / 2;
            break;
          }
          tot += weight * it->second;
        }
      }
    }
  } else
    return 0;
  return tot;
}

double HitFinderTools::calculateAreaTriangle(
    const JPetHit &hit, std::map<int, double> thrToTOT_side,
    std::map<int, double> trailingPoints, int thrNr, int thr) {
  double vcp = 0.0;
  double area = 0.0;
  int thrn_1, thrn, thrNrn_1, thrNrn;

  thrn = thr;
  if (thr == 80) {
    thrn_1 = 30;
    thrNrn_1 = 1;
    thrNrn = 2;
  } else if (thr == 190) {
    thrn_1 = 30;
    thrNrn_1 = 1;
    thrNrn = 3;
  } else if (thr == 300) {
    thrn_1 = 80;
    thrNrn_1 = 2;
    thrNrn = 4;
  }

  auto trailn_1 = trailingPoints.find(
      thrNrn_1); /// first [1,4] corresponds to the point (thr)
  auto trailn = trailingPoints.find(
      thrNrn); // second stores the time value tot +=
               // (trailSearch->second - leadSearch->second);

  auto TOTn_1 =
      thrToTOT_side.find(thrn_1); /// first [1,4] corresponds to the point (thr)
  auto TOTn =
      thrToTOT_side.find(thrn); // second stores the time value tot +=
                                // (trailSearch->second - leadSearch->second);
  // std::cout << "thn_1 and thn " << TOTn_1->first << " " << TOTn->first << "
  // tailn_1 and trailn " << trailn_1->second << " " << trailn->second << "
  // TOTn_1 and TOTn " << TOTn_1->second << "  " << TOTn->second << std::endl;

  int thn_1 = TOTn_1->first;
  int thn = TOTn->first;
  double tn_1 = TMath::Abs(trailn_1->second);
  double tn = TMath::Abs(trailn->second);
  double totn_1 = TOTn_1->second;
  double totn = TOTn->second;

  // std::cout << "thn_1 and thn " << thn_1 << " " << thn << " tailn_1 and
  // trailn " << tn_1 << " " << tn << " TOTn_1 and TOTn " << totn_1 << "  " <<
  // totn << std::endl;

  double num1 = ((thn) - (thn_1));
  double denom1 = ((tn) - (tn_1));
  double factor1 = num1 / denom1;
  double num2 = (((tn_1) * (totn_1)) - ((tn) * (totn_1)));
  double denom2 = (totn - totn_1);
  double factor2 = num2 / denom2;
  double vcp1 = factor1 * factor2;
  // std::cout << "num1 denom1 num2 denom2 " << num1 << " " << denom1 << " " <<
  // num2 << " " << denom2 << std::endl; std::cout << factor1 << " * " <<
  // factor2 << std::endl; std::cout << vcp1 << std::endl;

  vcp = vcp1 + thn_1;
  vcp = vcp - thn;
  vcp = vcp * totn;

  area = vcp / 2.;
  std::cout << "Calculated vcp " << vcp << " area " << area << std::endl;
  return area;
}
