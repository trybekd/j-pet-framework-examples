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
 *  @file OPSCleaner.h
 */

#ifndef OPSCLEANER_H 
#define OPSCLEANER_H 

#include <vector>
#include <map>
#include <string>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetHit/JPetHit.h>
#include <JPetEvent/JPetEvent.h>
#include "JPetOpsEvent.h"
#include "OPSAnalysisTools.h"

#ifdef __CINT__
#	define override
#endif

enum MCEventType{
  ANY,
  SIGNAL,
  POSSIBLE_SIGNAL,
  ALL_BCG,
  B2B_SCAT,
  B2B_PROMPT,
  RANDOM,
  OTHER
};

struct EventInfo{
  bool event_odd;
  MCEventType type;
  ops_analysis_tools::SmallestAngles angles_xy;
  ops_analysis_tools::SmallestAngles angles_3d;
  ops_analysis_tools::Kinematics3g kinematics;
  ops_analysis_tools::Operators operators;
  ops_analysis_tools::Distances distances;
  ops_analysis_tools::ScatterTests scatter_tests;
};

class OPSCleaner : public JPetUserTask{
public:
  OPSCleaner(const char * name);
  virtual ~OPSCleaner(){}
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;
protected:
  bool fIsMC = false;
  void bookHisto(TH1* h);
  TH1F* getHisto1D(std::string name, MCEventType type, bool selected);
  TH2F* getHisto2D(std::string name, MCEventType type, bool selected);
  void fillHistos(const JPetOpsEvent& event, const EventInfo& evt_info, MCEventType type, bool selected);
  MCEventType classifyEvent(const JPetTimeWindowMC* time_window_mc, const JPetEvent& event);
  
  const std::string fAngleSumCutKey = "OPSCleaner_angles_sum_cut_float";
  float fAngleSumCut;
  
  const double kSpeedOfLight = 29.9792458; // cm  / ns
  const int kNsteps = 10;
  
  int fEventCouters[10];

  std::map<MCEventType, std::string> fHistoSuffices = {
    {ANY, "_all"},
    {SIGNAL, "_sig"},
    {POSSIBLE_SIGNAL, "_pos_sig"},
    {ALL_BCG, "_bcg"},
    {B2B_SCAT, "_b2b_scat"},
    {B2B_PROMPT, "_b2b_prompt"},
    {RANDOM, "_random"},
    {OTHER, "_other"}
  };
  
};
#endif /*  !OPSCLEANER_H */










