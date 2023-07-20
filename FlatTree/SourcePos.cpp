//
// Created by dtr on 20.07.23.
//

#include "SourcePos.h"

#include <JPetGeantEventPack/JPetGeantEventPack.h>

SourcePos::SourcePos(const char* name) : JPetUserTask(name) {}

std::map<unsigned int, TVector3> SourcePos::idSourcePosMap = {};

SourcePos::~SourcePos() {}

bool SourcePos::init() { return true; }

bool SourcePos::exec()
{
  if (const auto& mcEventPack = dynamic_cast<JPetGeantEventPack* const>(fEvent))
  {
    auto pos = mcEventPack->GetEventInformation()->GetVtxPosition();
    for (int i = 0; i < mcEventPack->GetNumberOfHits(); ++i)
    {
      auto evNum = mcEventPack->GetHit(i)->GetEvtID();
      idSourcePosMap[evNum] = pos;
    }
  }
  return true;
}

bool SourcePos::terminate() { return true; }
