//
// Created by dtr on 20.07.23.
//

#ifndef J_PET_FRAMEWORK_EXAMPLES_SOURCEPOS_H
#define J_PET_FRAMEWORK_EXAMPLES_SOURCEPOS_H

#include <JPetUserTask/JPetUserTask.h>

class SourcePos : public JPetUserTask
{

public:
  SourcePos(const char* name);
  virtual ~SourcePos();
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;

  static std::map<unsigned int, TVector3> idSourcePosMap;
};

#endif // J_PET_FRAMEWORK_EXAMPLES_SOURCEPOS_H
