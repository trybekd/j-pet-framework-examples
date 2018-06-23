/**
<<<<<<< HEAD
 *  @copyright Copyright 2017 The J-PET Framework Authors. All rights reserved.
=======
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
>>>>>>> kdulski_EvCatTools
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
 *  @file EventCategorizerTools.h
 */

<<<<<<< HEAD
#ifndef _EVENTCATEGORIZERTOOLS_H_
#define _EVENTCATEGORIZERTOOLS_H_

#include <JPetStatistics/JPetStatistics.h>
#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>
#define kLightVelocity_cm_ns 29.9792458
#define kUndefinedValue 999.0

/**
 * @brief Tools for Event Categorization
 *
 * Lots of tools in constatnt developement.
 */
class EventCategorizerTools
{
public:
  static bool checkFor2Gamma(const JPetEvent& event, JPetStatistics& stats,
    bool saveHistos, double b2bSlotThetaDiff);
  static bool checkFor3Gamma(const JPetEvent& event, JPetStatistics& stats, bool saveHistos);
  static bool checkForPrompt(const JPetEvent& event, JPetStatistics& stats,
    bool saveHistos, double deexTOTCutMin, double deexTOTCutMax);
  static bool checkForScatter(const JPetEvent& event, JPetStatistics& stats,
    bool saveHistos, double scatterTOFTimeDiff);
  static double calculateTOT(const JPetHit& hit);
  static double calculateDistance(const JPetHit& hit1, const JPetHit& hit2);
  static double calculateScatteringTime(const JPetHit& hit1, const JPetHit& hit2);
  static double calculateScatteringAngle(const JPetHit& hit1, const JPetHit& hit2);
  static double calculateTOF(const JPetHit& firstHit, const JPetHit& latterHit);
  static TVector3 calculateAnnihilationPoint(const JPetHit& firstHit, const JPetHit& latterHit);
};

#endif /*  !EVENTCATEGORIZERTOOLS_H */
=======
#ifndef EVENTCATEGORIZERTOOLS_H
#define EVENTCATEGORIZERTOOLS_H
#include <vector>
#include <JPetHit/JPetHit.h>
#include <JPetStatistics/JPetStatistics.h>
#include <TRotation.h>
#define kLightVelocityCmS 29.979246

class EventCategorizerTools
{
public:
  
	static double calcTOT( JPetHit hit );
	static double calcScattAngle( JPetHit hit1, JPetHit hit2 );
	static double calcScattTime( JPetHit hit1, JPetHit hit2 );
	
	static double calcDistance( JPetHit hit1, JPetHit hit2 );
	static double calcAngle( JPetHit hit1, JPetHit hit2 );
	static double calcAngle2D( JPetHit hit1, JPetHit hit2 );
	static std::vector<double> calcAnglesFrom3Hit( JPetHit hit1, JPetHit hit2, JPetHit hit3 );
	static std::vector<double> calcAngles2DFrom3Hit( JPetHit hit1, JPetHit hit2, JPetHit hit3 );
	static double calcDistanceOfSurfaceAndZero( JPetHit hit1, JPetHit hit2, JPetHit hit3 );

	static TVector3 recoPosition( const JPetHit & hit1, const JPetHit & hit2);
	
	// Not officially introduced
	static TVector3 recoPosition3Hit( JPetHit hit1, JPetHit hit2, JPetHit hit3 );
	static TVector3 findIntersection( TVector3 hit1Pos, TVector3 hit2Pos, TVector3 hit3Pos, double t21, double t31 );
	static double findMinFromQuadraticFit( std::vector<double> arg, std::vector<double> val );
	static double findMinFromDerrivative( std::vector<double> arg, std::vector<double> val );
	static std::vector< std::vector<double> > findIntersectionPointsOfCircles( TVector3 hit1Pos, TVector3 hit2Pos, TVector3 hit3Pos, double R1, double R2, double R3, double R13, double R21, double R32 );
	// Not officially introduced
	
	static int checkIfScattered( JPetHit hit1, JPetHit hit2, double errorInterval );
	
	static TVector3 decayInto3PosReco_AlekVersion( JPetHit hit1, JPetHit hit2, JPetHit hit3 );
	static TVector3 decayInto3PosReco_AlekVersion_withAddition( JPetHit hit1, JPetHit hit2, JPetHit hit3 );
	
	static double calcDistFromCentres( TVector3 sol1, TVector3 gamma1, TVector3 gamma2, TVector3 gamma3 );
	static double normalizeTime( JPetHit hit1 );
	static double normalizeTimeToPoint( JPetHit hit1, TVector3 point );

	
};

#endif
>>>>>>> kdulski_EvCatTools
