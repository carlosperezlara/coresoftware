//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: PHG4GDMLWriteParamvol.hh 77528 2013-11-25 13:06:36Z gcosmo $
//
//
// class PHG4GDMLWriteParamvol
//
// Class description:
//
// GDML class for writing parameterised entities dimensions.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _PHG4GDMLWRITEPARAMVOL_INCLUDED_
#define _PHG4GDMLWRITEPARAMVOL_INCLUDED_

#include "PHG4GDMLWriteSetup.hh"

class G4Box;
class G4Trd;
class G4Trap;
class G4Tubs;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Ellipsoid;
class G4Para;
class G4Hype;
class G4Polycone;
class G4Polyhedra;
class G4VPhysicalVolume;

class PHG4GDMLWriteParamvol : public PHG4GDMLWriteSetup
{

 public:

   virtual void ParamvolWrite(xercesc::DOMElement*,
                              const G4VPhysicalVolume* const);
   virtual void ParamvolAlgorithmWrite(xercesc::DOMElement* paramvolElement,
                              const G4VPhysicalVolume* const paramvol);

 protected:

   PHG4GDMLWriteParamvol();
   virtual ~PHG4GDMLWriteParamvol();

   void Box_dimensionsWrite(xercesc::DOMElement*, const G4Box* const);
   void Trd_dimensionsWrite(xercesc::DOMElement*, const G4Trd* const);
   void Trap_dimensionsWrite(xercesc::DOMElement*, const G4Trap* const);
   void Tube_dimensionsWrite(xercesc::DOMElement*, const G4Tubs* const);
   void Cone_dimensionsWrite(xercesc::DOMElement*, const G4Cons* const);
   void Sphere_dimensionsWrite(xercesc::DOMElement*, const G4Sphere* const);
   void Orb_dimensionsWrite(xercesc::DOMElement*, const G4Orb* const);
   void Torus_dimensionsWrite(xercesc::DOMElement*, const G4Torus* const);
   void Ellipsoid_dimensionsWrite(xercesc::DOMElement*, const G4Ellipsoid* const);
   void Para_dimensionsWrite(xercesc::DOMElement*, const G4Para* const);
   void Hype_dimensionsWrite(xercesc::DOMElement*, const G4Hype* const);
   void Polycone_dimensionsWrite(xercesc::DOMElement*, const G4Polycone* const);
   void Polyhedra_dimensionsWrite(xercesc::DOMElement*, const G4Polyhedra* const);
   void ParametersWrite(xercesc::DOMElement*,
                        const G4VPhysicalVolume* const, const G4int&);

};

#endif
