

#ifndef TCL_GENINTERFACEPTS_H
#define TCL_GENINTERFACEPTS_H

#include <stdlib.h>
#include <string.h>

#include <tcl.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <EmbeddedBeamInterface.h>
#include <EmbeddedBeamInterfaceL.h>
#include <EmbeddedBeamInterfaceP.h>
#include <EmbeddedBeamInterfaceP2.h>
#include <EmbeddedBeamInterfaceAL.h>
#include <EmbeddedBeamInterfaceAL2.h>
#include <EmbeddedBeamContact.h>
#include <EmbeddedBeam.h>
#include <EmbeddedEPBeamInterface.h>
#include <EmbeddedBeamToe.h>
#include <EmbeddedBeamToeP.h>

const double PI = 3.14159265359;

int invIsoMapping(Vector nodesX, Vector nodesY, Vector nodesZ, double Px, double Py, double Pz, double & xi, double & eta, double & zeta, bool & inBounds);

#endif