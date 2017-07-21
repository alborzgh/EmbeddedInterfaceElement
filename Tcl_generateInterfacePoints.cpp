#include <stdlib.h>
#include <string.h>
#include <vector>

#include <Element.h>
#include <OPS_Stream.h>
#include <Domain.h>
#include <Node.h>
#include <ElementIter.h>
#include <ElasticBeam3d.h>
#include <CrdTransf.h>
#include <TclModelBuilder.h>
#include <elementAPI.h>

#include "Tcl_generateInterfacePoints.h"


#ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
extern PartitionedDomain theDomain;
#else
extern Domain theDomain;
#endif

int
TclCommand_GenerateInterfacePoints(ClientData clientData, Tcl_Interp *interp, int argc,
	TCL_Char **argv)
{
	int numArgsRemaining = argc - 1;
	int curArgPos = 1;
	int res = 0;

	int beamTag;
	double radius = 1.0;
	int nP = 100;
	int nL = 20;
	int shape = 1;

	bool radiusSet = false;
	bool nPset = false;
	bool nLset = false;
	bool readingEleTag = false;
	bool shapeSet = false;
	bool eleSetDefined = false;
	bool debugFlag = false;
	std::vector <int> eleTags;

	int transfTag = 1;
	CrdTransf *theTransf = OPS_GetCrdTransf(transfTag)->getCopy3d();

	if (numArgsRemaining < 2)
	{
		opserr << "Need beam tag and geometric transformation tag." << endln;
		return -1;
	}
	else
	{
		beamTag = strtol(argv[curArgPos], NULL, 10);
		numArgsRemaining--; curArgPos++;
		transfTag = strtol(argv[curArgPos], NULL, 10);
		numArgsRemaining--; curArgPos++;
	}


	while (numArgsRemaining > 0)
	{
		if ((strcmp(argv[curArgPos], "-shape") == 0) && numArgsRemaining >= 2)
		{
			readingEleTag = false;
			numArgsRemaining--; curArgPos++;
			if (strcmp(argv[curArgPos], "circle") == 0)
			{
				shapeSet = true;
				shape = 1;
				numArgsRemaining--; curArgPos++;
			}
			else if (strcmp(argv[curArgPos], "square") == 0)
			{
				shapeSet = true;
				shape = 2;
				numArgsRemaining--; curArgPos++;
				opserr << "The procedure for square shapes not implemented yet." << endln;
			}
			else if (strcmp(argv[curArgPos], "rectangle") == 0)
			{
				shapeSet = true;
				shape = 4;
				numArgsRemaining--; curArgPos++;
				opserr << "The procedure for rectangle shapes not implemented yet." << endln;
			}
			else if (strcmp(argv[curArgPos], "hexagon") == 0)
			{
				shapeSet = true;
				shape = 6;
				numArgsRemaining--; curArgPos++;
				opserr << "The procedure for hexagon shapes not implemented yet." << endln;
			}
			else {
				numArgsRemaining--; curArgPos++;
				opserr << "Unknown shape." << endln;
			}

		}
		else if ((strcmp(argv[curArgPos], "-radius") == 0) && numArgsRemaining >= 2)
		{
			readingEleTag = false;
			radiusSet = true;
			numArgsRemaining--; curArgPos++;
			radius = strtod(argv[curArgPos], NULL);
			numArgsRemaining--; curArgPos++;
		}
		else if ((strcmp(argv[curArgPos], "-nP") == 0) && numArgsRemaining >= 2)
		{
			readingEleTag = false;
			nPset = true;
			numArgsRemaining--; curArgPos++;
			nP = strtol(argv[curArgPos], NULL, 10);
			numArgsRemaining--; curArgPos++;
		}
		else if ((strcmp(argv[curArgPos], "-nL") == 0) && numArgsRemaining >= 2)
		{
			readingEleTag = false;
			nLset = true;
			numArgsRemaining--; curArgPos++;
			nL = strtol(argv[curArgPos], NULL, 10);
			numArgsRemaining--; curArgPos++;
		}
		else if ((strcmp(argv[curArgPos], "-ele") == 0) && numArgsRemaining >= 2)
		{
			readingEleTag = true;
			eleSetDefined = true;
			numArgsRemaining--; curArgPos++;
			eleTags.push_back(strtol(argv[curArgPos], NULL, 10));
			numArgsRemaining--; curArgPos++;
		}
		else if ((strcmp(argv[curArgPos], "-eleRange") == 0) && numArgsRemaining >= 3)
		{
			readingEleTag = false;
			eleSetDefined = true;
			numArgsRemaining--; curArgPos++;
			int startEle = strtol(argv[curArgPos], NULL, 10);
			numArgsRemaining--; curArgPos++;
			int endEle = strtol(argv[curArgPos], NULL, 10);
			numArgsRemaining--; curArgPos++;

			for (int ii = startEle; ii <= endEle; ii++)
				eleTags.push_back(ii);
		}
		else {
			if (readingEleTag)
			{
				eleTags.push_back(strtol(argv[curArgPos], NULL, 10));
				numArgsRemaining--; curArgPos++;
				continue;
			}
			else {
				opserr << "Unknown argument ..." << endln;
				res = -1;
				break;
			}
		}
	}



	Element *theElement;
	theElement = theDomain.getElement(beamTag);
	if (strcmp(theElement->getClassType(), "ElasticBeam3d") != 0)
	{
		opserr << "Beam element of type " << theElement->getClassType() << "not supported." << endln;
		return -1;
	}

	


	
	ElementIter& theElements = theDomain.getElements();
	int maxTag = 0;
	while ((theElement = theElements()) != 0)
	{
		if (!eleSetDefined)
			if (strcmp(theElement->getClassType(), "Brick") == 0)
				eleTags.push_back(theElement->getTag());
		if (theElement->getTag() > maxTag)
			maxTag = theElement->getTag();
	}
	int startTag = maxTag + 1;

	


	Node** solidNodesPtr;
	Node** beamNodePtr;

	beamNodePtr = theDomain.getElement(beamTag)->getNodePtrs();
	double L1x = beamNodePtr[0]->getCrds()(0);
	double L1y = beamNodePtr[0]->getCrds()(1);
	double L1z = beamNodePtr[0]->getCrds()(2);
	double L2x = beamNodePtr[1]->getCrds()(0);
	double L2y = beamNodePtr[1]->getCrds()(1);
	double L2z = beamNodePtr[1]->getCrds()(2);


	// // calculate the rotation information(the centerline is calculated using a
	// // rotation operation)
	// double ux = L2x - L1x;
	// double uy = L2y - L1y;
	// double uz = L2z - L1z;
	// double L = sqrt(pow(ux, 2.0) + pow(uy, 2) + pow(uz, 2));
	// ux /= L;
	// uy /= L;
	// uz /= L;
	// 
	// double sintheta, costheta, tantheta;
	// if (ux == 0)
	// {
	// 	sintheta = 0;
	// 	costheta = 1;
	// }
	// else
	// {
	// 	tantheta = uy / ux;
	// 	costheta = 1 / sqrt(1.0 + pow(tantheta, 2.0));
	// 	sintheta = tantheta * costheta;
	// }
	// 
	// double singamma = uz;
	// double cosgamma = sqrt(1 - pow(singamma, 2));
	// 
	// // create the pile points
	// Vector t(nP);
	// for (int ii = 0; ii < nP; ii++)
	// {
	// 	t(ii) = 0 + (double)ii / nP * 2.0 * PI;
	// }
	// 
	// Vector cX(nL*nP), cY(nL*nP), cZ(nL*nP);
	// for (int ii = 0; ii < nL; ii++)
	// 	for (int jj = 0; jj < nP; jj++)
	// 	{
	// 		cX(ii * nP + jj) = radius * cos(t(jj));
	// 		cY(ii * nP + jj) = radius * sin(t(jj));
	// 		cZ(ii * nP + jj) = 0 + (double)ii / (nL - 1) * L;
	// 	}
	// 
	// Vector cXr = cX*costheta*singamma - cY*sintheta + cZ*costheta*cosgamma + L1x;
	// Vector cYr = cX*sintheta*singamma + cY*costheta + cZ*sintheta*cosgamma + L1y;
	// Vector cZr = cZ*singamma - cX*cosgamma + L1z;


	// use the geometric transformation object to generate beam surface points
	if (theTransf->initialize(beamNodePtr[0], beamNodePtr[1]))
	{
		opserr << "generateInterfacePoints: Error initializing coordinate transformation";
		return -1;
	}
	Vector loc_x(3), loc_y(3), loc_z(3);
	theTransf->getLocalAxes(loc_x, loc_y, loc_z);

	// create the pile points
	Vector t(nP);
	for (int ii = 0; ii < nP; ii++)
	{
		t(ii) = 0 + (double)ii / nP * 2.0 * PI;
	}

	double L = sqrt(pow(L2x - L1x, 2.0) + pow(L2y - L1y, 2) + pow(L2z - L1z, 2));

	Vector cX(nL*nP), cY(nL*nP), cZ(nL*nP);
	Vector loc_rho(nL*nP), loc_theta(nL*nP);
	for (int ii = 0; ii < nL; ii++)
		for (int jj = 0; jj < nP; jj++)
		{
			cX(ii * nP + jj) = 0 + (double)ii / (nL - 1) * L;
			cY(ii * nP + jj) = radius * cos(t(jj));
			cZ(ii * nP + jj) = radius * sin(t(jj));

			loc_rho(ii * nP + jj) = cX(ii * nP + jj) / L;
			loc_theta(ii * nP + jj) = t(jj);
		}

	Vector cXr = cX*loc_x(0) + cY*loc_y(0) + cZ*loc_z(0) + L1x;
	Vector cYr = cX*loc_x(1) + cY*loc_y(1) + cZ*loc_z(1) + L1y;
	Vector cZr = cX*loc_x(2) + cY*loc_y(2) + cZ*loc_z(2) + L1z;


	if (debugFlag)
		for (int ii = 0; ii < nL; ii++)
			for (int jj = 0; jj < nP; jj++)
			{
				opserr << "point " << ii * nP + jj + 1 << " : " << cXr(ii * nP + jj) << " " << cYr(ii * nP + jj) << " " << cZr(ii * nP + jj) << endln;
			}
	

	Vector tempX(8), tempY(8), tempZ(8);
	int contactElemCount = 0, contactPointCount = 0;
	bool contactElemFlag = false;
	double xi, eta, zeta;
	bool inBounds = false;

	for (int ii = 0; ii < eleTags.size(); ii++)
	{
		theElement = theDomain.getElement(eleTags[ii]);
		solidNodesPtr = theElement->getNodePtrs();
		if (debugFlag) 
			opserr << endln << "Solid element " << eleTags[ii] << " coordiantes:" << endln;
		for (int jj = 0; jj < 8; jj++)
		{
			if (debugFlag)
				opserr << "    Node " << jj + 1 << " coordinates = " << solidNodesPtr[jj]->getCrds();
			tempX(jj) = solidNodesPtr[jj]->getCrds()(0);
			tempY(jj) = solidNodesPtr[jj]->getCrds()(1);
			tempZ(jj) = solidNodesPtr[jj]->getCrds()(2);
		}
		for (int jj = 0; jj < nP*nL; jj++)
		{
			invIsoMapping(tempX, tempY, tempZ, cXr(jj), cYr(jj), cZr(jj), xi, eta, zeta, inBounds);
			if (inBounds)
			{
				contactElemFlag = true;
				if (debugFlag) opserr << "Beam tag : " << beamTag << ", Solid tag : " << eleTags[ii] << ", Real Coordinates = (" << cXr(jj) << "," << cYr(jj) << "," << cZr(jj) << "), Iso Coordinates = (" << xi << "," << eta << "," << zeta << ")" << endln;


				maxTag++;
				theElement = new EmbeddedBeamInterface(maxTag, beamTag, eleTags[ii], transfTag, loc_rho(jj), loc_theta(jj), xi, eta, zeta, radius);

				// if one of the above worked
				if (theElement != 0) {
					if (theDomain.addElement(theElement) == false) {
						opserr << "WARNING could not add element with tag: " << theElement->getTag() << " and of type: "
							<< theElement->getClassType() << " to the Domain\n";
						delete theElement;
						return -1;
					}
				}

				contactPointCount++;
			}
		}
		if (contactElemFlag)
			contactElemCount++;
	}

	opserr << "Number of elements in contact: " << contactElemCount << ", number of contact points: " << contactPointCount << endln;
	opserr << "EmbeddedBeamInterface Elements " << startTag << " to " << maxTag << " were created." << endln;


	if (debugFlag)
	{
		if (radiusSet)
			opserr << "Radius is set to " << radius << endln;
		else
			opserr << "No radius is set!!! Assuming 1.0" << endln;
		if (nPset)
			opserr << "nP is set to " << nP << endln;
		else
			opserr << "nP not set!!! Assuming " << nP << endln;
		if (nLset)
			opserr << "nL is set to " << nL << endln;
		else
			opserr << "nL not set!!! Assuming " << nL << endln;
		if (shapeSet)
			opserr << "Shape is set to " << shape << endln;
		else
			opserr << "Default shape (circle) is assumed." << endln;
		if (eleTags.size() > 0)
			opserr << "Here are the elements: \n     { ";
		for (int ii = 0; ii < eleTags.size(); ii++)
			opserr << eleTags[ii] << " ";
		if (eleTags.size() > 0)
			opserr << "}" << endln;
	}


	return 0;

}

int 
invIsoMapping(Vector nodesX, Vector nodesY, Vector nodesZ, double Px, double Py, double Pz, double & xi, double & eta, double & zeta, bool & inBounds)
{


	// This function returns the parameters of the isoparametric element using
	// trilinear shape functions.nodesX, nodesY and nodesZ contain the
	// actual coordinates of the element and p_x, p_y and p_z are the actual
	// coordinates of the point for which the parameters of the isoparametric
	// shape functions are needed.Obviously any assumptions regarding the
	// trilinear shape functions apply here as well.

	inBounds = true;

	double x1 = nodesX(0); double y1 = nodesY(0); double z1 = nodesZ(0);
	double x2 = nodesX(1); double y2 = nodesY(1); double z2 = nodesZ(1);
	double x3 = nodesX(2); double y3 = nodesY(2); double z3 = nodesZ(2);
	double x4 = nodesX(3); double y4 = nodesY(3); double z4 = nodesZ(3);
	double x5 = nodesX(4); double y5 = nodesY(4); double z5 = nodesZ(4);
	double x6 = nodesX(5); double y6 = nodesY(5); double z6 = nodesZ(5);
	double x7 = nodesX(6); double y7 = nodesY(6); double z7 = nodesZ(6);
	double x8 = nodesX(7); double y8 = nodesY(7); double z8 = nodesZ(7);

	double a = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 - 8 * Px;
	double a_xi = -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8;
	double a_eta = -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8;
	double a_zeta = -x1 - x2 - x3 - x4 + x5 + x6 + x7 + x8;
	double a_xe = x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8;
	double a_xz = x1 - x2 - x3 + x4 - x5 + x6 + x7 - x8;
	double a_ez = x1 + x2 - x3 - x4 - x5 - x6 + x7 + x8;
	double a_xez = -x1 + x2 - x3 + x4 + x5 - x6 + x7 - x8;

	double b = y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 - 8 * Py;
	double b_xi = -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8;
	double b_eta = -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8;
	double b_zeta = -y1 - y2 - y3 - y4 + y5 + y6 + y7 + y8;
	double b_xe = y1 - y2 + y3 - y4 + y5 - y6 + y7 - y8;
	double b_xz = y1 - y2 - y3 + y4 - y5 + y6 + y7 - y8;
	double b_ez = y1 + y2 - y3 - y4 - y5 - y6 + y7 + y8;
	double b_xez = -y1 + y2 - y3 + y4 + y5 - y6 + y7 - y8;

	double c = z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 - 8 * Pz;
	double c_xi = -z1 + z2 + z3 - z4 - z5 + z6 + z7 - z8;
	double c_eta = -z1 - z2 + z3 + z4 - z5 - z6 + z7 + z8;
	double c_zeta = -z1 - z2 - z3 - z4 + z5 + z6 + z7 + z8;
	double c_xe = z1 - z2 + z3 - z4 + z5 - z6 + z7 - z8;
	double c_xz = z1 - z2 - z3 + z4 - z5 + z6 + z7 - z8;
	double c_ez = z1 + z2 - z3 - z4 - z5 - z6 + z7 + z8;
	double c_xez = -z1 + z2 - z3 + z4 + z5 - z6 + z7 - z8;

	Matrix R(3, 3);
	R(0, 0) = a_xi; R(0, 1) = a_eta; R(0, 2) = a_zeta;
	R(1, 0) = b_xi; R(1, 1) = b_eta; R(1, 2) = b_zeta;
	R(2, 0) = c_xi; R(2, 1) = c_eta; R(2, 2) = c_zeta;

	Vector res(3);
	res(0) = -a; res(1) = -b; res(2) = -c;

	Vector A(3);
	R.Solve(res, A);


	xi = A(0);
	eta = A(1);
	zeta = A(2);

	for (int ii = 1; ii <= 20; ii++)
	{
		res(0) = a + a_xi*xi + a_eta*eta + a_zeta*zeta + a_xe*xi*eta + a_xz*xi*zeta + a_ez*eta*zeta + a_xez*xi*eta*zeta;
		res(1) = b + b_xi*xi + b_eta*eta + b_zeta*zeta + b_xe*xi*eta + b_xz*xi*zeta + b_ez*eta*zeta + b_xez*xi*eta*zeta;
		res(2) = c + c_xi*xi + c_eta*eta + c_zeta*zeta + c_xe*xi*eta + c_xz*xi*zeta + c_ez*eta*zeta + c_xez*xi*eta*zeta;

		if (res.Norm() < 1.0e-10)
		{
			if ((xi < -1) || (xi > 1) || (eta < -1) || (eta > 1) || (zeta < -1) || (zeta > 1))
				inBounds = false;

			// opserr << "Calculated in " << ii - 1 << " steps." << endln;
			return 0;
		}

		R(0, 0) = a_xi + a_xe * eta + a_xz * zeta + a_xez * eta * zeta;
		R(0, 1) = a_eta + a_xe * xi + a_ez * zeta + a_xez * xi  * zeta;
		R(0, 2) = a_zeta + a_xz * xi + a_ez * eta + a_xez * xi  * eta;
		R(1, 0) = b_xi + b_xe * eta + b_xz * zeta + b_xez * eta * zeta;
		R(1, 1) = b_eta + b_xe * xi + b_ez * zeta + b_xez * xi  * zeta;
		R(1, 2) = b_zeta + b_xz * xi + b_ez * eta + b_xez * xi  * eta;
		R(2, 0) = c_xi + c_xe * eta + c_xz * zeta + c_xez * eta * zeta;
		R(2, 1) = c_eta + c_xe * xi + c_ez * zeta + c_xez * xi  * zeta;
		R(2, 2) = c_zeta + c_xz * xi + c_ez * eta + c_xez * xi  * eta;

		R.Solve(res, A);

		xi -= A(0);
		eta -= A(1);
		zeta -= A(2);

	}

	// opserr << "Did not converge in 20 steps." << endln;
	return -1;
}
