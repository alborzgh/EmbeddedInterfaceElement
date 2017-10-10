/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: Alborz Ghofrani, Diego Turello, Pedro Arduino, U.Washington 
// Created: May 2017
// Description: This file contains the class definition for EmbeddedBeam.

#ifndef EmbedBeam_h
#define EmbedBeam_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <ConstraintHandler.h>

#include <vector>
#include <set>
#include <map>

// number of dimensions
#define EB_NUM_DIM  3

class Node;
class NDMaterial;
class Response;
class CrdTransf;

class EmbeddedBeam : public Element
{
public:
    EmbeddedBeam(int tag);
    EmbeddedBeam(int tag, int beamTag, std::vector <int> solidTag, int crdTransfTag,
        std::vector <double>  beamRho, std::vector <double>  beamTheta, std::vector <double>  solidXi,
        std::vector <double>  solidEta, std::vector <double>  solidZeta, double radius, double area, 
        const std::vector<int>& connectedBeamTags = std::vector<int>());
    EmbeddedBeam();
    ~EmbeddedBeam();

    const char *getClassType(void) const { return "EmbeddedBeam"; };

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);

    const Vector &getResistingForce(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
        &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc,
        OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

protected:

private:
    // private attributes - a copy for each object of the class
    int EB_numNodes, EB_numDOF;

    ID externalNodes; // Tags of beam and solid nodes

    Node **theNodes;
    Element **theConnctedBeamElems;

    int *theSolidTag;
    int *theConnectedBeamTag;
    int *solidNodeTags;
    int theBeamTag;

    Element* theBeam;
    Vector m_beamDisp;

    SP_Constraint **theConstraints;

    Vector		m_InterfaceForces;	// force vector
    Matrix		m_InterfaceStiffness;	// stiffness matrix

    std::map<int, int> m_nodeMap;

    // iso-parametric coordinates
    Vector m_solid_xi;
    Vector m_solid_eta;
    Vector m_solid_zeta;
    Vector m_beam_rho;
    Vector m_beam_theta;

    bool m_isExternalBeamConnected;
    int m_isFirstUpdate;
    int  m_connectedNode;

    int  m_numSolidNodes, m_numEmbeddedPoints, m_numConnectedBeams;

    // shape functions
    Vector  m_Ns;
    double  m_Hb1, m_Hb2, m_Hb3, m_Hb4, m_Nb1, m_Nb2;
    double  m_dH1, m_dH2, m_dH3, m_dH4;

    double	m_beam_radius;	 // beam Radius
    double  m_beam_length;   // beam length
    double  m_area;          // interface element area

    CrdTransf* crdTransf;  // pointer to coordinate tranformation object

    // copied from BeamContact3D
    double mchi;                // twist rotation from end 1 to end 2
    Matrix mQa;                 // coordinate transform for node a
    Matrix mQb;                 // coordinate transform for node b
    Matrix mQc;
    Vector mc1;                 // tangent vector at project point c
    Vector m_Ba_rot_n, m_Bb_rot_n;
    Matrix mBphi, mBu, mHf;
    Matrix mA, mB, mB_inv, mKbb;

    void ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu);            // method to compute Bphi and Bu, used in ComputeB and update
    void ComputeHf(Matrix &Hf, double theta);                   // method to compute Hf
    void UpdateTransforms(void);                                // method to update Qa, Qb
    void ComputeQc();
    int	 updateShapeFuncs(double xi, double eta, double zeta, double rho);             // method to update shape functions

    Matrix ExponentialMap(Vector theta); // function returns the exponential map of a vector
    Matrix ComputeSkew(Vector theta);    // function returns skew matrix of given vector
    Vector CrossProduct(const Vector &V1, const Vector &V2); // cross product (does not exist in Vector Class!)
    Matrix Transpose(int dim1, int dim2, const Matrix &M);   // functions returns the tranpose of Matrix M (does not exist in Matrix Class!)

    Vector Geta1(void);                 // returns a1 = mQa(:,0)      
    Vector Getb1(void);                 // returns b1 = mQb(:,0)
    void   Setc1(Vector c1_vec);        // sets member vector c1
    Vector Getc1(void);                 // returns member vector c1

    // define some functions for recording results purposes
    Vector GetInteractionPtDisp();
    Vector GetInteractionPtForce();

};

#endif

