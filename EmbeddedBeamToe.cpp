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
// Description: This file contains the class definition for EmbeddedBeamToe.

#include <EmbeddedBeamToe.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CrdTransf.h>
#include <elementAPI.h>
#include <cmath>
#include <NodeIter.h>

static int num_EmbeddedBeamToe = 0;
static const double m_Pi = 3.14159265359;

void *
OPS_EmbeddedBeamToe(void)
{
    // TODO: The OPS_EmbedBeam() needs to be completed

    if (num_EmbeddedBeamToe == 0) {
        num_EmbeddedBeamToe++;
        opserr << "EmbeddedBeamToe element - Written: A.Ghofrani, D.Turello, P.Arduino, U.Washington\n";
    }

    Element *theElement = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "Want: EmbeddedBeamToe tag? \n";
        return 0;
    }

    int iData[1];
    int eleTag = 0;
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid integer data: element EmbeddedBeamToe" << endln;
        return 0;
    }

    eleTag = iData[0];

    theElement = new EmbeddedBeamToe(iData[0]);

    if (theElement == 0) {
        opserr << "WARNING could not create element of type EmbeddedBeamToe\n";
        return 0;
    }

    return theElement;
}


EmbeddedBeamToe::EmbeddedBeamToe(int tag) : Element(tag, ELE_TAG_EmbeddedBeamToe)
{

}

EmbeddedBeamToe::EmbeddedBeamToe(int tag, int beamTag, std::vector <int> solidTag, int crdTransfTag,
    std::vector <double>  beamRho, std::vector <double>  beamTheta, std::vector <double>  solidXi, std::vector <double>  solidEta,
    std::vector <double>  solidZeta, std::vector <double> radius, double beam_radius, double area) : Element(tag, ELE_TAG_EmbeddedBeamToe),
    m_beam_radius(beam_radius), m_area(area), theBeamTag(beamTag),
    m_Ba_rot_n(3), m_Bb_rot_n(3),
    mQa(3, 3), mQb(3, 3), mQc(3, 3), mc1(3),
    mBphi(3, 6), mBu(3, 6), mHf(3, 6), m_Ns(8)
{
    // get domain to have access to element tags and their nodes (this cannot be done in setDomain because number of nodes and dof's need to be known.)
    #ifdef _PARALLEL_PROCESSING
    #include <PartitionedDomain.h>
    extern PartitionedDomain theDomain;
    #else
    extern Domain theDomain;
    #endif

    // initialize information for solids in contact
    std::set <int> uniqueNodeTags;
    m_numEmbeddedPoints = solidTag.size();
    theSolidTag         = new int[m_numEmbeddedPoints];
    solidNodeTags       = new int[8 * m_numEmbeddedPoints];
    m_beam_rho          = m_beam_r = m_beam_theta = m_solid_xi = m_solid_eta = m_solid_zeta = Vector(m_numEmbeddedPoints);

    Element *theElement;
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        theSolidTag[ii]     = solidTag[ii];
        m_solid_xi(ii)      = solidXi[ii];
        m_solid_eta(ii)     = solidEta[ii];
        m_solid_zeta(ii)    = solidZeta[ii];
        m_beam_rho(ii)      = beamRho[ii];
        m_beam_theta(ii)    = beamTheta[ii];
        m_beam_r(ii)        = radius[ii];

        theElement = theDomain.getElement(solidTag[ii]);

        // opserr << "Point " << ii + 1 << " : element " << solidTag[ii] << " at (" << solidXi[ii] << "," << solidEta[ii] << "," << solidZeta[ii] << ") , beam: " << beamTag << " at (" << beamRho[ii] << "," << beamTheta[ii] << ")" << endln;

        for (int jj = 0; jj < 8; jj++)
        {
            uniqueNodeTags.insert(theElement->getNodePtrs()[jj]->getTag());
            solidNodeTags[ii * 8 + jj] = theElement->getNodePtrs()[jj]->getTag();
        }
    
    }

    // update the number of connected nodes and number of dof's and update the connectivity information
    m_numSolidNodes = (int)uniqueNodeTags.size();
    EBT_numNodes = m_numSolidNodes;
    EBT_numDOF   = m_numSolidNodes * 3;

    externalNodes = ID(EBT_numNodes);
    theNodes = new Node*[EBT_numNodes];

    int count = 0;
    for (std::set <int>::iterator it = uniqueNodeTags.begin(); it != uniqueNodeTags.end(); ++it)
    {
        m_nodeMap[*it] = count;
        externalNodes(count) = *it;
        count++;
    }

    // get the beam element from the domain to get its nodes
    theBeam = theDomain.getElement(beamTag);

    // get the sp-constraints to apply the displacement to the beam nodes
    theConstraints = new SP_Constraint*[6];

    // get memory for internal variables
    m_InterfaceForces    = Vector(EBT_numDOF);
    m_InterfaceStiffness = Matrix(EBT_numDOF, EBT_numDOF);
    mA     = Matrix(3*m_numSolidNodes, 6);
    mB     = Matrix(6, 6);
    mB_inv = Matrix(6, 6);
    mKbb   = Matrix(6, 6);

    // get the coordinate transformation object
    crdTransf = OPS_GetCrdTransf(crdTransfTag)->getCopy3d();

}

EmbeddedBeamToe::EmbeddedBeamToe()
    : Element(0, ELE_TAG_EmbeddedBeamToe)
{

}

EmbeddedBeamToe::~EmbeddedBeamToe()
{
    // TODO: There are some dybnamically allocated objects that need to be removed
}

int
EmbeddedBeamToe::getNumExternalNodes(void) const
{
    return EBT_numNodes;
}

const ID&
EmbeddedBeamToe::getExternalNodes(void)
{
    return externalNodes;
}

Node **
EmbeddedBeamToe::getNodePtrs(void)
{
    return theNodes;
}

int
EmbeddedBeamToe::getNumDOF(void)
{
    return EBT_numDOF;
}

int
EmbeddedBeamToe::revertToLastCommit(void)
{
    // TODO: Check if something needs to be done
    return 0;
}

int
EmbeddedBeamToe::revertToStart(void)
{
    // TODO: Check if something needs to be done
    return 0;
}


const Matrix&
EmbeddedBeamToe::getTangentStiff(void)
{
    return m_InterfaceStiffness;
}

const Matrix&
EmbeddedBeamToe::getInitialStiff(void)
{
    return this->getTangentStiff();
}

const Vector&
EmbeddedBeamToe::getResistingForce(void)
{
    m_InterfaceForces.Zero();

    Vector sDisp(3 * m_numSolidNodes);
    Vector bForces(6);

    // get the solid node displacements
    for (int ii = 0; ii < m_numSolidNodes; ii++)
    {
        sDisp(3 * ii) = theNodes[ii]->getTrialDisp()(0);
        sDisp(3 * ii + 1) = theNodes[ii]->getTrialDisp()(1);
        sDisp(3 * ii + 2) = theNodes[ii]->getTrialDisp()(2);
    }

    // get the beam external forces
    for (int ii = 0; ii < 6; ii++)
        bForces(ii) = theBeam->getResistingForce()(ii);

    m_InterfaceForces = (m_InterfaceStiffness * sDisp) - (mA * mB_inv * bForces);

    return m_InterfaceForces;
}

int
EmbeddedBeamToe::sendSelf(int commitTag, Channel &theChannel)
{
    // TODO: needs to be completed
    return 0;
}

int
EmbeddedBeamToe::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
    &thEBTroker)
{
    // TODO: needs to be completed
    return 0;
}

int
EmbeddedBeamToe::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return 0;
}

void
EmbeddedBeamToe::Print(OPS_Stream &s, int flag)
{
    // TODO: needs to be completed
    return;
}

Response*
EmbeddedBeamToe::setResponse(const char **argv, int argc,
    OPS_Stream &s)
{
    // TODO: needs to be updated
    if (strcmp(argv[0], "locationBeam") == 0 || strcmp(argv[0], "locBeam") == 0) {
        return new ElementResponse(this, 1, Vector(3));

    }
    else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "disp") == 0) {
        return new ElementResponse(this, 2, Vector(3*m_numEmbeddedPoints));

    }
    else if (strcmp(argv[0], "beamCL") == 0 || strcmp(argv[0], "beamCenterLine") == 0) {
        return new ElementResponse(this, 3, Vector(3));

    }
    else if (strcmp(argv[0], "c1") == 0 || strcmp(argv[0], "tangent") == 0) {
        return new ElementResponse(this, 4, Vector(3));

    }
    else if (strcmp(argv[0], "c2") == 0 || strcmp(argv[0], "perp2") == 0) {
        return new ElementResponse(this, 5, Vector(3));

    }
    else if (strcmp(argv[0], "c3") == 0 || strcmp(argv[0], "perp3") == 0) {
        return new ElementResponse(this, 6, Vector(3));

    }
    else if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "globalForce") == 0) {
        return new ElementResponse(this, 7, Vector(3));

    }
    else {
        opserr << "EmbeddedBeamToe Recorder, " << argv[0] << "is an unknown recorder request"
            << "  Element tag : " << this->getTag() << endln;
        return 0;
    }
}

int
EmbeddedBeamToe::getResponse(int responseID, Information &eleInformation)
{
    // TODO: needs to be updated
    if (responseID == 1) { // location

        return eleInformation.setVector(Vector(3));

    }
    else if (responseID == 2) { // location

        return eleInformation.setVector(GetInteractionPtDisp());
    }
    else if (responseID == 3) { // centerline

        return eleInformation.setVector(Vector(3));
    }
    else if (responseID == 4) { // c1
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 0);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 5) { // c2
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 1);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 6) { // c3
        Vector temp(3);
        for (int ii = 0; ii < 3; ii++)
            temp(ii) = mQc(ii, 2);
        return eleInformation.setVector(temp);
    }
    else if (responseID == 7) { // contact force
        Vector temp(3);
        return eleInformation.setVector(temp);
    }
    else {
        opserr << "EmbeddedBeamToe, tag = " << this->getTag()
            << " -- unknown request" << endln;
        return -1;
    }
}

int
EmbeddedBeamToe::setParameter(const char **argv, int argc, Parameter &param)
{
    return 0;
}

int
EmbeddedBeamToe::updateParameter(int parameterID, Information &info)
{
    return 0;
}

void
EmbeddedBeamToe::setDomain(Domain *theDomain)
{
    // get pointers to the element nodes
    for (int ii = 0; ii < m_numSolidNodes; ii++)
    {
        theNodes[ii] = theDomain->getNode(externalNodes(ii));
        if (theNodes[ii] == 0)
        {
            opserr << "Could not find node " << externalNodes(ii) << "." << endln;
            return;
        }
    }

    // set the sp constraints
    for (int ii = 0; ii < 1; ii++)
        for (int jj = 0; jj < 6; jj++)
        {
            // check if an existing SP_COostraint exists for that dof at the node
            bool found = false;
            SP_ConstraintIter &theExistingSPs = theDomain->getSPs();
            SP_Constraint *theExistingSP = 0;
            while ((found == false) && ((theExistingSP = theExistingSPs()) != 0))
            {
                int spNodeTag = theExistingSP->getNodeTag();
                int dof = theExistingSP->getDOF_Number();
                if (theBeam->getNodePtrs()[ii]->getTag() == spNodeTag && jj == dof)
                {
                    // set the constraint pointer
                    found = true;
                    theConstraints[ii * 6 + jj] = theExistingSP;
                }
            }
            if (!found)
            {
                // create a new sp constraint and add it to the domain
                theConstraints[ii * 6 + jj] = new SP_Constraint(theBeam->getNodePtrs()[ii]->getTag(), jj, 1.0, false);
                theDomain->addSP_Constraint(theConstraints[ii * 6 + jj]);
            }
        }
    
    // initialize the transformation
    if (crdTransf->initialize(theBeam->getNodePtrs()[0], theBeam->getNodePtrs()[1]))
    {
        opserr << "EmbeddedBeamToe::setDomain(): Error initializing coordinate transformation";
        return;
    }
    
    // check if the beam has zero length
    m_beam_length = crdTransf->getInitialLength();
    if (m_beam_length < 1.0e-12) {
        opserr << "FATAL ERROR EmbeddedBeamToe (tag: " << this->getTag() << ") : "
            << "Beam element has zero length." << endln;
        return;
    }
    
    // update local coordinate systems
    Vector initXAxis(3);
    Vector initYAxis(3);
    Vector initZAxis(3);
    crdTransf->getLocalAxes(initXAxis, initYAxis, initZAxis);
    // fill mQa
    for (int i = 0; i < 3; i++) {
        mQa(i, 0) = initXAxis(i);
        mQa(i, 1) = initYAxis(i);
        mQa(i, 2) = initZAxis(i);
    }
    // set mQb = mQa : beam column element requires zero initial twist
    // if mQa = mQb -> mchi = 0
    mQc = mQb = mQa;
    mchi = 0;

    Vector c1(3), c2(3), c3(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c1(ii) = mQc(ii, 0);
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }
    Setc1(c1);
    
    // calculate A and B
    mA.Zero();
    mB.Zero();
    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));
        mHf.Zero();
        ComputeHf(mHf, m_beam_r(ii), m_beam_theta(ii));


        Element * theElement = theDomain->getElement(theSolidTag[ii]);

        for (int jj = 0; jj < 8; jj++)
        {
            int nodeInA = m_nodeMap[theElement->getNodePtrs()[jj]->getTag()];

            // opserr << "Element " << theSolidTag[ii] << " - Node " << theElement->getNodePtrs()[jj]->getTag() << " which is " << nodeInA << " in the local Mat." << endln;

            for (int kk = 0; kk < 6; kk++)
            {
                mA(3 * nodeInA, kk)     += m_Ns(jj) * mHf(0, kk);
                mA(3 * nodeInA + 1, kk) += m_Ns(jj) * mHf(1, kk);
                mA(3 * nodeInA + 2, kk) += m_Ns(jj) * mHf(2, kk);
            }
        }

        // update the matrices that define kinematics of the points on the beam surface
        Matrix Hb(3, 6);

        ComputeBphiAndBu(mBphi, mBu);
        
        Hb = mBu - (m_beam_r(ii)*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;

        Matrix HbT = Transpose(3, 6, Hb);
        mB += HbT * mHf;
    }

    // update B inverse
    mB.Invert(mB_inv);

    // get the beam tangent matrix
    for (int ii = 0; ii < 6; ii++) 
        for (int jj = 0; jj < 6; jj++) 
            mKbb(ii,jj) = theBeam->getTangentStiff()(ii,jj);

    m_InterfaceStiffness.Zero();

    Matrix ABinv = mA * mB_inv;
    Matrix ABinvT(6, 3 * m_numSolidNodes);
    ABinvT = Transpose(3 * m_numSolidNodes, 6, ABinv);

    // update the element tangent stiffness matrix
    m_InterfaceStiffness = ABinv * mKbb * ABinvT;

    this->DomainComponent::setDomain(theDomain);
    return;
}

int
EmbeddedBeamToe::update(void)
{
    // TODO: the coordinate systems and kinematics need to be updated for larger deformations

    return 0;
}

int
EmbeddedBeamToe::commitState(void)
{
    int retVal = 0;

    // use the sp constraints to push the beam displacements to the beam nodes
    extern ConstraintHandler *theHandler;  // need to get access to the constraint handler

    // calculate beam displacements
    Vector bDisp(6);
    Vector sDisp(3 * m_numSolidNodes);
    for (int ii = 0; ii < m_numSolidNodes; ii++)
    {
        sDisp(3 * ii) = theNodes[ii]->getTrialDisp()(0);
        sDisp(3 * ii + 1) = theNodes[ii]->getTrialDisp()(1);
        sDisp(3 * ii + 2) = theNodes[ii]->getTrialDisp()(2);
    }

    Matrix ABinv = mA * mB_inv;
    Matrix ABinvT(6, 3 * m_numSolidNodes);
    ABinvT = Transpose(3 * m_numSolidNodes, 6, ABinv);

    bDisp = ABinvT * sDisp;

    // update the sp constraints
    for (int ii = 0; ii < 6; ii++)
        theConstraints[ii]->applyConstraint(bDisp(ii));

    // need to let the constraint handler know about the update
    if (!(theHandler == NULL))
        theHandler->applyLoad();

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
        opserr << "EmbeddedBeamToe::commitState() - failed in base class";
    }

    return retVal;
}

int EmbeddedBeamToe::updateShapeFuncs(double xi, double eta, double zeta, double rho)
{
    if ((xi < -1.0) || (xi > 1.0) || (eta < -1.0) || (eta > 1.0) || (zeta < -1.0) || (zeta > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    if ((rho < -1.0) || (rho > 1.0))
    {
        opserr << "Error in shape function." << endln;
        return -1;
    }

    double rho2 = rho * rho;
    double rho3 = rho * rho2;

    m_Ns(0) = -0.125 * (xi - 1) * (eta - 1) * (zeta - 1);
    m_Ns(1) =  0.125 * (xi + 1) * (eta - 1) * (zeta - 1);
    m_Ns(2) = -0.125 * (xi + 1) * (eta + 1) * (zeta - 1);
    m_Ns(3) =  0.125 * (xi - 1) * (eta + 1) * (zeta - 1);
    m_Ns(4) =  0.125 * (xi - 1) * (eta - 1) * (zeta + 1);
    m_Ns(5) = -0.125 * (xi + 1) * (eta - 1) * (zeta + 1);
    m_Ns(6) =  0.125 * (xi + 1) * (eta + 1) * (zeta + 1);
    m_Ns(7) = -0.125 * (xi - 1) * (eta + 1) * (zeta + 1);

    m_Hb1 = 0.125 * (4.0 - 6.0 * rho + 2.0 * rho3);
    m_Hb3 = 0.125 * (4.0 + 6.0 * rho - 2.0 * rho3);
    m_Hb2 = 0.125 * m_beam_length * (1.0 - rho - rho2 + rho3);
    m_Hb4 = 0.125 * m_beam_length * (-1.0 - rho + rho2 + rho3);

    m_Nb1 = 0.5 * (1 - rho);
    m_Nb2 = 0.5 * (1 + rho);

    m_dH1 = 0.75 * (-1.0 + rho2);
    m_dH3 = 0.75 * (1.0 - rho2);
    m_dH2 = 0.125 * m_beam_length * (-1.0 - 2.0 * rho + 3.0 * rho2);
    m_dH4 = 0.125 * m_beam_length * (-1.0 + 2.0 * rho + 3.0 * rho2);

    return 0;
}

Vector
EmbeddedBeamToe::CrossProduct(const Vector &V1, const Vector &V2)
{
    Vector V3(3);

    V3(0) = V1(1)*V2(2) - V1(2)*V2(1);
    V3(1) = V1(2)*V2(0) - V1(0)*V2(2);
    V3(2) = V1(0)*V2(1) - V1(1)*V2(0);

    return V3;
}

Vector
EmbeddedBeamToe::Geta1(void) 
{
    Vector a1(3);
    int i;

    for (i = 0; i<3; i++) {
        a1(i) = mQa(i, 0);
    }

    return a1;
}

Vector
EmbeddedBeamToe::Getb1(void) 
{
    Vector b1(3);
    int i;

    for (i = 0; i<3; i++) {
        b1(i) = mQb(i, 0);
    }

    return b1;
}

void
EmbeddedBeamToe::Setc1(Vector c1_vec) 
{
    mc1 = c1_vec;

    return;
}

Vector
EmbeddedBeamToe::Getc1(void) {
    return mc1;
}

Vector EmbeddedBeamToe::GetInteractionPtDisp()
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector c2(3), c3(3);
    Vector bDisp(6);
    Matrix Hb(3, 6);
    Vector ptDisp(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < 6; ii++)
        bDisp(ii)     = theBeam->getNodePtrs()[0]->getTrialDisp()(ii);


    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeBphiAndBu(mBphi, mBu);

        Hb = mBu - (m_beam_r(ii)*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        ptDisp = Hb * bDisp;

        for (int jj = 0; jj < 3; jj++)
            res(3 * ii + jj) = ptDisp(jj);
    }

    return res;
}

Vector EmbeddedBeamToe::GetInteractionPtForce()
{
    Vector res(3 * m_numEmbeddedPoints);
    Vector c2(3), c3(3);
    Vector bDisp(6);
    Matrix Hb(3, 6);
    Vector ptDisp(3);

    // update local coordinate system
    for (int ii = 0; ii < 3; ii++)
    {
        c2(ii) = mQc(ii, 1);
        c3(ii) = mQc(ii, 2);
    }

    for (int ii = 0; ii < 6; ii++)
        bDisp(ii) = theBeam->getNodePtrs()[0]->getTrialDisp()(ii);


    for (int ii = 0; ii < m_numEmbeddedPoints; ii++)
    {
        // update the interpolation functions for displacements and interaction forces
        updateShapeFuncs(m_solid_xi(ii), m_solid_eta(ii), m_solid_zeta(ii), m_beam_rho(ii));

        // update the matrices that define kinematics of the points on the beam surface
        ComputeBphiAndBu(mBphi, mBu);

        Hb = mBu - (m_beam_r(ii)*(cos(m_beam_theta(ii))*ComputeSkew(c2) + sin(m_beam_theta(ii))*ComputeSkew(c3))) * mBphi;
        ptDisp = Hb * bDisp;

        for (int jj = 0; jj < 3; jj++)
            res(3 * ii + jj) = ptDisp(jj);
    }

    return res;
}


Matrix
EmbeddedBeamToe::Transpose(int dim1, int dim2, const Matrix &M)
{
    // copied from transpose function in Brick.cpp

    Matrix Mtran(dim2, dim1);

    for (int i = 0; i < dim1; i++)
        for (int j = 0; j < dim2; j++)
            Mtran(j, i) = M(i, j);

    return Mtran;
}


Matrix
EmbeddedBeamToe::ComputeSkew(Vector th)
{
    Matrix skew_th(3, 3);

    skew_th(0, 0) = 0.0;
    skew_th(0, 1) = -th(2);
    skew_th(0, 2) = th(1);
    skew_th(1, 0) = th(2);
    skew_th(1, 1) = 0.0;
    skew_th(1, 2) = -th(0);
    skew_th(2, 0) = -th(1);
    skew_th(2, 1) = th(0);
    skew_th(2, 2) = 0.0;

    return skew_th;
}

void
EmbeddedBeamToe::ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu)
{
    int i, j;
    Matrix dummy1(3, 3);
    Matrix dummy2(3, 3);
    Matrix dummy3(3, 3);
    Matrix dummy4(3, 3);
    double L = m_beam_length / 2.0;

    Bphi.Zero();
    Bu.Zero();
    
    // Compute Bphi(0:2, 3:5)
    dummy1.Zero();
    dummy2.Zero();
    dummy3.Zero();
    dummy4.Zero();

    // dummy1 = N1 * Qc*(E1 dyadic E1)
    dummy1(0, 0) = m_Nb1*mQc(0, 0);
    dummy1(1, 0) = m_Nb1*mQc(1, 0);
    dummy1(2, 0) = m_Nb1*mQc(2, 0);
    // dummy1 += dH2 * Qc*P1
    dummy1(0, 1) = m_dH2*mQc(0, 1)/L;  // dH2 * mQc(0:2,1:2)
    dummy1(1, 1) = m_dH2*mQc(1, 1)/L;
    dummy1(2, 1) = m_dH2*mQc(2, 1)/L;
    dummy1(0, 2) = m_dH2*mQc(0, 2)/L;
    dummy1(1, 2) = m_dH2*mQc(1, 2)/L;
    dummy1(2, 2) = m_dH2*mQc(2, 2)/L;
    // dummy2 = Qa^T
    dummy2 = Transpose(3, 3, mQa);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ dH2 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bphi(i, 3 + j) = dummy3(i, j);

    // Reuse parts of dummy1 and dummy2 to calculate Bu(0:2,0:2)
    // dummy1 += H1 * Qc*P1
    dummy1(0, 1) = m_Hb1*mQc(0, 1);  // H1 * mQc(0:2,1:2)
    dummy1(1, 1) = m_Hb1*mQc(1, 1);
    dummy1(2, 1) = m_Hb1*mQc(2, 1);
    dummy1(0, 2) = m_Hb1*mQc(0, 2);
    dummy1(1, 2) = m_Hb1*mQc(1, 2);
    dummy1(2, 2) = m_Hb1*mQc(2, 2);
    // dummy3 = Qc * (N1*(E1 dyadic E1)+ H1 * P1) * Qa^T
    dummy3 = dummy1*dummy2;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bu(i, j) = dummy3(i, j);

    // Reuse dummy2 and Compute Bphi(0:2, 0:2) and Bu(0:2, 3:5)
    dummy1.Zero();
    dummy3.Zero();
    // dummy1 = Qc*E^R*P1 (E^R is the skew symmetric meatrix for E1 cross product => E1 x a = [E^R].a)
    dummy1(0, 1) =  mQc(0, 2);
    dummy1(0, 2) = -mQc(0, 1);
    dummy1(1, 1) =  mQc(1, 2);
    dummy1(1, 2) = -mQc(1, 1);
    dummy1(2, 1) =  mQc(2, 2);
    dummy1(2, 2) = -mQc(2, 1);
    // dummy3 = Qc*E^R*P1*Qa^T
    dummy3 = dummy1*dummy2;
    // Compute Bphi(0:2, 0:2)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bphi(i, j) = m_dH1 / L * dummy3(i, j);
    // Compute Bu(0:2, 3:5)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            Bu(i, 3 + j) = -m_Hb2 * dummy3(i, j);

    return;
}

void EmbeddedBeamToe::ComputeHf(Matrix & Hf, double radius, double theta)
{
    Hf.Zero();
    double piR2 = m_Pi * m_beam_radius * m_beam_radius;
    double oneOverPiR2 = 1.0 / piR2 * m_area;
    double oneOverPiR4 = m_Pi * oneOverPiR2 / piR2;
    for (int ii = 0; ii < 3; ii++)
    {
        for (int jj = 0; jj < 3; jj++)
            Hf(ii, jj)     = oneOverPiR2 * m_Nb1 * mQa(jj, ii);

        Hf(0, ii + 3)  =  4.0 * oneOverPiR4 * m_Nb1 * radius * (mQa(ii, 1) * sin(theta) - mQa(ii, 2) * cos(theta));
        Hf(1, ii + 3)  = -2.0 * oneOverPiR4 * mQa(ii, 0) * m_Nb1 * radius * sin(theta);
        Hf(2, ii + 3)  =  2.0 * oneOverPiR4 * mQa(ii, 0) * m_Nb1 * radius * cos(theta);
    }
    
    Hf = mQc * Hf;
    return;
}

void
EmbeddedBeamToe::UpdateTransforms(void)
{
    Vector temp_a(6);           // trial disp/rot vector at a(total disp/rot)
    Vector temp_b(6);           // trial disp/rot vector at a(total disp/rot) 
    Vector rot_a(3);            // incr. rot vector at a (from n->n+1)
    Vector rot_b(3);            // incr. rot vector at a (from n->n+1)
    Matrix Omega(3, 3);         // Matrix used for Exponential Map

                                               // Recalculate incremental rotations from step n to n+1
    temp_a = theNodes[m_numSolidNodes]->getTrialDisp();
    temp_b = theNodes[m_numSolidNodes + 1]->getTrialDisp();
    
    for (int ii = 0; ii < 3; ii++) {
        rot_a(ii)  = temp_a(ii + 3) - m_Ba_rot_n(ii);
        rot_b(ii)  = temp_b(ii + 3) - m_Bb_rot_n(ii);
    }

    // Perform exponential update of Qa
    //   calculate exponential map of current incremental rotations
    Omega = ExponentialMap(rot_a);
    //   calculate new Qa
    mQa = Omega*mQa;    

    // Perform exponential update of Qb
    //   calculate exponential map of current incremental rotations
    Omega = ExponentialMap(rot_b);
    //   calculate new Qb
    mQb = Omega*mQb;


    return;
}

void
EmbeddedBeamToe::ComputeQc()
{
    Vector c1(3);         // tangent vector at projection point, c
    Vector a1(3);         // tangent vector at a
    Vector b1(3);         // tangent vector at b
    Vector temp(3);       // dummy vector for use in calcs
    Matrix Qc_df(3, 3);   // Drill free transformation matrix for c
    Matrix Qc_chi(3, 3);  // Twist transformation matrix for c
    Matrix Qb_df(3, 3);   // Drill free transf. matrix from a to b

                                                // Fill tangent vectors
    a1 = Geta1();
    b1 = Getb1();
    c1 = Getc1();
    temp.Zero();

    // Calculate the drill free transformation from a to c, Qc_df
    temp = CrossProduct(a1, c1);
    Qc_df = ExponentialMap(temp);

    // Calculate the drill free transformation from a to b, Qb_df
    // for determination of twist angle mchi
    temp = CrossProduct(a1, b1);
    Qb_df = ExponentialMap(temp);
    Qb_df = Qb_df * mQa;

    // mchi = arcsin( b3 dot b2_df) = arcsin(mQb(:,2) dot Qb_df(:,1))
    // WATCH SIGN!!!!
    mchi = mQb(0, 2)*Qb_df(0, 1) + mQb(1, 2)*Qb_df(1, 1) + mQb(2, 2)*Qb_df(2, 1);
    mchi = -asin(mchi);

    // Calculate twist transformation from a to c, Qc_df
    // based upon linear scaling of twist angle: mxi * mchi * c1
    temp = m_Nb2*mchi*c1;
    Qc_chi = ExponentialMap(temp);

    mQc = (Qc_chi * Qc_df) * mQa;

    return;
}

Matrix
EmbeddedBeamToe::ExponentialMap(Vector th)
{
    double theta = th.Norm();             // vector norm
    Matrix sk_theta(3,3);    // skew of vector
    if (theta > 1.0e-10)
        sk_theta = ComputeSkew(th / theta);
    else
        sk_theta.Zero();
    Matrix sk_theta2 = sk_theta * sk_theta; // dyadic product of vector
    Matrix Q(3, 3);           // Exonential Map Vector  

    
    Q.Zero();
    
    Matrix meye1(3, 3);
    meye1(0, 0) = 1.0;
    meye1(1, 1) = 1.0;
    meye1(2, 2) = 1.0;

    Q = meye1 + sin(theta) * sk_theta + (1 - cos(theta)) * sk_theta2;


    return Q;
}