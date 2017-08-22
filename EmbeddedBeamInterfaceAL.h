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
// Description: This file contains the class definition for EmbeddedBeamInterfaceAL.

#ifndef EmbedBeamInterfaceAL_h
#define EmbedBeamInterfaceAL_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// number of nodes per element
#define EBIAL_NUM_NODE 10
// number of dimensions
#define EBIAL_NUM_DIM  3
// degrees of freedom per element
#define EBIAL_NUM_DOF  36

class Node;
class NDMaterial;
class Response;
class CrdTransf;

class EmbeddedBeamInterfaceAL : public Element
{
  public:
	EmbeddedBeamInterfaceAL(int tag);
	EmbeddedBeamInterfaceAL(int tag, int beamTag, int solidTag, int crdTransfTag, double beamRho, double beamTheta, double solidXi, double solidEta, double solidZeta, double radius, double area);
    EmbeddedBeamInterfaceAL();
    ~EmbeddedBeamInterfaceAL();

    const char *getClassType(void) const {return "EmbeddedBeamInterfaceAL";};

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
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
	
  protected:
    
  private:
    // private attributes - a copy for each object of the class

    ID externalNodes; // Tags of quad and beam nodes

	int theSolidTag, theBeamTag;

    Node *theNodes[10];
	
	static Vector		m_InterfaceForces;	// force vector
	static Matrix		m_InterfaceStiffness;	// stiffness matrix
	static const double m_Pi;

    // iso-parametric coordinates
    double m_solid_xi;
    double m_solid_eta;
    double m_solid_zeta;
    double m_beam_rho;
    double m_beam_theta;

    // shape functions
    Vector  m_Ns;
    double  m_Hb1, m_Hb2, m_Hb3, m_Hb4, m_Nb1, m_Nb2;
    double  m_dH1, m_dH2, m_dH3, m_dH4;

    double	m_beam_radius;	 // beam Radius
    double  m_beam_length;   // beam length
	double	m_ep;		     // penalty parameter
    double  m_area;          // interface element area
    Vector  m_lambda;        // Lagrange parameter

    CrdTransf* crdTransf;  // pointer to coordinate tranformation object

    double  m_Force;

    Vector m_Ba_rot_n, m_Bb_rot_n;
    Vector m_Ba_disp_n, m_Bb_disp_n;
    Vector m_Bcl_pos, m_Bcl_pos_n;
    Vector m_B_loc, m_S_disp;
    Vector m_Ba1, m_Bb1;
    Vector m_pos;
    Vector m_lambda_n;

    int	   updateShapeFuncs();

    // copied from BeamContact3D
    double mchi;                // twist rotation from end 1 to end 2
    Matrix mQa;                 // coordinate transform for node a
    Matrix mQb;                 // coordinate transform for node b
    Matrix mQc;
    Vector mc1;                 // tangent vector at project point c
    Matrix mBphi, mBu;

    void ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu);            // method to compute Bphi and Bu, used in ComputeB and update
    void UpdateTransforms(void);         // method to update Qa, Qb
    void ComputeQc();                    // method to compute Qc from Qa and Qb

    Matrix ExponentialMap(Vector theta); // function returns the exponential map of a vector
    Matrix ComputeSkew(Vector theta);    // function returns skew matrix of given vector
    Vector CrossProduct(const Vector &V1, const Vector &V2); // cross product (does not exist in Vector Class!)
    Matrix Transpose(int dim1, int dim2, const Matrix &M);   // functions returns the tranpose of Matrix M (does not exist in Matrix Class!)

    Vector Geta1(void);                 // returns a1 = mQa(:,0)      
    Vector Getb1(void);                 // returns b1 = mQb(:,0)
    void   Setc1(Vector c1_vec);        // sets member vector c1
    Vector Getc1(void);                 // returns member vector c1
	
};

#endif

