# EmbeddedInterfaceElement
Implementation of the embedded beam in solid interface element

## Add **EmbeddedBeamInterface** element to OpenSees:
On windows, you can add EmbeddedBeamInterface.h and EmbeddedBeamInterface.cpp files to the element project in the Visual Studio solution file under UWelements, add a unique tag to the %OpenSees%/SRC/classTags.h (like "#define ELE_TAG_EmbeddedBeamInterface  501") for this element. In the TclElementCommand.cpp do the following:
* First add this line to the beginning of the file where you see a lot of extern commands (around line 130):
extern void *OPS_EmbeddedBeamInterface(void);
* And add this "else if" expression to the  lines further down where you see all these "else if" statements (around line 930):
'''c++
else if ((strcmp(argv[1], "EmbeddedBeamInterface") == 0)) {  
	void *theEle = OPS_EmbeddedBeamInterface();
	if (theEle != 0)
		theElement = (Element *)theEle;
	else {
		opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
	    return TCL_ERROR;
	}
}'''
	
This should be enough to introduce this element to OpenSees. On linux, besides every change you applied to the files described above, you need to modify %OpenSees%/SRC/element/UWelements/Makefile to include EmbeddedBeamInterface.o, and in %OpenSees%/SRC/ change the Makefile to include "$(FE)/element/UWelements/EmbeddedBeamInterface.o" in the ELE_LIBS variable. You should be able to compile OpenSees applying these modifications!

## Add generateInterfacePoints command to OpenSees:
* In the TclModelBuilder.cpp file add the following lines:
  1. Add these somwhere around line 440:
  '''
  // Added by Alborz Ghofrani - UW
		int
		TclCommand_GenerateInterfacePoints(ClientData clientData,
			Tcl_Interp *interp,
			int argc,
			TCL_Char **argv);
	'''
  2. Add these somewhere around line 660 (Tcl_CreateCommand functions):
	'''
    // Added by Alborz Ghofrani - UW
		Tcl_CreateCommand(interp, "generateInterfacePoints",
			TclCommand_GenerateInterfacePoints,
			(ClientData)NULL, NULL);
	'''
  3. Add this somewhere around line 780 (Tcl_DeleteCommand functions):
	'''
    Tcl_DeleteCommand(theInterp, "generateInterfacePoints"); // Added by Alborz Ghofrani - UW
	'''
  4. Add these to the end of the file:
	'''
    // Added by Alborz Ghofrani - UW
		extern int
		TclCommand_GenerateInterfacePoints(ClientData clientData, Tcl_Interp *interp, int argc,
			TCL_Char **argv);
	'''
* Add files Tcl_generateInterfacePoints.h and .cpp to the element project under UWelements.
	

## If OpenSees compiles without any problems, you should be able to run the test script inside Examples folder and see 24 elements being generated.
	
Alborz
