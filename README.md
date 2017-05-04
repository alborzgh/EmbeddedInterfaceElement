# EmbeddedInterfaceElement
Implementation of the embedded beam in solid interface element

In order to add this element to OpenSees, you can add these two files to the element project in the Visual Studio solution file, add a unique tag to the %OpenSees%/SRC/classTags.h (like "#define ELE_TAG_EmbeddedBeamInterface  501") for this element and you should be able to compile OpenSees in Windows.

On linux, you need to add the unique tag for this element to the "classTags.h" file, inside %OpenSees%/SRC/element/UWelements/ modify the Makefile to include EmbeddedBeamInterface.o, and in %OpenSees%/SRC/ change the Makefile to include "$(FE)/element/UWelements/EmbeddedBeamInterface.o" in the ELE_LIBS variable. You should be able to compile OpenSees applying these modifications!

Alborz
