#include "hello_maya.h"
#include "customDeformer.h"
#include <maya/MFnPlugin.h>
#include <maya/MArgDatabase.h>

// define EXPORT for exporting dll functions
#define EXPORT _declspec(dllexport)

// Maya Plugin creator function
void* helloMaya::creator()
{
	return new helloMaya;
}

// Plugin doIt function
MStatus helloMaya::doIt(const MArgList& argList)
{
	MStatus status;
	MGlobal::displayInfo("Hello World!");

	MArgDatabase argData(syntax(), argList);
	
	MString commandArgValue0 = "NoName";
	MString commandArgValue1 = "-1";
	try {
        status = argData.getCommandArgument(0, commandArgValue0);
    }
    catch (const MStatus& status) {
        commandArgValue0 = "NoNameProvided";
    }

	try {
		status = argData.getCommandArgument(1, commandArgValue1);
	}
	catch (const MStatus& status) {
		commandArgValue1 = "-9999";
	}

	MString name_text = "\"Name: " + commandArgValue0 + "\"";
	MString id_text = "\"ID: " + commandArgValue1 + "\"";

	MString command = "if (`window -exists \"myWindow1\"`)\n";
	command += "{ deleteUI -window \"myWindow1\";}\n";

	command += "window -title \"Hello Maya\" -w 200 -h 100 myWindow1;\n";

	command += "columnLayout -adjustableColumn false -rowSpacing 5;\n";

	command += "text -al center -label " + name_text + ";\n";
	command += "text -al center -label " + id_text + ";\n";

	command += "button -al center -label \"OK\" -width 50;\n";
	command += "showWindow myWindow1;";

	MString result;

	status = MGlobal::executeCommand(command, result);

	return status;
}

// Add command call syntax
MSyntax helloMaya::newSyntax() {
	MSyntax syntax;
	syntax.addArg(MSyntax::kString);
	syntax.addArg(MSyntax::kString);
	return syntax;
}

// Initialize Maya Plugin upon loading
EXPORT MStatus initializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj, "CIS660", "1.0", "Any");
	status = plugin.registerCommand("helloMaya", helloMaya::creator, helloMaya::newSyntax);
	if (!status)
		status.perror("registerCommand failed");


	status = plugin.registerNode("customDeformer", customDeformer::id, customDeformer::creator, customDeformer::initialize, MPxNode::kDeformerNode);
	if (!status)
		status.perror("register customDeformer node failed");
	return status;

}

// Cleanup Plugin upon unloading
EXPORT MStatus uninitializePlugin(MObject obj)
{

	MStatus status;
	MFnPlugin plugin(obj);
	status = plugin.deregisterCommand("helloMaya");
	if (!status)
		status.perror("deregisterCommand failed");

	status = plugin.deregisterNode(customDeformer::id);

	return status;
}