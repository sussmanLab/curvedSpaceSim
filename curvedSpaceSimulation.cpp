#include "std_include.h"
#include <tclap/CmdLine.h>

using namespace std; // bad form
using namespace TCLAP;
int main(int argc, char*argv[])
    {

    //First, we set up a basic command line parser with some message and version
    CmdLine cmd("simulations in curved space!",' ',"V0.0");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programBranchSwitchArg("z","programBranchSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<string> meshSwitchArg("m","meshSwitch","filename of the mesh you want to load",false,"../exampleMeshes/torus_isotropic_remesh.off","string",cmd);
    //parse the arguments
    cmd.parse( argc, argv );
    //define variables that correspond to the command line parameters
    int programBranch = programBranchSwitchArg.getValue();
    string meshName = meshSwitchArg.getValue();

    if(!fileExists(meshName))
        {
        cout << meshName << " not found." << endl;
        return 0;
        }
    cout << "program executing branch " << programBranch << " and loading mesh from file" << meshName << endl;

    return 0;
    };
