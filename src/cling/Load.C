// Load script that links all the required libraries.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@alumni.iu.edu
// -----------------------------------------------------------------------------

void Load()
{
    TString lib_ext   = gSystem->GetSoExt();

    //----------------------------------------------------------------------
    // Core library

    TString main_dir  = gSystem->Getenv("TRIANGLETOOLS");

    // Load the main library files
    TString lib  = main_dir + "/lib/libTRIANGLETOOLS." + lib_ext;

    // Headers
    TString core  = main_dir + "/src"; 
    TString plot  = main_dir + "/src/plotting"; 
    TString cube  = main_dir + "/src/cubature"; 

    if (!gSystem->AccessPathName(lib.Data()))
    {
        Int_t lib_loaded = gSystem->Load(lib.Data());
        if (lib_loaded < 0) Fatal("Load()", "Library not loaded sucessfully!");

        gInterpreter->AddIncludePath( core.Data());
        gInterpreter->AddIncludePath( plot.Data());
        gInterpreter->AddIncludePath( cube.Data());
        gInterpreter->AddIncludePath( main_dir.Data());
    }
    else
    {
        Warning("Load()", "triangleTools library not found! Looked in: %s", lib.Data());
    }
}