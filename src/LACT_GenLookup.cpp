/*

            Main Program Generate Lookup Tables

*/
#include "LACT_Lookup.h"
#include "LACT_RunPara.h"
int main(int argc, char** argv)
{
    LACT_RUNPARA* lact_runpara = new LACT_RUNPARA();
    lact_runpara->ProcessCommandLine(argc, argv, "true");
    LACTLookup* lact_lookup = new LACTLookup(lact_runpara->GetLookupName());

    for( auto input : lact_runpara->input_file)
    {
        TFile* inputfile = TFile::Open(input.c_str(), "read");
        TTree* input_tree = (TTree*)inputfile->Get("JointRecTree");
        lact_lookup->GetData(input_tree);
    }
}

