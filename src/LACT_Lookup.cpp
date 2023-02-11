#include "LACT_Lookup.h"
#include "LACT_TableCalculatorData.h"
LACTLookup::LACTLookup(std::string name)
{
    lookupfile = TFile::Open(name.c_str(), "recreate");

}
void LACTLookup::GetData(TTree* RecTree)
{
    LactRecTree = RecTree;
    lactrec = new LACTRecEvent();
    RecTree->SetBranchAddress("LactRecEvent", &lactrec);

}

void LACTLookup::InitLookupTableData()
{

    TableCalculatorData[MRSW] = new LACT_TableCalculatorData();
    TableCalculatorData[MRSW]->fEnergy = false;
    TableCalculatorData[MRSW]->fFillVariable = "width";

    TableCalculatorData[MRSL] = new LACT_TableCalculatorData();
    TableCalculatorData[MRSL]->fEnergy = false;
    TableCalculatorData[MRSL]->fFillVariable = "length";

    TableCalculatorData[EREC] = new LACT_TableCalculatorData();
    TableCalculatorData[EREC]->fEnergy = true;
    TableCalculatorData[EREC]->fFillVariable = "size/Energy";
    int zenithbin = 0;

}


void LACTLookup::terminate()
{
    for(auto Data: TableCalculatorData)
        Data->terminate(lookupfile);
}

void LACTLookup::Loop()
{
    if( LactRecTree)
    {
        for( int i = 0; i < LactRecTree->GetEntries(); i++)
        {
            LactRecTree->GetEntry(i);


        }
    }
}