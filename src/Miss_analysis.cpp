#include "TFile.h"
#include "TProfile2D.h"
#include "LACTRecEvent.h"
#include "LACT_RunPara.h"
#include "TMath.h"
#include "TTree.h"




int main( int argc, char** argv)
{
    std::string out_file = "dst.root";
    LACT_RUNPARA* lactrun = new LACT_RUNPARA();
    lactrun->ProcessCommandLine(argc, argv);

    TFile* out_root = TFile::Open(lactrun->out_file.c_str(), "recreate");

    TProfile2D* miss = new TProfile2D("miss1", "Miss Versus Rp and log10Size", 50, 2, 7, 75, 0, 750);
    TProfile2D* miss2 = new TProfile2D("miss2", "Miss Versus Rp and MCenergy", 40, -1, 3, 75, 0, 500);
    
    for( auto input : lactrun->input_file)
    {
        TFile* input_root = TFile::Open(input.c_str(), "read");
        LACTRecEvent* lactrec = new LACTRecEvent();
        if(input_root->IsZombie())
        {
            continue;
        }
        TTree* lactrectree = (TTree*) input_root->Get("LactRectree");
        lactrectree->SetBranchAddress("LactRecEvent", &lactrec);
        for( int i = 0; i < lactrectree->GetEntries(); i++)
        {
            lactrectree->GetEntry(i);
            if( lactrec->direction_error > 0)
            {
                for( int i = 0 ; i < lactrec->GetNtel(); i++)
                {
                    miss->Fill(log10(lactrec->GetTelSize(i)), lactrec->miss[i] * TMath::RadToDeg(), lactrec->weight);
                    miss2->Fill(log10(lactrec->MCenergy), lactrec->miss[i] * TMath::RadToDeg(), lactrec->weight);
                }
            }
        }

        input_root->Close();
    }
    out_root->cd();
    miss->Write();
    miss2->Write();
    out_root->Close();

}