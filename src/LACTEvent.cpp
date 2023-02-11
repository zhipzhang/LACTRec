#include "LACTEvent.h"
#include "LACT_TelData.h"
#include "LACT_Telconfig.h"
#include "TRandom2.h"
#include <iostream>
#include <vector>

TRandom2* rd2  = new TRandom2();
void LACTEvent::Init( TTree* event_tree)
{
    int nEntries = event_tree->GetEntries();
    std::cout << "There are " << nEntries << " Entries" << std::endl;
    event_tree->SetBranchAddress("runNumber", &Runnumber);
    event_tree->SetBranchAddress("eventNumber", &Eventnumber);

    event_tree->SetBranchAddress("ntel ", &ntel);
    event_tree->SetBranchAddress("MCe0", &MCenergy);
    event_tree->SetBranchAddress("MCxcore", &MCxcore);
    event_tree->SetBranchAddress("MCycore", &MCycore);
    event_tree->SetBranchAddress("MCze", &MCze);
    event_tree->SetBranchAddress("MCaz", &MCaz);
    event_tree->SetBranchAddress("MCprim", &particle_id);
    event_tree->SetBranchAddress("hmax", &hmax);
    event_tree->SetBranchAddress("xmax", &xmax);
    event_tree->SetBranchAddress("flag", &flag);
    event_tree->SetBranchAddress("tel_az", &Point_az);
    event_tree->SetBranchAddress("tel_al", &Point_al);
    event_tree->SetBranchAddress("weight", &weight);
    event_tree->SetBranchAddress("ntrig", &ntrig);

    event_tree->SetBranchAddress("Pe", pe_list);
    event_tree->SetBranchAddress("Paz", tel_az);
    event_tree->SetBranchAddress("Pal", tel_al);
    event_tree->SetBranchAddress("ltrig_list", &Trig_list);

}

bool LACTEvent::GetData()
{
    MCal = 90 -MCze;
    if( flag == 0)
    {
        return false;
    }
    if (flag == 1) 
    {
        for( int i = 0; i < ntrig; i++)
        {
            int itel = Trig_list[i] - 1;
            tel_data.push_back( new LACT_TelData(npix[itel]));
            tel_data[i]->SetTelid(itel);
            for( int j = 0; j < npix[itel]; j++)
            {
                tel_data[i]->GetPe()[j] = rd2->PoissonD(8.91) - 8.91 + pe_list[i][j];
            }
            tel_data[i]->SetPointDirection(tel_az[itel], tel_al[itel]);

        }
    }
    
    return true;
}

void LACTEvent::Reset()
{
    for( auto data : tel_data)
    {
        delete data;
    }
    tel_data.clear();
}

LACTEvent::LACTEvent()
{

}