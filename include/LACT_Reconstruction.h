#ifndef _LACT_REC_H
#define _LACT_REC_H

#include "Limits_defined.h"
#include "LACT_Telconfig.h"
#include "LACT_TelData.h"
#include "LACTEvent.h"
#include "LACTRecEvent.h"
#include "LACT_TableCalculator.h"
#include "LACT_TableCalculatorData.h"
#include "LACT_RunPara.h"
#include <vector>
#include <map>
#include "TTree.h"
#include "TFile.h"

using namespace std;
class LACT_Reconstruction
{
    private:
        enum E_value {MRSW, MRSL, EREC};
        vector<LACT_TelConfig*> tel_config;

        int num_tel;
        bool UseSimpleSteroRec = true;                     // True Use Simple Direction Rec

        int num_only;
        bool known[LACT_MAXTEL]{0};

        int min_tel;
        int min_good_images;
        int min_size;
        double max_dist;
        double max_energy;
        double min_energy;
        int nums_err;

        int tail1;
        int tail2;

        bool havelookup;
        std::string lookup_file;
        TFile* LookupTableFile;

        bool DrawMode;
        std::vector<int> DrawEvents;

        //LACTTableCalculator* TableCalculator;
        //std::map<int, LACT_TableCalculatorData *>TableData;

    public:
        LACT_Reconstruction();
        ~LACT_Reconstruction();
        void SetLookupFile(std::string& filename)
        {
            havelookup = true;
            lookup_file = filename;
        }
        void SetOnlyTel(int i)
        {
            if( i > 0)
            {
                known[i] = 1;
                num_only++;
            }

        }
        void SetDrawMode()
        {
            DrawMode = true;
        }
        bool OnlyDraw()
        {
            return DrawMode;
        }
        void Draw(TTree*, LACTEvent*);
        void SetEventPix(LACTEvent* );
        void Draw_Events(LACTEvent*, int);
        void display(LACTRecEvent*, LACT_TelData*,int ievent, int i);
        void GetCommandConfig(LACT_RUNPARA* );
        void GetTelConfig(TTree* config_tree);
        void ComputePixNeighbor();
        void ComputeMoments(LACTEvent* event, LACTRecEvent* rec);
        void ConvertToRad(LACTRecEvent* rec);
        int CheckImageQuality(LACTEvent* event, LACTRecEvent* rec);
        bool SimpleSteroRec(LACTEvent* event, LACTRecEvent* rec);
        void SetMcData(LACTEvent* event, LACTRecEvent* rec);
        bool DirectionRec(LACTEvent* event, LACTRecEvent* rec);
        void FillTelRp(LACTRecEvent* rec);
        void EventRec(LACTEvent* event, LACTRecEvent* rec);
        void FillLookupTable(LACTRecEvent* );
        void InitLookupTableData();


};

























#endif