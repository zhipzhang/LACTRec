#include "LACT_RunPara.h"
#include <cstdlib>
#include <cstring>
#include "straux.h"


#ifdef IHEP
std::string prefix = "root://eos01.ihep.ac.cn/";
#else
std::string prefix = "";
#endif
void LACT_RUNPARA::ProcessCommandLine(int argc, char **argv)
{
    while( argc > 2)
        {
                if((strcmp(argv[1], "--lookup")) == 0 && lookup_file.empty())
                {
                        havelookup = true;
                        lookup_file = prefix +argv[2];
                        argc -= 2;
                        argv += 2;
                        continue;
                }
                if( (strcmp(argv[1], "--auto-lookup")) == 0 && !havelookup && lookup_file.empty())
                {
                        writelookup = true;
                        lookup_file = "lookup.out";
                        argc --;
                        argv ++;
                        continue;
                }
                if( (strcmp(argv[1], "--Draw")) == 0)
                {
                        DrawMode = true;
                        special_entries = atoi(argv[2]);
                        argc -= 2;
                        argv += 2;
                        continue;
                }
                if( (strcmp(argv[1], "--outfile")) == 0)
                {
                        out_file =  prefix +argv[2];
                        argc -= 2;
                        argv += 2;
                        continue;
                }
                if(strcmp(argv[1], "--only-telescopes") ==0 || 
                        strcmp(argv[1],"--only-telescope") == 0)
                {
                        SelectTel = true;
                        char word[20];
                        int ipos;
                        while (getword(argv[2], &ipos, word, sizeof(word)-1, ',','\n') > 0 )
                        {  
                                int tel_idx = atoi(word);
                                if(tel_idx > 0)
                                {
                                        Only_Telescope.push_back(tel_idx -1);
                                        printf("Only Telescope %d\n" , tel_idx);
                                }
                        }
                        argc -= 2;
                        argv += 2;
                        continue;
                }
                if(strcmp(argv[1], "--draw") == 0)
                {
                        DrawMode = true;
                        char word[20];
                        int ipos;
                        while( getword(argv[2], &ipos, word, sizeof(word) - 1, ',', '\n') > 0)
                        {
                                int event_idx = atoi(word);
                                if(event_idx >= 0)
                                {
                                        Draw_Events.push_back(event_idx);
                                }
                        }
                        argc -= 2;
                        argv += 2;
                        continue;

                }
                else 
                {
                        break;
                }

        }
        while( argc > 1)
        {
            input_file.push_back(prefix +argv[1]);
            argc --;
            argv ++;
        }

}

bool LACT_RUNPARA::ProcessCommandLine(int argc, char** argv, bool lookup_flag)
{
        if(lookup_flag)
        {
                writelookup = true;
        }
        else 
        {
                return false;
        }
        if ( writelookup)
        {
                while( argc > 2)
                {
                        if( strcmp(argv[1], "--out_file") == 0)
                        {
                                lookup_file = argv[2];
                                argc -= 2;
                                argv += 2;
                                continue;
                        }
                        else 
                        {
                                break;
                        }
                }

                while( argc > 1)
                {
                        input_file.push_back(prefix + argv[1]);
                        argc--;
                        argv++;
                }
        }
        return true;
}

LACT_RUNPARA::LACT_RUNPARA()
{
        havelookup = writelookup = DrawMode = SelectTel = false;
        lookup_file = "";
        out_file = "dst.root";

}