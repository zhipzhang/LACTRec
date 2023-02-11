#include "LACT_Reconstruction.h"
#include "LACTEvent.h"
#include "LACTRecEvent.h"
#include "LACT_TableCalculatorData.h"
#include "LACT_Telconfig.h"
#include "Limits_defined.h"
#include "LACT_Utilities.h"
#include "LACT_TelData.h"
#include "RtypesCore.h"
#include "TCanvas.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TString.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>

//default setting
LACT_Reconstruction::LACT_Reconstruction()
{
    min_tel = 4;
    min_good_images = 2;
    min_size = 200;
    max_dist = 4.;
    tail1 = 5;
    tail2 = 10;

    havelookup = false;
    num_only = 0; 
    
    DrawMode = false;
    lookup_file = "";

}

void LACT_Reconstruction::GetCommandConfig(LACT_RUNPARA * LACT_Runpara)
{

        if( LACT_Runpara->havelookup)
        {
                SetLookupFile(LACT_Runpara->lookup_file);
        }
        if( LACT_Runpara->GetOnlyNum() > 0)
        {
                for( auto i: LACT_Runpara->GetOnlyTels())
                {
                        SetOnlyTel(i);
                }

        }
        if( LACT_Runpara->DrawMode)
        {
            SetDrawMode();
            DrawEvents = LACT_Runpara->Draw_Events;
        }
}
void LACT_Reconstruction::GetTelConfig(TTree *config_tree)
{
    int n = config_tree->GetEntries();
    if( n > LACT_MAXTEL)
    {
        std::cout << " There are Too mamy Telescopes ! Program Exited" << std::endl;
        exit(EXIT_FAILURE);
    }
    num_tel = n;

    float *x_pix, *y_pix, *pix_size;
    int *pixel_shape;
    int npix;
    float focal_length, effective_focal_length;
    float pos[3];

    x_pix = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    y_pix = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    pix_size = (float*) malloc(LACT_MAXPIXELS * sizeof(float));
    pixel_shape = (int *) malloc(LACT_MAXPIXELS * sizeof(int));

    config_tree->SetBranchAddress("Tel_pos", pos);
    config_tree->SetBranchAddress("Npix", &npix);
    config_tree->SetBranchAddress("Focal_length", &focal_length);
    config_tree->SetBranchAddress("Pixel_Shape", pixel_shape);
    config_tree->SetBranchAddress("Pixel_Size", pix_size);
    config_tree->SetBranchAddress("X_Pix", x_pix);
    config_tree->SetBranchAddress("Y_Pix", y_pix);


    for( int i = 0 ; i < n; i++)
    {
        config_tree->GetEntry(i);
        tel_config.push_back(new LACT_TelConfig());

        if(tel_config[i]->SetNpix(npix) < 0)
        {
            std::cout << "There are too many pixels !" << std::endl;
            exit(EXIT_FAILURE);
        }
        if(tel_config[i]->SetTelPos(3, pos) < 0)
        {
            std::cout << "Error Setting the Tel_Pos" << std::endl;
        }
        tel_config[i]->SetPix_Shape_Size(pix_size[0], pixel_shape[0]);
        if( tel_config[i]->SetPixPos(x_pix, y_pix) < 0)
        {
            std::cout << "Error Setting pixel Position" << std::endl;
            exit(EXIT_FAILURE);
        }
        tel_config[i]->SetFocalLength(focal_length);
    }
    for( int i = 0; i < n; i++)
    {
        tel_config[i]->ComputeNeighbor();
    }

    free(x_pix);
    free(y_pix);
    free(pix_size);
    free(pixel_shape);
}

void LACT_Reconstruction::ComputeMoments(LACTEvent *event, LACTRecEvent *rec)
{

    
    for( auto itel = event->GetTelData().begin(); itel != event->GetTelData().end(); ++itel)
    {
        int tel_id = (*itel)->GetTelid();
        double sx = 0, sxx = 0, sxy = 0, sy = 0, syy  = 0, sA = 0. ,sx3 = 0.;
        double beta, cb, sb;
        double a = 0, b = 0;
        if( num_only > 0)
        {
            if(!known[tel_id])
                continue;
        }
        if( (*itel)->GetFlag() )
        {
            rec->AddTel();
            rec->tel_id.push_back(tel_id);
            rec->SetTelAz((*itel)->GetTelAzimuth());
            rec->SetTelAl((*itel)->GetTelAltitude());

        }
        else 
        {
            continue;
        }
        double dx, dy;
        Utilities::angles_to_offset(event->GetAzimuth() * TMath::DegToRad(), event->GetAltitude() * TMath::DegToRad(), (*itel)->GetTelAzimuth() * TMath::DegToRad()
                , (*itel)->GetTelAzimuth() * TMath::DegToRad(), tel_config[tel_id]->GetFocalLength(),  &dx, &dy);
        
        rec->SetTrueSource(dx, dy);

        for( int j = 0; j < (*itel)->GetImagePixnum(); j++)
        {
            int ipix = (*itel)->GetImagePixId(j);
            double x = tel_config[tel_id]->GetPixX(ipix) - dx;
            double y = tel_config[tel_id]->GetPixY(ipix) - dy;
            double A = (*itel)->GetPe(ipix);
            sA  += A;
            sx  += (A * x);
            sxx += (A * x) * x;
            sxy += (A * x) * y;
            sy  += (A * y);
            syy += (A * y) * y;
        }
        sx /= sA;
        sy /= sA;
        sxx = sxx/sA - sx * sx;
        sxy = sxy/sA - sx * sy;
        syy = syy/sA - sy * sy;
        if (fabs(sxy) > 1e-3 * fabs(sxx) && fabs(sxy) > 1e-3 * fabs(syy) )
        {
            double p1 = syy - sxx, p2 = sxy*sxy;
            double q, r1, r2;
            if ( p2 > 1e-8*(p1*p1) )
                q = p1 + sqrt(p1*p1+4.*p2);
            else
                q = 2.*p2;
            b = 0.5 * q/sxy; // solve  the equation
            a = sy - b*sx;
            if ( (r1 = syy + 2.*p2/q) > 0. )
                    rec->SetLength(sqrt(r1));
            else
                    rec->SetLength(0);
            if ( (r2 = sxx - 2.*p2/q) > 0. )
                    rec->SetWidth(sqrt(r2));
            else
                    rec->SetWidth(0.);
        }
        else 
        {
            if ( fabs(syy) < 1e-8*fabs(sxx) )
                syy = 0.;
            else if ( fabs(sxx) < 1e-8*fabs(syy) )
                sxx = 0.;
            if( sxx > syy && syy >= 0.)
            {
                rec->SetLength(sqrt(sxx));
                rec->SetWidth(sqrt(syy));
                b = 0.;
                a = sy;
            }
            else if(syy >= 0. && sxx >= 0.)
            {
                rec->SetLength(sqrt(syy));
                rec->SetWidth(sqrt(sxx));
                b = 100000.;
                a = sy - b * sx;
            }
        }
        beta = atan(b);
        cb   = cos(beta);
        sb   = sin(beta);

        rec->SetSize(sA);
        rec->SetImageSource(sx + dx, sy + dy);
        sx = sx3 = 0.;
        for( int k = 0; k < (*itel)->GetImagePixnum(); k++)
        {
            int ipix = (*itel)->GetImagePixId(k);
            double x = tel_config[tel_id]->GetPixX(ipix);
            double y = tel_config[tel_id]->GetPixY(ipix);
            double A = (*itel)->GetPe(ipix);
            double xp;
            xp   = cb * (x - sx) + sb *( y - sy);
            sxx += (A * xp) * xp;
            sx3 += ((A * xp) * xp) * xp;
        }
        if( sx3 / pow(sxx, 1.5) < 0)
        {
            beta = beta * TMath::RadToDeg() + 180;
        }
        else 
        {
            beta = beta * TMath::RadToDeg();
        }
        rec->SetAlpha(beta * TMath::DegToRad());
    }
}


void LACT_Reconstruction::ConvertToRad(LACTRecEvent *rec)
{
    for( int i = 0; i < rec->GetNtel(); i++)
    {
        int tel_id = rec->GetTelid(i);
        rec->ConvertRad(i, tel_config[tel_id]->GetFocalLength());
    }

}

int LACT_Reconstruction::CheckImageQuality(LACTEvent* event, LACTRecEvent* rec)
{
    if(rec->GetNtel() < min_tel)
    {
        return -1;
    }
    int pass = 0;
    int good_images = 0;
    for( int i = 0; i < rec->GetNtel(); i++)
    {
        if( rec->GetTelSize(i) > min_size)
        {
            pass++;
            rec->SetPassSize(i);
            if(rec->GetTelDist(i)  * TMath::RadToDeg() < max_dist)
            {
                good_images++;
                rec->SetGoodImage(i);
            }
        }
        
    }
    if( pass < min_tel || good_images < min_good_images)
    {
        return -1;
    }
    rec->SetNPassSize(pass);
    rec->SetNGood(good_images);
    return 1;
}

bool LACT_Reconstruction::SimpleSteroRec(LACTEvent* event, LACTRecEvent *rec)
{
    int ntel = rec->GetNtel();
    double *x_ref  = new double[ntel];
    double *y_ref  = new double[ntel];
    double *alpha_ref = new double[ntel];
    double *xt     = new double[ntel];
    double *yt     = new double[ntel];

    for( int i = 0; i < ntel; i++)
    {
        x_ref[i] = y_ref[i] = alpha_ref[i] = 0.;
    }
    
    double w ; //weight we will use
    double sum_xs, sum_xs2, sum_ys, sum_ys2, sum_w;
    double xs, ys, angs; // intersect point and the angle between two lines
    double amp_red;
    double rec_az, rec_alt;
    double trans[3][3];
    double xh, yh, zh;
    double xc, yc;

    for( int i = 0; i < ntel; i++)
    {
        if(rec->GetStatus(i) >= 1)
        {
            int tel_id = rec->GetTelid(i);
            Utilities::cam_to_ref(rec->GetTelImageX(i), rec->GetTelImageY(i), rec->GetTelAlpha(i), rec->GetPointAz() * TMath::DegToRad()
                                , rec->GetPointAl() * TMath::DegToRad(), 0, rec->GetTelaz(i) * TMath::DegToRad(), rec->GetTelal(i) * TMath::DegToRad(), 1.0, &x_ref[i], &y_ref[i], &alpha_ref[i]);
        }
    }

    sum_xs = sum_ys = sum_w = sum_xs2 = sum_ys2 = 0.;
    for(int i = 0; i < ntel; i++)
    {
        if( rec->GetStatus(i) >= 1)
        {
            for( int j = 0; j < i; j++)
            {
                if( rec->GetStatus(j) >= 1)
                {
                    if(Utilities::intersect_lines(x_ref[i], y_ref[i], alpha_ref[i], x_ref[j], y_ref[j], alpha_ref[j], &xs, &ys, &angs) != 1)
                    {
                        continue;
                    }
                }
                else 
                {
                    continue;
                }
                amp_red  = (rec->GetTelSize(i) * rec->GetTelSize(j))/(rec->GetTelSize(i) + rec->GetTelSize(j));
                w        =  pow(amp_red * sin(angs) * (1 - rec->GetTelLength(i)/rec->GetTelWidth(i)) * (1 - rec->GetTelLength(j)/rec->GetTelWidth(j)), 2);
                sum_w   +=  w;
                sum_xs  +=  xs * w;
                sum_xs2 +=  xs * xs * w;
                sum_ys  +=  ys * w;
                sum_ys2 +=  ys * ys * w;

            }
        }
        else 
        {
            continue;
        }
    }
    if( fabs(sum_w) < 1e-10)
    {
        delete [] x_ref;
        delete [] y_ref;
        delete [] alpha_ref;
        delete [] xt;
        delete [] yt;
        nums_err++;
        return false;
    }
    sum_xs /= sum_w;
    sum_ys /= sum_w;
    rec->SetCameraSource(sum_xs, sum_ys);
    Utilities::offset_to_angles(sum_xs, sum_ys, rec->GetPointAz() * TMath::DegToRad(), rec->GetPointAl() * TMath::DegToRad(), 1.0, &rec_az, &rec_alt);
    rec_az -= (2.*TMath::Pi()) * floor(rec_az / (2 * TMath::Pi()));
    rec->SetRecDirection(rec_az * TMath::RadToDeg(), rec_alt * TMath::RadToDeg());

    Utilities::get_shower_trans_matrix(rec_az, rec_alt, trans);

    for( int i = 0; i < ntel; i++)
    {
        int tel_id = rec->GetTelid(i);
        xt[i]  = trans[0][0] * tel_config[tel_id]->GetTelpos(0) + 
                 trans[0][1] * tel_config[tel_id]->GetTelpos(1) +
                 trans[0][2] * tel_config[tel_id]->GetTelpos(2);
        yt[i]  = trans[1][0] * tel_config[tel_id]->GetTelpos(0) +
                 trans[1][1] * tel_config[tel_id]->GetTelpos(1) +
                 trans[1][2] * tel_config[tel_id]->GetTelpos(2);
    }
    sum_xs = sum_ys = sum_w = sum_xs2 = sum_ys2 = 0.;
    for( int i = 0; i < ntel; i++)
    {
        if(rec->GetStatus(i) >= 1)
        {
            for( int j = 0; j < i; j++)
            {
                if(rec->GetStatus(j) >= 1 )
                {
                    if( Utilities::intersect_lines(xt[i], yt[i], rec->GetTelAlpha(i), xt[j], yt[j], rec->GetTelAlpha(j), &xs, &ys, &angs) != 1)
                    {
                        continue;
                    }
                }
                else
                {
                    continue;
                }
                amp_red  = (rec->GetTelSize(i) * rec->GetTelSize(j))/(rec->GetTelSize(i) + rec->GetTelSize(j));
                w        =  pow(amp_red * sin(angs) * (1 - rec->GetTelLength(i)/rec->GetTelWidth(i)) * (1 - rec->GetTelLength(j)/rec->GetTelWidth(j)), 2);
                sum_w   +=  w;
                sum_xs  +=  xs * w;
                sum_xs2 +=  xs * xs * w;
                sum_ys  +=  ys * w;
                sum_ys2 +=  ys * ys * w;
            }
        }
        else 
        {
            continue;
        }
    }
    if( sum_w == 0)
    {
        return false;
    }
    xs = sum_xs / sum_w;
    ys = sum_ys / sum_w;
    xh = trans[0][0] * xs +
        trans[1][0] * ys;
    yh = trans[0][1] * xs +
        trans[1][1] * ys;
    zh = trans[0][2] * xs +
        trans[1][2] * ys;
    
    xc = xh - trans[2][0]*zh/trans[2][2];
    yc = yh - trans[2][1]*zh/trans[2][2];
    rec->SetRecCore(xc, yc);

    delete [] x_ref;
    delete [] y_ref;
    delete [] alpha_ref;
    delete [] xt;
    delete [] yt;
    return true;
}

void LACT_Reconstruction::SetMcData(LACTEvent *event, LACTRecEvent *rec)
{
    rec->GetMCData(event);
}
void LACT_Reconstruction::SetEventPix(LACTEvent * event)
{
    for(int i = 0; i < num_tel; i++)
    {
        event->SetTelpix(i, tel_config[i]->GetNpix());
    }

}

void LACT_Reconstruction::EventRec(LACTEvent *event, LACTRecEvent *rec)
{
    SetMcData(event, rec);
    for( auto itel_data = event->GetTelData().begin(); itel_data != event->GetTelData().end(); ++itel_data)
    {
        int tel_id = (*itel_data)->GetTelid();
        (*itel_data)->ImageClean(tail1, tail2,  tel_config[tel_id]->GetNeighbor());
    }
    ComputeMoments(event, rec);
    ConvertToRad(rec);
    if(CheckImageQuality(event, rec) > 0)
    {
        if(DirectionRec(event, rec))
        {
            double direction_error = Utilities::angle_between(event->GetAzimuth() * TMath::DegToRad(), event->GetAltitude() * TMath::DegToRad(), rec->GetRecAz() * TMath::DegToRad(), rec->GetRecAlt() * TMath::DegToRad());
            rec->SetDirectionError(direction_error * TMath::RadToDeg());
            FillTelRp(rec);
        }
        else 
        {
            nums_err ++;
        }
    }

}

void LACT_Reconstruction::FillTelRp(LACTRecEvent *rec)
{
    for( int i = 0; i < rec->GetNtel(); i++)
    {
        int tel_id = rec->GetTelid(i);
        double rp = Utilities::line_point_distance(rec->GetMCCoreX(), rec->GetMCCoreY(), 0, cos(rec->GetMCaz() * TMath::DegToRad()) * cos(rec->GetMCal() * TMath::DegToRad()), 
                                -sin(rec->GetMCaz() * TMath::DegToRad()) * cos(rec->GetMCal() * TMath::DegToRad()), sin(rec->GetMCal() * TMath::DegToRad()), tel_config[tel_id]->GetTelpos(0),tel_config[tel_id]->GetTelpos(1), tel_config[tel_id]->GetTelpos(2));
        rec->SetTelRp(rp);
        double rec_rp = Utilities::line_point_distance(rec->GetRecCoreX(), rec->GetRecCoreY(), 0, cos(rec->GetRecAz() * TMath::DegToRad()) * sin(rec->GetRecAlt() * TMath::DegToRad()),
                         -sin(rec->GetRecAz() * TMath::DegToRad()) * cos(rec->GetRecAlt() * TMath::DegToRad()), sin(rec->GetRecAlt() * TMath::DegToRad()), tel_config[tel_id]->GetTelpos(0), tel_config[tel_id]->GetTelpos(1), tel_config[tel_id]->GetTelpos(2));
        rec->SetTelRecRp(rec_rp);
    }
}
void LACT_Reconstruction::ComputePixNeighbor()
{
    for( auto itel_config = tel_config.begin(); itel_config != tel_config.end(); ++itel_config)
    {
        (*itel_config)->ComputeNeighbor();
    }
}

bool LACT_Reconstruction::DirectionRec(LACTEvent *event, LACTRecEvent *rec)
{
    if(UseSimpleSteroRec)
    {
         return SimpleSteroRec(event, rec);
    }
    else 
    {
        return false;
    }
}

/*void LACT_Reconstruction::InitLookupTableData()
{
    TableData[MRSW] = new LACT_TableCalculatorData();
    TableData[MRSW]->fEnergy = false;
    TableData[MRSW]->fFillVariable = "width";

    TableData[MRSL] = new LACT_TableCalculatorData();
    TableData[MRSL]->fEnergy = false;
    TableData[MRSL]->fFillVariable = "length";

    TableData[EREC] = new LACT_TableCalculatorData();
    TableData[EREC]->fEnergy = true;
    TableData[EREC]->fFillVariable = true;

}
*/
void LACT_Reconstruction::Draw(TTree* eventTree, LACTEvent* event)
{
    if(DrawMode && !DrawEvents.empty())
    {
        for( auto ievent : DrawEvents)
        {
            eventTree->GetEntry(ievent);
            SetEventPix(event);
            event->GetData();
            Draw_Events(event, ievent);
        }

    }
}

void LACT_Reconstruction::Draw_Events(LACTEvent* event, int ievent)
{
    LACTRecEvent* rec = new LACTRecEvent();
    EventRec(event, rec);
    for(int i = 0; i < rec->GetNtel(); i++)
    {

        int tel_id = rec->GetTelid(i);
        if(rec->good_image[i] >= 1)
        {
            int index = event->GetTelIndex(tel_id);
            display(rec, event->GetTelData(index), ievent, i);
        }

    }
}

void LACT_Reconstruction::display(LACTRecEvent *rec, LACT_TelData* iteldata, int ievent, int i)
{
    int tel_id = iteldata->GetTelid();
    TCanvas* camera_image =  new TCanvas(Form("Event %d of Camera %d", ievent, tel_id), Form("LACT IMAGE"), 1600, 1600);
    TH2Poly* camera       =  new TH2Poly(Form("Event %d Camera %d", ievent, tel_id),"camera", -6, 6, -6, 6);
    for( int i = 0; i < iteldata->GetImagePixnum(); i++)
    {
        int ipix = iteldata->GetImagePixId(i);
        double binsize = tel_config[tel_id]->GetPixSize() / tel_config[tel_id]->GetFocalLength() * TMath::RadToDeg();
        double x = tel_config[tel_id]->GetPixX(ipix) / tel_config[tel_id]->GetFocalLength() * TMath::RadToDeg();
        double y = tel_config[tel_id]->GetPixY(ipix) / tel_config[tel_id]->GetFocalLength() * TMath::RadToDeg();
        double bin_x[4] = {x - 0.5 * binsize, x + 0.5 * binsize, x + 0.5 * binsize, x - 0.5 * binsize};
        double bin_y[4] = {y - 0.5 * binsize, y - 0.5 * binsize, y + 0.5 * binsize, y + 0.5 * binsize};
        camera->AddBin(4, bin_x, bin_y);
        camera->Fill(x, y, iteldata->GetPe(ipix));
    }
    camera->SetStats(0);
    camera->Draw("colz");

    TGraph* g1 = new TGraph();
    g1->SetPoint(g1->GetN(), rec->GetTrueCameraY(i) , rec->GetTrueCameraY(i));
    g1->SetMarkerStyle(4);
    g1->SetMarkerSize(4);
    g1->SetMarkerColor(2);
    TGraph* g2 = new TGraph();
    g2->SetPoint(g1->GetN(), rec->rec_camerax, rec->rec_cameray);
    g2->SetMarkerStyle(5);
    g2->SetMarkerColor(2);
    g2->SetMarkerSize(4);

    TMultiGraph* mg = new TMultiGraph();
    mg->Add(g1);
    mg->Add(g2);
    camera_image->cd();
    mg->Draw("p");

    TEllipse* ellipse = new TEllipse(rec->GetTelImageX(i) * TMath::RadToDeg(), rec->GetTelImageY(i) * TMath::RadToDeg(),
                                    rec->GetTelLength(i) * TMath::RadToDeg(), rec->GetTelWidth(i) * TMath::RadToDeg(), 0, 360, rec->GetTelAlpha(i) * TMath::RadToDeg());
    ellipse->SetLineWidth(2);
    ellipse->SetLineColor(2);
    ellipse->SetFillStyle(0);
    ellipse->Draw();

    TEllipse* el2 = new TEllipse(0, 0, 5.0, 5.0, 0, 360);
    el2->SetLineWidth(2);
    el2->SetLineColor(kBlack);
    el2->SetFillStyle(0);
    el2->Draw();
    
    TLine * line = new TLine(rec->GetTelImageX(i) * TMath::RadToDeg()  - 5 * TMath::Cos(rec->GetTelAlpha(i)), rec->GetTelImageY(i) * TMath::RadToDeg() - 5 * TMath::Sin(rec->GetTelAlpha(i)), rec->GetTelImageX(i) * TMath::RadToDeg() + 5 * TMath::Cos(rec->GetTelAlpha(i)), rec->GetTelImageY(i) * TMath::RadToDeg() + 5 * TMath::Sin(rec->GetTelAlpha(i)));
    line->SetLineWidth(2);
    line->SetLineColor(4);
    line->Draw();

    TPaveText *pavet = new TPaveText(-6, 6.3, 6, 7.6);
    pavet->SetFillStyle(0);
    pavet->Draw("same");
    std::cout << "Begin Save Event" <<ievent << " Tel " << tel_id << std::endl; 
    std::cout << "Camera Size is " << rec->GetTelSize(i) << std::endl;
    camera_image->SaveAs(Form("Event%d_camera%d.png", ievent, tel_id));
    delete camera_image;
    delete camera;
}