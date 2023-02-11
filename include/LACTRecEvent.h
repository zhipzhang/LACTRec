/*
    Class Which will be Written to the disk
*/

#ifndef _LACT_RECEVENT_H
#define _LACT_RECEVENT_H

#include "TObject.h"
#include <vector>
#include <map>
#include <cmath>
#include "LACTEvent.h"

class LACTRecEvent : public TObject
{
    public:
        double particle_id;
        double rec_energy;

        double MCenergy;
        double MCxcore;
        double MCycore;
        double rec_x;
        double rec_y;

        double MCaz;
        double MCal;
        double weight;
        double rec_azimuth;
        double rec_altitude;
        double rec_camerax;
        double rec_cameray;
        double direction_error;
        double Point_Az;                             // Ideal Point Position
        double Point_Al;

        int ntel = 0;
        std::vector<int>  tel_id;
        std::vector< double> length;                 // All Expressed in Rad
        std::vector< double> width;
        std::vector< double> size;
        std::vector< double> dist;
        std::vector< double> image_x;
        std::vector< double> image_y;
        std::vector< double> true_x;
        std::vector< double> true_y;
        std::vector< double> alpha;
        std::vector< double> tel_az;              //! Considered Pointing Error
        std::vector< double> tel_al;              //! Considered Pointing Error
        std::vector< double> rp; 

        std::vector< int>  good_image;
        int ngood_images;
        int npass_size;
        
        double mrsw;
        double mrsl;

        LACTRecEvent();
        virtual ~LACTRecEvent();
        void Reset();
        double GetRecAz()
        {
            return rec_azimuth;
        }
        double GetRecAlt()
        {
            return rec_altitude;
        }
        void SetDirectionError( double error)
        {
            direction_error = error;
        }
        
        void SetTrueSource( double x, double y) 
        {
            true_x.push_back(x);
            true_y.push_back(y);
        }
        void SetLength(double l)
        {
            length.push_back(l);
        }
        void SetWidth( double w)
        {
            width.push_back(w);
        }
        void SetSize(double s)
        {
            size.push_back(s);
        }
        void SetCameraSource(double x, double y)
        {
            rec_camerax = x;
            rec_cameray = y;
        }
        void SetImageSource(double x, double y)
        {
            image_x.push_back(x);
            image_y.push_back(y);
        }
        void SetAlpha(double a)
        {
            alpha.push_back(a);
        }
        void AddTel()
        {
            ntel++;
            good_image.push_back(0);
        }
        int GetNtel()
        {
            return ntel;
        }
        int GetTelid(int i)
        {
            return tel_id[i];
        }
        void ConvertRad(int i, double focal)
        {
            length[i] = length[i] / focal;
            width[i]  = width[i]  /focal;
            image_x[i] = image_x[i] / focal;
            image_y[i] = image_y[i] / focal;
            true_x[i]  = true_x[i]  / focal;
            true_y[i]  = true_y[i]  / focal;
        }
        double  GetTelSize(int i )
        {
            return size[i];
        }             // Be careful i is not the Tel_id
        double GetTelDist(int i)
        {
            return sqrt(pow(image_x[i], 2) + pow(image_y[i],2));
        }
        void SetPassSize(int i)
        {
            good_image[i] = 1;
        }
        void SetGoodImage(int i)
        {
            good_image[i] = 2;
        }
        void SetNGood(int n)
        {
            ngood_images = n;
        }
        void SetNPassSize(int n)
        {
            npass_size = n;
        }
        int GetPassSize()
        {
            return npass_size;
        }
        int GetStatus(int i)
        {
            return good_image[i];               // 0 : mean < min_size 1: > min_size and dist > max_dist 2: good images
        }
        double GetTelImageX(int i)
        {
            return image_x[i];
        }
        double GetTelImageY(int i)
        {
            return image_y[i];
        }
        double GetTelAlpha(int i)
        {
            return alpha[i];
        }
        double GetPointAz()
        {
            return Point_Az;
        }
        double GetPointAl()
        {
            return Point_Al;
        }
        double GetTelaz(int i)
        {
            return tel_az[i];
        }
        double GetTelal(int i)
        {
            return tel_al[i];
        }
        void SetTelAz(double az)
        {
            tel_az.push_back(az);
        }
        void SetTelAl(double al)
        {
            tel_al.push_back(al);
        }
        double GetTelLength(int i)
        {
            return length[i];
        } 
        double GetTelWidth(int i)
        {
            return width[i];
        }
        void SetRecDirection(double az, double al)
        {
            rec_azimuth = az;
            rec_altitude = al;
        }
        void SetRecCore(double x, double y)
        {
            rec_x = x;
            rec_y = y;
        }
        double  GetTrueCamereX(int i)
        {
            return true_x[i];
        }
        double GetTrueCameraY(int i)
        {
            return true_y[i];
        }
        void GetMCData(LACTEvent*) ;
        void Display(LACTEvent*);

        ClassDef(LACTRecEvent, 1)


};



























#endif