#include "LACTRecEvent.h"
#include "LACTEvent.h"
#include "Rtypes.h"

ClassImp(LACTRecEvent)
void LACTRecEvent::GetMCData(LACTEvent *event)
{
    MCenergy = event->GetMCenergy();
    MCxcore  = event->GetMCxcore();
    MCycore  = event->GetMCycore();
    MCal     = event->GetMCal();
    MCaz     = event->GetMCaz();
    particle_id = event->GetParticleid();
    weight   = event->GetWeight();
    Point_Al  = event->GetPointAl();
    Point_Az  = event->GetPointAz();
    

}

LACTRecEvent::LACTRecEvent()
{
    particle_id = -1;
    rec_energy = MCenergy = MCxcore = MCycore = rec_x = rec_y = -9999;
    MCal = MCaz = weight = rec_altitude = rec_azimuth = rec_camerax = rec_cameray = direction_error = Point_Al = Point_Al = -999;
    ngood_images = npass_size = 0;
    mrsl = mrsw = - 9999;
}

void LACTRecEvent::Reset()
{
    particle_id = -1;
    rec_energy = MCenergy = MCxcore = MCycore = rec_x = rec_y = -9999;
    MCal = MCaz = weight = rec_altitude = rec_azimuth = rec_camerax = rec_cameray = direction_error = Point_Al = Point_Al = -999;
    ngood_images = npass_size = 0;
    ntel = 0;
    mrsl = mrsw = -9999;
    tel_id.clear();
    length.clear();
    width.clear();
    size.clear();
    dist.clear();
    image_x.clear();
    image_y.clear();
    true_x.clear();
    true_y.clear();
    alpha.clear();
    tel_az.clear();
    tel_al.clear();
    rp.clear();
    rec_rp.clear();
}

LACTRecEvent::~LACTRecEvent()
{
    tel_id.clear();
    length.clear();
    width.clear();
    size.clear();
    dist.clear();
    image_x.clear();
    image_y.clear();
    true_x.clear();
    true_y.clear();
    alpha.clear();
    tel_az.clear();
    tel_al.clear();
    rp.clear();
    rec_rp.clear();
}