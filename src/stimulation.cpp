#include "stimulation.h"
#include <cmath>
#include <iostream>

BiphasicStim::BiphasicStim(double amplitude, double pcl, double stim_shift, double pulseDuration)
: amplitude(amplitude), pcl(pcl), stim_shift(stim_shift), pulseDuration(pulseDuration)
{	
}

double BiphasicStim::i_stim(double t) const
{
    const double fmt = std::fmod(t, pcl) - std::round(stim_shift);
    if (fmt >= 0 && fmt < pulseDuration) {
        return 2 * amplitude / M_PI * std::atan(std::tan((2 * M_PI * fmt) / (2 * pulseDuration)));
    } else {
		return 0;
	}
}

double StimulationNone::i_stim(double t) const
{
	return 0;
}
BiphasicStim_CaSR_Protocol::BiphasicStim_CaSR_Protocol(double amplitude, double pcl_start, double pcl_end,
								double growth_time, double pcl_end_duration, double stim_shift, double pulseDuration)
: amplitude(amplitude), pcl_start(pcl_start), pcl_end(pcl_end), growth_time(growth_time), pcl_end_duration(pcl_end_duration),
	stim_shift(stim_shift), pulseDuration(pulseDuration)
{
}
/*
double BiphasicStim_CaSR_Protocol::i_stim(double t) const
{
	double pcl;
	if (t < growth_time) {
		//pcl = pcl_start;
		pcl = ( (growth_time - t) / growth_time * pcl_start + t / growth_time * pcl_end);
	} else if (t < growth_time + pcl_end_duration) {
		pcl = pcl_end;
	} else {
		//reset pcl
		pcl = pcl_start;
	}

	
    const double fmt = std::fmod(t, pcl) - std::round(stim_shift);
    if (fmt >= 0 && fmt < pulseDuration) {
        return 2 * amplitude / M_PI * std::atan(std::tan((2 * M_PI * fmt) / (2 * pulseDuration)));
    } else {
		return 0;
	}
}
*/


double BiphasicStim_CaSR_Protocol::i_stim(double t) const
{
	t -= stim_shift;
	const double N = growth_time;
	double pcl;
	
	double fmt = 0;
	for (int i = 0; i < N; i++) {
		const double stim_moment_growth = (double)i/(N-1) * ((N - (double)(i+1)/2) * pcl_start + (double)(i-1)/2*pcl_end);
	    const double fmt_x = t - std::round(stim_moment_growth);
		if (fmt_x < pulseDuration) {
			fmt = fmt_x;
			break;
		}
	}
	const double Ggrowth_time = N/2 * (pcl_start + pcl_end);
	
	if (fmt == 0 && t < Ggrowth_time + pcl_end_duration) {
		pcl = pcl_end;
	} else {
		//reset pcl
		pcl = pcl_start;
	}

	if (fmt == 0)
		fmt = std::fmod(t, pcl);
    if (fmt >= 0 && fmt < pulseDuration) {
        return 2 * amplitude / M_PI * std::atan(std::tan((2 * M_PI * fmt) / (2 * pulseDuration)));
    } else {
		return 0;
	}
}
