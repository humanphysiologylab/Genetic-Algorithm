#ifndef STIMULATION_H
#define STIMULATION_H

class StimulationBase
{
public:
	virtual double i_stim(double t) const= 0;
};

class StimulationNone:
	public StimulationBase
{
public:
	double i_stim(double t) const;
};
class BiphasicStim:
	public StimulationBase
{
	double amplitude, pcl, stim_shift, pulseDuration;
public:
	BiphasicStim(double amplitude, double pcl, double stim_shift, double pulseDuration);
	double i_stim(double t) const;
};

class BiphasicStim_CaSR_Protocol:
	public StimulationBase
{
	double amplitude, pcl_start, pcl_end, growth_time, pcl_end_duration, stim_shift, pulseDuration;
public:
	BiphasicStim_CaSR_Protocol(double amplitude, double pcl_start, double pcl_end, double growth_time, double pcl_end_duration, double stim_shift, double pulseDuration);

	double i_stim(double t) const;
};

#endif
