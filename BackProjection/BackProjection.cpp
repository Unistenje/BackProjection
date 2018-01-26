#include <iostream>
#include "stdafx.h"
#include "armadillo"
#include <math.h> 
#include <string>

#define PI 3.14159265

const int xResCam = 1024;
const int yResCam = 768;
const int xResPro = 1024;
const int yResPro = 768;

const int startKut = 80;
const int endKut = 120;

const int offProjector = 500;
const int offCameraX = 260;
const int offCameraY = 260;

const double cuttoffXmin = 0.0;
const double cuttoffYmin = 0.0;
const double cuttoffZmin = 0.0;
const double cuttoffXmax = 0.0;
const double cuttoffYmax = 0.0;
const double cuttoffZmax = 0.0;


using namespace arma;
using namespace std;

int main(int argc, char *argv[], char *envp[])
{
	if (argc > 2) {
		cout << "Application expects 2 arguments";
		return 0;
	}

	for (int d = startKut; d < endKut; d++) 
	{
		if (d != 90) 
		{
			double k = tan(d*PI / 180.0);
			double a = -1.0;
			double b = -a / k;

			double minKut = 90.0;
			double maxKut = 0.0;

			vector<double> sviMid;
			double sum = 0.0;

			for (int i = offProjector; i < (xResPro + yResPro) - offProjector; i++)
			{
				double u = 0.0;
				double v = 0.0;
				if (i <= 768)
				{
					v = 768.0 - i;
				}
				else
				{
					u = i - 768.0;
				}

				double c = (-a*u) - (b*v);

				//TODO: za kut "d" i pravac "l = [a b c]"	
				double *results;
				results = ProjektivneMatrice(a, b, c);
				double min = *results;
				if (min < minKut && min != 0.0) 
				{
					minKut = min;
				}
				double max = *(results + 1);
				if (maxKut > maxKut && max != 0.0)
				{
					maxKut = max;
				}
				double avg = *(results + 2);
				sum += avg;
				double mid = *(results + 3);
				sviMid.push_back(mid);

				int gjgh = 0;
			}
			if (sviMid.size() != 0)
			{
				int mid = sviMid.size() / 2;
				double midValue = sviMid.at(mid);
				double avg = sum / sviMid.size();
			}
			else 
			{
				double midValue = 0.0;
				double avg = 0.0;
			}
		}		
	}
	
    return 0;
}

double * ProjektivneMatrice(double a, double b, double c) 
{
	mat PCamera;
	mat PProjector;
	mat l;

	l << a << endr
		<< b << endr
		<< c << endr;

	PCamera << 1786.603948201665 << 0.0 << 573.7549315111892 << 0.0 << endr
		<< 0.0 << 1785.706892739146 << 350.2829444929243 << 0.0 << endr
		<< 0.0 << 0.0 << 1.0 << 0.0 << endr;

	PProjector << 2825.995772353162 << 274.0176300938570 << 1000.840054395188 << -242440.0087389195 << endr
		<< -190.6590424264215 << 3186.957811683789 << 434.8773172000017 << 37778.08800992880 << endr
		<< -0.132385180956547 << 0.229897452819599 << 0.964168722294061 << 64.524629745266440 << endr;

	mat rav = PProjector.t()*l;

	mat M = PCamera(span(0, 2), span(0, 2));
	
	mat p4 = PCamera(span(0, 2), span(3, 3));

	mat Minv = M.i();

	mat MinvP4 = -(Minv)*p4;

	mat jedan;
	jedan = 1.0;

	mat nula;
	nula = 0.0;

	mat C = join_cols(MinvP4, jedan);

	///C.print("C:");

	mat skCam;

	double minKut = 90.0;
	double maxKut = 0.0;

	///int maxBrojUListi = (xResCam-2*offCameraX)*(yResCam-1*offCameraY);

	vector<double> sviKutevi;
	double sum = 0.0;
	static double results[4];
	bool noCuttoff = cuttoffXmin == cuttoffXmax;

	for (int i = offCameraX; i < xResCam - offCameraX; i++)
	{
		for (int j = offCameraY; j < yResCam - offCameraY; j++)
		{			
			skCam << i+1 <<endr
				<< j+1 <<endr
				<< 1 << endr;

			mat temp1 = M.i()*skCam;
			mat D = join_cols(temp1, nula);
			mat L = (C*D.t()) - (D*C.t());
			mat pom = L*rav;
			mat sjec;
			sjec << pom(0, 0) / pom(3, 0) << endr
				<< pom(1, 0) / pom(3, 0) << endr
				<< pom(2, 0) / pom(3, 0) << endr;
			double x = sjec(0, 0);
			double y = sjec(1, 0);
			double z = sjec(2, 0);
			if (noCuttoff || (x <= cuttoffXmax && y <= cuttoffYmax && z <= cuttoffZmax && x >= cuttoffXmin && y >= cuttoffYmin && z >= cuttoffZmin))
			{
				double gore = (rav(0, 0)*L(3, 0) + rav(1, 0)*L(3, 1) + rav(2, 0)*L(3, 2));
				double dolje = (double)sqrt((double)pow(rav(0, 0), 2) + (double)pow(rav(1, 0), 2) + (double)pow(rav(2, 0), 2))
					*(double)sqrt((double)pow(L(3, 0), 2) + (double)pow(L(3, 1), 2) + (double)pow(L(3, 2), 2));
				double alpha = (acos(gore / dolje)) * 180.0 / PI;
				if (alpha > 90)
				{
					alpha = alpha - 90;
				}
				else if (alpha <= 90)
				{
					alpha = 90 - alpha;
				}

				///Priprema podataka za povratak

				sviKutevi.push_back(alpha);
				sum += alpha;
				if (alpha < minKut)
				{
					minKut = alpha;
				}
				if (alpha > maxKut)
				{
					maxKut = alpha;
				}
			}

		}
	}
	if (sviKutevi.size() != 0) 
	{
		int mid = sviKutevi.size() / 2;
		double midValue = sviKutevi.at(mid);
		double avg = sum / sviKutevi.size();
		results[0] = minKut;
		results[1] = maxKut;
		results[2] = avg;
		results[3] = midValue;
	}
	else 
	{
		results[0] = 0.0;
		results[1] = 0.0;
		results[2] = 0.0;
		results[3] = 0.0;
	}	
	return results;
}

