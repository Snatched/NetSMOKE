#include <iostream>
#include "ref_ligs.h"
#include "OpenSMOKEpp"

#include "maps/Maps_CHEMKIN"

OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>*				thermodynamicMapXML;

using namespace std;




int main()
{
	double c_sample, h_sample, o_sample, others_sample, totallig, LIGC, LIGO, LIGH;
	int molepres, masspres;
	FILE * pFile;
	pFile = fopen("Ligin_properties.info", "w");

	c_sample = 0.607; // input variables
	h_sample = 0.059;
	o_sample = 0.334;
	others_sample = 0; // all the other constituents of the sample such as H20, N, S
	totallig = c_sample + o_sample + h_sample + others_sample;

	if (totallig != 1)
	{
		printf("%s \n\n", "//////////////////// ERROR: Input sum different from 1.000 ////////////////////");
		fprintf(pFile, "%s\n%s\n", "INPUT DATA", "//////////////////// ERROR: Input sum different from 1.000 ////////////////////");
		system("pause");
		return 0;
	}

	molepres = 1;  // set=1 to obtain results in molar fractions
	masspres = 0;  // set=1 to obtain results in mass fractions

	ref_ligs Lig_sample;
	Lig_sample.set_values(c_sample, h_sample, o_sample, others_sample);

	if (molepres == masspres)
	{
		printf("%s \n\n", "//////////////////// ERROR: masspres or molepres are both set to 0 or 1. ////////////////////");
		fprintf(pFile, "%s\n%s%f\n%s%f\n%s%f\n%s%f\n\n%s\n%s\n", "INPUT DATA", "C: ", c_sample, "H: ", h_sample, "O: ", o_sample, "Other components: ", others_sample, "OUTPUT DATA", "//////////////////// ERROR: masspres or molepres are both set to 0 or 1. ////////////////////");
		system("pause");
		return 0;
	}

	if (molepres == 1)
	{
		LIGC = Lig_sample.ligc_mole();
		LIGO = Lig_sample.ligo_mole();
		LIGH = Lig_sample.ligh_mole();
		fprintf(pFile, "%s\n%s%f\n%s%f\n%s%f\n%s%f\n\n%s\n%s\n%s%f\n%s%f\n%s%f\n", "INPUT DATA", "C: ", c_sample, "H: ", h_sample, "O: ", o_sample, "Other components: ", others_sample, "OUTPUT DATA", "Molar fractions", "LIGC: ", LIGC, "LIGH: ", LIGH, "LIGO: ", LIGO);
	}

	if (masspres == 1)
	{
		LIGC = Lig_sample.ligc_mass();
		LIGO = Lig_sample.ligo_mass();
		LIGH = Lig_sample.ligh_mass();
		fprintf(pFile, "%s\n%s%f\n%s%f\n%s%f\n%s%f\n\n%s\n%s\n%s%f\n%s%f\n%s%f\n", "INPUT DATA", "C: ", c_sample, "H: ", h_sample, "O: ", o_sample, "Other components: ", others_sample, "OUTPUT DATA", "Mass fractions", "LIGC: ", LIGC, "LIGH: ", LIGH, "LIGO: ", LIGO);
	}


	printf("%s %f \n%s %f \n%s %f \n\n", "LIG-C: ", LIGC, "LIG-H: ", LIGH, "LIG-O: ", LIGO);
	system("pause");

	return 0;
}



