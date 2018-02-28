
namespace OpenSMOKE
{

	void BiomassInterpreter::SetValues(double C_fraction, double H_fraction, double O_fraction, double others_fraction) {
		C = C_fraction;
		H = H_fraction;
		O = O_fraction;
		others_comp = others_fraction;
		CheckValues();
	}

	void BiomassInterpreter::CheckValues(){

		double totallig = c_sample + o_sample + h_sample + others_sample;
		
		if (totallig != 1)
		{
			printf("%s \n\n", "//////////////////// ERROR: Input sum different from 1.000 ////////////////////");
			fprintf(pFile, "%s\n%s\n", "INPUT DATA", "//////////////////// ERROR: Input sum different from 1.000 ////////////////////");
			OpenSMOKE::FatalErrorMessage();
		}
	}

}