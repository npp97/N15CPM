#include <R.h>

/* A typical C trick to get readable names for the parameters is #define
   This method is simple and efficient, but there are, of course,
   other possibilities that use dynamic variables.
*/

static double parms[32];

#define C_NRec_2_NH4		  	parms[0]
#define Vmax_NH4_2_NRec			parms[1]
#define Km_NH4_2_NRec		  	parms[2]
#define Vmax_NH4_2_NLab			parms[3]
#define Km_NH4_2_NLab			  parms[4]
#define Vmax_NLab_2_NH4			parms[5]
#define Km_NLab_2_NH4			  parms[6]
#define Vmax_NO3_2_NRec			parms[7] 
#define Km_NO3_2_NRec			  parms[8]
#define C_NRec_2_NO3			  parms[9] 
#define K_NH4ads_2_NH4			parms[10] 
#define K_NH4_2_NH4ads			parms[11]
#define Vmax_NO3_2_NH4			parms[12] 
#define Km_NO3_2_NH4			  parms[13]
#define K_NH4_2_NO2nit			parms[14]
#define K_NO2nit_2_NO3			parms[15]
#define K_NO3_2_NO3sto			parms[16]
#define K_NO3sto_2_NO3			parms[17]
#define K_NO3_2_NO2den			parms[18]
#define K_NRec_2_NO2org			parms[19]
#define K_NO2org_2_N2Oorg		parms[20]
#define K_NO2nit_2_N2Onit		parms[21]
#define K_NO2den_2_N2Oden		parms[22]
#define K_NO2denNRec_2_N2Ocod	parms[23]
#define K_N2Oorg_2_N2			  parms[24]
#define K_N2Ocod_2_N2			  parms[25]
#define K_N2Oden_2_N2			  parms[26]
#define K_N2Onit_2_N2			  parms[27]
#define K_N2Oorg_2_N2OorgE		parms[28]
#define K_N2Ocod_2_N2OcodE		parms[29]
#define K_N2Oden_2_N2OdenE		parms[30]
#define K_N2Onit_2_N2OnitE		parms[31]


/* initializer  */
void initmod(void(*odeparms)(int *, double *))
{
	int N = 32;
	odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs(int *neq, double *t, double *y, double *ydot,
	double *yout, int *ip)
{
	//----------Full
	double NLab_2_NH4, NH4_2_NLab, NRec_2_NH4, NRec_2_NO3, NO3_2_NRec, NH4ads_2_NH4, NO3_2_NH4, NO3_2_NO3sto, NO3sto_2_NO3;
	double NH4_2_NRec, NH4_2_NH4ads;
	double NLab, NH4, NO3, NH4ads, NRec, NO3sto;
	double dNLab, dNH4, dNO3, dNH4ads, dNRec, dNO3sto;

	double NO2,NO2org,NO2nit,NO2den;
	double N2Oorg,N2Ocod,N2Onit,N2Oden;
	double N2O,EN2Oorg,EN2Ocod,EN2Onit,EN2Oden;
	double N2;

	double NO3_2_NO2den, NH4_2_NO2nit,NO2nit_2_NO3,NRec_2_NO2org;
	double NO2org_2_N2Oorg,NO2nit_2_N2Onit,NO2den_2_N2Oden,NO2denNRec_2_N2Ocod;

	double dNO2org,dNO2nit,dNO2den,dNO2;
	double dN2Oorg,dN2Oden,dN2Onit,dN2Ocod;
	double dEN2Oorg,dEN2Onit,dEN2Oden,dEN2Ocod,dEN2O;
	double dRN2Oorg,dRN2Onit,dRN2Oden,dRN2Ocod,dN2;

//-----------N15

	double NLab_2_NH4_15, NH4_2_NLab_15, NRec_2_NH4_15, NRec_2_NO3_15, NO3_2_NRec_15, NH4ads_2_NH4_15, NO3_2_NH4_15, NO3_2_NO3sto_15, NO3sto_2_NO3_15;
	double NH4_2_NRec_15, NH4_2_NH4ads_15;
	double NLab_15, NH4_15, NO3_15, NH4ads_15, NRec_15, NO3sto_15;
	double dNLab_15, dNH4_15, dNO3_15, dNH4ads_15, dNRec_15, dNO3sto_15;

	double NO2_15,NO2org_15,NO2nit_15,NO2den_15;
	double N2Oorg_15,N2Ocod_15,N2Onit_15,N2Oden_15;
	double N2O_15,EN2Oorg_15,EN2Ocod_15,EN2Onit_15,EN2Oden_15;
	double N2_15;

	double NO3_2_NO2den_15, NH4_2_NO2nit_15,NO2nit_2_NO3_15,NRec_2_NO2org_15;
	double NO2org_2_N2Oorg_15,NO2nit_2_N2Onit_15,NO2den_2_N2Oden_15,NO2denNRec_2_N2Ocod_15;

	double dNO2org_15,dNO2nit_15,dNO2den_15,dNO2_15;
	double dN2Oorg_15,dN2Oden_15,dN2Onit_15,dN2Ocod_15,dN2O_15;
	double dEN2Oorg_15,dEN2Onit_15,dEN2Oden_15,dEN2Ocod_15,dEN2O_15;
	double dRN2Oorg_15,dRN2Onit_15,dRN2Oden_15,dRN2Ocod_15,dN2_15;


	//

	NLab = y[0];
	NLab_15 = y[1];

	NH4 = y[2];
	NH4_15 = y[3];

	NO3 = y[4];
	NO3_15 = y[5];

	NH4ads = y[6];
	NH4ads_15 = y[7];

	NRec = y[8];
	NRec_15 = y[9];

	NO3sto = y[10];
	NO3sto_15 = y[11];

	NO2 = y[12];
	NO2_15 = y[13];

	NO2org = y[14];
	NO2den = y[15];
	NO2nit = y[16];

	NO2org_15 = y[17];
	NO2den_15 = y[18];
	NO2nit_15 = y[19];

/*	NO2org = 0.07*NO2;
	NO2den = 0.50*NO2;
	NO2nit = 0.43*NO2;

	NO2org_15 = 0.07*NO2_15;
	NO2den_15 = 0.50*NO2_15;
	NO2nit_15 = 0.43*NO2_15;
*/
	N2O = y[20];       //EN2O
	N2O_15 = y[21];    //EN2O_15

	N2Oorg = y[22];
	N2Ocod = y[23];
	N2Onit = y[24];
	N2Oden = y[25];

	N2Oorg_15 = y[26];
	N2Ocod_15 = y[27];
	N2Onit_15 = y[28];
	N2Oden_15 = y[29];
	
	EN2Oorg = y[30];
	EN2Ocod = y[31];
	EN2Onit = y[32];
	EN2Oden = y[33];
	
	EN2Oorg_15 = y[34];
	EN2Ocod_15 = y[35];
	EN2Onit_15 = y[36];
	EN2Oden_15 = y[37];

	N2 = y[38];
	N2_15 = y[39];
/*	N2Oorg = 1.980*N2O;
	N2Ocod = 0.087*N2O;
	N2Onit = 0.230*N2O;
    N2Oden = 1.940*N2O;

	N2Oorg_15 = 1.980*N2O_15;
	N2Ocod_15 = 0.087*N2O_15;
	N2Onit_15 = 0.230*N2O_15;
    N2Oden_15 = 1.940*N2O_15;
*/
   // NH3 = 0;
   // N2O = y[12];
   // N2 = 0;
  //----full  
//NLab
	NLab_2_NH4 = Vmax_NLab_2_NH4 * NLab / (Km_NLab_2_NH4 + NLab);
	NH4_2_NLab = Vmax_NH4_2_NLab * NH4 / (Km_NH4_2_NLab + NH4);

	dNLab = NH4_2_NLab - NLab_2_NH4; //1

//NRec

	NRec_2_NH4 = C_NRec_2_NH4;
	NRec_2_NO3 = C_NRec_2_NO3;
	NO3_2_NRec = Vmax_NO3_2_NRec * NO3 / (Km_NO3_2_NRec + NO3);
	NH4_2_NRec = Vmax_NH4_2_NRec * NH4 / (Km_NH4_2_NRec + NH4);
	NRec_2_NO2org = K_NRec_2_NO2org * NRec;

//NH4ads
	NH4ads_2_NH4 = K_NH4ads_2_NH4 * NH4ads;
	NH4_2_NH4ads = K_NH4_2_NH4ads * NH4;
	dNH4ads = NH4_2_NH4ads - NH4ads_2_NH4;//4

//DNRA
	NO3_2_NH4 = NO3 * Vmax_NO3_2_NH4 / (Km_NO3_2_NH4 + NO3);
	NH4_2_NO2nit = K_NH4_2_NO2nit * NH4;

//NH4
	dNH4 = NLab_2_NH4 - NH4_2_NLab + NRec_2_NH4 - NH4_2_NRec + NH4ads_2_NH4 - NH4_2_NH4ads + NO3_2_NH4 - NH4_2_NO2nit; //3

//NO3sto
	NO3_2_NO3sto = K_NO3_2_NO3sto * NO3;
	NO3sto_2_NO3 = K_NO3sto_2_NO3 * NO3sto;
	dNO3sto = NO3_2_NO3sto - NO3sto_2_NO3;

//NO3
	NO3_2_NO2den = K_NO3_2_NO2den * NO3;
	NO2nit_2_NO3 = K_NO2nit_2_NO3 * NO2nit;

	dNO3 = NRec_2_NO3 - NO3_2_NRec  - NO3_2_NH4 + NO2nit_2_NO3 - NO3_2_NO3sto + NO3sto_2_NO3 - NO3_2_NO2den;

/* N2O */
	NO2org_2_N2Oorg = K_NO2org_2_N2Oorg * NO2org * NO2org;
	NO2nit_2_N2Onit = K_NO2nit_2_N2Onit * NO2nit * NO2nit;
	NO2den_2_N2Oden	= K_NO2den_2_N2Oden * NO2den * NO2den;
	NO2denNRec_2_N2Ocod = K_NO2denNRec_2_N2Ocod * NO2den * NRec;

//dNRec
	dNRec = NO3_2_NRec + NH4_2_NRec - NRec_2_NH4 - NRec_2_NO3 - NRec_2_NO2org - K_NO2denNRec_2_N2Ocod * NRec; //2

//NO2org

	dNO2org = NRec_2_NO2org - NO2org_2_N2Oorg;

//NO2den
	dNO2den = NO3_2_NO2den - NO2den_2_N2Oden - K_NO2denNRec_2_N2Ocod * NO2den;

//NO2nit
	dNO2nit = NH4_2_NO2nit - NO2nit_2_N2Onit - NO2nit_2_NO3;
//NO2
    dNO2 = dNO2org + dNO2den + dNO2nit;

//N2Oorg
    dEN2Oorg = K_N2Oorg_2_N2OorgE * N2Oorg;
    dEN2Onit = K_N2Onit_2_N2OnitE * N2Onit;
    dEN2Oden = K_N2Oden_2_N2OdenE * N2Oden;
    dEN2Ocod = K_N2Ocod_2_N2OcodE * N2Ocod;

    dEN2O = dEN2Oorg + dEN2Onit + dEN2Oden + dEN2Ocod;

    dRN2Oorg = K_N2Oorg_2_N2 * N2Oorg;
    dRN2Onit = K_N2Onit_2_N2 * N2Onit;
    dRN2Oden = K_N2Oden_2_N2 * N2Oden;
    dRN2Ocod = K_N2Ocod_2_N2 * N2Ocod;

    dN2 = dRN2Oorg + dRN2Onit + dRN2Oden + dRN2Ocod;

    dN2Oorg = NO2org_2_N2Oorg - dEN2Oorg - dRN2Oorg;
    dN2Onit = NO2nit_2_N2Onit - dEN2Onit - dRN2Onit;
    dN2Oden = NO2den_2_N2Oden - dEN2Oden - dRN2Oden;
    dN2Ocod = NO2denNRec_2_N2Ocod - dEN2Ocod - dRN2Ocod;


//---N15
//NLab-15
	NLab_2_NH4_15 = Vmax_NLab_2_NH4 * NLab_15 / (Km_NLab_2_NH4 + NLab_15);;
	NH4_2_NLab_15 = Vmax_NH4_2_NLab * NH4_15 / (Km_NH4_2_NLab + NH4_15);

	dNLab_15 = NH4_2_NLab_15 - NLab_2_NH4_15; //1

//NRec-15

	NRec_2_NH4_15 = C_NRec_2_NH4 * NRec_15 / NRec;
	NRec_2_NO3_15 = C_NRec_2_NO3 * NRec_15 / NRec;
	NO3_2_NRec_15 = Vmax_NO3_2_NRec * NO3_15 / (Km_NO3_2_NRec + NO3_15);
	NH4_2_NRec_15 = Vmax_NH4_2_NRec * NH4_15 / (Km_NH4_2_NRec + NH4_15);
	NRec_2_NO2org_15 = K_NRec_2_NO2org * NRec_15;

//	dNRec_15 = NO3_2_NRec_15 + NH4_2_NRec_15 - NRec_2_NH4_15 - NRec_2_NO3_15; //2

//NH4ads-15
	NH4ads_2_NH4_15 = K_NH4ads_2_NH4 * NH4ads_15;
	NH4_2_NH4ads_15 = K_NH4_2_NH4ads * NH4_15;

	dNH4ads_15 = NH4_2_NH4ads_15 - NH4ads_2_NH4_15;//4

//DNRA
	NO3_2_NH4_15 = Vmax_NO3_2_NH4 * NO3_15 / (Km_NO3_2_NH4 + NO3_15);
	NH4_2_NO2nit_15 = K_NH4_2_NO2nit * NH4_15;

//NH4-15
	dNH4_15 = NLab_2_NH4_15 - NH4_2_NLab_15 + NRec_2_NH4_15 - NH4_2_NRec_15 + NH4ads_2_NH4_15 - NH4_2_NH4ads_15 + NO3_2_NH4_15 - NH4_2_NO2nit_15; //3

//NO3sto-15
	NO3_2_NO3sto_15 = K_NO3_2_NO3sto * NO3_15;
	NO3sto_2_NO3_15 = K_NO3sto_2_NO3 * NO3sto_15;
	dNO3sto_15 = NO3_2_NO3sto_15 - NO3sto_2_NO3_15;
//NO3-15
	NO3_2_NO2den_15 = K_NO3_2_NO2den * NO3_15;
	NO2nit_2_NO3_15 = K_NO2nit_2_NO3 * NO2nit_15;

	dNO3_15 = NRec_2_NO3_15 - NO3_2_NRec_15  - NO3_2_NH4_15 + NO2nit_2_NO3_15 - NO3_2_NO3sto_15 + NO3sto_2_NO3_15 - NO3_2_NO2den_15;


/* N2O */
	NO2org_2_N2Oorg_15 = K_NO2org_2_N2Oorg * NO2org_15 * NO2org_15;
	NO2nit_2_N2Onit_15 = K_NO2nit_2_N2Onit * NO2nit_15 * NO2nit_15;
	NO2den_2_N2Oden_15	= K_NO2den_2_N2Oden * NO2den_15 * NO2den_15;
	NO2denNRec_2_N2Ocod_15 = K_NO2denNRec_2_N2Ocod * NO2den_15 * NRec_15;

//dNRec
	dNRec_15 = NO3_2_NRec_15 + NH4_2_NRec_15 - NRec_2_NH4_15 - NRec_2_NO3_15 - NRec_2_NO2org_15 - K_NO2denNRec_2_N2Ocod * NRec_15; //2

//NO2org
	dNO2org_15 = NRec_2_NO2org_15 - NO2org_2_N2Oorg_15;

//NO2den
	dNO2den_15 = NO3_2_NO2den_15 - NO2den_2_N2Oden_15 - K_NO2denNRec_2_N2Ocod * NO2den_15;

//NO2nit
	dNO2nit_15 = NH4_2_NO2nit_15 - NO2nit_2_N2Onit_15 - NO2nit_2_NO3_15;

//NO2
    dNO2_15 = dNO2org_15 + dNO2den_15 + dNO2nit_15;

//N2Oorg
    dEN2Oorg_15 = K_N2Oorg_2_N2OorgE * N2Oorg_15;
    dEN2Onit_15 = K_N2Onit_2_N2OnitE * N2Onit_15;
    dEN2Oden_15 = K_N2Oden_2_N2OdenE * N2Oden_15;
    dEN2Ocod_15 = K_N2Ocod_2_N2OcodE * N2Ocod_15;

	dEN2O_15 = dEN2Oorg_15 + dEN2Onit_15 + dEN2Oden_15 + dEN2Ocod_15;

	dRN2Oorg_15 = K_N2Oorg_2_N2 * N2Oorg_15;
	dRN2Onit_15 = K_N2Onit_2_N2 * N2Onit_15;
	dRN2Oden_15 = K_N2Oden_2_N2 * N2Oden_15;
	dRN2Ocod_15 = K_N2Ocod_2_N2 * N2Ocod_15;

	dN2_15 = dRN2Oorg_15 + dRN2Onit_15 + dRN2Oden_15 + dRN2Ocod_15;

	dN2Oorg_15 = NO2org_2_N2Oorg_15 - dEN2Oorg_15 - dRN2Oorg_15;
	dN2Onit_15 = NO2nit_2_N2Onit_15 - dEN2Onit_15 - dRN2Onit_15;
	dN2Oden_15 = NO2den_2_N2Oden_15 - dEN2Oden_15 - dRN2Oden_15;
	dN2Ocod_15 = NO2denNRec_2_N2Ocod_15 -dEN2Ocod_15 - dRN2Ocod_15;

//Output;
	ydot[0] = dNLab;
	ydot[1] = dNLab_15;
	ydot[2] = dNH4;
	ydot[3] = dNH4_15;
	ydot[4] = dNO3;
	ydot[5] = dNO3_15;

	ydot[6] = dNH4ads;
	ydot[7] = dNH4ads_15;
	ydot[8] = dNRec;
	ydot[9] = dNRec_15;
	ydot[10] = dNO3sto;
	ydot[11] = dNO3sto_15;

	ydot[12] = dNO2;
	ydot[13] = dNO2_15;

	ydot[14] = dNO2org;
	ydot[15] = dNO2den;
	ydot[16] = dNO2nit;
	ydot[17] = dNO2org_15;
	ydot[18] = dNO2den_15;
	ydot[19] = dNO2nit_15;

	ydot[20] = dEN2O;
	ydot[21] = dEN2O_15;

	ydot[22] = dN2Oorg;
	ydot[23] = dN2Ocod;
	ydot[24] = dN2Onit;
	ydot[25] = dN2Oden;
	
	ydot[26] = dN2Oorg_15;
	ydot[27] = dN2Ocod_15;
	ydot[28] = dN2Onit_15;
	ydot[29] = dN2Oden_15;
	
	ydot[30] = dEN2Oorg;
	ydot[31] = dEN2Ocod;
	ydot[32] = dEN2Onit;
	ydot[33] = dEN2Oden;

	ydot[34] = dEN2Oorg_15;
	ydot[35] = dEN2Ocod_15;
	ydot[36] = dEN2Onit_15;
	ydot[37] = dEN2Oden_15;

	ydot[38] = dN2;
	ydot[39] = dN2_15;

	yout[0] = NH4;
	yout[1] = NO3;
	yout[2] = NH4_15;
	yout[3] = NO3_15;
	yout[4] = NO2;
	yout[5] = NO2_15;
	yout[6] = N2O;
	yout[7] = N2O_15;
	yout[8] = N2;
	yout[9] = N2_15;

}
