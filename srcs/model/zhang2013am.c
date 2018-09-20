#include <R.h>

/* A typical C trick to get readable names for the parameters is #define
   This method is simple and efficient, but there are, of course,
   other possibilities that use dynamic variables.
*/

static double parms[18];

#define C_NRec_2_NH4    parms[0]
#define Vmax_NH4_2_NRec parms[1]
#define Km_NH4_2_NRec   parms[2]
#define Vmax_NH4_2_NLab parms[3]
#define Km_NH4_2_NLab   parms[4]
#define Vmax_NLab_2_NH4 parms[5]
#define Km_NLab_2_NH4   parms[6]
#define Vmax_NO3_2_NRec parms[7] 
#define Km_NO3_2_NRec   parms[8]
#define C_NRec_2_NO3    parms[9] 
#define K_NH4ads_2_NH4  parms[10] 
#define K_NH4_2_NH4ads  parms[11]
#define Vmax_NO3_2_NH4  parms[12] 
#define Km_NO3_2_NH4    parms[13]
#define Vmax_NH4_2_NO3  parms[14]
#define Km_NH4_2_NO3    parms[15]
#define K_NO3_2_NO3sto  parms[16]
#define K_NO3sto_2_NO3  parms[17]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N = 18;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  //----------Full
  double NLab_2_NH4,NH4_2_NLab,NRec_2_NH4,NRec_2_NO3,NO3_2_NRec,NH4ads_2_NH4,NO3_2_NH4,NH4_2_NO3,NO3_2_NO3sto,NO3sto_2_NO3;
  double NH4_2_NRec,NH4_2_NH4ads;
  double NLab,NH4,NO3,NH4ads,NRec,NO3sto;
  double dNLab,dNH4,dNO3,dNH4ads,dNRec,dNO3sto;  
 
  //-----------N15
  double NLab_2_NH4_15,NH4_2_NLab_15,NRec_2_NH4_15,NRec_2_NO3_15,NO3_2_NRec_15,NH4ads_2_NH4_15,NO3_2_NH4_15,NH4_2_NO3_15;
  double NH4_2_NRec_15,NH4_2_NH4ads_15,NO3_2_NO3sto_15,NO3sto_2_NO3_15;
  double NLab_15,NH4_15,NO3_15,NH4ads_15,NRec_15,NO3sto_15;
  double dNLab_15,dNH4_15,dNO3_15,dNH4ads_15,dNRec_15,dNO3sto_15;  
  //

  NLab = y[0];
  NLab_15 = y[1];
  
  NH4  = y[2];
  NH4_15 = y[3];
  
  NO3 = y[4];
  NO3_15 = y[5];
  
  NH4ads =y[6];
  NH4ads_15 =y[7];
  
  NRec = y[8];
  NRec_15 = y[9];
  
  NO3sto = y[10];
  NO3sto_15 = y[11];
    
  //N2O = y[10];
 // NH3 = 0;
 // N2O = y[12];
 // N2 = 0;
//----full  
  NLab_2_NH4 = Vmax_NLab_2_NH4*NLab/(Km_NLab_2_NH4 + NLab);
  NH4_2_NLab = Vmax_NH4_2_NLab*NH4/(Km_NH4_2_NLab + NH4);
  
  dNLab = NH4_2_NLab - NLab_2_NH4; //1
  
  NRec_2_NH4 = C_NRec_2_NH4;
  NRec_2_NO3 = C_NRec_2_NO3;
  NO3_2_NRec = Vmax_NO3_2_NRec*NO3/(Km_NO3_2_NRec+NO3);
  NH4_2_NRec = Vmax_NH4_2_NRec*NH4/(Km_NH4_2_NRec+NH4);
  dNRec = NO3_2_NRec + NH4_2_NRec - NRec_2_NH4 - NRec_2_NO3; //2
  
  NH4ads_2_NH4 = K_NH4ads_2_NH4 * NH4ads;
  NH4_2_NH4ads = K_NH4_2_NH4ads * NH4;
  
  NO3_2_NH4 =NO3* Vmax_NO3_2_NH4/(Km_NO3_2_NH4 + NO3);
  NH4_2_NO3 = Vmax_NH4_2_NO3 * NH4/(Km_NH4_2_NO3+NH4);
  
  dNH4 = NLab_2_NH4 - NH4_2_NLab + NRec_2_NH4 - NH4_2_NRec + NH4ads_2_NH4 - NH4_2_NH4ads +NO3_2_NH4 - NH4_2_NO3; //3
  
  dNH4ads = NH4_2_NH4ads - NH4ads_2_NH4;//4
  
  NO3_2_NO3sto = K_NO3_2_NO3sto*NO3;
  NO3sto_2_NO3 = K_NO3sto_2_NO3*NO3sto;
  
  dNO3 = NRec_2_NO3 - NO3_2_NRec + NH4_2_NO3 - NO3_2_NH4 - NO3_2_NO3sto + NO3sto_2_NO3;//5
  dNO3sto = NO3_2_NO3sto - NO3sto_2_NO3;
  
  //---N15
  NLab_2_NH4_15 = Vmax_NLab_2_NH4*NLab_15/(Km_NLab_2_NH4 + NLab_15);;
  NH4_2_NLab_15 = Vmax_NH4_2_NLab*NH4_15/(Km_NH4_2_NLab + NH4_15);
  
  dNLab_15 = NH4_2_NLab_15 - NLab_2_NH4_15; //1
  
  NRec_2_NH4_15 = C_NRec_2_NH4*NRec_15/NRec;
  NRec_2_NO3_15 = C_NRec_2_NO3*NRec_15/NRec;
  NO3_2_NRec_15 = Vmax_NO3_2_NRec*NO3_15/(Km_NO3_2_NRec+NO3_15);
  NH4_2_NRec_15 = Vmax_NH4_2_NRec*NH4_15/(Km_NH4_2_NRec+NH4_15);
  
  dNRec_15 = NO3_2_NRec_15 + NH4_2_NRec_15 - NRec_2_NH4_15 - NRec_2_NO3_15; //2
  
  NH4ads_2_NH4_15 = K_NH4ads_2_NH4 * NH4ads_15;
  NH4_2_NH4ads_15 = K_NH4_2_NH4ads * NH4_15;
  
  NO3_2_NO3sto_15 = K_NO3_2_NO3sto * NO3_15;
  NO3sto_2_NO3_15 = K_NO3sto_2_NO3 * NO3sto_15;
  
  NO3_2_NH4_15 = Vmax_NO3_2_NH4*NO3_15/(Km_NO3_2_NH4 + NO3_15);
  NH4_2_NO3_15 = Vmax_NH4_2_NO3 * NH4_15/(Km_NH4_2_NO3+NH4_15);
  
  dNH4_15 = NLab_2_NH4_15 - NH4_2_NLab_15 + NRec_2_NH4_15 - NH4_2_NRec_15 + NH4ads_2_NH4_15 - NH4_2_NH4ads_15 +NO3_2_NH4_15 - NH4_2_NO3_15; //3
  
  dNH4ads_15 = NH4_2_NH4ads_15 - NH4ads_2_NH4_15;//4
  
  dNO3_15 = NRec_2_NO3_15 + NO3sto_2_NO3_15 - NO3_2_NO3sto_15 - NO3_2_NRec_15 + NH4_2_NO3_15 - NO3_2_NH4_15;//5
  dNO3sto_15 = NO3_2_NO3sto_15 - NO3sto_2_NO3_15;
 
  //	 double dNLab,dNH4,dNO3,dNH4ads,dNRec,dNO2nit,dNO2den,dNO2org,dN2O,dN2Oorg,dN2Ocod,dN2Oden,dN2Onit,dN2; 	
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
  ydot[10]= dNO3sto;
  ydot[11]= dNO3sto_15;
  
  yout[0] = NH4;
  yout[1] = NO3;
  yout[2] = NH4_15;
  yout[3] = NO3_15;
//  yout[4] = N2O;
}
