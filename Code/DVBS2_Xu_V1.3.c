
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dvbs2.h"

// Gaussian Random bit generator
double gaussianRandom(void);
int DVBS2DECODER(int *, double *);

#define Ndata 7200
#define Ncode 16200
double a;
double In_PR1[16200];
double In_PR0[16200];
double In_PdtR2[9000];

int main(void){
	double SNR = 0.2;    //Gives 0 BER DVBS2
	int Npar = Ncode-Ndata;  //parity check bit length
	int iter;

	int c[Ncode] = {0};
	memset(c,0,Ncode*(sizeof(int)));

	double z[Ncode] = {0};

	int max_sim = 100; //Maximum number of simulations

	printf("\n\n\n\nCLocks per sec = %d\n\n", (int) CLOCKS_PER_SEC);
	int i,j;
	int ret;

	for(i=0;i<Ncode;i++){
		In_PR1[i]=1;
		In_PR0[i]=1;
		if(i<Npar)
		In_PdtR2[i]=1;
	}

	int sim=0;
	double r_xb[Ndata] = {0};
	int xb[Ndata];

	int ** Addr;
	FILE *fp1 = fopen("Addr.txt","r");
	Addr = malloc(Ndata*sizeof(int*));
	for (i = 0; i< Ndata; i++){
	  Addr[i] = malloc(8*sizeof(int));
	}
	for (i = 0; i < Ndata; i++) {
	  for (j = 0; j < 8; j++) {
	    ret = fscanf(fp1, "%d", &Addr[i][j]);
	  }
	}
	fclose(fp1);
	
	clock_t startx = clock();
	while(sim<max_sim){
		srand((unsigned) time(NULL));
		memset(r_xb, 0, sizeof r_xb);

		for (i = 0; i < Ndata; i++){
		  r_xb[i]=gaussianRandom();
		}
		
		for(i=0;i<Ndata;i++){
		  if(r_xb[i]>0){xb[i]=1;}
		  else if(r_xb[i]<=0){xb[i]=-1;}
		}
		
		int x[Ndata];
		
		for(i=0;i<Ndata;i++){
		  if(xb[i]>0){x[i]=1;}
		  else if(xb[i]<=0){x[i]=0;}
		}
		
		
		int par[Npar];
		for (i=0;i<Npar;i++){
		  par[i]=0;
		}
		
		for (i=0;i<360*5;i++){
		  for (j=0;j<8;j++){
		    par[Addr[i][j]-1] = (par[Addr[i][j]-1]+x[i])%2;
		  }
		}
		for (i=360*5;i<360*20;i++){
		  for (j=0;j<3;j++){
		    par[Addr[i][j]-1] = (par[Addr[i][j]-1]+x[i])%2;
		  }
		}
		
		int Xx[Ncode];
		for (i=0;i<Npar;i++){
		  if(i>0){
		    par[i]=(par[i]+par[i-1])%2;
		  }
		  Xx[i+Ndata]=par[i];
		}
		
		for (i=0;i<Ndata;i++){
		  Xx[i]=x[i];
		}
	
		
		memset(z,0,Ncode*(sizeof(double)));
		for(i=0;i<Ncode;i++){
		  z[i] = 2 * Xx[i] -1;
		}
		
		double noise[Ncode];
		for (i = 0; i < Ncode; i++){
		  noise[i]=gaussianRandom();
		}
		
		double Npower=0;
		for (i=0;i<Ncode;i++)Npower = Npower + noise[i]*noise[i];
		
		double Spower=0;
		for (i=0;i<Ncode;i++)Spower = Spower + z[i]*z[i];
		
		for (i=0;i<Ncode;i++){
		  noise[i]=noise[i]*pow(Spower/Npower,0.5)*pow(10,-SNR/20);
		}
		
		Npower = 0;
		for (i=0;i<Ncode;i++) Npower = Npower + noise[i]*noise[i];
		
		SNR = Spower/Npower;
		a=sqrt(SNR);	//a = sqare root of (2*Rate*SNR), where Rate = 0.5
		SNR = 10*log10(SNR);
		
		for(i=0;i<Ncode;i++){
		  z[i]=z[i]+noise[i];
		}
		
		//Reinvented Gallager’s Decoder
		iter = DVBS2DECODER(c, z);
		
		sim++;
	}
	
	clock_t endx = clock();

	printf("Number of iterations, sim = %d\t max_sim  = %d\n", sim, max_sim);
	printf("CLOCKS_PER_SEC = %f\n", (double) CLOCKS_PER_SEC);
	double run_time = ((double) (endx-startx)) / (double) CLOCKS_PER_SEC;
	
	printf("BPSK total run time - %f (ms), %f (s) elapsed.\n",run_time * 1000 , run_time);
	printf("Average time elapsed per simulation = %f (ms)\n",run_time/(double) sim * 1000);
	printf("Average time elapsed per iteration of decoding = %f (ms)\n",run_time/(double) sim * 1000/iter);
	printf("Data rate = %f (Msps)\n", (double) Ndata / run_time * (double) sim /1e6);
	
	int x_hat[Ndata];
	for(i=0;i<Ndata;i++){
	  x_hat[i]=c[i];
	}
	
	int sign_z[Ndata];
	
	for (i=0;i<Ndata;i++)sign_z[i] = (z[i]>0) ? 1:(-1);
	
	double BER = 0;
	for (i=0;i<Ndata;i++)BER = BER + abs(xb[i]-sign_z[i])/2;
	
	double BER_RAW;
	double BERDVBS2;

	BER_RAW = BER/Ndata; //Raw BER of transmitted codeword
		
	double BER_x_hat=0;   //The number of error bits in the decoded codeword
	for(i=0;i<Ndata;i++) BER_x_hat = BER_x_hat + abs(xb[i]-2*x_hat[i]+1)/2;
	BERDVBS2 = BER_x_hat/Ndata; // Bit error rate of the decoded codeword
	
	printf("SNR = %g\n",SNR);
	printf("BER_RAW = %f\n",BER_RAW);
	printf("BERDVBs2 = %f \n",BERDVBS2);   //Displaying BER
	printf("Number of decoder iterations = %d\n",iter);
	printf("\n");
	
	return 0;


}

double gaussianRandom(void){
	double v1, v2, s;

	do {
		v1 = 2 * ((double) rand() / RAND_MAX)-1;
		v2 = 2 * ((double) rand() / RAND_MAX)-1;
		s = v1 * v1 + v2 * v2;
	} while (s>=1||s==0);
	s = sqrt((-2*log(s))/s);
	return v1*s;
}


int DVBS2DECODER(int *c, double *z){

	unsigned int Npar = Ncode-Ndata;
	unsigned int i;
	unsigned int no_nonzero = 48599;
	int iter=0;
	unsigned int success = 0;
	unsigned int max_iter=150;

	double sumQ;
	double sumq;
	double F1[Ncode],F0[Ncode];
	double Q1[no_nonzero],Q0[no_nonzero],dtQ[no_nonzero],q1[Ncode],q0[Ncode];
	double R1[no_nonzero],R0[no_nonzero],Updated_dtQ[no_nonzero],Prod_dtQ[Npar];
	double prod_R1[no_nonzero],prod_R0[no_nonzero],Prod_R1[Ncode],Prod_R0[Ncode];
	unsigned int ss[Npar];
	unsigned int s[Npar];
	unsigned int fail;
	for (i=0;i<Ncode;i++){
		F1[i]=1/(1+exp((-2)*a*z[i]));
		F0[i]=1-F1[i];
	}

	//Initialization
	memset(Q0,0x00,sizeof(Q0));
	memset(Q1,0x00,sizeof(Q1));
	for(i=0;i<no_nonzero;i++){
		Q1[i]=F1[cw_colindex[i]];
		Q0[i]=F0[cw_colindex[i]];
	}

	while ((success ==0) & (iter<max_iter)){
		iter++;
		fail = 0;

		//Horizontal step
		memset(R0,0x00,sizeof(R0));
		memset(R1,0x00,sizeof(R1));
		memset(q0,0x00,sizeof(q0));
		memset(q1,0x00,sizeof(q1));
		memset(dtQ,0x00,sizeof(dtQ));
	
		memcpy(Prod_R1,In_PR1,sizeof(In_PR1));
		memcpy(Prod_R0,In_PR0,sizeof(In_PR0));
		memcpy(Prod_dtQ,In_PdtR2,sizeof(In_PdtR2));

		for (i=0;i<no_nonzero;i++)
		{
			dtQ[i]=Q0[i]-Q1[i];
			Prod_dtQ[cw_rowindex[i]]=Prod_dtQ[cw_rowindex[i]]*dtQ[i];
		}
		
		//Horizontal step
		for (i=0;i<no_nonzero;i++)
		{	
			Updated_dtQ[i]=Prod_dtQ[cw_rowindex[i]]/dtQ[i];
			R0[i]=0.5*(1+Updated_dtQ[i]);
			R1[i]=1-R0[i];	//0.5-0.5*Updated_dtQ

			Prod_R1[cw_colindex[i]]=Prod_R1[cw_colindex[i]]*R1[i];
			Prod_R0[cw_colindex[i]]=Prod_R0[cw_colindex[i]]*R0[i];
		}

		//Vertical step
		for (i=0;i<no_nonzero;i++)
		{
			sumQ=0;
			prod_R1[i]=Prod_R1[cw_colindex[i]]/R1[i];
			Q1[i]=F1[cw_colindex[i]]*prod_R1[i];
			prod_R0[i]=Prod_R0[cw_colindex[i]]/R0[i];
			Q0[i]=F0[cw_colindex[i]]*prod_R0[i];	
			sumQ=Q1[i]+Q0[i];
			Q1[i]=Q1[i]/sumQ;
			Q0[i]=1-Q1[i];
			if(i<Ncode)
			{
				sumq=0;
				q0[i]=F0[i]*Prod_R0[i];
				q1[i]=F1[i]*Prod_R1[i];
				sumq=q1[i]+q0[i];
				q1[i]=q1[i]/sumq;
				q0[i]=1-q1[i];
				c[i]=0;
				if(q0[i]<=q1[i]){c[i]=1;}
			}
			ss[cw_rowindex[i]]=ss[cw_rowindex[i]]+(unsigned int)c[cw_colindex[i]];
		}

		for(i=0;i<Npar;i++){
			s[i]=ss[i]%2;
			ss[i]=0;
			if(s[i]==1){fail++;}
		}
		if (fail ==0){success = 1;}
	}
	return iter;
}
