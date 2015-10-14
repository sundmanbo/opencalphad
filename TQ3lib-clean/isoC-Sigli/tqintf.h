//define PARALLEL 0 no loop with parallelization
//define PARALLEL 1 declared loops with parallelization
#define PARALLEL 1
#define NCPU 12
#define MAXEL 20
#define MAXPH 400
#define PHFIXED 2
#define PHENTERED 0
#define PHSUS -3
#define GRID 0
#define NOGRID -1
#define TCtoTK 273.15



#include "octqc.h"
#include <string>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <sstream>
#include <omp.h>
#include <string> 


extern"C"
{
    void c_tqini(int, void *);                                                  // initiates the OC package
    void c_tqrfil(char *, void *);                                              // read all elements from a TDB file
    //void c_tqgcom(int *, char[MAXEL][24], void **);                           // get system component names. At present the elements
    void c_tqrpfil(char *, int, char **, void *);                               // read TDB file with selection of elements
    //void c_tqgnp(int *, void **);                                             // get total number of phases and composition sets
    void c_tqgpi(int *, char *, void *);                                        // get index of phase phasename
	void c_tqgpn(int, char *, void *);                                          // get name of phase+compset tuple with index phcsx
    void c_tqgetv(char *, int, int, int *, double *, void *);                   // get equilibrium results using state variables
    void c_tqsetc(char *, int, int, double, int *, void *);                     // set condition
    void c_tqce(char *, int, int, double *, void *);                            // calculate quilibrium with possible target
    //void c_tqgnp(int, gtp_equilibrium_data **);                               // get total number of phases and composition sets
    void examine_gtp_equilibrium_data(void *);                                  //
    //void c_getG(int, void *);
    //void c_calcg(int, int, int, int, void *);
    void c_tqgphc1(int, int * , int *, int *, double *, double *, double *,
                                                                        void *);
    void c_tqsphc1(int, double *, double *, void *);
    void c_tqcph1(int, int, int *, double *, double *, double *, double *, double *, void *);
	void c_reset_condition_T( void *);
	void c_Change_Status_Phase(char *, int,double,void *);
	void c_List_Conditions(void *);
	void c_checktdb(char *);
	void c_newEquilibrium(char *,int *);
	void c_selecteq(int ,void *);
    void c_copy_equilibrium(void *,char *,void *);
	void c_set_status_globaldata();
	int c_errors_number();
	void c_new_gtp();
}

extern"C" int  c_ntup;                                                          //
extern"C" int  c_nel;                                                           // number of elements
extern"C" int  c_maxc;                                                          //
extern"C" char *c_cnam[MAXEL];                                                     // character array with all element names
extern"C" double c_gval[24];
extern"C" int c_noofcs(int);

using namespace std;
void Get_Ceq(const int iceq,void *ceq){
	c_selecteq(iceq,ceq);
	//cout << "-> Adress of ceq-Storage: [" << ceq << "]" <<endl;
}
void Initialize(void *ceq)
{
   int n = 0;
	void *ceq2 = NULL;
    //===============
    c_tqini(n, ceq);
    //===============

   //cout << "-> Adress of ceq-Storage: [" << ceq << "]" <<endl;
   
   
};

int Create_New_Ceq_and_Return_ID(const string &Ceq_Name){
	int ieq;
	char *filename = strcpy((char*)malloc(Ceq_Name.length()+1), Ceq_Name.c_str());
	c_newEquilibrium(filename,&ieq);
	return ieq;
}
void Get_Ceq_pointer(const int &ieq, void *ceq){
	c_selecteq(ieq,&ceq);
	
}


void GetAllElementsFromDatabase(string tdbfilename){
	 char *filename = strcpy((char*)malloc(tdbfilename.length()+1), tdbfilename.c_str());
	 c_checktdb(filename);
}

void ReadDatabase(string tdbfilename, void *ceq)
{
    char *filename = strcpy((char*)malloc(tdbfilename.length()+1), tdbfilename.c_str());

    //======================
    c_tqrfil(filename, ceq);
    //======================

    /*cout << "-> Element Data: [";
    for(int i = 0; i < c_nel; i++)
    {
       cout << c_cnam[i];
        if(i < c_nel-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<endl;
	*/
};

void ReadDatabaseLimited(string tdbfilename, vector<string> elnames, void *ceq)
{
    char *filename = strcpy((char*)malloc(tdbfilename.length()+1), tdbfilename.c_str());
    char *selel[elnames.size()];
    for(size_t i = 0; i < elnames.size(); i++)
    {
        char *tempchar
             = strcpy((char*)malloc(elnames[i].length()+1), elnames[i].c_str());
        selel[i] = tempchar;
    }

    //==============================================
    c_tqrpfil(filename, elnames.size(), selel, ceq);
    //==============================================
/*
    cout << "-> Element Data: [";
    for(int i = 0; i < c_nel; i++)
    {
        cout << c_cnam[i];
        if(i < c_nel-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;
*/

};

void ReadPhases(vector<string> &phnames, void *ceq)
{
    phnames.clear();
	phnames.resize(c_ntup);
	

    for(int i = 1; i < c_ntup+1; i++)
    {
        char phn[24];

        //==========================
        c_tqgpn(i, phn, ceq);
        //==========================
		
        
		int index;
		c_tqgpi(&index,phn,ceq);
		phnames[index-1]=phn;
    }
/*
    cout << "-> Phase Data: [";
    for(size_t i = 0; i < phnames.size(); i++)
    {
        cout << i<< " "<<phnames[i];
        if(i < phnames.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << " [" << &ceq << "]" <<
    endl;
*/	
};
void ResetTemperature(void *ceq){
	c_reset_condition_T(ceq);
}
void Change_Phase_Status(string name,int nystat,double val,void *ceq){
//nystat=0 :Entered
//nystat=2 :Fixed
	char *phasename = strcpy((char*)malloc(name.length()+1), name.c_str());
	c_Change_Status_Phase(phasename,nystat,val,ceq);
	
}
void SetTemperature(const double &T, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "T";
  //  if (T < 1.0) T = 1.0;

    //=========================================
    c_tqsetc(par, n1, n2, T, &cnum, ceq);
    //=========================================

   // cout << "-> Set Temperature to: [" << T << "]" << " [" << &ceq << "]" <<
   // endl;
};

void SetPressure(const double &P, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "P";
   // if (P < 1.0) P = 1.0;

    //=========================================
    c_tqsetc(par, n1, n2, P, &cnum, ceq);
    //=========================================

//    cout << "-> Set Pressure to: [" << P << "]" << " [" << &ceq << "]" <<
//    endl;
};

void SetMoles(double N, void *ceq)
{
    int cnum;
    int n1 = 0;
    int n2 = 0;
    char par[60] = "N";

    //=========================================
    c_tqsetc(par, n1, n2, N, &cnum, ceq);
    //=========================================

 //   cout << "-> Set Moles to: [" << N << "]" << " [" << &ceq << "]" <<
 //   endl;
};

void SetComposition(vector<double>& X, void *ceq, const int &i_ref,string &compo_unit)
{
    int cnum;

    int n2 = 0;

    
	
    char par[60];
    strcpy(par,compo_unit.c_str());
    
    for (size_t i = 0; i < c_nel; i++)
    {
       if (X[i] < 1.0e-8) X[i] = 1.0e-8;                                       // Check and fix, if composition is below treshold

        if(not (i == i_ref))
        {                                                                       // Set and print composition, if element 'i' is not the reference/(last) element
            //==================================================
            c_tqsetc(par, i+1, n2, X[i], &cnum, ceq);
            //==================================================

 //           cout << "-> Set Composition of " << c_cnam[i] << " to: [" <<
 //                        X[i] << "]" << " [" << &ceq << "]" <<
 //           endl;
        }
        else
        {                                                                       // Print composition, if element 'i' is the reference/(last) element
           double X_ref = 1;
            for(size_t j = 0; j < i; j++)
            {
                X_ref -= X[j];
            }

//            cout << "-> Set Composition of " << c_cnam[i] << " to: [" <<
//                         X_ref << "]" << " [" << &ceq << "]" <<
//            endl;
        }
    }
};

void SetConstituents(int phidx, vector<double> y, void *ceq)
{
    int stable1 = phidx;
    double extra[MAXPH];
    double yfr[y.size()];
    for(size_t i = 0; i < y.size(); i++)
    {
        yfr[i] = y[i];
    }

    //===============================
    c_tqsphc1(stable1,yfr,extra,ceq);
    //===============================

    cout << "-> Set Constituents to: [";
    for(int i = 0; i < y.size(); i++)
    {
        cout << i << ": " << yfr[i];
        if(i < y.size()-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;
};


void SelectSinglePhase(int PhIdx, void *ceq)
{
    //
};
void List_Conditions(void *ceq){
	c_List_Conditions(ceq);
}
void CalculateEquilibrium(void *ceq, int n1, int &i_error)
{
	
	
	
	i_error=0;
    char target[60] = " ";
   
    int n2 = 0;
    double val;
	int iter=0;
	do{
		if (not (i_error==0)) {
			cout<<" !!!! Equilibrium not converged & trying again iter="<<iter<<endl;
//			c_List_Conditions(ceq);
			
		}
		
		iter+=1;
		
		//======================================
		c_tqce(target, n1, n2, &val, ceq);
		//======================================
		i_error=c_errors_number();
	}while((not(i_error==0))and(iter<1));
   
};

void GetGibbsData(int phidx, void *ceq)
{
    int n2 = 2;
    int n3;
    double gtp[6];
    double dgdy[100];
    double d2gdydt[100];
    double d2gdydp[100];
    double d2gdy2[100];

    //=================================================================
    c_tqcph1(phidx, n2, &n3, gtp, dgdy, d2gdydt, d2gdydp, d2gdy2, ceq);
    //=================================================================

    cout << "-> Read Gibbs Data G: [";
    for(int i = 0; i < 6; i++)
    {
        cout << gtp[i];
        if(i < 5)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data dGdY: [";
    for(int i = 0; i < n3; i++)
    {
        cout << dgdy[i];
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data d2GdYdT: [";
    for(int i = 0; i < n3; i++)
    {
        cout << d2gdydt[i];
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    cout << "-> Read Gibbs Data d2GdYdP: [";
    for(int i = 0; i < n3; i++)
    {
        cout << d2gdydp[i];
        if(i < n3-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;

    int kk=n2*(n2+1)/2;

    cout << "-> Read Gibbs Data d2GdY2: [";
    for(int i = 0; i < kk; i++)
    {
        cout << d2gdy2[i];
        if(i < kk-1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;
};

void ReadPhaseFractions(vector<string> phnames, vector<double>& phfract,
                                                                      void *ceq)
{
    double npf[MAXPH];
    char statevar[60] = "NP";
    int n1 =  -1;//-1
    int n2 =  0;
    int n3 = MAXPH;//sizeof(npf) / sizeof(npf[0]);

    //========================================
    c_tqgetv(statevar, n1, n2, &n3, npf, ceq);
    //========================================
	
    for(int i = 0; i < phnames.size(); i++){
	/*
	char phn[24];
	c_tqgpn(i+1, phn, ceq);
	size_t index=0;
		for (size_t j=0;j<phnames.size();j++){
			
			if  (phnames[j]==phn){
				index=j;
				break;
			}
		}
		*/
	phfract[i]=npf[i];
    //phfract[index]=npf[i];
	//cout<<i<<" " <<phnames[i]<<" : "<<  phfract[i] <<endl;
	}

    
};
double ReadTemperature(void *ceq)
{
    double npf[1];
    char statevar[60] = "T";
    int n1 = 0;
    int n2 = 0;
    int n3 = 1;
	double TK;

    //========================================
    c_tqgetv(statevar, n1, n2, &n3, npf, ceq);
    //========================================

    
    TK=npf[0];
	return(TK);
};
void ReadConstituentFractions(vector<string> phnames, vector<double> phfract,
                                     vector< vector<double> > &elfract, void *ceq)
{
   
    
	
    double pxf[10*MAXPH];
    for (int i = 1; i < c_ntup+1; i++)
    {
		char phn[24];
		c_tqgpn(i, phn, ceq);
		size_t index=0;
		for (size_t j=0;j<phnames.size();j++){
			
			if  (phnames[j]==phn){
				index=j;
				break;
			}
		}
        if (phfract[index] > 1e-10)
        {
            char statevar[60] =  "X";

            int n2 = -1;                                                        //composition of stable phase n2 = -1 means all fractions
            int n4 = sizeof(pxf)/sizeof(pxf[0]);

            //=======================================
            c_tqgetv(statevar, i, n2, &n4, pxf, ceq);
            //=======================================
			
            for (int k = 0; k < n4; k++)
            {
			
                elfract[index][k]=pxf[k];
				
            }
			
			
           // cout << "-> Constituent Fractions for " << phnames[i-1] <<" [";
			//cout << "-> Constituent Fractions for " << phnames[index]<<" [";
            for (int k = 0; k < n4; k++)
            {
                //cout << c_cnam[k] << ": " << elfract[index][k];
                if(k < n4-1)
                {
                   // cout << ", ";
                }
            }
            //cout << "]" << " [" << &ceq << "]" <<endl;
        }
    }
};

void ListExtConstituentFractions(int phidx, vector<string> phnames, void *ceq)
{
    int stable1 = phidx;
    int nlat;
    int nlatc[MAXPH];
    int conlista[MAXPH];
    double yfr[MAXPH];
    double sites[MAXPH];
    double extra[MAXPH];

    //======================================================================
    c_tqgphc1(stable1, &nlat, nlatc, conlista, yfr, sites, extra, ceq);
    //======================================================================

    cout << "-> Extended Constituent Fractions for " << phnames[stable1-1]
         << " [" << extra[0] << " moles of atoms/formula unit]";
    int consti = 0;
    for(int i = 0; i < nlat; i++)
    {
        cout << " [";
        for(int j = 0; j < nlatc[i]; j++)
        {
            cout << "Const. " << consti << ": " << yfr[consti];
            if(j < nlatc[i]-1)
            {
                cout << ", ";
            }
            consti += 1;
        }
        cout << "]_(" << sites[i] << ")";
    }
    cout << endl;
};

std::string IntToString ( int number )
{
	std::string mystr;
	std::stringstream out;
	out << number;
	mystr = out.str();
  return mystr;
}
// Write the results of a given equilibrium
// el_reduced_names: vector of names elements with non zero composition
// phnames: vector of names phases that can appear for these elements
// phfract: atomic fraction of these phases after equilibrium
// elfract[i][j]: atomic composition of element i in phase j
// ceqh: pointer for the given equilibrium calculation
// mode: 1 write only atomic fractions of phases after equilibrium
// mode: 1 write atomic fractions + compositions of phases after equilibrium

void Write_Results_Equilibrium(const vector<string> &el_reduced_names, const vector<string> &phnames, vector<double> &phfract, 
						  vector< vector<double> > &elfract, void *ceqh,const int &mode){
	
	//-------------------------------List Results-------------------------------
	
	ReadPhaseFractions(phnames, phfract, &ceqh);                                 // Read the amount of stable phases
	
	if (mode >1)  ReadConstituentFractions(phnames, phfract, elfract, &ceqh);                  // Read the composition of each stable phase
	
	double TC=ReadTemperature(&ceqh)-TCtoTK;
	
	cout<<endl;
	cout<<" ===================================== "<<endl;
	cout<<"    New Equilibrium at : "<<TC<<" C"<<endl;
	List_Conditions(&ceqh);
	cout<<" --------------------------------------- "<<endl;	
	for (size_t i=0; i<phnames.size(); i++){
		if (phfract[i]>0){
			if (mode >1) cout<<" --------------------------------------- "<<endl;
			cout<<"         "<<phnames[i]<<" fat%= "<<phfract[i]*100<<endl;
			if (mode >1) {
				
				cout<<" --------------------------------------- "<<endl;
				
				
				for (size_t j=0; j<el_reduced_names.size();j++){
					if (elfract[i][j]>1e-10) cout<<"        "<<el_reduced_names[j]<<" = "<<elfract[i][j]*100<<" (at%)"<<endl;
				}
			}
			
		}
	}
	
}

// ***************************************************************************************************************
// find all the transitions temperatures for a given alloy composition and accuracy
// if you want to run the program with parallelization 
// you need to declare 	bool parallel =true; 
// you need to uncomment: #pragma omp parallel for
// if you want to run the program without parallelization PARALLEL is set to 0 (see top of this file)
// if no parallelization the standart equilibrium pointer is used and we do not enter new equilmibria to save timme
// if parallelization (here on 10 equilibria) we need to enter 10 new equilibria
// this is performed with the 3 commands:
// Ceq_Name=root+IntToString(i); in order to have a different name for each equilibrium
// iceq=Create_New_Ceq_and_Return_ID(Ceq_Name); iceq is the index in the equilibrium vector eqlista of OC3  
// Store_Equilibria.push_back(iceq); all the indexes are stored in the vector Store_Equilibria
// here you scan the temperature and we create a vector of the different temperatures that will be used in the parallel calculation
void Find_Transitions( const double &TK_start,const int &nstep,const double &step_TK,vector<double> &W, const vector<string> &phnames,vector<double> &Transitions,const vector<string> &el_reduced_names,const bool first_iteration,  const bool last_iteration, vector<int> &Store_Equilibria,vector< string > &Phase_transitions_mixture, void *ceq,const double required_accuracy_on_TK){				  
	
	int iceq=0;
	vector<double> phfract;
	phfract.resize(phnames.size(),0.);
	
	vector< vector<double> > elfract;                                           // Array including all equilibrium compositions
	elfract.resize(phnames.size(),vector<double>(el_reduced_names.size(),0.));
	
	
	
	
	
	vector<double> TKCE;
	vector< vector<double> > CeqFract;
	
	TKCE.resize(0);
	CeqFract.resize(0);
	double TK_end=TK_start+(nstep-1)*step_TK;
	double TK=TK_start;
	int nstep_total=nstep;
	if (not first_iteration) nstep_total+=1;
	for (int i=0; i<nstep_total;i++){
		TKCE.push_back(TK);//here you scan the temperature and we create a vector of the different temperatures that will be used in the parallel calculation
		TK+=step_TK;
	}
	CeqFract.resize(TKCE.size(),vector<double>(phnames.size(),0.));
   
 
	size_t max_number_of_phase=0;
// the three lines below trigger parallelism for the nex for {....} loop if PARALLEL is not 0
    omp_set_num_threads(NCPU);
//	cout<<"number of threads detected:"<<omp_get_num_procs()<<endl;
#if PARALLEL>0
#pragma omp parallel for default(none), schedule(dynamic), shared(TKCE,W,Store_Equilibria, el_reduced_names,phnames,phfract,CeqFract)
#endif	    
	for (int i=0; i<TKCE.size();i++){
		
		void *ceqi= NULL;
		if ((PARALLEL>0)) {
		
			c_selecteq(Store_Equilibria[i], &ceqi);// retrieve the pointer with index stored in Store_Equilibria
			
		}else{
	//		ceqi=ceq;// if no parallelization use STANDART EQUILIBRIUM
			c_selecteq(1, &ceqi);
		}
		for (int k=0;k<phnames.size();k++) Change_Phase_Status(phnames[k],PHENTERED,0.0,&ceqi);
		Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceqi);// 
 //		Change_Phase_Status("FCC_A1",PHENTERED,0.1,&ceqi);//
//		cout<<"T="<<TKCE[i]<<endl;
		SetTemperature(TKCE[i], &ceqi); // set temperature for specific equilibrium
//		List_Conditions(&ceqi);
        int i_error=0;
//		List_Conditions(&ceqi);
		if (PARALLEL>0){
			CalculateEquilibrium(&ceqi,NOGRID,i_error); // perform an equilibrium calculation
		}else{
			CalculateEquilibrium(&ceqi,NOGRID,i_error); // perform an equilibrium calculation
		}
		
		if (not(i_error==0)){
			/*cout<<" equilibrium calculation not converged in transition subroutine for the following conditions"<<endl;
			cout<<" TK="<<TKCE[i]<<" "<< ReadTemperature(&ceqi)<<endl;
			cout<<" composition:"<<endl;
			for (size_t i=0;i<el_reduced_names.size();i++) {
				cout<<el_reduced_names[i]<<" (w%): "<<W[i]<<endl;
			}
			*/
			
			exit(EXIT_FAILURE);
		}
		ReadPhaseFractions(phnames, phfract, &ceqi);// get the phase fraction of all phases
		for (size_t j=0; j<phnames.size(); j++){
			if (phfract[j]>0) CeqFract[i][j]=phfract[j];
		}
                                                                      
	}
/*	
	for (int i=0; i<TKCE.size();i++){
		cout<<i<<" ["<<TK_start<<","<<TK_end<<"]  ---->"<<TKCE[i]<<endl;
	}
*/	
	
	
	// analyse the results of ech equilibrium stored in CeqFract[i][j] is the index of the equilibrium J the index of the phase
	for (int i=0; i<TKCE.size()-1;i++){
		/*cout<<i<<" ["<<TK_start<<","<<TK_end<<"]  ---->"<<TKCE[i]<<endl;
		
		for (size_t j=0; j<phnames.size(); j++){
			if (CeqFract[i][j]>0) cout<<" "<<phnames[j];
		}
		cout<<endl;
		*/
		for (size_t j=0; j<phnames.size(); j++){
			
			if ((!(CeqFract[i][j]<1e-8)&&(CeqFract[i+1][j]<1e-8))||(!(CeqFract[i][j]>1e-8)&&(CeqFract[i+1][j]>1e-8)))
			{
  //            a transition has been detected			
  //			cout<<"********transition at: "<<TKCE[i]<<endl;
 //				cout<<"phase:"<<phnames[j]<<" "<<CeqFract[i][j]<<" "<<CeqFract[i+1][j]<<endl;
				int i_value;	
				if (! last_iteration) {
					Transitions.push_back(TKCE[i]);
				}else
				{
					if (Transitions.size()==0) {
						Transitions.push_back(TKCE[i]);
						Phase_transitions_mixture.push_back("");
						for (size_t k=0; k<phnames.size(); k++){
							if (CeqFract[i][k]>0) {
								Phase_transitions_mixture.back()+=phnames[k];
								Phase_transitions_mixture.back()+=" + ";
							}
						}
					}
					
					Transitions.push_back(TKCE[i+1]);
					Phase_transitions_mixture.push_back("");
					for (size_t k=0; k<phnames.size(); k++){
						if (CeqFract[i+1][k]>0) {
							Phase_transitions_mixture.back()+=phnames[k];
							Phase_transitions_mixture.back()+=" + ";
						}
					}
				}
				
				j=phnames.size()+1;// exit the loop 
			}
		}
	}
}

// ***************************************************************************************************************
// find all the transitions temperatures for a given alloy composition
// step_TK: first interval of temperature used
// n_step : number of steps set to NSTEP
// between TK_start and TK_end=TK_start+(n_step-1)*step_TK;
// required_accuracy_on_TK: self explanatory
// W : weight composition of elements
// phnames: vector of names phases that can appear for these elements
// el_reduced_names: vector of names elements with non zero composition
// ceq: pointer for the given equilibrium calculation used to pass the standart equilibrium in non parallel computation
// parallelization option is in Find_Transitions
// see Find_Transitions for comments on parallelization

void Global_Find_Transitions(double &TK_start,const int &n_step,double &TK_end,const double required_accuracy_on_TK, vector<double> &W, const vector<string> &phnames,const vector<string> &el_reduced_names, void *ceq, const int &i_ref, const string &compo_unit){				  
	string mycompo_unit=compo_unit;
	double TK_end_ini, TK_start_ini;
	TK_start_ini=TK_start;
	TK_end_ini=TK_end;
//	cout<<"sntep="<<n_step<<endl;
	double step_TK=(TK_end-TK_start)/(double)(n_step-1);
	double old_step_TK=TK_end-TK_start;
	int number_of_loops=(int)( log10( fabs(step_TK)/required_accuracy_on_TK)/log10(n_step)+1);

	cout<<"number of loops"<<number_of_loops<<endl;
	
	vector<double> phfract;
	vector<double> Transitions1;
	vector<double> Transitions0;
	vector<int> Store_Equilibria;
	Store_Equilibria.resize(0);
	phfract.resize(phnames.size(),0.);
	string root;
	Transitions1.push_back(TK_start);
	
	
	//c_no_ph_creation();
	root= "CEQ_";
	string Ceq_Name=root;
	
	if (PARALLEL>0) {
	int iceq;
		for (int i=0; i<n_step+1;i++){
			void *ceqi= NULL;
			Ceq_Name=root+IntToString(i);//in order to have a different name for each equilibrium
			iceq=Create_New_Ceq_and_Return_ID(Ceq_Name);// iceq is the index in the equilibrium vector eqlista of OC3 
//			cout<<Ceq_Name<<" "<<iceq<<endl;	
			Store_Equilibria.push_back(iceq);//all the indexes are stored in the vector Store_Equilibria
			c_selecteq(iceq, &ceqi);
//			Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceqi);// 
			SetPressure(1.0E5, &ceqi);// Set Pressure when ceqi is created (for the first loop of Global_Find_Transitions)
			SetMoles(1.0, &ceqi); // Set Number of moles when ceqi is created
			SetComposition(W, &ceqi,i_ref,mycompo_unit);// Set the composition  when ceqi is created
			SetTemperature(1200, &ceqi);
			c_set_status_globaldata();
		//---------------------Compute Equilibrium----------------------------
			int i_error=0;
            for (int k=0;k<phnames.size();k++) Change_Phase_Status(phnames[k],PHENTERED,0.0,&ceqi);
			Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceqi);// 
			CalculateEquilibrium(&ceqi,NOGRID,i_error);  // option GRID
			
		}
	}
	bool first_iteration=true;
	vector< string > Phase_transitions_mixture;
	for (size_t k=0; k<number_of_loops;k++){
		
		cout<<"         loop n:"<<k+1<<" increment of T="<< step_TK<<endl;
		
		if (k>0) first_iteration=false;
		Transitions0.resize(0);
		for (size_t i=0;i<Transitions1.size();i++) {
			Transitions0.push_back(Transitions1[i]);
//			cout<<i<<" "<<Transitions1[i]<<endl;
		}
		Transitions1.resize(0);
		bool last_iteration=false;
		if (k==number_of_loops-1) last_iteration=true;
		
		
		for (size_t i=0; i<Transitions0.size();i++){
//			cout<<"treating transition : "<<Transitions0[i]<<endl;
			TK_start=Transitions0[i];
			
			Find_Transitions(TK_start,n_step,step_TK,W,phnames,Transitions1,el_reduced_names,first_iteration,last_iteration,Store_Equilibria,Phase_transitions_mixture, ceq,required_accuracy_on_TK);
		}
		
		double old_step_TK=step_TK;
		step_TK=step_TK/n_step;
	}
	cout<<endl;
	cout<<endl;
	cout<<"===================================================="<<endl;
	cout<<" TQ Parallel: ";
	if (PARALLEL==0) {
		cout<<"N0";
	}
	else{
		cout<<"Yes";
		cout<<" / number of threads: "<<NCPU;
	}
	
	cout<<endl;
	cout<<"=========================================================="<<endl;
	cout<<" Here are the transition temperatures that have been found "<<endl;
	cout<<" in the temperature range ["<<TK_end_ini-TCtoTK<<","<<TK_start_ini-TCtoTK<<"] C"<<endl;
	cout<<" for the following composition:                    "<<endl;
	cout<<endl;
	for (size_t i=0;i<el_reduced_names.size();i++) {
		cout<<"      "<<el_reduced_names[i]<<" ("<<mycompo_unit<<"): "<<W[i]<<endl;
	}
	cout<<" -------------------------------------------------------- "<<endl;
	for (size_t i=0;i<Transitions1.size();i++) {
		cout<<" "<<i<<" "<<Transitions1[i]-TCtoTK<<" C  "<<Phase_transitions_mixture[i]<<endl;
	}
}
//************************************************************************************************************************************************************************


//************************************************************************************************************************************************************************

void Random_Equilibrium_Loop(double &TK_min,double &TK_max, vector<double> &W_ini, const vector<string> &phnames,const vector<string> &el_reduced_names, const void *ceq, const int i_ref, const string &compo_unit, const int &total_number_of_loops){	
	
	string mycompo_unit=compo_unit;	
	vector<  vector< double> > LISTCOMPO;
	vector< double> LISTTK;
	vector<double> Wrand;
	Wrand.resize(W_ini.size(),0.);
	int total_number_of_errors=0;
	string root= "CEQ_rand_";
	string Ceq_Name=root;
	vector<int> Store_Equilibria;
	Store_Equilibria.resize(0);
	int number_of_loops=NCPU*20;
	int jter_print=0;
//	if (PARALLEL==0) number_of_loops=1;
	
	LISTCOMPO.resize(number_of_loops,vector<double>(W_ini.size(),0.));
	LISTTK.resize(number_of_loops,0.);
	omp_set_num_threads(NCPU);
	cout<<"number of threads detected:"<<omp_get_num_procs()<<endl;
	string parallel_mode=" TQ Parallel: ";
	
		if (PARALLEL==0) {
			parallel_mode+="N0";
		}
		else{
			parallel_mode+="Yes";
			parallel_mode+=" / number of threads: "+IntToString(NCPU)+"  ";
		}
	
	if (PARALLEL>0) {
		int iceq;
		for (int i=0; i<number_of_loops;i++){
			void *ceqi= NULL;
			Ceq_Name=root+IntToString(i);//in order to have a different name for each equilibrium
			iceq=Create_New_Ceq_and_Return_ID(Ceq_Name);// iceq is the index in the equilibrium vector eqlista of OC3 
		//cout<<Ceq_Name<<" "<<iceq<<endl;	
			Store_Equilibria.push_back(iceq);//all the indexes are stored in the vector Store_Equilibria
			c_selecteq(iceq, &ceqi);
			for (int k=0;k<phnames.size();k++) Change_Phase_Status(phnames[k],PHENTERED,0.0,&ceqi);
 		    Change_Phase_Status("LIQUID",PHENTERED,0.5,&ceqi);// 
			Change_Phase_Status("LIQUID",PHENTERED,0.5,&ceqi);// 
			SetPressure(1.0E5, &ceqi);// Set Pressure when ceqi is created (for the first loop of Global_Find_Transitions)
			SetMoles(1.0, &ceqi); // Set Number of moles when ceqi is created
			SetComposition(W_ini, &ceqi,i_ref,mycompo_unit);// Set the composition  when ceqi is created
			SetTemperature(2000., &ceqi);
			
		//---------------------Compute Equilibrium----------------------------
			int i_error=0;
//			for (int k=0;k<phnames.size();k++) Change_Phase_Status(phnames[k],PHENTERED,0.0,&ceqi);
//			Change_Phase_Status("LIQUID",PHENTERED,0.5,&ceqi);// 
//			Change_Phase_Status("FCC_A1",PHENTERED,0.5,&ceqi);// 
			CalculateEquilibrium(&ceqi,NOGRID,i_error);  // option GRID
		}
	}
	
	int iter=0;
	do{
		iter+=1;
		int i_error=0;

		for (int k=0; k<number_of_loops;k++){
			
			void *ceqi= NULL;
			if ((PARALLEL>0)) {
				c_selecteq(Store_Equilibria[k], &ceqi);// retrieve the pointer with index stored in Store_Equilibria
				
			}else{
		//		ceqi=ceq;// if no parallelization use STANDART EQUILIBRIUM
				c_selecteq(1, &ceqi);
			}
			
			
	 //		Change_Phase_Status("FCC_A1",PHENTERED,0.1,&ceqi);//
	//		cout<<"T="<<TKCE[i]<<endl;
			Wrand[i_ref]=1.0;
				
			for (int i=0;i<W_ini.size();i++){
				double xrand01=((double)rand()/(double)RAND_MAX);
				
				if (not(i==i_ref)) {
					Wrand[i]=xrand01*W_ini[i];
					Wrand[i_ref]-=Wrand[i];
				}
			}
			if (Wrand[i_ref]<0){
				cout<<" reference element with negative concentration"<<endl;
				exit(EXIT_FAILURE);
			}
			
			double xrand01=((double)rand()/(double)RAND_MAX);
			
			double TK=TK_min+(TK_max-TK_min)*xrand01;
			LISTTK[k]=TK;
			for (int i=0;i<W_ini.size();i++){
				LISTCOMPO[k][i]=Wrand[i];
			}
			
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceqi);
			Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceqi);// 
 // 			Change_Phase_Status("FCC_A1",PHENTERED,0.5,&ceqi);// 
			SetComposition(Wrand, &ceqi, i_ref,mycompo_unit);// Set the composition  when ceqi is created
			SetTemperature(TK, &ceqi); // set temperature for specific equilibrium
			
	}
/*
#if PARALLEL>0
#pragma omp parallel for 
#endif
	for (int k=0; k<number_of_loops;k++){
			
			void *ceqi= NULL;
			if ((PARALLEL>0)) {
				c_selecteq(Store_Equilibria[k], &ceqi);// retrieve the pointer with index stored in Store_Equilibria
				
			}else{
		//		ceqi=ceq;// if no parallelization use STANDART EQUILIBRIUM
				c_selecteq(1, &ceqi);
			}
			
	
			int i_error=0;	
			CalculateEquilibrium(&ceqi,NOGRID,i_error);
			if (not(i_error==0)){
				cout<<" equilibrium calculation not converged in transition subroutine for the following conditions"<<endl;
				cout<<" TK="<<LISTTK[k]<<endl;
				cout<<" composition:"<<endl;
				for (size_t i=0;i<el_reduced_names.size();i++) {
					cout<<el_reduced_names[i]<<" (w%): "<<LISTCOMPO[k][i]<<endl;
				}
				
			
			}
}
	for (int k=0; k<number_of_loops;k++){
				void *ceqi= NULL;
			if ((PARALLEL>0)) {
				c_selecteq(Store_Equilibria[k], &ceqi);// retrieve the pointer with index stored in Store_Equilibria
				
			}else{
		//		ceqi=ceq;// if no parallelization use STANDART EQUILIBRIUM
				c_selecteq(1, &ceqi);
			}
			
			SetTemperature(LISTTK[k], &ceqi); // set temperature for specific equilibrium
//			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceqi);
//  		    Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceqi);// 
	}

*/


#if PARALLEL>0
#pragma omp parallel for default(none),schedule(dynamic), private (k,i_error), shared(LISTTK,LISTCOMPO,Store_Equilibria, total_number_of_errors, jter_print,number_of_loops,el_reduced_names)
#endif	    
	for (int k=0; k<number_of_loops;k++){
			
			void *ceqi= NULL;
			if ((PARALLEL>0)) {
				c_selecteq(Store_Equilibria[k], &ceqi);// retrieve the pointer with index stored in Store_Equilibria
				
			}else{
		//		ceqi=ceq;// if no parallelization use STANDART EQUILIBRIUM
				c_selecteq(1, &ceqi);
			}
			
	
			int i_error=0;
			
			
			
			CalculateEquilibrium(&ceqi,NOGRID,i_error); // perform an equilibrium calculation
			
			
			if (not(i_error==0)){
				/*cout<<" equilibrium calculation not converged in transition subroutine for the following conditions"<<endl;
				cout<<" TK="<<LISTTK[k]<<endl;
				cout<<" composition:"<<endl;
				for (int i=0;i<el_reduced_names.size();i++) {
					cout<<el_reduced_names[i]<<" (w%): "<<LISTCOMPO[k][i]<<endl;
				}
				*/
				SetTemperature(LISTTK[k]-1, &ceqi); // set temperature for specific equilibrium
				CalculateEquilibrium(&ceqi,GRID,i_error); // perform an equilibrium calculation
				if (i_error==0) {
					//cout<<" case fixed"<<endl;
				}else{
					total_number_of_errors+=1;
				}
                               
				
//				exit(EXIT_FAILURE);
			}
																	  
		}
		
		jter_print+=number_of_loops;
		if (jter_print>199) {
			cout<<parallel_mode<<" ====>total number of random tests:"<<iter*number_of_loops<<"  number of errors : "<<total_number_of_errors<<endl;
                        jter_print=0;
		}

	}while((iter*number_of_loops)<total_number_of_loops);
	cout<<parallel_mode<<" ====>total number of random tests:"<<iter*number_of_loops<<"  number of errors : "<<total_number_of_errors<<endl;	
}		  

	

