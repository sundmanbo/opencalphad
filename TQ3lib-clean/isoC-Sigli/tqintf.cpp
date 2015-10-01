#include "tqintf.h"
#include <iostream>
#include <algorithm> 
#include <string> 
#include <sys/time.h>
#include <time.h>

#include<omp.h>

using namespace std;



int main(int argc, char **argv)
{
	bool general_case=true;

	string TDBFILE = "";
 
 
    
	vector<string> eldatabase;
	TDBFILE="ctec";
	
	//------------------- read all the elements from a given tdb database--------------------------
 if (general_case){
	cout<<" Enter the name of the TDB file (*.tdb) : ";
	cin>>TDBFILE;
	cout<<endl;
	}
	TDBFILE+=".tdb";
	
	GetAllElementsFromDatabase(TDBFILE);
	
	cout<<" the following elements are in the database:"<<endl;
	for (size_t i=0; i<c_nel; i++){
		cout<<c_cnam[i]<<" / ";
		eldatabase.push_back(c_cnam[i]);
	}
	cout<<endl;
	
	
	
	//-------------- give the name of the refence element (Solvant for example)-------------------
	string el_ref="AL";
if (general_case){
	cout << " give the name of the reference element (solvant for example) : ";
	cin>>el_ref;
	cout<<endl;
	transform(el_ref.begin(), el_ref.end(), el_ref.begin(), ::toupper);
}	
	
	
	// vector of composition for all the elements of the database
	vector<double> Compo_all_el;
	Compo_all_el.resize(eldatabase.size(),0.);
/*	
   Ag 0
   AL 1
   Cr 2
   cu3
   fe4
   li5
   mg6
   mn7
   sc8
   si9
   ti10
   v11
   zn12
   zr13
   
   */
	
	


if (general_case){
	//read composition for all the elements of the database
	Read_Composition(eldatabase,Compo_all_el,el_ref);
}else{
     Compo_all_el[0]=0.7e-2;//Ag
	Compo_all_el[3]=2e-2;//Cu
	Compo_all_el[4]=0.08e-2;//Fe
 	Compo_all_el[5]=2e-2;//Li
	Compo_all_el[6]=2e-2;//Mg
 	Compo_all_el[7]=0.35e-2;//Mn
 	Compo_all_el[8]=0.1e-2;//Sc
	Compo_all_el[9]=0.04e-2;//Si
	Compo_all_el[10]=0.005e-2;//Ti
	Compo_all_el[12]=0.6e-2;//Zn
	Compo_all_el[13]=0.14e-2;//Zr
	
	Compo_all_el[1]=1.;//ref
	for (int i=0; i<Compo_all_el.size();i++){
		if (!(i==1)) Compo_all_el[1]-=Compo_all_el[i];
	}
}
	
	//compo_ref=Compo_all_el[1];
	
	// reduced vector of elements containing only elements with a non zero compo
	vector<double> W;
	vector<string> el_reduced_names;
	el_reduced_names.resize(0);	// Array including selected elements
	W.resize(0);

	for (size_t i=0; i<eldatabase.size();++i){
		if (Compo_all_el[i]>1e-10){
			el_reduced_names.push_back(eldatabase[i]);
			W.push_back(Compo_all_el[i]);
		}
	}
	c_set_status_globaldata();
	void *ceq =0;                                                              // Pointer to the OpenCalphad storage
	Initialize(&ceq);                                                          // Initialize OpenCalphad and allocate memory to the first equilibrium
	
	
    ReadDatabaseLimited(TDBFILE, el_reduced_names, &ceq);                       // Define TDB-file and read only selected elements (non zero composition)

    vector<string> phnames;                                                     // Array including all phase names
	ReadPhases(phnames, &ceq);                                                  // Read Phases data in tdb file
    vector<double> phfract;                                                     // Array including all phase fractions
	phfract.resize(phnames.size(),0.);
    vector< vector<double> > elfract;                                           // Array including all equilibrium compositions
	elfract.resize(phnames.size(),vector<double>(el_reduced_names.size(),0.));

	//-----------------------Initialize conditions-----------------------
	                                                  
	SetPressure(1.0E5, &ceq);                                                     // Set Pressure
	SetMoles(1.0, &ceq);                                                          // Set Number of moles
	SetComposition(W, &ceq);                                                    // Set Composition of the system
	double TK=1500;                                                             //Sest temperature 
	double TC;
	SetTemperature(TK, &ceq);
	c_set_status_globaldata();
	//---------------------Compute Equilibrium----------------------------
	int i_error=0;
	List_Conditions(&ceq);
	Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceq);// 
	
	CalculateEquilibrium(&ceq,NOGRID,i_error);  // option GRID
	
	if (i_error==0){
		Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,2);
	}
	
	
	
	/*	
	//---------------------Compute the liquidus of an alloy
	// not always converging 
	ResetTemperature(&ceq);	                                                   //remove condition on temperature
	Change_Phase_Status("LIQUID",PHFIXED,0.999,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...
	CalculateEquilibrium(&ceq, 0);      // this is why the previous equilibrium was at high T to have liquid present
	TK=ReadTemperature(&ceq);
	TC=TK-TCtoTK;
	cout<<" ----> Liquidus is: "<<TC<<" C"<<endl;
	Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,2);
		
	//---------------------Compute the solidus of an alloy
	// not always converging 
	Change_Phase_Status("LIQUID",PHFIXED,0.001,&ceq);
	CalculateEquilibrium(&ceq);
	TK=ReadTemperature(&ceq);
	TC=TK-273.15;
	cout<<" ----> Solidus is: "<<TC<<" C"<<endl;
	CalculateEquilibrium(&ceq);
	Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,2);


	Change_Phase_Status("LIQUID",PHENTERED,0.000001,&ceq);
	SetTemperature(TK, &ceq);
	CalculateEquilibrium(&ceq);
	
	Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,2);
	*/

	
	struct timeval start1, end1;

    long seconds, useconds;    
	double elapsed_time;

	gettimeofday(&start1, NULL);// get the present time
	
	double TK_start=1500;
	double step_TK=-250;//number of steps is predefined and equal to 10
	
	double required_accuracy_on_TK=1e-2;
	
	
	// find all the transitions temperatures for a given alloy composition
	Global_Find_Transitions(TK_start,NSTEP,step_TK,required_accuracy_on_TK, W, phnames,el_reduced_names,ceq);
	
	gettimeofday(&end1, NULL);
	
	seconds  = end1.tv_sec  - start1.tv_sec;
    useconds = end1.tv_usec - start1.tv_usec;

    elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;
	
	cout<<"elapsed time (s)= "<<elapsed_time<<endl;
    return 0;
}



