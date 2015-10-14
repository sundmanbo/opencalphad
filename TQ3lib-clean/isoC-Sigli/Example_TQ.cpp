#include "tqintf.h"
#include <iostream> 
#include <fstream> 
#include <cstdlib>

#include <algorithm> 
#include <string> 
#include <sys/time.h>
#include <time.h>
#include <omp.h>

using namespace std;



int main(int argc, char **argv)
{
	vector<string> eldatabase;
	string el_ref;
	bool compo_in_percent=false;
	string compo_unit="W%";
	string temp_unit="C";
	vector<double> Compo_all_el;
    vector<double> Compo_all_el_old;
	int i_ref=0;
	vector<double> W;
	vector<string> el_reduced_names;
	el_reduced_names.resize(0);	// Array including selected elements
	W.resize(0);
	int i_error=0;
	void *ceq =0;  // Pointer to the OpenCalphad storage
	double TK=2000;
	double TC=2000;
	vector<string> phnames; // Array including all phase names
	vector<double> phfract; // Array including all phase fractions
	vector< vector<double> > elfract;                                           // Array including all equilibrium compositions
	struct timeval start1, end1;

    long seconds, useconds;    
	double elapsed_time;

	char command[255];
	omp_set_num_threads(NCPU);
	cout <<" name of the input file :"<<argv[1] << endl ;
	ifstream inputfile;
	inputfile.open(argv[1]);
	string TDBFILE = "";
	char charname[255];
	char comment[255];
	char myline[1024];
	int line_number=1;
	
	while (!inputfile.eof()){
		inputfile.getline(myline,1024,'\n');
		string strmyline(myline);
		//cout<<strmyline<<endl;
		int i = strmyline.find("<");
		if ((i<0) or (i>strmyline.size())) {
			cout<<" error in command line < not found in line:"<<line_number<<endl;
			exit(EXIT_FAILURE);
		}
		
	
		// TDB_FILE_NAME 
		// DEFINE_REF_ELEMENT
                // DEFINE_UNIT_COMPO_INPUT W W% W X%
		// DEFINE_COMPOSITION	
		string strcommand=strmyline.substr(0,i);
		strmyline.erase(0,i+1);
		cout<<" command: "<<strcommand<<endl;
		// *************************************************************************************************************
		if(strcommand=="TDB_FILE_NAME"){
			
			
			
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}	
			
			string name =strmyline.substr(0,i);
			cout<<name<<endl;
			TDBFILE=name+".TDB";
			
			ifstream f(TDBFILE.c_str());
			if (not(f.good())){
				cout<<"tdb file "<<TDBFILE<<" not found"<<endl;
				exit(EXIT_FAILURE);
			}
			f.close();
			
			GetAllElementsFromDatabase(TDBFILE);
		
			cout<<" the following elements are in the database:"<<endl;
			for (size_t i=0; i<c_nel; i++){
				cout<<c_cnam[i]<<" / ";
				eldatabase.push_back(c_cnam[i]);
			}
			cout<<endl;
			Compo_all_el.resize(eldatabase.size(),0.);
			Compo_all_el_old.resize(eldatabase.size(),0.);
		}
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_REF_ELEMENT"){
			
			
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}	
			
			string name =strmyline.substr(0,i);
			cout<<name<<endl;
			el_ref=name;
			transform(el_ref.begin(), el_ref.end(), el_ref.begin(), ::toupper);// to have it in CAPITAL LETTERS
			bool found_el_ref=false;
			for (size_t i=0; i<eldatabase.size();i++){
				if (el_ref==eldatabase[i]){
					found_el_ref=true;
					break;
				}
			}
			if (not found_el_ref){
				cout<<"reference element not fond in database"<<endl;
				exit(EXIT_FAILURE);
			}
			cout << " name of the reference element (solvant for example) : "<<el_ref<<endl;;
		}
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_UNIT_COMPO_INPUT"){
			
			
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}	
			string strcharname =strmyline.substr(0,i);
			
			i = strcharname.find("%");
			if (i==1) compo_in_percent=true;
			strcharname.erase(1,1);
			compo_unit=strcharname;
			if (not((compo_unit=="W")or(compo_unit=="X"))){
				cout<<"problem detected in composition input units"<<endl;
				exit(EXIT_FAILURE);			
			}
         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_UNIT_TEMP_INPUT"){
			
			
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}	
			string strcharname =strmyline.substr(0,i);
			
			
			temp_unit=strcharname;
			if (not((temp_unit=="C")or(temp_unit=="K"))){
				cout<<"problem detected in temperature input units"<<endl;
				exit(EXIT_FAILURE);			
			}
         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_COMPOSITION"){
			for (int i=0;i<Compo_all_el_old.size();i++) Compo_all_el_old[i]=Compo_all_el[i];
			double compo_ref=1.0;
			
			
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}	
			string strcharname =strmyline.substr(0,i);
			i=strcharname.find("=");
			double factor=1.0;
			if (compo_in_percent) factor=0.01;
			
			while (not((i<0) or (i>strcharname.size()) )){	
				string element_name=strcharname.substr(0,i);
				int j=strcharname.find("/");
				string strcompo=strcharname.substr(i+1,j-i-1);
				cout<<element_name<<":"<<strcompo<<endl;
				transform(element_name.begin(), element_name.end(),element_name.begin(), ::toupper);// to have it in CAPITAL LETTERS
				strcharname.erase(0,j+1);
				if (j<0) strcharname="";
				i=strcharname.find("=");
				bool found_el_ref=false;
				for (size_t k=0; k<eldatabase.size();k++){
					if (element_name==eldatabase[k]){
						found_el_ref=true;
						Compo_all_el[k]=atof(strcompo.c_str())*factor;
						compo_ref-=Compo_all_el[k];
						break;
					}
				}
				if (not found_el_ref){
					cout<<"error in composition definition in line:"<<line_number<<endl;
					cout<<"element "<<element_name<<" not present in the database"<<endl;
				}

				
			}
			
			for (size_t i=0; i<eldatabase.size();++i){
				if (eldatabase[i]==el_ref) {
					Compo_all_el[i]=compo_ref;
					break;
				}
			}
			
			for (size_t i=0; i<eldatabase.size();++i){
				if (Compo_all_el[i]>1e-10){
					el_reduced_names.push_back(eldatabase[i]);
					W.push_back(Compo_all_el[i]);
				}
			}

			for (size_t i=0; i<el_reduced_names.size();++i){
				if (el_reduced_names[i]==el_ref) i_ref=i;
			}

			c_set_status_globaldata();
			                                                           
			Initialize(&ceq);                                                          // Initialize OpenCalphad and allocate memory to the first equilibrium


			ReadDatabaseLimited(TDBFILE, el_reduced_names, &ceq);                       // Define TDB-file and read only selected elements (non zero composition)
                                                     
			ReadPhases(phnames, &ceq);                                                  // Read Phases data in tdb file
			                                                     
			phfract.resize(phnames.size(),0.);
			
			elfract.resize(phnames.size(),vector<double>(el_reduced_names.size(),0.));
			SetPressure(1.0E5, &ceq);                                                     // Set Pressure
			SetMoles(1.0, &ceq);                                                          // Set Number of moles
			SetComposition(W, &ceq,i_ref, compo_unit);                                                    // Set Composition of the system
			TK=2000;                                                  //Sest temperature 
			SetTemperature(TK, &ceq);
			c_set_status_globaldata();
			//---------------------Compute Equilibrium----------------------------
			
		//	List_Conditions(&ceq);
		//	Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceq);// 

			CalculateEquilibrium(&ceq,NOGRID,i_error);  // option GRID

		}
		//************************************************************************************************************	
		// end of if(strcommand=="DEFINE_COMPOSITION")	
		//************************************************************************************************************	
		else if(strcommand=="LIQUIDUS"){
			
			SetTemperature(2000, &ceq);
			CalculateEquilibrium(&ceq,NOGRID,i_error);  // option GRID
			ResetTemperature(&ceq);	                                                   //remove condition on temperature
			Change_Phase_Status("LIQUID",PHFIXED,0.9999,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...
			CalculateEquilibrium(&ceq,NOGRID,i_error);  // option GRID    // this is why the previous equilibrium was at high T to have liquid present
			TK=ReadTemperature(&ceq);
			TC=TK-TCtoTK;
			cout<<" ----> Liquidus is: "<<TC<<" C"<<endl;
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceq);//
			//Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,2);
		}
		// *************************************************************************************************************
		else if(strcommand=="SOLIDUS"){
									
			SetTemperature(2000, &ceq);
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceq);// 	
			CalculateEquilibrium(&ceq,NOGRID,i_error);  // option GRID
			ResetTemperature(&ceq);	                                                   //remove condition on temperature
			Change_Phase_Status("LIQUID",PHFIXED,0.00001,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...
			CalculateEquilibrium(&ceq,NOGRID,i_error);  // option GRID    // this is why the previous equilibrium was at high T to have liquid present
			TK=ReadTemperature(&ceq);
			TC=TK-TCtoTK;
			cout<<" ----> solidus is: "<<TC<<" C"<<endl;
			Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,1);
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status("LIQUID",PHENTERED,1.0,&ceq);//
		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_TRANSITION_TEMPERATURES"){
			
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTK_start =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			
			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTK_end =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			
			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string straccuracy =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strnstep=strmyline.substr(0,i);
			
			double TK_start=atof(strTK_start.c_str());
			double TK_end=atof(strTK_end.c_str());
			double required_accuracy_on_TK=atof(straccuracy.c_str());
			int nstep=atoi(strnstep.c_str());
			
			cout<<TK_start<<" "<<TK_end<<" "<<required_accuracy_on_TK<<" "<<nstep<<endl;
			
			if (temp_unit=="C"){
				TK_start+=TCtoTK;
				TK_end+=TCtoTK;
			}
			
			gettimeofday(&start1, NULL);// get the present time

			// ************************************************************************************************************************************************************************
			Global_Find_Transitions(TK_start,nstep,TK_end,required_accuracy_on_TK, W, phnames,el_reduced_names,ceq,i_ref, compo_unit);// find all the transitions temperatures for a given alloy composition
			// ************************************************************************************************************************************************************************
			gettimeofday(&end1, NULL);
	
			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;


			cout<<endl;
			cout<<"elapsed time (s)= "<<elapsed_time<<endl;
		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_RANDOM_EQUILIBRIA"){
			
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTK_min =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			
			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTK_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			
			
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strnloops=strmyline.substr(0,i);
			
			
			double TK_max=atof(strTK_min.c_str());
			double TK_min=atof(strTK_max.c_str());
			int nloops=atoi(strnloops.c_str());
			
			if (temp_unit=="C"){
				TK_max+=TCtoTK;
				TK_min+=TCtoTK;
			}
			
			Random_Equilibrium_Loop(TK_min,TK_max, W, phnames,el_reduced_names, ceq,i_ref,compo_unit,nloops);
			
		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_EQUILIBRIUM"){
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strT=strmyline.substr(0,i);
			
			strmyline.erase(0,i+1);
		
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strni=strmyline.substr(0,i);
			
			double TK=atof(strT.c_str());
			int idetail=atoi(strni.c_str());
			
			if (temp_unit=="C"){
				TK+=TCtoTK;
			}
			SetTemperature(TK, &ceq);
			CalculateEquilibrium(&ceq,NOGRID,i_error); 
			if (not(i_error==0)) {
				cout<<" equilibrium failed to converge at line:"<<line_number<<endl;
			}
			else{
				Write_Results_Equilibrium(el_reduced_names,phnames,phfract,elfract,ceq,idetail);	
			}	
				
		}
		// *************************************************************************************************************
		else{
			cout<<"command line not recognized at line:"<<line_number<<endl;
			
		}
			
		line_number+=1;
	}
	

	
    return 0;
}




