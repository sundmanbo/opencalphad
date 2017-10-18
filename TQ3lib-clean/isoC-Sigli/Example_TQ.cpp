#define CTEC 0// should be 0 for non CTEC users

# if CTEC<1
#include "tqintf.h"
#endif
#if CTEC>0
#include  "CTEC.h"
#endif
#include <sys/time.h>







using namespace std;
int main(int argc, char **argv)
{
	int ncpu=1;
	vector<int> Store_Equilibria;
	Store_Equilibria.resize(0);
	vector<string> Store_Equilibria_compo_unit;
	vector<string> Suspended_phase_list;
	Suspended_phase_list.resize(0);

	Store_Equilibria.resize(0);
	Store_Equilibria_compo_unit.resize(0);
	vector<string> eldatabase;
	string el_ref;
	bool compo_in_percent=false;
	string compo_unit="W";
	string temp_unit="C";
	vector<double> Compo_all_el;
    vector<double> Compo_all_el_old;
	int i_ref=0;
	vector<double> W;
	vector<double> MU;
	vector<string> el_reduced_names;
	el_reduced_names.resize(0);	// Array including selected elements
	W.resize(0);

	int i_error=0;
	void *ceq =0;  // Pointer to the OpenCalphad storage
	double TK=2000;
	double TC=2000;
	double TK_Liquidus=10;
	vector<string> phnames; // Array including all phase names
	vector<double> phfract; // Array including all phase fractions
	vector< vector<double> > elfract;                                           // Array including all equilibrium compositions
	struct timeval start1, end1;
	ofstream file;
    long seconds, useconds;
	double elapsed_time;


	char command[255];
	omp_set_num_threads(ncpu);
	//cout <<" name of the input file :"<<argv[1] << endl ;
	ifstream inputfile;
	inputfile.open(argv[1]);
	string TDBFILE = "";
	char charname[255];
	char comment[255];
	char myline[1024];
	int line_number=1;
	bool data_base_already_read=false;
	vector<string> month;
	month.resize(13,"");
	month[1]="January";
	month[2]="February";
	month[3]="March";
	month[4]="April";
	month[5]="May";
	month[6]="June";
	month[7]="July";
	month[8]="August";
	month[9]="September";
	month[10]="October";
	month[11]="November";
	month[12]="December";







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

		// *************************************************************************************************************
		if(strcommand=="TDB_FILE_NAME"){
			cout<<setw(50)<<strcommand<<" ";


			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string name =strmyline.substr(0,i);
			cout<<name<<endl;
			cout<<endl;
			TDBFILE=name;
			file<<" name of the thermodynamic data base: "<<TDBFILE<<endl;
			file<<endl;
			bool CTEC_activated =false;


			ifstream f(TDBFILE.c_str());
			if (not(f.good())){
				cout<<"tdb file "<<TDBFILE<<" not found"<<endl;
				exit(EXIT_FAILURE);
			}

			f.close();
			#if CTEC<1
			{
				GetAllElementsFromDatabase(TDBFILE);
			}
			#endif
			//************************ applicable for CTEC only (WITH CTEC=1)
		    #if CTEC>0
			{
				CTECGetAllElementsFromDatabase(TDBFILE);
			}
			#endif
			//end *************************** applicable for CTEC only

			cout<<" the following elements are in the database:"<<endl;
			cout<<" ";
			file<<" elements in the database: ";
			for (size_t i=0; i<c_nel; i++){
				string mystr(c_cnam[i]);
				All_Capital_Letters(mystr);
				cout<<mystr<<" / ";
				eldatabase.push_back(mystr);
				file<<mystr<<" ";
			}
			cout<<endl;
			file<<endl;
			Compo_all_el.resize(eldatabase.size(),0.);
			Compo_all_el_old.resize(eldatabase.size(),0.);
		}
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_REF_ELEMENT"){
			cout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string name =strmyline.substr(0,i);

			el_ref=name;
			All_Capital_Letters(el_ref);
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
			cout <<el_ref<<endl;;
			file<<" name of the reference element: "<<name<<endl;
		}
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_UNIT_COMPO_INPUT"){
			cout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			cout<<strcharname<<endl;
			file<<" unit used for composition input: "<<strcharname<<endl;
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
		else if(strcommand=="DEFINE_NCPU"){
			cout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			cout<<strcharname<<endl;
			ncpu=atoi(strcharname.c_str());
			omp_set_num_threads(ncpu);
			file<<" number of threads used for some of the parallel TQ calculations: "<<ncpu<<endl;
         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_UNIT_TEMP_INPUT"){
			cout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			cout<<strcharname<<endl;

			temp_unit=strcharname;
			if (not((temp_unit=="C")or(temp_unit=="K"))){
				cout<<"problem detected in temperature input units"<<endl;
				exit(EXIT_FAILURE);
			}
			file<<" units used for input temperatures: "<<temp_unit<<endl;
         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_OUTPUT_FILE_NAME"){


			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			    time_t now = time(0);
		    tm *ltm = localtime(&now);

			file.open (strcharname.c_str());
			file<<endl;
			file<<"*************************************************************************************************************"<<endl;
			file<<endl;
			file<<"                                     "<<OCVERSION<<endl;
			file<<"                Computation performed on: "<<ltm->tm_mday<<" "<<month[ 1 + ltm->tm_mon]<<" "<<1900 + ltm->tm_year<<" , "<< 1 + ltm->tm_hour << "h:"<< 1 + ltm->tm_min << "mn:"<< 1 + ltm->tm_sec<<"s"<< endl;
			file<<endl;
			file<<"*************************************************************************************************************"<<endl;
			file<<endl;
			// current date/time based on current system


			cout<<endl;
			cout<<"*************************************************************************************************************"<<endl;
			cout<<endl;
			cout<<"                                     "<<OCVERSION<<endl;
			cout<<"                 Computation performed on: "<<ltm->tm_mday<<" "<<month[ 1 + ltm->tm_mon]<<" "<<1900 + ltm->tm_year<<" , "<< 1 + ltm->tm_hour << "h:"<< 1 + ltm->tm_min << "mn:"<< 1 + ltm->tm_sec<<"s"<< endl;
			cout<<endl;
			cout<<"*************************************************************************************************************"<<endl;
			cout<<endl;

         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_COMPOSITION"){
			cout<<setw(50)<<strcommand<<" ";
			for (int i=0;i<Compo_all_el_old.size();i++) Compo_all_el_old[i]=Compo_all_el[i];
			double compo_ref=1.0;
			if (data_base_already_read) for (size_t k=0; k<el_reduced_names.size();k++) W[k]=0.0;

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string strcharname =strmyline.substr(0,i);
			cout<<strcharname<<endl;
			i=strcharname.find("=");
			double factor=1.0;
			if (compo_in_percent) factor=0.01;

			while (not((i<0) or (i>strcharname.size()) )){
				string element_name=strcharname.substr(0,i);
				int j=strcharname.find("/");
				string strcompo=strcharname.substr(i+1,j-i-1);
				All_Capital_Letters(element_name);
				strcharname.erase(0,j+1);
				if (j<0) strcharname="";
				i=strcharname.find("=");
				bool found_el=false;
				if (not data_base_already_read){
					for (size_t k=0; k<eldatabase.size();k++){
						if (element_name==eldatabase[k]){
							found_el=true;
							if (not(element_name==el_ref)) {
								Compo_all_el[k]=atof(strcompo.c_str())*factor;
								compo_ref-=Compo_all_el[k];
							}else{
								cout<<"you are not supposed to give the componsition for the reference element"<<endl;
							}
							break;
						}
					}
					if (not found_el){
						cout<<"error in composition definition in line:"<<line_number<<endl;
						cout<<"element "<<element_name<<" not present in the database"<<endl;
					}
				}
				else{
					for (size_t k=0; k<el_reduced_names.size();k++){
						if (element_name==el_reduced_names[k]){
							found_el=true;
							if (not(element_name==el_ref)) {
								W[k]=atof(strcompo.c_str())*factor;
								compo_ref-=W[k];
							}
							break;
						}
					}
					if (not found_el){
						cout<<" you are not allowed to add a new element afer a first DEFINE_COMPOSITION. New element :"<<element_name<<endl;
						exit(EXIT_FAILURE);
					}
				}


			}
			if (not data_base_already_read){
				file<<"*************************************************************************************************************"<<endl;
				for (size_t i=0; i<eldatabase.size();++i){
					if (eldatabase[i]==el_ref) {
						Compo_all_el[i]=compo_ref;
						break;
					}
				}
			}
			else{
				W[i_ref]=compo_ref;
			}





			if (not data_base_already_read){
				file<<" first composition analyzed (which determines the set of elements to be read in the database):"<<endl;
				bool firts_el=true;
				file<<"[composition]"<<endl;
				file<<" ";
				for (size_t i=0; i<eldatabase.size();++i){
					if (Compo_all_el[i]>1e-10){
						if (not firts_el) file<<" / ";
						el_reduced_names.push_back(eldatabase[i]);
						W.push_back(Compo_all_el[i]);
						file<<eldatabase[i]<<"="<<Compo_all_el[i];

						firts_el=false;
					}
				}
				file<<endl;file<<endl;

				for (size_t i=0; i<el_reduced_names.size();++i){
					if (el_reduced_names[i]==el_ref) i_ref=i;
				}
				c_set_status_globaldata();

				Initialize(&ceq);                                                          // Initialize OpenCalphad and allocate memory to the first equilibrium
				#if CTEC<1
				{
					ReadDatabaseLimited(TDBFILE, el_reduced_names, &ceq);                       // Define TDB-file and read only selected elements (non zero composition)
				}
				#endif
				//************************ applicable for CTEC only (WITH CTEC=1)
				#if CTEC>0
				{
					CTECReadDatabaseLimited(TDBFILE, el_reduced_names, &ceq);                       // Define TDB-file and read only selected elements (non zero composition)

				}
				#endif
				//end *************************** applicable for CTEC only




				data_base_already_read=true;
				ReadPhases(phnames, &ceq);                                                  // Read Phases data in tdb file
				cout<<" list of possible phases in the system :"<<endl;
				file<<" list of possible phases for this set of elements:"<<endl;
				int length=0;
				for (int i=0;i<phnames.size();i++){
					cout<<" "<<phnames[i];
					file<<" "<<phnames[i];length=length+1+phnames[i].length();
					if (length>100){
						length=0;
						cout<<endl;
						file<<endl;
					}
				}
				cout<<endl;
				file<<endl;file<<endl;
				phfract.resize(phnames.size(),0.);

				elfract.resize(phnames.size(),vector<double>(el_reduced_names.size(),0.));

				SetPressure(1e5, &ceq);                                                     // Set Pressure

				SetMoles(1.0, &ceq);                                                          // Set Number of moles
				MU.resize(W.size(),0.);
			}
			else{
				file<<"=================================================================="<<endl;
				file<<" new composition:"<<endl;;
				file<<"[composition]"<<endl;
				bool firts_el=true;
				for (size_t i=0; i<el_reduced_names.size();++i){
					if (not (i==i_ref)){
						if (not firts_el) file<<" / ";
						file<<el_reduced_names[i]<<"="<<W[i];
						firts_el=false;
					}
				}
				file<<endl;file<<endl;
			}

			SetComposition(W, &ceq,i_ref, compo_unit);                                                    // Set Composition of the system
			TK=2000;                                                  //Set temperature
			SetTemperature(TK, &ceq);

			//---------------------Compute Equilibrium----------------------------

		//	List_Conditions(&ceq);




			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);


		}

		//************************************************************************************************************
		// end of if(strcommand=="DEFINE_COMPOSITION")
		//************************************************************************************************************
		else if(strcommand=="LIQUIDUS"){
			cout<<setw(50)<<strcommand<<endl;
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strLIQUID =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strLIQUID);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strSOLSOL =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strSOLSOL);
			cout<<strLIQUID<<" "<<strSOLSOL<<endl;
			double TK=1300.;
			SetTemperature(TK, &ceq);

			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);

			ResetTemperature(&ceq);	 //remove condition on temperature


			Change_Phase_Status(strLIQUID,PHFIXED,0.9999,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...

			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);  // option GRID    // this is why the previous equilibrium was at high T to have liquid present


			if (not (i_error==0)){
				cout<<"*";
				double targeted_fraction=1.0-1e-4;
				double temperature_accuracy=1e-3;
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
				find_TK_for_a_given_Liquid_fraction(TK, i_error,strLIQUID,strSOLSOL,targeted_fraction, temperature_accuracy, ceq, phnames, Suspended_phase_list);
			}else{
				TK=ReadTemperature(&ceq);
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
			}

			if (i_error==0){
				TK_Liquidus=TK;
				TC=TK-TCtoTK;

				cout<<" ----> liquidus is: "<<TC<<" C"<<endl;
				file<<" The Liquidus is: "<<TC<<" C"<<endl;
				cout<<endl;
			}else{
				cout<<" liquidus not converged"<<endl;
				file<<" liquidus not converged"<<endl;
			}

		}
		// *************************************************************************************************************
		else if(strcommand=="SOLIDUS"){
			cout<<setw(50)<<strcommand<<endl;
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strLIQUID =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strLIQUID);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strSOLSOL =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strSOLSOL);
			cout<<strLIQUID<<" "<<strSOLSOL<<endl;


			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
			double TK=2000;
			SetTemperature(TK, &ceq);
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);

			ResetTemperature(&ceq);	 			//remove condition on temperature

			Change_Phase_Status(strLIQUID,PHFIXED,0.00001,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...

			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);  // option GRID    // this is why the previous equilibrium was at high T to have liquid present
			if (not i_error==0){
				cout<<"*";
				double targeted_fraction=1e-4;
				double temperature_accuracy=1e-3;
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
				find_TK_for_a_given_Liquid_fraction(TK, i_error,strLIQUID,strSOLSOL,targeted_fraction, temperature_accuracy, ceq, phnames,Suspended_phase_list);
			}else{
				TK=ReadTemperature(&ceq);
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
			}

			if (i_error==0){

				TC=TK-TCtoTK;

				cout<<" ----> solidus is: "<<TC<<" C"<<endl;
				file<<" The Solidus is: "<<TC<<" C"<<endl;
				cout<<endl;
			}else{
				cout<<" solidus not converged"<<endl;
				file<<" solidus not converged"<<endl;
			}


		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_TRANSITION_TEMPERATURES"){
			cout<<setw(50)<<strcommand<<" ";
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
			cout<<strTK_start<<" / "<<strTK_end<<" / "<<straccuracy <<" / "<<strnstep<<endl;
			double TK_start=atof(strTK_start.c_str());
			double TK_end=atof(strTK_end.c_str());
			double required_accuracy_on_TK=atof(straccuracy.c_str());
			int nstep=atoi(strnstep.c_str());


			if (temp_unit=="C"){
				TK_start+=TCtoTK;
				TK_end+=TCtoTK;
			}

			gettimeofday(&start1, NULL);// get the present time

			// ************************************************************************************************************************************************************************
			Global_Find_Transitions(file,TK_start,nstep,TK_end,required_accuracy_on_TK, W, phnames,el_reduced_names,ceq,i_ref, compo_unit,ncpu, Store_Equilibria, Store_Equilibria_compo_unit,Suspended_phase_list);// find all the transitions temperatures for a given alloy composition
			// ************************************************************************************************************************************************************************
			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;
			cout<<"Store_Equilibria.size()="<<Store_Equilibria.size()<<endl;
			cout<<" elapsed time for the transition temperature routine (s)= "<<elapsed_time<<endl;
			cout<<endl;
			cout<<endl;

			TK_Liquidus=TK_start;
		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_RANDOM_EQUILIBRIA"){
			cout<<setw(50)<<strcommand<<" ";
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

			cout<<strTK_min<<" / "<<strTK_max<<" / "<<strnloops<<endl;

			double TK_max=atof(strTK_min.c_str());
			double TK_min=atof(strTK_max.c_str());
			int nloops=atoi(strnloops.c_str());

			if (temp_unit=="C"){
				TK_max+=TCtoTK;
				TK_min+=TCtoTK;
			}

			Random_Equilibrium_Loop(TK_min,TK_max, W, phnames,el_reduced_names, ceq,i_ref,compo_unit,nloops,ncpu,Store_Equilibria, Store_Equilibria_compo_unit,Suspended_phase_list);

		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_EQUILIBRIUM"){
			cout<<setw(50)<<strcommand<<" ";
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
			cout<<strT<<" / "<<strni<<endl;
			if (temp_unit=="C"){
				TK+=TCtoTK;
			}
			SetTemperature(TK, &ceq);

			CalculateEquilibrium(&ceq,GRID,i_error,Suspended_phase_list);
			if (not(i_error==0)) {
				cout<<" equilibrium failed to converge at line:"<<line_number<<endl;
				file<<" equilibrium failed to converge at line:"<<line_number<<endl;
			}
			else{
				Write_Results_Equilibrium(file,el_reduced_names,phnames,phfract,elfract,ceq,idetail,compo_unit,MU);
			}
			cout<<endl;
		}
		// *************************************************************************************************************
		else if(strcommand=="SUSPEND_A_PHASE"){
		    cout<<setw(50)<<strcommand<<" ";

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strphasename=strmyline.substr(0,i);
			All_Capital_Letters(strphasename);
			bool phase_found=false;
			for (int i=0;i<phnames.size() and not phase_found;i++){
				if (phnames[i]==strphasename){
					phase_found=true;
				}
			}
			if (phase_found ){
				Suspended_phase_list.push_back(strphasename);
				//Change_Phase_Status(strphasename,PHSUS,0.0,&ceq);//
				cout<<strphasename<<" suspended"<<endl;
				file<<" "<<strphasename<<" suspended"<<endl;
			}else
			{
				cout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
			}


		}
		// *************************************************************************************************************
		else if(strcommand=="SUSPEND_ALL_PHASES_BUT_ONE"){
			cout<<setw(50)<<strcommand<<" ";
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strphasename=strmyline.substr(0,i);
			All_Capital_Letters(strphasename);
			bool phase_found=false;
			for (int i=0;i<phnames.size() and not phase_found;i++){
				if (phnames[i]==strphasename){
					phase_found=true;
				}
			}
			if (phase_found ){

				for (int i=0;i<phnames.size();i++) {
					if (not(strphasename==phnames[i]))Suspended_phase_list.push_back(phnames[i]);
				}
				cout<<" all phases have been suspended but :"<<strphasename<<endl;
				file<<" all phases have been suspended but :"<<strphasename<<endl;
			}else
			{
				cout<<"error in line "<<line_number<<" : phase does not exist"<<endl;

			}


		}
		// *************************************************************************************************************
		else if(strcommand=="RESTORE_A_PHASE"){
			cout<<setw(50)<<strcommand<<" ";
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strphasename=strmyline.substr(0,i);
			All_Capital_Letters(strphasename);
			bool phase_found=false;
			for (int i=0;i<phnames.size() and not phase_found;i++){
				if (phnames[i]==strphasename){
					phase_found=true;
				}
			}

			if (phase_found ){
				phase_found=false;
				int i_found=0;
				for (int i=0;i<Suspended_phase_list.size()and not phase_found;i++){
					if (strphasename==Suspended_phase_list[i]) {
						i_found=i;
						phase_found=true;
					}
				}
				if  (phase_found ) {
					Suspended_phase_list.erase(Suspended_phase_list.begin() +i_found);
					Change_Phase_Status(strphasename,PHENTERED,0.,&ceq);//
					cout<<strphasename <<" reactivated"<<endl;
					file<<" "<<strphasename <<" reactivated"<<endl;
				}
			}else
			{
				cout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
			}


		}
		// *************************************************************************************************************
		else if(strcommand=="RESTORE_ALL_PHASES"){
			cout<<setw(50)<<strcommand<<" ";
			Suspended_phase_list.resize(0);
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.,&ceq);
			Change_Phase_Status("LIQUID",PHENTERED,1.,&ceq);
			cout<< "all phases have been reactivated"<<endl;
			file<< "all phases have been reactivated"<<endl;
		}

		// *************************************************************************************************************
		else if(strcommand=="SCHEIL_SOLIDIFICATION"){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			cout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strLIQUID =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strLIQUID);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strSOLSOL =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strSOLSOL);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strtarget_delta_f_liq =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_min =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, fourth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);



			cout<<strLIQUID<<" / "<<strSOLSOL<<" / "<<strtarget_delta_f_liq<<" / "<<strdelta_T_min<<" / "<<strdelta_T_max<<endl;

			double target_delta_f_liq=atof(strtarget_delta_f_liq.c_str());
			double delta_T_min=atof(strdelta_T_min.c_str());
			double delta_T_max=atof(strdelta_T_max.c_str());

			if (TK_Liquidus<30){
				cout<<" you need to have a valid liquidus temperature to start a sheill calculation"<<endl;
				cout<<" liquidus="<<TK_Liquidus<<endl;
				exit(EXIT_FAILURE);
			}

			gettimeofday(&start1, NULL);// get the present time

			scheil_solidif(strLIQUID,strSOLSOL,file,el_reduced_names,phnames,ceq, W, target_delta_f_liq,delta_T_min,delta_T_max, TK_Liquidus,i_ref,compo_unit,Suspended_phase_list);

			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;

			cout<<" elapsed time for the scheil solidification routine (s)= "<<elapsed_time<<endl;
			cout<<endl;
			cout<<endl;

		}
		//***************************************************************************************************************
		else if((strcommand=="DIFF_SOLIDIFICATION")and (CTEC>0)){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			cout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, very first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strGradientFileOut =strmyline.substr(0,i);
			strmyline.erase(0,i+1);


			 i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strLIQUID =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strLIQUID);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strSOLSOL =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strSOLSOL);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strtarget_delta_f_liq =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, fourth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_min =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, fifth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, sixth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strhalf_sdas =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, seventh / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strNbincrement =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, heigth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdim =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTpoint =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strFl_end_hot_tearing =strmyline.substr(0,i);
			strmyline.erase(0,i+1);


			cout<<strLIQUID<<"/"<<strSOLSOL<<"/"<<strtarget_delta_f_liq<<"/"<<strdelta_T_min<<"/"<<strdelta_T_max<<"/"<<strhalf_sdas<<"/"<<strNbincrement<<"/"<<strdim<<"/"<<strTpoint<<endl;

			double target_delta_f_liq=atof(strtarget_delta_f_liq.c_str());
			double delta_T_min=atof(strdelta_T_min.c_str());
			double delta_T_max=atof(strdelta_T_max.c_str());
			double half_sdas=atof(strhalf_sdas.c_str());
			int Nb_increment=atoi(strNbincrement.c_str());
			double dim=atof(strdim.c_str());
			double Tpoint=atof(strTpoint.c_str());
			double Fl_end_hot_tearing=atof(strFl_end_hot_tearing.c_str());

			if (TK_Liquidus<30){
				cout<<" you need to have a valid liquidus temperature to start a sheill calculation"<<endl;
				cout<<" liquidus="<<TK_Liquidus<<endl;
				exit(EXIT_FAILURE);
			}

			gettimeofday(&start1, NULL);// get the present time
			#if CTEC>0
			back_diff_solidif(strGradientFileOut,strLIQUID,strSOLSOL,file,el_reduced_names,phnames,ceq, W, target_delta_f_liq,delta_T_min,delta_T_max, TK_Liquidus,i_ref,compo_unit,half_sdas,dim,Tpoint,Nb_increment,Store_Equilibria,Store_Equilibria_compo_unit,Suspended_phase_list,Fl_end_hot_tearing);
			#endif
			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;

			cout<<" elapsed time for the back-diffusion solidification routine (s)= "<<elapsed_time<<endl;
			cout<<endl;
			cout<<endl;

		}
		// *************************************************************************************************************
		else if((strcommand=="HOMOGENIZING")and (CTEC>0)){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			vector < double > TC1;
			vector < double > TC2;
			vector < double > segments_time_h;
			cout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, very first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strGradientFileIn =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			 i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strLIQUID =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strLIQUID);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strSOLSOL =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			All_Capital_Letters(strSOLSOL);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strprintresultevery_s =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, fourth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strprintgradientevery_h =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, fifth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			cout<<strLIQUID<<"/"<<strSOLSOL<<"/"<<strprintresultevery_s<<"/"<<strprintgradientevery_h<<"/"<<strdelta_T_max<<endl;
			double printresultevery_s=atof(strprintresultevery_s.c_str());
			double printgradientevery_h=atof(strprintgradientevery_h.c_str());
			double delta_T_max=atof(strdelta_T_max.c_str());
			bool end_of_thermal_cycle_detected=false;
			cout<<setw(50)<<"thermal cycle:";
			while (not end_of_thermal_cycle_detected){
				i = strmyline.find("/");
				if ((i<0) or (i>strmyline.size())) {
					cout<<" error in command line, sixth / not found in line:"<<line_number<<endl;
					exit(EXIT_FAILURE);
				}
				string strT1 =strmyline.substr(0,i);
				strmyline.erase(0,i+1);

				i = strmyline.find("/");
				if ((i<0) or (i>strmyline.size())) {
					cout<<" error in command line, seventh / not found in line:"<<line_number<<endl;
					exit(EXIT_FAILURE);
				}
				string strT2 =strmyline.substr(0,i);
				strmyline.erase(0,i+1);


				i = strmyline.find("/");

				if ((i<0) or (i>strmyline.size())) {
					i = strmyline.find(">");

					if ((i<0) or (i>strmyline.size())) {
						cout<<" error in command line, > not found in line:"<<line_number<<endl;
						exit(EXIT_FAILURE);
					}
					end_of_thermal_cycle_detected=true;
				}
				string strtime_h =strmyline.substr(0,i);
				strmyline.erase(0,i+1);



				cout<<"/"<<strT1<<"/"<<strT2<<"/"<<strtime_h;


				double valueT1=atof(strT1.c_str());
				int valueT2=atof(strT2.c_str());
				double valuestime_h=atof(strtime_h.c_str());

				TC1.push_back(valueT1);
				TC2.push_back(valueT2);
				segments_time_h.push_back(valuestime_h);
			}
			cout<<endl;

			string strGradientFileOut = strGradientFileIn.substr(0,strGradientFileIn.length()-4); // remove .txt
			strGradientFileOut+="andhomo.txt";
			SetTemperature(1000, &ceq);
			CalculateEquilibrium(&ceq,GRID,i_error,Suspended_phase_list);
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			gettimeofday(&start1, NULL);// get the present time
			#if CTEC>0
			homo(strGradientFileIn,strGradientFileOut,printresultevery_s,printgradientevery_h,strLIQUID,strSOLSOL,file,el_reduced_names,phnames,delta_T_max,i_ref,Store_Equilibria,Store_Equilibria_compo_unit,TC1, TC2,segments_time_h,Suspended_phase_list);
			#endif
			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;

			cout<<" elapsed time for the back-diffusion solidification routine (s)= "<<elapsed_time<<endl;
			cout<<endl;
			cout<<endl;

		}
		else if((strcommand=="PROPERTIES")and (CTEC>0)){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			vector < double > TC1;
			vector < double > TC2;
			vector < double > segments_time_h;
			cout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strFile =strmyline.substr(0,i);
			strmyline.erase(0,i+1);


			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strSOLSOL =strmyline.substr(0,i);
			All_Capital_Letters(strSOLSOL);
			cout<<endl;
			#if CTEC>0
			compute_properties(strFile,strSOLSOL,el_reduced_names, W, compo_unit,i_ref, &ceq, phnames);
			#endif

		}
		else if((strcommand=="FIX_A_PHASE")){

			cout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strPhase =strmyline.substr(0,i);
			strmyline.erase(0,i+1);



			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strvalue =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			double value=atof(strvalue.c_str());;

			Change_Phase_Status(strPhase,PHFIXED,value,&ceq);//
			cout<<endl;
		}
		else if((strcommand=="COMPUTE_EQUILIBRIUM_WITH_TEMPERATURE_CHANGED_BY")){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max

			cout<<setw(50)<<strcommand<<" ";


			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strvalue =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				cout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strni=strmyline.substr(0,i);
			int idetail=atoi(strni.c_str());

			double value=atof(strvalue.c_str());;
			cout<<value;
			TK=ReadTemperature(&ceq)+value;
			SetTemperature(TK, &ceq);
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			Write_Results_Equilibrium(file,el_reduced_names,phnames,phfract,elfract,ceq,idetail,compo_unit,MU);
			cout<<endl;

		}

		// *************************************************************************************************************
		else if(strcommand==""){
		}
		// *************************************************************************************************************
		else{
			cout<<setw(50)<<strcommand<<" ";
			cout<<"command line not recognized at line:"<<line_number<<endl;

		}

		line_number+=1;
	}


	file.close();
    return 0;
}
