#define CTEC 0// should be 0 for non CTEC users

# if CTEC<1
#include "ocasiintf.h"
#define OCVERSION "Open CalPhad Software Interface July 2016"
#endif
#if CTEC>0
#include "CTEC.h"
#define OCVERSION "OC_Prophase July 2016"
#endif
#include <sys/time.h>

/*
	CTEC = 0  Compilation Open Source
	CTEC = 1  Compilation Constellium Unix Server
	CTEC = 2  Compilation Constellium Windows 7 PC
*/





using namespace std;
int main(int argc, char **argv)
{
	int ncpu=1;
	vector<int> Store_Equilibria;
	Store_Equilibria.resize(0);
	vector<string> Store_Equilibria_compo_unit;
	vector<string> Suspended_phase_list;
	Suspended_phase_list.resize(0);
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
	size_t i_compo=0;
	string strcompo="Compo";
	string strcomponb="";

	size_t i_eq=0;
	string strEqui="Equilibrium";
	string strEquinb="";

	string myequi="";
	string element_file="elements.txt";


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
	int line_number=1;
	bool data_base_already_read=false;
	bool right_to_use=true;
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
	string strLIQUID="LIQUID";
	string strSOLSOL="FCC_A1";

	#if CTEC==2
	{
		right_to_use=InitInstance();
	}
	#endif





	while (!inputfile.eof() ){
		string strmyline;
		std::getline(inputfile, strmyline);
		cout<<strmyline<<endl;
		int i = strmyline.find("<");
		if ((i<0) or (i>strmyline.size())) {
			sout<<" error in command line < not found in line:"<<line_number<<endl;
			exit(EXIT_FAILURE);
		}


		// TDB_FILE_NAME
		// DEFINE_REF_ELEMENT
                // DEFINE_UNIT_COMPO_INPUT W W% W X%
		// DEFINE_COMPOSITION
		string strcommand=strmyline.substr(0,i);
		strmyline.erase(0,i+1);

		// *************************************************************************************************************
		if(strcommand=="ELEMENT_FILE_NAME"){
			sout<<setw(50)<<strcommand<<" ";


			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			 element_file =strmyline.substr(0,i);
			sout<<element_file<<endl;
			sout<<endl;


			ifstream f(element_file.c_str());
			if (not(f.good())){
				sout<<"element file "<<element_file<<" not found"<<endl;
				exit(EXIT_FAILURE);
			}

			f.close();
		}
		// *************************************************************************************************************
		else if(strcommand=="TDB_FILE_NAME"){
			sout<<setw(50)<<strcommand<<" ";


			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string name =strmyline.substr(0,i);
			sout<<name<<endl;
			sout<<endl;
			TDBFILE=name;
			file<<" name of the thermodynamic data base: "<<TDBFILE<<endl;
			file<<endl;
			bool CTEC_activated =false;


			ifstream f(TDBFILE.c_str());
			if (not(f.good())){
				sout<<"tdb file "<<TDBFILE<<" not found"<<endl;
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

			sout<<" the following elements are in the database:"<<endl;
			sout<<" ";
			file<<" elements in the database: ";
			for (size_t i=0; i<c_nel; i++){
				string mystr(c_cnam[i]);
				All_Capital_Letters(mystr);
				sout<<mystr<<" / ";
				eldatabase.push_back(mystr);
				file<<TAB<<mystr;
			}
			sout<<endl;
			file<<endl;
			Compo_all_el.resize(eldatabase.size(),0.);
			Compo_all_el_old.resize(eldatabase.size(),0.);
		}
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_REF_ELEMENT"){
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
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
				sout<<"reference element not fond in database"<<endl;
				exit(EXIT_FAILURE);
			}
			sout <<el_ref<<endl;;
			file<<" name of the reference element: "<<TAB<<name<<endl;
		}
		else if(strcommand=="DEFINE_LIQUID_NAME"){
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string name =strmyline.substr(0,i);

			strLIQUID=name;
			All_Capital_Letters(strLIQUID);

			sout <<strLIQUID<<endl;;
			file<<" name of the liquid phase: "<<strLIQUID<<endl;
		}
		else if(strcommand=="DEFINE_SOLSOL_NAME"){
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string name =strmyline.substr(0,i);

			strSOLSOL=name;
			All_Capital_Letters(strSOLSOL);

			sout <<strSOLSOL<<endl;;
			file<<" name of the main solid solution phase: "<<strSOLSOL<<endl;
		}
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_UNIT_COMPO_INPUT"){
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			sout<<strcharname<<endl;
			file<<" unit used for composition input: "<<strcharname<<endl;
			i = strcharname.find("%");
			if (i==1) compo_in_percent=true;
			strcharname.erase(1,1);
			compo_unit=strcharname;
			if (not((compo_unit=="W")or(compo_unit=="X"))){
				sout<<"problem detected in composition input units"<<endl;
				exit(EXIT_FAILURE);
			}
         }
		 // *************************************************************************************************************
		else if(strcommand=="DEFINE_NCPU"){
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			sout<<strcharname<<endl;
			ncpu=atoi(strcharname.c_str());
			omp_set_num_threads(ncpu);
			file<<" number of threads used for some of the parallel TQ calculations: "<<ncpu<<endl;
         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_UNIT_TEMP_INPUT"){
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			sout<<strcharname<<endl;

			temp_unit=strcharname;
			if (not((temp_unit=="C")or(temp_unit=="K"))){
				sout<<"problem detected in temperature input units"<<endl;
				exit(EXIT_FAILURE);
			}
			file<<" units used for input temperatures: "<<temp_unit<<endl;
         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_OUTPUT_FILE_NAME"){


			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strcharname =strmyline.substr(0,i);
			    time_t now = time(0);
		    tm *ltm = localtime(&now);

			file.open (strcharname.c_str());
			if (not right_to_use){
				file<<"Licence not Valid"<<endl;
				sout<<"Licence not Valid"<<endl;
				file.close();
				exit(EXIT_FAILURE);
			}
			file<<endl;
			file<<"*************************************************************************************************************"<<endl;
			file<<endl;
			file<<"                                     "<<OCVERSION<<endl;
			file<<"                Computation performed on: "<<ltm->tm_mday<<" "<<month[ 1 + ltm->tm_mon]<<" "<<1900 + ltm->tm_year<<" , "<<  ltm->tm_hour << "h:"<<  ltm->tm_min << "mn:"<<  ltm->tm_sec<<"s"<< endl;
			file<<endl;
			file<<"*************************************************************************************************************"<<endl;
			file<<endl;
			// current date/time based on current system


			sout<<endl;
			sout<<"*************************************************************************************************************"<<endl;
			sout<<endl;
			sout<<"                                     "<<OCVERSION<<endl;
			sout<<"                 Computation performed on: "<<ltm->tm_mday<<" "<<month[ 1 + ltm->tm_mon]<<" "<<1900 + ltm->tm_year<<" , "<< ltm->tm_hour << "h:"<<  ltm->tm_min << "mn:"<< ltm->tm_sec<<"s"<< endl;
			sout<<endl;
			sout<<"*************************************************************************************************************"<<endl;
			sout<<endl;

         }
		// *************************************************************************************************************
		else if(strcommand=="DEFINE_COMPOSITION"){
			sout<<setw(50)<<strcommand<<" ";
			for (int i=0;i<Compo_all_el_old.size();i++) Compo_all_el_old[i]=Compo_all_el[i];
			double compo_ref=1.0;
			if (data_base_already_read) for (size_t k=0; k<el_reduced_names.size();k++) W[k]=0.0;
			if (not data_base_already_read){
				Initialize(&ceq);       				// Initialize OpenCalphad and allocate memory to the first equilibrium
				//sout<<"initialization of first equi"<<endl;
			}

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}

			string strcharname =strmyline.substr(0,i);
			sout<<strcharname<<endl;

			double factor=1.0;
			if (compo_in_percent) factor=0.01;

			string strcharname_tampon=strcharname;

			bool re_read_compo=true;
			size_t iter=0;
			while (re_read_compo){
				iter+=1;
				re_read_compo=false;
				strcharname=strcharname_tampon;
				i=strcharname.find("=");
				sout<<strcharname<<endl;
				//sout<<"i="<<i<<endl;
				while (not((i<0) or (i>strcharname.size()))){
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
									sout<<"you are not supposed to give the composition for the reference element"<<endl;
								}
								break;
							}
						}
						if (not found_el){
							sout<<"error in composition definition in line:"<<line_number<<endl;
							sout<<"element "<<element_name<<" not present in the database"<<endl;
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

							re_read_compo=true;
							sout<<" -----------> New element detected :"<<element_name<<endl;
							i=-1;
							// c_new_gtp doest not reinitialize all
							//exit(EXIT_FAILURE);
						}

					}
					if ((iter==1) and (data_base_already_read) ) re_read_compo=true;// even if the compo is the same re-initialize everything

				}

				for (size_t k=0; ((k<el_reduced_names.size()) and (not re_read_compo));k++){
					//sout<<el_reduced_names[k]<<" : "<<W[k]<<endl;
					if ((W[k]<1e-15) and (not k==i_ref)){
						re_read_compo=true;
						sout<<" -----------> element detected with zero composition :"<<el_reduced_names[k]<<endl;
					}
				}
				if (not re_read_compo){
					if (Suspended_phase_list.size()>0){
						re_read_compo=true;
						sout<<" 	all possible phases will be allowed again with new composition" <<endl;
					}
				}
				
				if (re_read_compo){
					//sout<<"re_read_compo"<<endl;
					c_new_gtp();
					//sout<<"c_new_gtp"<<endl;
					//in gtp3E ant the end of new_gtp call init_gtp(intv,dblv) must be commented
					Initialize(&ceq);
					
					//sout<<"initialize(&ceq)"<<endl;
					data_base_already_read=false;
					el_reduced_names.resize(0);
					Suspended_phase_list.resize(0);
					W.resize(0);
					Store_Equilibria.resize(0);
					Store_Equilibria_compo_unit.resize(0);
					for (size_t i=0; i<eldatabase.size();++i) Compo_all_el[i]=0.;
					compo_ref=1.0;
					//sout<<"reinitialize done"<<endl;
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




			i_compo+=1;
			i_eq=0;
			char titre_base[1024] ;
			sprintf(titre_base, "%04i", i_compo) ;


			strcomponb=strcompo+"("+titre_base+")";

			if (not data_base_already_read){
				//file<<" first composition analyzed (which determines the set of elements to be read in the database):"<<endl;
				bool firts_el=true;

				file<<"["<<strcomponb<<"]"<<TAB;

				for (size_t i=0; i<eldatabase.size();++i){
					if (Compo_all_el[i]>1e-10){

						el_reduced_names.push_back(eldatabase[i]);
						W.push_back(Compo_all_el[i]);
						file<<eldatabase[i]<<TAB<<Compo_all_el[i]/factor<<TAB;

						firts_el=false;
					}
				}
				file<<endl;

				for (size_t i=0; i<el_reduced_names.size();++i){
					if (el_reduced_names[i]==el_ref) i_ref=i;
				}
//				c_set_status_globaldata();

//				Initialize(&ceq);                                                          // Initialize OpenCalphad and allocate memory to the first equilibrium
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
				sout<<"reading phases"<<endl;
				ReadPhases(phnames, &ceq);                                                  // Read Phases data in tdb file
				sout<<" list of possible phases in the system :"<<endl;
				file<<"["<<strcomponb<<".List_of_possible_phases]"<<TAB;
				int length=0;
				for (int i=0;i<phnames.size();i++){
					sout<<" "<<phnames[i];
					file<<phnames[i]<<TAB;
					length=length+1+phnames[i].length();
					if (length>100){
						length=0;
						sout<<endl;
					}
				}
				sout<<endl;
				file<<endl;
				phfract.resize(phnames.size(),0.);

				elfract.resize(phnames.size(),vector<double>(el_reduced_names.size(),0.));

				SetPressure(1e5, &ceq);                                                     // Set Pressure

				SetMoles(1.0, &ceq);                                                          // Set Number of moles
				MU.resize(W.size(),0.);
			}
			else{
				file<<"=================================================================="<<endl;
				file<<"["<<strcomponb<<"]"<<TAB;
				bool firts_el=true;
				for (size_t i=0; i<el_reduced_names.size();++i){
					if (not (i==i_ref)){
						if (not firts_el) file<<" / ";
						file<<el_reduced_names[i]<<"="<<W[i];
						firts_el=false;
					}
				}
				file<<endl;
			}

			SetComposition(W, &ceq,i_ref, compo_unit);                                                    // Set Composition of the system
			TK=1200;                                                  //Set temperature
			SetTemperature(TK, &ceq);

			//---------------------Compute Equilibrium----------------------------
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);








		}

		//************************************************************************************************************
		// end of if(strcommand=="DEFINE_COMPOSITION")
		//************************************************************************************************************
		else if(strcommand=="LIQUIDUS"){
			sout<<setw(50)<<strcommand<<endl;
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, first > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			
			double TK=1300.;
			SetTemperature(TK, &ceq);
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
			
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);

			//ResetTemperature(&ceq);	 //remove condition on temperature


			//Change_Phase_Status(strLIQUID,PHFIXED,0.9999,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...

			//CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);  // option GRID    // this is why the previous equilibrium was at high T to have liquid present


			//if (not (i_error==0))
			{
				sout<<"*";
				double targeted_fraction=1.0-1e-4;
				double temperature_accuracy=1e-3;
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
				find_TK_for_a_given_Liquid_fraction(TK, i_error,strLIQUID,strSOLSOL,targeted_fraction, temperature_accuracy, ceq, phnames, Suspended_phase_list);
			}
			/*
			else{
				TK=ReadTemperature(&ceq);
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
				CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			}
			*/
			file<<"["<<strcomponb<<".Liquidus]"<<TAB;
			if (i_error==0){
				TK_Liquidus=TK;
				TC=TK-TCtoTK;

				sout<<" ----> liquidus is: "<<TC<<" C"<<endl;

				file<<TC<<endl;
				sout<<endl;
			}else{
				sout<<" liquidus not converged"<<endl;

				file<<"-1000"<<endl;
				sout<<endl;
			}

		}
		// *************************************************************************************************************
		else if(strcommand=="SOLIDUS"){
			sout<<setw(50)<<strcommand<<endl;
			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, first > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}



			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
			double TK=2000;
			SetTemperature(TK, &ceq);
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			/*
			ResetTemperature(&ceq);	 			//remove condition on temperature

			Change_Phase_Status(strLIQUID,PHFIXED,0.00001,&ceq);// 				   // ask the liquid phase to have an atomic fraction of 0.99...

			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);  // option GRID    // this is why the previous equilibrium was at high T to have liquid present
			if (not i_error==0)
			*/{
				sout<<"*";
				double targeted_fraction=1e-4;
				double temperature_accuracy=1e-3;
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
				find_TK_for_a_given_Liquid_fraction(TK, i_error,strLIQUID,strSOLSOL,targeted_fraction, temperature_accuracy, ceq, phnames,Suspended_phase_list);
			}
			/*
			else{
				TK=ReadTemperature(&ceq);
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);//
				CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			}
			*/
			file<<"["<<strcomponb<<".Solidus]"<<TAB;
			if (i_error==0){

				TC=TK-TCtoTK;

				sout<<" ----> solidus is: "<<TC<<" C"<<endl;
				file<<TC<<endl;
				sout<<endl;
			}else{
				sout<<" solidus not converged"<<endl;
				file<<"-1000"<<endl;
				sout<<endl;
			}


		}
		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_TRANSITION_TEMPERATURES"){
			sout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTK_start =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTK_end =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string straccuracy =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strnstep=strmyline.substr(0,i);
			sout<<strTK_start<<" / "<<strTK_end<<" / "<<straccuracy <<" / "<<strnstep<<endl;
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
			Global_Find_Transitions(strLIQUID,strSOLSOL,file,TK_start,nstep,TK_end,required_accuracy_on_TK, W, phnames,el_reduced_names,ceq,i_ref, compo_unit,ncpu, Store_Equilibria, Store_Equilibria_compo_unit,Suspended_phase_list,strcomponb);// find all the transitions temperatures for a given alloy composition
			// ************************************************************************************************************************************************************************
			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;
			sout<<"Store_Equilibria.size()="<<Store_Equilibria.size()<<endl;
			sout<<" elapsed time for the transition temperature routine (s)= "<<elapsed_time<<endl;
			sout<<endl;
			sout<<endl;

			TK_Liquidus=TK_start;
		}

		// *************************************************************************************************************
		else if(strcommand=="COMPUTE_EQUILIBRIUM"){
			sout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strT=strmyline.substr(0,i);

			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strni=strmyline.substr(0,i);

			double TK=atof(strT.c_str());
			int idetail=atoi(strni.c_str());
			sout<<strT<<" / "<<strni<<endl;
			if (temp_unit=="C"){
				TK+=TCtoTK;
			}
			SetTemperature(TK, &ceq);
			/*for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
			Change_Phase_Status(strSOLSOL,PHENTERED,0.5,&ceq);//
			Change_Phase_Status(strLIQUID,PHENTERED,0.5,&ceq);//

			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			if (not(i_error==0)) {
				for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
				Change_Phase_Status(strSOLSOL,PHENTERED,1.0,&ceq);//
				//Change_Phase_Status(strLIQUID,PHENTERED,0.5,&ceq);
				CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
				if (not(i_error==0)) {
					for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.0,&ceq);
					//Change_Phase_Status(strSOLSOL,PHENTERED,1.0,&ceq);//
					Change_Phase_Status(strLIQUID,PHENTERED,1.0,&ceq);
					CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
				}

			}*/
			Safer_CalculateEquilibrium (ceq,NOGRID,i_error,Suspended_phase_list,strLIQUID,strSOLSOL,phnames);
			char titre_base[1024] ;
			i_eq+=1;
			sprintf(titre_base, "%04i", i_eq) ;


			strEquinb=strEqui+"("+titre_base+")";
			myequi="["+strcomponb+"."+strEquinb;


			if (not(i_error==0)) {

				file<<myequi<<".Status]"<<TAB<<"Failed"<<endl;

			}
			else{
				file<<myequi<<".Status]"<<TAB<<"Good"<<endl;
				Write_Results_Equilibrium(file,el_reduced_names,phnames,phfract,elfract,ceq,idetail,compo_unit,MU,temp_unit,myequi);
			}

			sout<<endl;


		}

		// *************************************************************************************************************
		else if(strcommand=="SUSPEND_PHASES"){
		    sout<<setw(50)<<strcommand<<" ";
			file<<"["<<strcomponb<<".Phase_suspended]";
			i = strmyline.find("/");


			while(not(((i<0) or (i>strmyline.size()))))
			{
				string strphasename=strmyline.substr(0,i);
				All_Capital_Letters(strphasename);
				bool phase_found=false;
				for (int j=0;j<phnames.size() and not phase_found;j++){
					if (phnames[j]==strphasename){
						phase_found=true;
					}
				}
				if (phase_found ){
					Suspended_phase_list.push_back(strphasename);
					//Change_Phase_Status(strphasename,PHSUS,0.0,&ceq);//
					sout<<strphasename<<" suspended"<<endl;
					file<<TAB<<strphasename;
				}else
				{
					sout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
				}


				strmyline.erase(0,i+1);
				i = strmyline.find("/");
			}

			{
				i = strmyline.find(">");
				if ((i<0) or (i>strmyline.size())) {
					sout<<" error in command line, > not found in line:"<<line_number<<endl;
					exit(EXIT_FAILURE);
				}
				string strphasename=strmyline.substr(0,i);
				All_Capital_Letters(strphasename);
				bool phase_found=false;
				for (int j=0;j<phnames.size() and not phase_found;j++){
					if (phnames[j]==strphasename){
						phase_found=true;
					}
				}
				if (phase_found ){
					Suspended_phase_list.push_back(strphasename);
					//Change_Phase_Status(strphasename,PHSUS,0.0,&ceq);//
					sout<<strphasename<<" suspended"<<endl;
					file<<TAB<<strphasename;
				}else
				{
					sout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
				}
			}
			file<<endl;
		}
		// *************************************************************************************************************
		else if(strcommand=="SUSPEND_ALL_PHASES_BUT_ONE"){
			sout<<setw(50)<<strcommand<<" ";
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
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
				sout<<" all phases have been suspended but :"<<strphasename<<endl;
				file<<"["<<strcomponb<<".All_phase_suspended_but]"<<TAB<<strphasename<<endl;
			}else
			{
				sout<<"error in line "<<line_number<<" : phase does not exist"<<endl;

			}


		}
		// *************************************************************************************************************
		else if(strcommand=="RESTORE_PHASES"){
			sout<<setw(50)<<strcommand<<" ";
			i = strmyline.find(">");
			file<<"["<<strcomponb<<".Phase_restored]";
			while(not(((i<0) or (i>strmyline.size()))))
			{
				string strphasename=strmyline.substr(0,i);
				All_Capital_Letters(strphasename);
				bool phase_found=false;
				for (int j=0;j<phnames.size() and not phase_found;j++){
					if (phnames[j]==strphasename){
						phase_found=true;
					}
				}

				if (phase_found ){
					phase_found=false;
					int i_found=0;
					for (int j=0;j<Suspended_phase_list.size()and not phase_found;j++){
						if (strphasename==Suspended_phase_list[j]) {
							i_found=j;
							phase_found=true;
						}
					}
					if  (phase_found ) {
						Suspended_phase_list.erase(Suspended_phase_list.begin() +i_found);
						Change_Phase_Status(strphasename,PHENTERED,0.,&ceq);//
						sout<<strphasename <<" reactivated";
						file<<TAB<<strphasename<<endl;
					}
				}
				else
				{
					sout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
				}


				strmyline.erase(0,i+1);
				i = strmyline.find("/");
			}

			{
				i = strmyline.find(">");
				if ((i<0) or (i>strmyline.size())) {
					sout<<" error in command line, > not found in line:"<<line_number<<endl;
					exit(EXIT_FAILURE);
				}
				string strphasename=strmyline.substr(0,i);
				All_Capital_Letters(strphasename);
				bool phase_found=false;
				for (int j=0;j<phnames.size() and not phase_found;j++){
					if (phnames[j]==strphasename){
						phase_found=true;
					}
				}
				if (phase_found ){
					phase_found=false;
					int i_found=0;
					for (int j=0;j<Suspended_phase_list.size()and not phase_found;j++){
						if (strphasename==Suspended_phase_list[j]) {
							i_found=j;
							phase_found=true;
						}
					}
					if  (phase_found ) {
						Suspended_phase_list.erase(Suspended_phase_list.begin() +i_found);
						Change_Phase_Status(strphasename,PHENTERED,0.,&ceq);//
						sout<<strphasename <<" reactivated"<<endl;
						file<<TAB<<strphasename;
					}
				}
				else
				{
					sout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
				}
			}
			file<<endl;






			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
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
					sout<<strphasename <<" reactivated"<<endl;
					file<<"["<<strcomponb<<".Phase_restored]"<<TAB<<strphasename<<endl;
				}
			}
			else
			{
				sout<<"error in line "<<line_number<<" : phase does not exist"<<endl;
			}


		}
		// *************************************************************************************************************
		else if(strcommand=="RESTORE_ALL_PHASES"){
			sout<<setw(50)<<strcommand<<" ";
			Suspended_phase_list.resize(0);
			for (int i=0;i<phnames.size();i++) Change_Phase_Status(phnames[i],PHENTERED,0.,&ceq);
			Change_Phase_Status(strLIQUID,PHENTERED,1.,&ceq);
			sout<< "all phases have been reactivated"<<endl;
			file<<"["<<strcomponb<<".All_hase_restored]"<<endl;
		}

		// *************************************************************************************************************
		else if(strcommand=="SCHEIL_SOLIDIFICATION"){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, very first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strGradientFileOut =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, second / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strtarget_delta_f_liq =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_min =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, fourth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);



			sout<<strGradientFileOut<<" / "<<strtarget_delta_f_liq<<" / "<<strdelta_T_min<<" / "<<strdelta_T_max<<endl;

			double target_delta_f_liq=atof(strtarget_delta_f_liq.c_str());
			double delta_T_min=atof(strdelta_T_min.c_str());
			double delta_T_max=atof(strdelta_T_max.c_str());

			if (TK_Liquidus<30){
				sout<<" you need to have a valid liquidus temperature to start a sheill calculation"<<endl;
				sout<<" liquidus="<<TK_Liquidus<<endl;
				exit(EXIT_FAILURE);
			}

			gettimeofday(&start1, NULL);// get the present time

			scheil_solidif(strGradientFileOut,strLIQUID,strSOLSOL,file,el_reduced_names,phnames,ceq, W, target_delta_f_liq,delta_T_min,delta_T_max, TK_Liquidus,i_ref,compo_unit,Suspended_phase_list,strcomponb);

			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;

			sout<<" elapsed time for the scheil solidification routine (s)= "<<elapsed_time<<endl;
			sout<<endl;
			sout<<endl;

		}
		//***************************************************************************************************************
		else if((strcommand=="DIFF_SOLIDIFICATION")and (CTEC>0)){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			sout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, very first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strGradientFileOut =strmyline.substr(0,i);
			strmyline.erase(0,i+1);




			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strtarget_delta_f_liq =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, fourth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_min =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, fifth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, sixth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strhalf_sdas =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, seventh / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strNbincrement =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, heigth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdim =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strTpoint =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strFl_end_hot_tearing =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			
			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strFastSolidification =strmyline.substr(0,i);
			strmyline.erase(0,i+1);


			sout<<strtarget_delta_f_liq<<"/"<<strdelta_T_min<<"/"<<strdelta_T_max<<"/"<<strhalf_sdas<<"/"<<strNbincrement<<"/"<<strdim<<"/"<<strTpoint<<"/"<<strFastSolidification<<endl;

			double target_delta_f_liq=atof(strtarget_delta_f_liq.c_str());
			double delta_T_min=atof(strdelta_T_min.c_str());
			double delta_T_max=atof(strdelta_T_max.c_str());
			double half_sdas=atof(strhalf_sdas.c_str());
			int Nb_increment=atoi(strNbincrement.c_str());
			double dim=atof(strdim.c_str());
			double Tpoint=atof(strTpoint.c_str());
			double Fl_end_hot_tearing=atof(strFl_end_hot_tearing.c_str());
			int iFastSolidification=atoi(strFastSolidification.c_str());
			bool FastSolidification=false;
			if (iFastSolidification==1) FastSolidification=true;

			if (TK_Liquidus<30){
				sout<<" you need to have a valid liquidus temperature to start a sheill calculation"<<endl;
				sout<<" liquidus="<<TK_Liquidus<<endl;
				exit(EXIT_FAILURE);
			}

			gettimeofday(&start1, NULL);// get the present time
			#if CTEC>0
			back_diff_solidif(strGradientFileOut,strLIQUID,strSOLSOL,file,el_reduced_names,phnames,ceq, W, target_delta_f_liq,delta_T_min,delta_T_max, TK_Liquidus,i_ref,compo_unit,half_sdas,dim,Tpoint,Nb_increment,Store_Equilibria,Store_Equilibria_compo_unit,Suspended_phase_list,Fl_end_hot_tearing,strcomponb,element_file, FastSolidification);
			#endif
			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;
			file<<"["<<strcomponb<<".Back_diff_solidif.CPU(s)"<<TAB<<elapsed_time<<endl;
			sout<<" elapsed time for the back-diffusion solidification routine (s)= "<<elapsed_time<<endl;
			sout<<endl;
			sout<<endl;

		}
		// *************************************************************************************************************
		else if((strcommand=="HOMOGENIZING")and (CTEC>0)){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			vector < double > TC1;
			vector < double > TC2;
			vector < double > segments_time_h;
			sout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, very first / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strGradientFileIn =strmyline.substr(0,i);
			strmyline.erase(0,i+1);


			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, third / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strprintresultevery_s =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, fourth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strprintgradientevery_h =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, fifth / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strdelta_T_max =strmyline.substr(0,i);
			strmyline.erase(0,i+1);
			sout<<strprintresultevery_s<<"/"<<strprintgradientevery_h<<"/"<<strdelta_T_max<<endl;
			double printresultevery_s=atof(strprintresultevery_s.c_str());
			double printgradientevery_h=atof(strprintgradientevery_h.c_str());
			double delta_T_max=atof(strdelta_T_max.c_str());
			bool end_of_thermal_cycle_detected=false;
			sout<<setw(50)<<"thermal cycle:";
			while (not end_of_thermal_cycle_detected){
				i = strmyline.find("/");
				if ((i<0) or (i>strmyline.size())) {
					sout<<" error in command line, sixth / not found in line:"<<line_number<<endl;
					exit(EXIT_FAILURE);
				}
				string strT1 =strmyline.substr(0,i);
				strmyline.erase(0,i+1);

				i = strmyline.find("/");
				if ((i<0) or (i>strmyline.size())) {
					sout<<" error in command line, seventh / not found in line:"<<line_number<<endl;
					exit(EXIT_FAILURE);
				}
				string strT2 =strmyline.substr(0,i);
				strmyline.erase(0,i+1);


				i = strmyline.find("/");

				if ((i<0) or (i>strmyline.size())) {
					i = strmyline.find(">");

					if ((i<0) or (i>strmyline.size())) {
						sout<<" error in command line, > not found in line:"<<line_number<<endl;
						exit(EXIT_FAILURE);
					}
					end_of_thermal_cycle_detected=true;
				}
				string strtime_h =strmyline.substr(0,i);
				strmyline.erase(0,i+1);



				sout<<"/"<<strT1<<"/"<<strT2<<"/"<<strtime_h;


				double valueT1=atof(strT1.c_str());
				int valueT2=atof(strT2.c_str());
				double valuestime_h=atof(strtime_h.c_str());

				TC1.push_back(valueT1);
				TC2.push_back(valueT2);
				segments_time_h.push_back(valuestime_h);
			}
			sout<<endl;

			string strGradientFileOut = strGradientFileIn.substr(0,strGradientFileIn.length()-4); // remove .txt
			strGradientFileOut+="andhomo.txt";
			SetTemperature(1000, &ceq);
			//CalculateEquilibrium(&ceq,GRID,i_error,Suspended_phase_list);
			CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			gettimeofday(&start1, NULL);// get the present time
			#if CTEC>0
			homo(strGradientFileIn,strGradientFileOut,printresultevery_s,printgradientevery_h,strLIQUID,strSOLSOL,file,el_reduced_names,phnames,delta_T_max,i_ref,Store_Equilibria,Store_Equilibria_compo_unit,TC1, TC2,segments_time_h,Suspended_phase_list,strcomponb,element_file);
			#endif
			gettimeofday(&end1, NULL);

			seconds  = end1.tv_sec  - start1.tv_sec;
			useconds = end1.tv_usec - start1.tv_usec;

			elapsed_time = ((double)(((seconds) * 1000 + useconds/1000.0) + 0.5))/1000.;

			sout<<" elapsed time for the back-diffusion solidification routine (s)= "<<elapsed_time<<endl;
			sout<<endl;
			sout<<endl;

		}
		else if((strcommand=="PROPERTIES")and (CTEC>0)){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max
			vector < double > TC1;
			vector < double > TC2;
			vector < double > segments_time_h;
			sout<<setw(50)<<strcommand<<" ";

			int i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strFile =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			file<<myequi<<".Physical_properties.outputfile]"<<TAB;
			file<<strFile<<endl;

			sout<<endl;
			#if CTEC>0
			compute_properties(strFile,strSOLSOL,el_reduced_names, W, compo_unit,i_ref, &ceq, phnames,file,strcomponb,element_file);
			#endif

		}
		else if((strcommand=="FIX_A_PHASE")){

			sout<<setw(50)<<strcommand<<" ";
			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strPhase =strmyline.substr(0,i);
			strmyline.erase(0,i+1);



			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strvalue =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			double value=atof(strvalue.c_str());;

			Change_Phase_Status(strPhase,PHFIXED,value,&ceq);//
			sout<<endl;
		}
		else if((strcommand=="COMPUTE_EQUILIBRIUM_WITH_TEMPERATURE_CHANGED_BY")){
			//parameter target_delta_f_liq
			//paramter delta_T_min
			//paramter delta_T_max

			sout<<setw(50)<<strcommand<<" ";


			int i = strmyline.find("/");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, / not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strvalue =strmyline.substr(0,i);
			strmyline.erase(0,i+1);

			i = strmyline.find(">");
			if ((i<0) or (i>strmyline.size())) {
				sout<<" error in command line, > not found in line:"<<line_number<<endl;
				exit(EXIT_FAILURE);
			}
			string strni=strmyline.substr(0,i);
			int idetail=atoi(strni.c_str());

			double value=atof(strvalue.c_str());;
			sout<<value;
			i_eq+=1;

			char titre_base[1024] ;
			sprintf(titre_base, "%04i", i_eq) ;


			strEquinb=strEqui+"("+titre_base+")";
			myequi="["+strcomponb+"."+strEquinb;


			TK=ReadTemperature(&ceq)+value;
			SetTemperature(TK, &ceq);
			//CalculateEquilibrium(&ceq,NOGRID,i_error,Suspended_phase_list);
			Safer_CalculateEquilibrium (ceq,NOGRID,i_error,Suspended_phase_list,strLIQUID,strSOLSOL,phnames);


			if (not(i_error==0)) {

				file<<myequi<<".Status]"<<TAB<<"Failed"<<endl;
			}
			else{
				file<<myequi<<".Status]"<<TAB<<"Good"<<endl;
				Write_Results_Equilibrium(file,el_reduced_names,phnames,phfract,elfract,ceq,idetail,compo_unit,MU,temp_unit,myequi);
			}

			sout<<endl;



		}
		else if((strcommand=="")){
		}
		// *************************************************************************************************************
		else
		{
			sout<<setw(50)<<strcommand<<" ";
			sout<<"command line not recognized at line:"<<line_number<<endl;

		}

		line_number+=1;
		//sout<<" line:"<<line_number<<endl;
	}
	
			
			return -1;
	

	file.close();
    return 1;
}
