//Johnrob Y. Bantang, Natinal Institute of Physics
//Created: 03 October 2002
// Makes new C files
// usage: newC filename
//Modifications:
// >> 21 Jan 2003, Johnrob
//   included the constant AUTHOR and AFFILIATION for portability

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>
#include <time.h>
#include <string.h>

const char *const AUTHOR= "Johnrob Y. Bantang";
const char *const AFFILIATION= "National Institute of Physics";
const char *const EXTENSION= ".cpp";

int main(int argc,char **argv){
	if(argc!=2){
		cout<<"usage: newC filename"<<endl;
		exit(0);
	}
	char *fname= new char[strlen(argv[1])+5];
	ostrstream out(fname,strlen(argv[1])+5);
	out<<argv[1]<<EXTENSION;
	time_t date=time(NULL);
	
	ofstream file(fname,ios::nocreate,0);
	//opens normal file that **already exists**;
	if(file){
		for(int i=1;i<5;i++)
			cout<<"WARNING! file aready exists!"<<endl;
		cout<<endl<<"you can type"<<endl<<endl;
		cout<<"\t\"head "<<fname<<"\""<<endl;
		cout<<"in command line to *view version*"<<endl<<endl;
		cout<<"please enter 1 to OVERWRITE this file"<<endl;
		cout<<"type anything to cancel"<<endl;
		for(int i=1;i<5;i++)
			cout<<"WARNING! file aready exists!"<<endl;
		int n;
		cin>>n;
		if(n!=1){
			cout<<"*no* file is created... exiting..."<<endl<<endl;
			exit(0);
		}
		file.close();
		file.open(fname);
		if(!file)
			cout<<"**cannot create new file!!**"<<endl;
		cout<<"OLD FILE: "<<fname<<" *overwritten!*"<<endl;
	}
	if(!file){
		file.open(fname);
		if(!file){
			cout<<"**cannot create new file!!**"<<endl;
			exit(0);
		}
		cout<<endl<<"NEW FILE created: "<<fname<<endl<<endl;
	}
	
	cout<<"\tcreating contents for the new C++ file: "<<endl<<endl;
	//creating headers...
	file<<"//filename: \""<<fname<<"\""<<endl;
	file<<"//"<<AUTHOR<<", "<<AFFILIATION<<endl;
	file<<"//Created: "<<asctime( localtime(&date) );
		//writes the time and date today; endl already in asctime();
	file<<"//"<<endl;
	file<<"//Comments:"<<endl;
	file<<"// >>"<<endl;
	file<<"//"<<endl;
	file<<"//This file is generated using the \"newC generator\"..."<<endl;
	file<<"//Modifications:"<<endl;
	file<<"// >>"<<endl<<endl;
	file<<"#include <iostream.h>"<<endl;
	file<<"#include <math.h>"<<endl<<endl;
	//starting the main body...
	file<<"int main(int argc,char **argv){"<<endl;
	file<<"//\tif(argc!= ...){"<<endl;
	file<<"//\t\tcout<<\"usage: "<<argv[1]<<".exe ... \"<<endl;"<<endl;
	file<<"//\t\texit(1);"<<endl;
	file<<"//\t}"<<endl;
	file<<"\t//write the main body here"<<endl;
	file<<"return(0);"<<endl<<"}"<<endl;
	
	delete fname;
  file.close();
  cout<<endl<<"\tCREATION SUCCESSFUL!"<<endl<<endl;
return 0;
}
