// dummy C routine that returns a legal filename 

#include <stdio.h>
#include <string.h>
#include "tinyfiledialogs.h"

char const * tinyopen(
		      int const typ)
{
  char const * lFilterPatterns1[1] = {"*.TDB"};
  char const * lFilterPatterns2[1] = {"*.ocu"};
  char const * lFilterPatterns3[1] = {"*.OCM"};
  char const * lFilterPatterns4[1] = {"*.OCD"};
  char const * lFilterPatterns5[1] = {"*.plt"};
  char const * lFilterPatterns6[1] = {"*.PDB"};
  char const * p2;
  //printf("start of tinydummy \n");
  //printf("input value of typ: %i \n",typ);
  //printf("now copy the string \n");
  //strcpy(filename,"C:\\User\\Bosse\\Document\\Software\\openfile\\test.TDB");
  if(typ<0)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Output file name",
					       "",
					       1,
					       lFilterPatterns1,
					       NULL,
					       0);
    }
  else if(typ==1)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       1,
					       lFilterPatterns1,
					       NULL,
					       0);
    }
  else if(typ==2)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       1,
					       lFilterPatterns2,
					       NULL,
					       0);
      //p2="C:\\User\\Bosse\\Document\\Software\\openfile\\test.UNF";
    }
  else if(typ==3)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       1,
					       lFilterPatterns3,
					       NULL,
					       0);
      //p2="C:\\User\\Bosse\\Document\\Software\\openfile\\test.UNF";
    }
  else if(typ==4)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       1,
					       lFilterPatterns4,
					       NULL,
					       0);
      //p2="C:\\User\\Bosse\\Document\\Software\\openfile\\test.UNF";
    }
  else if(typ==5)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       1,
					       lFilterPatterns5,
					       NULL,
					       0);
      //p2="C:\\User\\Bosse\\Document\\Software\\openfile\\test.UNF";
    }
  else if(typ==6)
    {
      //lTheOpenFileName = tinyfd_openFileDialog(
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       1,
					       lFilterPatterns6,
					       NULL,
					       0);
      //p2="C:\\User\\Bosse\\Document\\Software\\openfile\\test.UNF";
    }
  else
    {
      //no default extension
      p2 = tinyfd_openFileDialog(
					       "Input file name",
					       "",
					       0,
					       NULL,
					       NULL,
					       0);
      //p2="C:\\User\\Bosse\\Document\\Software\\openfile\\test.DAT";
    }
  //if (! lTheOpenFileName)
  if (! p2)
    {
      tinyfd_messageBox(
			"Error",
			"Open file name is NULL",
			"ok",
			"error",
			1);
      return NULL ;
    }
  //printf("return name: %s \n",p2);
  //printf("end of tinydummy \n");
  //return lTheOpenFileName;
  return p2;
}
