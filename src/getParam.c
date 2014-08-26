/*  		
		4P v1.0 - Parallel Processing of Polymorphism Panels
    Copyright (C) 2013  Andrea Benazzo, Alex Panziera and Giorgio Bertorelle 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

int getParamChar(char *param, const char *flag, int argc, char *argv[])
{
/*this function read input file name from argv vector (parameters passed via command line)
specified with the "-f" flag and also check if there are more than one input file specified*/
int i;
int ok=1;
for (i=0;i<argc;i++){
	if (strncmp(argv[i],flag,2)==0){
		if (ok==1){
			if ((strcmp(flag,"-h")!=0)&(strcmp(flag,"-v")!=0)&(strcmp(flag,"-y")!=0)){
				strcpy(param,argv[(i+1)]);
				ok=0;
			}
			else{
				strcpy(param,"null");
				ok=0;
			}
		}
		else{
			perror("More than 1 input parameter of the same type!");
			return 1;
		}
	}
}
return ok;
}

int checkSumstatPar(char name[], int *par1, float *par2, char *type, const int npar)
/* This function open the file "sumstat.par" and check the presence
of specific parameters inside */
{
FILE *parfile;
char buf[100];
char *buf_split;
int ok=0;
int i;
parfile=fopen("./sumstat.par","r");
if( parfile==NULL ){
	printf("ERROR: opening input file!\n");
       	perror("Impossible to open sumstat.par\n");
       	return 0;
}
while ( fgets(buf,sizeof(buf),parfile) !=NULL ){
	choppy(buf);
	if ((buf[0]=='/')&&(buf[1]='/')){continue;}
	buf_split=strtok(buf,"#");
	i=0;
	while (buf_split != NULL){
		i++;
		if (strncmp(buf_split,name,strlen(name))==0){ok=1;}
		else{
			if ((ok==1)&&(npar==2)){
				if (i==2){*par1=atoi(buf_split);}
				if (i==3){strncpy(type,buf_split,strlen(type));return 0;}
			}
			else if ((ok==1)&&(npar==3)){
				if (i==2){*par1=atoi(buf_split);}
				if (i==3){*par2=atof(buf_split);}
				if (i==4){strncpy(type,buf_split,strlen(type));return 0;}
			}
		}//else
	buf_split=strtok(NULL,"#");
	}//2while
}//1while
fclose(parfile);
}//end
