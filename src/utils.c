/*  		
		4P v1.0 - Parallel Processing of Polymorphism Panels
    Copyright (C) 2013  Andrea Benazzo, Alex Panziera and Giorgio Bertorelle 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "struct.h"

struct ind_gen;

void logo()
/*This function print the logo to the standard output*/
{
printf("----------------------------------------------------------------------------------------\n");
printf("|                            4P v1.0 - command line                     	       |\n");
printf("|                                                                                      |\n");
printf("|        Software for Parallel Processing of Polymorphism Panels created by:           |\n");
printf("|               Andrea Benazzo, University of Ferrara, Italy                           |\n");
printf("|                Alex Panziera, University of Ferrara, Italy		               |\n");
printf("|             Giorgio Bertorelle, University of Ferrara, Italy		               |\n");
printf("|                                                                                      |\n");
printf("----------------------------------------------------------------------------------------\n");
}


void help()
{
printf("4P HELP PAGE:\n");
printf("\t-f <file>\tGenotypes input file name\n");
printf("\t-i <number>\tInput file type (0:plink,1:arlequin,2:vcf)\n");
//printf("\t-g <number>\tInput type of data (genotipic=0/haplotipic=1)\n");
printf("\t-m <file>\tMap input file name (only for plink files)\n");
printf("\t-n <number>\tNumber of individuals\n");
printf("\t-s <number>\tNumber of SNPs\n");
printf("\t-p <file>\tPopulations input file name (only for vcf files)\n");
printf("\t-a <file>\tAncestral allele input file name (only for plink and vcf files)\n");
printf("\t-t <number>\tNumber of cores used for computations\n");
printf("\t-v\t\tOutput genotype informations to the standard output\n");
printf("\t-h\t\tPrint this help\n");
printf("\nUsage: 4p -f <inputfile> -m <mapfile> -i <inputtype> -n <nind> -s <nsnp> -p <populationfile> -a <ancallfile> -t <ncores>\n");
}


void choppy(char *s)
{
	s[strcspn(s,"\n")]='\0';
}


int uniquePop(struct ind_gen *matgen, int n, int *pos, int *u_ptr)
{
/*THIS FUNCTION NEEDS INPUT CONTROLS*/
/*This function find the number of different populations among all individuals
entered in the input matrix. The pointer pos return the positions of the unique
populations and the pointer u_ptr return the number of unique pop. The function
returns 0 if it is all right*/
int i,j,k=0;
int ck,ck1;
char *prev;

prev=(char*)malloc(sizeof(matgen[0].pop));
pos[0]=0;
*(u_ptr)=1;
strcpy(prev,matgen[0].pop);
for (i=1;i<n;i++){
	if (ck1=strcmp(matgen[i].pop,prev)!=0){
		ck=1;
		for (j=0;j<*u_ptr;j++){
			if (strcmp(matgen[pos[j]].pop,matgen[i].pop)==0){ck=0;} 
		}
	if (ck==1){k++;pos[k]=i;(*u_ptr)++;strcpy(prev,matgen[i].pop);}
	}
}
free(prev);
return 0;
}

char flipBase(char base)
/*this function return the allele state after a mutation event (transition)*/
{
char flip;
if (base=='A'){flip='G';}
if (base=='G'){flip='A';}
if (base=='T'){flip='C';}
if (base=='C'){flip='T';}
return flip;
}

int outGen(struct ind_gen *matgen,long int row, long int col)
/*Output the genotype informations to the standard output (first 50 snps are shown)*/
{
long int i,j,nomit;
printf("Reading the loaded genotypes:\n");
for (i=0;i<row;i++){
	printf("POPID: %s\t",matgen[i].pop);
	printf("NAME: %s\t",matgen[i].name_ind);
	printf("PAT_ID: %s\t",matgen[i].paternal_id);
	printf("MAT_ID: %s\t",matgen[i].maternal_id);
	printf("SEX: %d\t",matgen[i].sex);
	printf("PHEN: %d\t",matgen[i].phen);
	printf("PHASE: %d\t",matgen[i].phase);
	printf("DP: %d\t",matgen[i].dp);
	printf("FT: %s\t",matgen[i].ft);
	printf("GL: %s\t",matgen[i].gl);
	printf("GLE: %s\t",matgen[i].gle);
	printf("PL: %s\t",matgen[i].pl);
	printf("GP: %s\t",matgen[i].gp);
	printf("GQ: %s\t",matgen[i].gq);
	printf("HQ: %s\t",matgen[i].hq);
	printf("PS: %s\t",matgen[i].ps);
	printf("PQ: %s\t",matgen[i].pq);
	printf("EC: %s\t",matgen[i].ec);
	printf("MQ: %s\t",matgen[i].mq);
	printf("\n");
	nomit=col-50;
	if (nomit<1){
		for (j=0;j<col;j++){
			if (j==0){printf("GEN1:");}
			printf("%c",matgen[i].gen1[j]);
		}
		printf("\n");
		for (j=0;j<col;j++){
			if (j==0){printf("GEN2:");}
			printf("%c",matgen[i].gen2[j]);
		}
		printf("\n---------------------------------------------------------------------------------------------------------\n");
	}
	else{
		for (j=0;j<50;j++){
			if (j==0){printf("GEN1:");}
			printf("%c",matgen[i].gen1[j]);
		}
		printf(" - %ld genotypes omitted\n",nomit);
		for (j=0;j<50;j++){
			if (j==0){printf("GEN2:");}
			printf("%c",matgen[i].gen2[j]);
		}
		printf(" - %ld genotypes omitted",nomit);
		printf("\n---------------------------------------------------------------------------------------------------------\n");
	}
}
return 0;
}


int outMap(struct ind_map *matmap, long int col)
/*This function print to to standard output the informations loaded in the map structure*/
{
long int i;
printf("SNP information loaded from mapfile:\n");
for (i=0;i<col;i++){
	printf("SNP nr: %ld\tCr: %d\trs: %s\tdist(Morgans): %d\tbp(position): %ld\tancAll: %c\trefAll: %c\taltAll: %s\tquality: %d\tfilter: %s\tinfo: %s\n",(i+1),matmap[i].crom,matmap[i].rs,matmap[i].gen_dist,matmap[i].bp_pos,matmap[i].anc,matmap[i].ref,matmap[i].alt,matmap[i].qual,matmap[i].filt,matmap[i].info);
}
return 0;
}

char getAll(char ref, char *alt, int x)
/*This function take the xth allele present in the matgen.alt field for a SNP
or take the reference if x is 0. If no valid conversion could be performed, a 'N' value is returned. */
{
char alt1[20];
char *s, *sps;
unsigned int i;

strcpy(alt1,alt);
if (x<0){printf("\nERROR: getAll: negative position detected: %d\n",x);return 'n';}
else{
	if (x==0){return ref;}
	else{
		if (x==1){return alt1[0];}
		else{
			s=strtok_r(alt1,",",&sps);//printf("\nfirst token: %s\n",s);
			if (s==NULL){printf("\nERROR: getAll: impossible to take the %dth allele from alt: %s, only one present.\n",x,alt1);return 'n';}
			i=2;
			while (i==x){
					s=strtok_r(NULL,",",&sps);//printf("\nsecond token: %s\n",s);
					if ((s==NULL) && (i<x)){printf("\nERROR: getAll: impossible to take the %dth allele from alt: %s, position out of bounds.\n",x,alt1);return 'n';}
					if (i==x){return s[0];}
					i++;
					}
			}
		}
	}
}

int updateAllInfo(struct ind_gen *matgen, struct ind_map *matmap, long int col, long int row, int inputType){
/*This function check the genotypes information stored in a ind_gen structure and update the ref/alt fields in the ind_map structure.*/
long int i,j,max;
short int k,pos;
char ref;
long int *tot;//A,G,T,C counter
const char base[5]="AGTC\0";
tot=malloc(4*sizeof(long int));
if (tot==NULL){printf("ERROR: updateAllInfo_mp: memory allocation failure.\n");return 1;}
if (inputType==1){
	for (i=0;i<col;i++){
		for (k=0;k<4;k++){tot[k]=0;}
		matmap[i].anc='A';
		for (j=0;j<row;j++){
			if (matgen[j].gen1[i]=='A'){tot[0]++;}
			else if (matgen[j].gen1[i]=='G'){tot[1]++;}
			if (matgen[j].gen2[i]=='A'){tot[0]++;}
			else if (matgen[j].gen2[i]=='G'){tot[1]++;}
		}
		if (tot[0]>tot[1]){matmap[i].ref='A';strcpy(matmap[i].alt,"G\0");}
		else{matmap[i].ref='G';strcpy(matmap[i].alt,"A\0");}
	}
}
else{
	for (i=0;i<col;i++){
		for (k=0;k<4;k++){tot[k]=0;}
		for (j=0;j<row;j++){
			if (matgen[j].gen1[i]=='A'){tot[0]++;}
			else if (matgen[j].gen1[i]=='G'){tot[1]++;}
			else if (matgen[j].gen1[i]=='T'){tot[2]++;}
			else if (matgen[j].gen1[i]=='C'){tot[3]++;}
			if (matgen[j].gen2[i]=='A'){tot[0]++;}
			else if (matgen[j].gen2[i]=='G'){tot[1]++;}
			else if (matgen[j].gen2[i]=='T'){tot[2]++;}
			else if (matgen[j].gen2[i]=='C'){tot[3]++;}	
		}
		max=0;
		for (k=0;k<4;k++){
			if (tot[k]>max){
				max=tot[k];
				ref=base[k];
			}
		}
		matmap[i].ref=ref;
		//printf("\nA: %ld G: %ld T: %ld C: %ld",tot[0],tot[1],tot[2],tot[3]);
		pos=0;
		for (k=0;k<4;k++){
			if ((base[k]!=ref)&(tot[k]>0)){	
				matmap[i].alt[pos]=base[k];
				pos++;
			}
		}
		if (pos>0){matmap[i].alt[pos]='\0';}
		else{strcpy(matmap[i].alt,"N\0");}
	}
}
free(tot);
return 0;
}

int updateAllInfo_mp(struct ind_gen *matgen, struct ind_map *matmap, long int col, long int row, int inputType, int nt){
/*This function check the genotypes information stored in a ind_gen structure and update the ref/alt fields in the ind_map structure - OpenMP accellerated.*/
long int i,j,max;
short int k,pos;
char ref;
long int *tot;//A,G,T,C counter
const char base[5]="AGTC\0";
int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/
tot=malloc(4*sizeof(long int));
if (tot==NULL){printf("ERROR: updateAllInfo_mp: memory allocation failure.\n");return 1;}
if (inputType==1){
	for (k=0;k<4;k++){tot[k]=0;}
	#pragma omp parallel for private(i,k,j) firstprivate(tot)
	for (i=0;i<col;i++){
		for (k=0;k<4;k++){tot[k]=0;}
		matmap[i].anc='A';
		for (j=0;j<row;j++){
			if (matgen[j].gen1[i]=='A'){tot[0]++;}
			else if (matgen[j].gen1[i]=='G'){tot[1]++;}
			if (matgen[j].gen2[i]=='A'){tot[0]++;}
			else if (matgen[j].gen2[i]=='G'){tot[1]++;}
		}
		if (tot[0]>tot[1]){matmap[i].ref='A';strcpy(matmap[i].alt,"G\0");}
		else{matmap[i].ref='G';strcpy(matmap[i].alt,"A\0");}
	}
}
else{
	for (k=0;k<4;k++){tot[k]=0;}
	max=0;pos=0;
	ref='N';
	#pragma omp parallel for private(i,k,j,max,ref,pos) firstprivate(tot)
	for (i=0;i<col;i++){
		for (k=0;k<4;k++){tot[k]=0;}
		for (j=0;j<row;j++){
			if (matgen[j].gen1[i]=='A'){tot[0]++;}
			else if (matgen[j].gen1[i]=='G'){tot[1]++;}
			else if (matgen[j].gen1[i]=='T'){tot[2]++;}
			else if (matgen[j].gen1[i]=='C'){tot[3]++;}
			if (matgen[j].gen2[i]=='A'){tot[0]++;}
			else if (matgen[j].gen2[i]=='G'){tot[1]++;}
			else if (matgen[j].gen2[i]=='T'){tot[2]++;}
			else if (matgen[j].gen2[i]=='C'){tot[3]++;}	
		}
		max=0;
		for (k=0;k<4;k++){
			if (tot[k]>max){
				max=tot[k];
				ref=base[k];
			}
		}
		matmap[i].ref=ref;
		//printf("\nA: %ld G: %ld T: %ld C: %ld",tot[0],tot[1],tot[2],tot[3]);
		pos=0;
		for (k=0;k<4;k++){
			if ((base[k]!=ref)&(tot[k]>0)){	
				matmap[i].alt[pos]=base[k];
				pos++;
			}
		}
		if (pos>0){matmap[i].alt[pos]='\0';}
		else{strcpy(matmap[i].alt,"N\0");}
	}
}
free(tot);
return 0;
}

void structDeall(struct ind_gen *matgen, struct ind_map *matmap, long int row){
long int i;
for (i=0;i<row;i++){
	free(matgen[i].gen1);
	free(matgen[i].gen2);
	}
free(matgen);
free(matmap);
}

