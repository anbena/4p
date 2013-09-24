/*		4P v1.0 - Parallel Processing of Polymorphism Panels
    Copyright (C) 2013  Andrea Benazzo, Alex Panziera and Giorgio Bertorelle

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "getParam.h"
#include "genStat.h"
#include "struct.h"
#include "utils.h"
#include "parsFiles.h"

int main(int argc, char *argv[]){
	/*variables for internal storage of the data*/
	struct ind_gen *matgen;
	struct ind_map *matmap;
	/*variables useful for reading example input file*/
	FILE *i_f, *i_m, *i_p, *i_a;//pointer to the input file
	FILE *outf;//pointer for output file creation
	/*main variables*/
	long int i,j,k;
	char ifilename[100],imapfilename[100],ipopfilename[100],iancfilename[100];
	char buf_input[100];
	char h[5],v[5];//var per help
	int chk;
	//int ck_gen;//variable for ploidy
	long int row;
	long int col;
	int inputType=10;//initialized to 10 but only 0,1,2 are permitted
	/*pointer for population distances computation*/
	double *outv;
	/*pointers for finding unique pops*/
	int *pos;//pointer for find unique populations positions
	int upop;//number of unique populations
	int *u_ptr;
	int *n;
	/*variables useful for check sumstat file*/
	int *par1;
	float *par2;
	char type[4];
	//variables for mp parallelization
	int nt=0;

	/*memory allocation for checking sumstat.par file*/
	par1=malloc(sizeof(int));
	par2=malloc(sizeof(float));

	/*PRINT LOGO*/
	logo();

	/* reading command line parameters and initialising main variables */
        chk=getParamChar(h,"-h",argc,argv);
	if (chk==0) {
		help();
		return 0;
	}
	
	chk=getParamChar(buf_input,"-i",argc,argv);
	if ((chk!=0)||(atoi(buf_input)<0)||(atoi(buf_input)>2)) {
		printf("ERROR: incorrect input file type(%d)!(0:plink,1:arlequin,2:vcf)\n",atoi(buf_input));
		return 0;
	}
	inputType=atoi(buf_input);
	printf("-i, Genotypes input file type:\t[%d]\n",inputType);

	chk=getParamChar(ifilename,"-f",argc,argv);
	if (chk!=0) {
		printf("ERROR: when reading genotypes input file name!\n");
		return 0;
	}
	i_f=fopen(ifilename,"r");
	printf("-f, Genotypes input file name:\t[%s]\n",ifilename);
	if( i_f==NULL ){
        	perror("Error opening genotypes input file.");
        	return 0;
	}
	
	chk=getParamChar(imapfilename,"-m",argc,argv);
	if (chk!=0 & inputType==0) {
		printf("ERROR: when reading map input file name!\n");
		return 0;
	}else{
		if (inputType==0){
			i_m=fopen(imapfilename,"r");
			printf("-m, Map input file name:\t[%s]\n",imapfilename);
			if( i_m==NULL ){
        			perror("Error opening map input file.");
        			return 0;
			}
		}else{
			printf("-m, Map input file name:\t[none]\n");
		}
		
	}

	chk=getParamChar(ipopfilename,"-p",argc,argv);
	if (chk!=0 & inputType==2) {
		printf("ERROR: when reading populations input file name! Note: a population file MUST be specified with vcf files.\n");
		return 0;
	}else{
		if (inputType==2){
			i_p=fopen(ipopfilename,"r");
			printf("-p, Population input file name:\t[%s]\n",ipopfilename);
			if( i_p==NULL ){
        			perror("Error opening populations input file.");
        			return 0;
			}
		}else{
			printf("-p, Population input file name:\t[none]\n");
		}
		
	}
	
	chk=getParamChar(iancfilename,"-a",argc,argv);
	if (chk==0 & inputType!=1) {
		i_a=fopen(iancfilename,"r");
		printf("-a, Ancestral allele file name:\t[%s]\n",iancfilename);
		if( i_a==NULL ){
        		perror("Error opening ancestral allele input file.");
        		return 0;
		}
		
	}
	else {
		if (chk==0 & inputType==1){
		printf("ERROR: when reading ancestral allele file name! Note: this file is only used for plink or vcf files.\n");
		return 0;
		}
		if (chk!=0){
			i_a=NULL;
			printf("-a, Ancestral allele file name:\t[none]\n");
		}
	}

	chk=getParamChar(buf_input,"-n",argc,argv);
	if ((chk!=0)||(atoi(buf_input)<1)) {
		printf("ERROR: when reading number of individuals!\n");
		return 0;
	}
	row=atol(buf_input);
	printf("-n, Number of individuals:\t[%ld]\n",row);

	chk=getParamChar(buf_input,"-s",argc,argv);
	if ((chk!=0)||(atoi(buf_input)<1)) {
		printf("ERROR: when reading input number of SNPs!\n");
		return 0;
	}
	col=atol(buf_input);
	printf("-s, Number of SNPs:\t\t[%ld]\n",col);

	/*chk=getParamChar(buf_input,"-g",argc,argv);
	if ((chk!=0)||(atoi(buf_input)>=1)||(atoi(buf_input)<0)) {
		printf("ERROR: when reading input type of data (0: genotipic/1: haplotipic)(only genotipic data are supported now)!\n");
		return 0;
	}
	ck_gen=atoi(buf_input);
	printf("-g, Ploidy (0:geno/1:haplo):\t%d\n",ck_gen);*/
	
	nt=omp_get_num_procs();
	chk=getParamChar(buf_input,"-t",argc,argv);
	if ((chk!=0)||(atoi(buf_input)<1)||(atoi(buf_input)>nt)) {
		//nt=omp_get_num_procs();
		printf("-t, Number of threads:\t\t[None specified, using %d threads]\n",nt);
	}
	else {
		nt=atoi(buf_input);
		if (nt>col){nt=(int)col;}
		printf("-t, Number of threads:\t\t[%d]\n",nt);
	}

	/*Output the genotypes iformation to the standard output*/
        chk=getParamChar(v,"-v",argc,argv);
	if (chk==0){printf("-v, Output genotype info:\t[yes]\n");}
	else{printf("-v, Output genotype info:\t[no]\n");}

	/* Defininig structure vectors to allocate inidividuals informations*/
	printf("Allocating internal data structures\t\t");
	matgen=malloc(row*sizeof(struct ind_gen));
	if (matgen==NULL){
			perror("ERROR: error allocating memory for individuals storage.\n");
			return 0;
			}
	matmap=malloc(col*sizeof(struct ind_map));
	if (matmap==NULL){
			perror("ERROR: error allocating memory for SNPs informations storage.\n");
			return 0;
			}
	printf("[done]\n");

	printf("Allocating genotypes\t\t\t\t");
	/*memory allocation to load genotypes and initialise it to 'N'*/
	for (i=0;i<row;i++){
		matgen[i].gen1=malloc((col+1)*sizeof(char));
		if (matgen[i].gen1==NULL){perror("ERROR: error allocating memory for genotypes storage in crom1.\n");return 0;}
		matgen[i].gen2=malloc((col+1)*sizeof(char));
		if (matgen[i].gen2==NULL){perror("ERROR: error allocating memory for genotypes storage in crom2.\n");return 0;}
		for (k=0;k<col;k++){
			matgen[i].gen1[k]='N';
			matgen[i].gen2[k]='N';
			}
		}
	printf("[done]\n");

	/*Initialise individuals - BRING IT OUTSIDE THE MAIN CODE*/
	printf("Initialising individuals\t\t\t");
	for (i=0;i<row;i++){
		strcpy(matgen[i].pop,"NULL");
		strcpy(matgen[i].name_ind,"NULL");
		strcpy(matgen[i].paternal_id,"NULL");
		strcpy(matgen[i].maternal_id,"NULL");
		matgen[i].sex=0;
		matgen[i].phen=0;
		matgen[i].phase=0;
		matgen[i].dp=-1;
		strcpy(matgen[i].ft,"NULL");
		strcpy(matgen[i].gl,"NULL");
		strcpy(matgen[i].gle,"NULL");
		strcpy(matgen[i].pl,"NULL");
		strcpy(matgen[i].gp,"NULL");
		strcpy(matgen[i].gq,"NULL");
		strcpy(matgen[i].hq,"NULL");
		strcpy(matgen[i].ps,"NULL");
		strcpy(matgen[i].pq,"NULL");
		strcpy(matgen[i].ec,"NULL");
		strcpy(matgen[i].mq,"NULL");
		}
		//outGen(matgen,row,col);
	printf("[done]\n");

	/*Initialise markers - BRING IT OUTSIDE THE MAIN CODE*/
	printf("Initialising genetic loci\t\t\t");
	for (i=0;i<col;i++){
		matmap[i].crom = 0;// short unsigned int crom;
		strncpy(matmap[i].rs , "-",1);// char rs[20];
		//strcpy(matmap[i].rs , "NULL");// char rs[20];
		//if(snprintf(matmap[i].rs,20,"%ld",i)<0){printf("ERROR: Init: failure initialising loci.\n");return 1;}//set the rs to the snp position
		matmap[i].gen_dist = 0;// short int gen_dist;
		matmap[i].bp_pos = i;//long int bp_pos;
		matmap[i].anc = 'N';// char anc[1];
		matmap[i].ref = 'N';// char ref[1];
		strcpy(matmap[i].alt , "NULL");// char alt[20];//A,T,C,G,N
		matmap[i].qual = 0;// int qual;
		strcpy(matmap[i].filt , "NULL");//char filt[10];
		strcpy(matmap[i].info , "NULL");//char info[100];
	}
	printf("[done]\n");

	/*LOAD INDIVIDUALS FROM PLINK PED/MAP FILE*/
	if (inputType==0){	
		if (loadPedSnp(matgen, col, row, i_f)!=0){perror("ERROR: Error attempting to load plink ped input file!\n");return 0;}
		fclose(i_f);
		if (loadMapSnp(matmap, col, i_m)!=0){perror("ERROR: Error attempting to load plink map input file!\n");return 0;}
		fclose(i_m);
		if (i_a!=NULL){ 
			if (loadAncAll(matmap, col, i_a)!=0){perror("ERROR: Error attempting to load ancestral allele informations from input file!\n");return 0;}
			fclose(i_a);
		}
		if (updateAllInfo(matgen, matmap, col, row, inputType)!=0){perror("ERROR: updateAllInfo: Error updating alleles state information!\n");return 0;}
	}
	
	/*LOAD INDIVIDUALS FROM ARLEQUIN ARP FILE*/
	if (inputType==1){
		if (loadArlSnp(matgen, col, row, i_f)!=0){perror("ERROR: Error attempting to load arlequin input file!\n");return 0;}
		fclose(i_f);
		if (nt<2){if (updateAllInfo(matgen, matmap, col, row, inputType)!=0){perror("ERROR: updateAllInfo: Error updating alleles state information!\n");return 0;}}
		else{if (updateAllInfo_mp(matgen, matmap, col, row, inputType, nt)!=0){perror("ERROR: updateAllInfo_mp: Error updating alleles state information!\n");return 0;}}
	}

	/*LOAD INDIVIDUALS FROM VARIANT CALL FORMAT FILE (VCF 4.1)*/
	if (inputType==2){
		if (loadVcfSnp(matgen, matmap, col, row, i_f)!=0){perror("ERROR: Error attempting to load vcf input file!\n");return 0;}
		fclose(i_f);
		if (loadPop(matgen, row, i_p)!=0){printf("ERROR: reading informations from populations input file!\n");return 0;}
		fclose(i_p);
		if (i_a!=NULL){ 
			if (loadAncAll(matmap, col, i_a)!=0){perror("ERROR: Error attempting to load ancestral allele informations from input file!\n");return 0;}
			fclose(i_a);
		}
	}

	/*PRINT GENOTYPE INFORMATIONS FROM PEDFILE*/
	if(strcmp(v,"null")==0){outGen(matgen,row,col);}
	/*PRINT MAP INFORMATIONS*/
	if(strcmp(v,"null")==0){outMap(matmap,col);}

	/*Find unique populations - pos: int vector with the position of the first individual of each different pop; upop: number of populations (int) */
	pos=(int*)malloc(row*sizeof(int));
	upop=0;
	u_ptr=&upop;
	if (uniquePop(matgen, row, pos,u_ptr)==0){
		n=malloc(upop*sizeof(int));
		printf("Summary of the %d different populations found:\n",upop);
		for (i=0;i<upop;i++){
			chk=0;
			printf("Population name [%s]",matgen[pos[i]].pop);
			for (j=0;j<row;j++){
				if (strcmp(matgen[j].pop,matgen[pos[i]].pop)==0){chk++;}
			}
			n[i]=chk;
			printf(" composed by %d individuals.\n",chk);
		}
	}
	else{perror("ERROR: impossible to indentify unique populations");return 0;}
	
	/*CHECK PARAMETERS IN SUMSTAT FILE AND BASE FREQUENCIES COMPUTATION*/
	*par1=-1;
	*par2=-1;
	strcpy(type,"NUL");
	if (checkSumstatPar("BASEFREQ", par1, par2, type, 2)==0){
		if (*par1==1){
			if (upop==1){
				if (computeBaseFreq(matgen,matmap,row,col,"whole",type,nt)!=0){
					printf("ERROR: impossible to compute whole population base frequencies!\n");
					return 0;
				}
			}
			else if(upop>1){
				for (i=0;i<upop;i++){
					if (computeBaseFreq(matgen,matmap,row,col,matgen[pos[i]].pop,type,nt)!=0){
						printf("ERROR: impossible to compute population specific base frequencies!\n");
						return 0;
					}
				}
				if (computeBaseFreq(matgen,matmap,row,col,"whole",type,nt)!=0){
					printf("ERROR: impossible to compute whole population base frequencies!\n");
					return 0;
				}
			}
		}
		else{printf("Skip Base frequency Computation\n");}
	}
	else{perror("ERROR: reading the BASEFREQ parameters from sumstat.par");}
	

	/*CHECK PARAMETERS IN SUMSTAT FILE AND OBSERVED & EXPECTED HETEROZIGOSITY COMPUTATION*/
	*par1=-1;
	*par2=-1;
	strcpy(type,"NUL");
	if (checkSumstatPar("HET", par1, par2, type, 3)==0){
		if (*par1==1){
			if (upop==1){
				if (hetOss(matgen,matmap,"whole",col,row,type, nt)!=0){
					printf("ERROR: impossible to compute whole population observed heterozigosity!\n");
					return 0;
				}
			}
			else if (upop>1){
				for (i=0;i<upop;i++){
					if (hetOss(matgen,matmap,matgen[pos[i]].pop,col,row,type,nt)!=0){
						printf("ERROR: impossible to compute population specific observed heterozigosity!\n");
						return 0;
					}
				}
				if (hetOss(matgen,matmap,"whole",col,row,type,nt)!=0){
					printf("ERROR: impossible to compute whole population observed heterozigosity!\n");
					return 0;
				}
			}
		}else{printf("Skip observed heterozigosity computation\n");}		
		if (*par2==1){
			if (upop==1){
				if (hetExp(matgen,matmap,"whole",col,row,type,nt)!=0){
					printf("ERROR: impossible to compute whole population expected heterozigosity!\n");
					return 0;
				}
			}
			else if (upop>1){
				for (i=0;i<upop;i++){
					if (hetExp(matgen,matmap,matgen[pos[i]].pop,col,row,type,nt)!=0){
						printf("ERROR: impossible to compute population specific expected heterozigosity!\n");
						return 0;
					}
				}
			}
			if (hetExp(matgen,matmap,"whole",col,row,type,nt)!=0){
				printf("ERROR: impossible to compute whole population expected heterozigosity!\n");
				return 0;
			}
		}else{printf("Skip expected heterozigosity computation\n");}			
	}
	else{perror("ERROR: reading the HET parameters from sumstat.par");}	

	/*AFS (Allele frequency spectrum) COMPUTATION - output all 2D possible combinations*/
	*par1=-1;
	*par2=-1;
	strcpy(type,"NUL");
	if (checkSumstatPar("AFS", par1, par2, type, 3)==0){
		if (*par1==1){
			if (*par2==1){
				for (i=0;i<upop;i++){
					if (afs(matgen, matmap, pos[i], -1, -1, col, row, atoi(type), 1, nt)!=0){
						printf("ERROR: impossible to compute unfolded AFS in population [%s]!\n",matgen[pos[i]].pop);
						return 0;
					}	
				}
			}
			if (*par2==2){
				if (upop<2){printf("WARNING: 2D AFS selected but only %i populations are present in input.\n",upop);}
				for (i=0;i<upop;i++){
					for (k=i+1;k<upop;k++){
						if (afs(matgen, matmap, pos[i], pos[k], -1, col, row, atoi(type), 1, nt)!=0){
							printf("ERROR: impossible to compute unfolded AFS between population [%s] and [%s]!\n",matgen[pos[i]].pop,matgen[pos[k]].pop);
							return 0;
						}
					}
				}
			}
			if (*par2==3){
				if (upop<3){printf("WARNING: 3D AFS selected but only %i populations are present in input.\n",upop);}
				for (i=0;i<upop;i++){
					for (k=i+1;k<upop;k++){
						for (j=k+1;j<upop;j++){
							if (afs(matgen, matmap, pos[i], pos[k], pos[j], col, row, atoi(type), 1, nt)!=0){
								printf("ERROR: impossible to compute unfolded AFS between population [%s], [%s] and [%s]!\n",matgen[pos[i]].pop,matgen[pos[k]].pop,matgen[pos[j]].pop);
								return 0;
							}
						}
					}
				}
			}
		}
		else{printf("Skip AFS Computation\n");}
	}else{perror("ERROR: reading the AFS parameters from sumstat.par");}

	/* POPULATION DISTANCES COMPUTATION - NEI GST73, NEI GST83, HED G'ST05, JOSTD08, W&C84 FST - */
	*par1=-1;
	*par2=-1;
	strcpy(type,"NUL");
	if (checkSumstatPar("DIST", par1, par2, type, 2)==0){
		if (*par1==1){
			outv=malloc(5*sizeof(double));//vector with gstNei1973, gstNei1983, g'stHedrick2005, JostD2008 and FstW&C1984
			if (outv==NULL){perror("ERROR: DIST: memory allocation failure.");}
			for (i=0;i<upop;i++){
				for (k=i+1;k<upop;k++){
					for(j=0;j<5;j++){outv[j]=-5.0;}
					if (nt<2){
						if (dist(matgen, matmap, matgen[pos[i]].pop, matgen[pos[k]].pop, col, row, outv, type)!=0){perror("ERROR: dist: impossible to compute between population distances!");return 0;}}
					else{
						if (dist_mp(matgen, matmap, matgen[pos[i]].pop, matgen[pos[k]].pop, col, row, outv, type, nt)!=0){perror("ERROR: dist_mp: impossible to compute between population distances!");return 0;}
						}
				}
			}
			free(outv);
		}else{printf("Skip population distances computation\n");}	
	}else{perror("ERROR: reading the DIST parameters from sumstat.par");}


	/*INDIVIDUAL SIMILARITY/DISSIMILARITY MATRIX COMPUTATION*/
	*par1=-1;
	*par2=-1;
	strcpy(type,"NUL");
	if (checkSumstatPar("DALL", par1, par2, type, 2)==0){
		if (*par1==1){
			outv=malloc(sizeof(double));
			if (outv==NULL){printf("ERROR: malloc: impossible to allocate outv variable.\n");return 0;}
			/*header printing*/
			if (type[0]=='0'){
				printf("-> Start computing similarity matrix between pairs of individuals using %ld SNPs...\n", col);
				outf=fopen("IND_SIM_MATRIX.txt","w");
				if(outf==NULL ){
        				perror("Error opening similarity output file.");
        				return 0;
				}		
			}
			else if (type[0]=='1'){
				printf("-> Start computing dissimilarity matrix between pairs of individuals using %ld SNPs...\n", col);
				outf=fopen("IND_DIS_MATRIX.txt","w");
				if(outf==NULL ){
        				perror("Error opening dissimilarity output file.");
        				return 0;
				}
			}
			else{printf("ERROR: SIM/DIS function: wrong DALL type spercified in the sumstat.par\n"); return 0;}
			/*print header*/
			for (i=0;i<row;i++){fprintf(outf,"%s\t",matgen[i].name_ind);}
			fprintf(outf,"\n");
			/*header end*/
			k=0;
			for (i=0;i<row;i++){
				//printf("\r- %1.0f%% analysed...",(float)i/row*100);
				////if (i%((int)(row/10))==0){printf("[%ld%%]\n",k*10);k++;}
				printf("\r");
				float ratio = (i+1)/(float)row;
   				int c = ratio * 20;				
				printf("%3d%% [", (int)(ratio*100) );
				for (k=0; k<c; k++){printf("=");}
				for (k=c; k<20; k++){printf(" ");}
				fprintf(outf,"%s\t",matgen[i].name_ind);
				for (j=0;j<=i;j++){
					if (type[0]=='0'){
						if (nt<2){
							if (shAll(matgen, i, j, col, outv, 1)==0){fprintf(outf,"%1.8f\t",*outv);}
							else{printf("ERROR: shAll: error during similarity/dissimilarity index computation.\n");fclose(outf);return 1;}
						}
						else{
							if (shAll_mp(matgen, i, j, col, outv, 1, nt)==0){fprintf(outf,"%1.8f\t",*outv);}
							else{printf("ERROR: shAll: error during similarity/dissimilarity index computation.\n");fclose(outf);return 1;}
						}
					}
					if (type[0]=='1'){
						if (nt<2){
							if (shAll(matgen, i, j, col, outv, 0)==0){fprintf(outf,"%1.8f\t",*outv);}
							else{printf("ERROR: shAll: error during similarity/dissimilarity index computation.\n");fclose(outf);return 1;}
						}
						else{
							if (shAll_mp(matgen, i, j, col, outv, 0, nt)==0){fprintf(outf,"%1.8f\t",*outv);}
							else{printf("ERROR: shAll: error during similarity/dissimilarity index computation.\n");fclose(outf);return 1;}
						}
					}
				}
				fprintf(outf,"\n");
			}
			free(outv);
			fclose(outf);
			printf("][done]\n");
		}else{printf("Skip individual distances computation\n");}
	}else{perror("ERROR: reading the DALL parameters from sumstat.par");}
		
	
	/*memory deallocation*/
	free(par1);
	free(par2);
	free(pos);
	free(n);
	structDeall(matgen, matmap, row);
	/*memory deallocation end*/
	return 0;
}








