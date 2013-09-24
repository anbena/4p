/*  		
		4P v1.0 - Parallel Processing of Polymorphism Panels
    Copyright (C) 2013  Andrea Benazzo, Alex Panziera and Giorgio Bertorelle 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "struct.h"
#include "utils.h"

struct ind_gen;
struct ind_map;

int computeBaseFreq(struct ind_gen *matgen, struct ind_map *matmap, long int row, long int col, char *namepop, char *type, int nt)
/*This function compute base frequencies for each SNP and compute the minor allele */
/* matgen is the structure with genotypes */
/* row is the number or individuals */
/* namepop is the name of the population where perform the computation; whole meand pooled populations */ 
/* v_maf is the vector to export minor base freq */
/*p is the vectors with the unique populations positions */
{
	long int i,j,ntmiss;
	int k,pos,ok=0;
	char iname[100];
	//const char base[]="ACTG";
	char base[5];
	float **tot;/*position corrispondence: A C T G 0 n°ind*/
	tot=(float**) malloc(6*sizeof(float*));
	if (tot==NULL){perror("ERROR: computeBaseFreq: main memory allocation failure.");}
	for (i=0;i<7;i++){
		tot[i]=(float *) malloc(col*sizeof(float));
		if (tot[i]==NULL){perror("ERROR: computeBaseFreq: locus specific memory allocation failure.");}
	}
	char b[col];
	float v_maf[col];
	double tota,totc,tott,totg,totm=0;
	double totsd[5]={0,0,0,0,0};
	float min;
	int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
	if (nt>iCPU){nt=iCPU;}
	omp_set_num_threads(nt);/* Now set the number of threads*/
	base[0]='A';base[1]='C';base[2]='T';base[3]='G';base[4]='\0';	
	FILE *freq_f;
	//array init
	for(i = 0; i < 6; i++) {
		for(j = 0; j < col; j++){
		tot[i][j] = 0;
		if (i==0){b[j]='N';v_maf[j]=0.0;}
		}
	}
	ntmiss=0;
	printf("-> Start computing base frequency in pop [%s] using %u SNPs\t",namepop,(unsigned)strlen(matgen[0].gen1));
		#pragma omp parallel for private(i,j,k,min,pos,ok) reduction(+:tota,totc,totg,tott,totm,ntmiss)
		for (i=0;i<col;i++){
				for (j=0;j<row;j++){
					if (strncmp(namepop,"whole",5)==0){
						if (matgen[j].gen1[i]=='A'){tot[0][i]++;tot[5][i]++;}
						if (matgen[j].gen1[i]=='C'){tot[1][i]++;tot[5][i]++;}
						if (matgen[j].gen1[i]=='T'){tot[2][i]++;tot[5][i]++;}
						if (matgen[j].gen1[i]=='G'){tot[3][i]++;tot[5][i]++;}
						if ((matgen[j].gen1[i]=='0')||(matgen[j].gen1[i]=='N')){tot[4][i]++;}
						if ((matgen[j].gen1[i]!='A')&&(matgen[j].gen1[i]!='C')&&(matgen[j].gen1[i]!='T')&&(matgen[j].gen1[i]!='G')&&(matgen[j].gen1[i]!='0')&&(matgen[j].gen1[i]!='N')){printf("1Corrupted SNP is: %c\n",matgen[j].gen1[i]);perror("ERROR: Base not recognised");exit(1);}
						if (matgen[j].gen2[i]=='A'){tot[0][i]++;tot[5][i]++;}
						if (matgen[j].gen2[i]=='C'){tot[1][i]++;tot[5][i]++;}
						if (matgen[j].gen2[i]=='T'){tot[2][i]++;tot[5][i]++;}
						if (matgen[j].gen2[i]=='G'){tot[3][i]++;tot[5][i]++;}
						if ((matgen[j].gen2[i]=='0')||(matgen[j].gen2[i]=='N')){tot[4][i]++;}
						if ((matgen[j].gen2[i]!='A')&&(matgen[j].gen2[i]!='C')&&(matgen[j].gen2[i]!='T')&&(matgen[j].gen2[i]!='G')&&(matgen[j].gen2[i]!='0')&&(matgen[j].gen2[i]!='N')){printf("2Corrupted SNP is: %c\n",matgen[j].gen2[i]);perror("ERROR: Base not recognised");exit(1);}
					}
					else if ( strncmp(namepop,"whole",strlen(namepop))!=0){
						if (strcmp(matgen[j].pop,namepop)==0){
							if (matgen[j].gen1[i]=='A'){tot[0][i]++;tot[5][i]++;}
							if (matgen[j].gen1[i]=='C'){tot[1][i]++;tot[5][i]++;}
							if (matgen[j].gen1[i]=='T'){tot[2][i]++;tot[5][i]++;}
							if (matgen[j].gen1[i]=='G'){tot[3][i]++;tot[5][i]++;}
							if ((matgen[j].gen1[i]=='0')||(matgen[j].gen1[i]=='N')){tot[4][i]++;}
							if ((matgen[j].gen1[i]!='A')&&(matgen[j].gen1[i]!='C')&&(matgen[j].gen1[i]!='T')&&(matgen[j].gen1[i]!='G')&&(matgen[j].gen1[i]!='0')&&(matgen[j].gen1[i]!='N')){
								printf("3Corrupted SNP is: %c\n",matgen[j].gen1[i]);
								perror("ERROR: Base not recognised");exit(1);
							}
							if (matgen[j].gen2[i]=='A'){tot[0][i]++;tot[5][i]++;}
							if (matgen[j].gen2[i]=='C'){tot[1][i]++;tot[5][i]++;}
							if (matgen[j].gen2[i]=='T'){tot[2][i]++;tot[5][i]++;}
							if (matgen[j].gen2[i]=='G'){tot[3][i]++;tot[5][i]++;}
							if ((matgen[j].gen2[i]=='0')||(matgen[j].gen2[i]=='N')){tot[4][i]++;}
							if ((matgen[j].gen2[i]!='A')&&(matgen[j].gen2[i]!='C')&&(matgen[j].gen2[i]!='T')&&(matgen[j].gen2[i]!='G')&&(matgen[j].gen2[i]!='0')&&(matgen[j].gen2[i]!='N')){
								printf("4Corrupted SNP is: %c\n",matgen[j].gen2[i]);
								perror("ERROR: Base not recognised");exit(1);
							}
						}
					}
					/*find MAF for each SNP*/
					min=2.0;
					pos=-1;
					ok=0;
					for (k=0;k<4;k++){
						if (tot[k][i]>0){ok++;}//the ok variable check for variability, if only one genotypes ok=1, if more ok>1
					}
					for (k=0;k<4;k++){//excluding miss data
						if (  ((tot[k][i]/tot[5][i])<min) && ((tot[k][i]/tot[5][i])!=0) ){min=(tot[k][i]/tot[5][i]);pos=k;}
					}
					if (pos>-1){
						//if(out=='T'){fprintf(maf_f,"%c\t%1.5f\t%1.5f\n",base[pos],min,(tot[4]/(j*2*1.0)));}
						if (ok>1){v_maf[i]=min;b[i]=base[pos];}
						else{v_maf[i]=0.0;b[i]='N';}
					}
					else{
						if ((tot[0][i]+tot[1][i]+tot[2][i]+tot[3][i])==0){v_maf[i]=0.0;}
						else{printf("ERROR: impossible to compute MAF!\n");perror("in function computeBaseFreq, impossible to compute MAF\n");exit(1);}
					}
				}
				if (tot[5][i]!=0){
					tota=tota+((double)tot[0][i]/tot[5][i]);
					totc=totc+((double)tot[1][i]/tot[5][i]);
					tott=tott+((double)tot[2][i]/tot[5][i]);
					totg=totg+((double)tot[3][i]/tot[5][i]);
				}
				else{
					ntmiss++;
				}
				totm=totm+((double)tot[4][i]/(tot[0][i]+tot[1][i]+tot[2][i]+tot[3][i]+tot[4][i]));
		}
	printf("[done]");
	strcpy(iname,"SNP_BASE_FREQUENCY_");
	strcat(iname,namepop);
	strcat(iname,".txt");
	if (type[0]=='1'){//print on a file
		freq_f=fopen(iname,"w");
		printf(" - Writing frequencies to file\t");
		fprintf(freq_f,"POS\tRS\tA\tC\tT\tG\tMISS\tMAFALL\tMAF\n");
		for (i=0;i<col;i++){
			if (tot[5][i]!=0){
				fprintf(freq_f,"%ld\t%s\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%c\t%1.5f\n",matmap[i].bp_pos,matmap[i].rs,tot[0][i]/(tot[5][i]),tot[1][i]/(tot[5][i]),tot[2][i]/(tot[5][i]),tot[3][i]/(tot[5][i]),tot[4][i]/(tot[0][i]+tot[1][i]+tot[2][i]+tot[3][i]+tot[4][i]),b[i],v_maf[i]);
			}
			else{
				fprintf(freq_f,"%ld\t%s\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%c\t%1.5f\n",matmap[i].bp_pos,matmap[i].rs,0.0,0.0,0.0,0.0,tot[4][i]/(tot[0][i]+tot[1][i]+tot[2][i]+tot[3][i]+tot[4][i]),b[i],v_maf[i]);
			}
		}
		fclose(freq_f);
	}
	else if (type[0]=='0'){//print only summary
		freq_f=fopen(iname,"w");
		printf(" - Writing frequencies summary to file\t");
		fprintf(freq_f,"mean(A)\tmean(C)\tmean(T)\tmean(G)\tmean(MISS)\tsd(A)\tsd(C)\tsd(T)\tsd(G)\tsd(MISS)\n");
		//tota=tota/(2*row*col);totc=totc/(2*row*col);tott=tott/(2*row*col);totg=totg/(2*row*col);totm=totm/(2*row*col);
		for (i=0;i<col;i++){
			if (tot[5][i]!=0){
				totsd[0]=totsd[0]+pow(((tot[0][i]/(tot[5][i]))-tota/(col-ntmiss)),2);
				totsd[1]=totsd[1]+pow(((tot[1][i]/(tot[5][i]))-totc/(col-ntmiss)),2);
				totsd[2]=totsd[2]+pow(((tot[2][i]/(tot[5][i]))-tott/(col-ntmiss)),2);
				totsd[3]=totsd[3]+pow(((tot[3][i]/(tot[5][i]))-totg/(col-ntmiss)),2);
			}
			totsd[4]=totsd[4]+pow(((tot[4][i]/(tot[0][i]+tot[1][i]+tot[2][i]+tot[3][i]+tot[4][i]))-totm/col),2);
		}
		totsd[0]=sqrt(totsd[0]/(col-ntmiss-1));
		totsd[1]=sqrt(totsd[1]/(col-ntmiss-1));
		totsd[2]=sqrt(totsd[2]/(col-ntmiss-1));
		totsd[3]=sqrt(totsd[3]/(col-ntmiss-1));
		totsd[4]=sqrt(totsd[4]/(col-1));
		fprintf(freq_f,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",tota/(col-ntmiss),totc/(col-ntmiss),tott/(col-ntmiss),totg/(col-ntmiss),totm/col,totsd[0],totsd[1],totsd[2],totsd[3],totsd[4]);
		fclose(freq_f);
	}
	for (i=0;i<6;i++){free(tot[i]);}
	free(tot);
	printf("[done]\n");
	return 0;
}


int hetOss(struct ind_gen *matgen, struct ind_map *matmap, char *namepop, long int col, long int row, char *type, int nt)
/*This function compute the observed heterozigosity for each snp */
{
long int j,i,c,ntmiss;
long int tot=0;
float v_ho[col];
float totsd,totv=0;
FILE *o_oh;
char iname[100];
int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/
for (j=0;j<col;j++){v_ho[j]=0;}
ntmiss=0;
printf("-> Start computing observed heterozigosity in pop [%s] using %u SNPs\t",namepop,(unsigned)strlen(matgen[0].gen1));
#pragma omp parallel for private(i,j,tot,c) reduction(+:totv,ntmiss) schedule(dynamic)
for (j=0;j<col;j++){
	tot=0;
	c=0;
	for (i=0;i<row;i++){
		if (strncmp(namepop,"whole",5)==0){
			if ((matgen[i].gen1[j]!='N')&&(matgen[i].gen1[j]!='0')&&(matgen[i].gen2[j]!='N')&&(matgen[i].gen2[j]!='0')){
				c++;
				if (matgen[i].gen1[j]!=matgen[i].gen2[j]){tot++;}
			}
		}
		else{
			if (strcmp(matgen[i].pop,namepop)==0){
				if ((matgen[i].gen1[j]!='N')&&(matgen[i].gen1[j]!='0')&&(matgen[i].gen2[j]!='N')&&(matgen[i].gen2[j]!='0')){
					c++;
					if (matgen[i].gen1[j]!=matgen[i].gen2[j]){tot++;}
				}
			}
		}
	}
	if (c!=0){
		v_ho[j]=((float) tot / (float) c);
		totv=totv+v_ho[j];
	}
	else{
		v_ho[j]=-9.0;//missing het
		ntmiss++;
	}
}
printf("[done]");
strcpy(iname,"HET_OBS_");
strcat(iname,namepop);
strcat(iname,".txt");
if (type[0]=='1'){//print on a file
	printf(" - Writing observed het. to file\t");
	o_oh=fopen(iname,"w");
	fprintf(o_oh,"POS\tRS\tH_OBS\n");
	for (j=0;j<col;j++){
		if (v_ho[j]!=-9){fprintf(o_oh,"%ld\t%s\t%1.5f\n",matmap[j].bp_pos,matmap[j].rs,v_ho[j]);}
		else{fprintf(o_oh,"%ld\t%s\tNA\n",matmap[j].bp_pos,matmap[j].rs);}
	}
	fclose(o_oh);
	printf("[done]\n");
}
else if (type[0]=='0'){//print only summary
	o_oh=fopen(iname,"w");
	printf(" - Writing observed het. summary to file\t");
	fprintf(o_oh,"%s_mean(H_OBS)\t%s_sd(H_OBS)\n",namepop,namepop);
	totv=totv/(col-ntmiss);
	totsd=0;
	for (j=0;j<col;j++){
		if (v_ho[j]!=-9){totsd=totsd+(pow((v_ho[j]-totv),2));}
	}
	fprintf(o_oh,"%1.5f\t%1.5f\n",totv,sqrt(totsd/(col-ntmiss-1)));
	fclose(o_oh);
	printf("[done]\n");
}
return 0;
}

int hetExp(struct ind_gen *matgen, struct ind_map *matmap, char *namepop, long int col, long int row, char *type, int nt){
/*This functionn compute the expected heterozigosity for multiallelic markers with the sample size correction of Nei1987*/
long int j,m,n1,ntmiss;
long int p1[4]={0,0,0,0};
float v_he[col];
float totv=0;
float totsd=0;
FILE *o_eh;
char iname[100];
int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/
for (j=0;j<col;j++){v_he[j]=0;}
ntmiss=0;
printf("-> Start computing expected heterozigosity in pop [%s] using %u SNPs\t",namepop,(unsigned)strlen(matgen[0].gen1));
#pragma omp parallel for private(m,j,n1,p1) reduction(+:totv,ntmiss) schedule(dynamic)
for (j=0;j<col;j++){	
	n1=0;
	//initialize allelic frequency vector
	for (m=0;m<4;m++){p1[m]=0.0;}
	/*allelic frequencies computation*/
	for (m=0;m<row;m++) {
		if (strncmp(namepop,"whole",5)==0) {
			if ((matgen[m].gen1[j]!='N')&&(matgen[m].gen1[j]!='0')) {
				n1++;
				if (matgen[m].gen1[j]=='A'){p1[0]++;}
				else if (matgen[m].gen1[j]=='T'){p1[1]++;}
				else if (matgen[m].gen1[j]=='C'){p1[2]++;}
				else if (matgen[m].gen1[j]=='G'){p1[3]++;}
			}
			if ((matgen[m].gen2[j]!='N')&&(matgen[m].gen2[j]!='0')) {
				n1++;
				if (matgen[m].gen2[j]=='A'){p1[0]++;}
				else if (matgen[m].gen2[j]=='T'){p1[1]++;}
				else if (matgen[m].gen2[j]=='C'){p1[2]++;}
				else if (matgen[m].gen2[j]=='G'){p1[3]++;}
			}
		}
		if (strcmp(matgen[m].pop,namepop)==0) {
			if ((matgen[m].gen1[j]!='N')&&(matgen[m].gen1[j]!='0')) {
				n1++;
				if (matgen[m].gen1[j]=='A'){p1[0]++;}
				else if (matgen[m].gen1[j]=='T'){p1[1]++;}
				else if (matgen[m].gen1[j]=='C'){p1[2]++;}
				else if (matgen[m].gen1[j]=='G'){p1[3]++;}
			}
			if ((matgen[m].gen2[j]!='N')&&(matgen[m].gen2[j]!='0')) {
				n1++;
				if (matgen[m].gen2[j]=='A'){p1[0]++;}
				else if (matgen[m].gen2[j]=='T'){p1[1]++;}
				else if (matgen[m].gen2[j]=='C'){p1[2]++;}
				else if (matgen[m].gen2[j]=='G'){p1[3]++;}
			}
		}
	}
	if (n1!=0){
		/*Unbiased expected heterozygosity (Nei 1987)*/	
		v_he[j]=(1.0-(pow(((float)p1[0]/(float)n1),2)+pow(((float)p1[1]/(float)n1),2)+pow(((float)p1[2]/(float)n1),2)+pow(((float)p1[3]/(float)n1),2)))*(n1/((float)n1-1));
		totv=totv+v_he[j];
	}
	else{
		v_he[j]=-9;//missing het
		ntmiss++;
	}
}
printf("[done]");
strcpy(iname,"HET_EXP_");
strcat(iname,namepop);
strcat(iname,".txt");
if (type[0]=='1'){//print on a file
	printf(" - Writing expected het. to file\t");
	o_eh=fopen(iname,"w");
	fprintf(o_eh,"POS\tRS\tH_EXP\n");
	for (j=0;j<col;j++){
		if (v_he[j]!=-9){fprintf(o_eh,"%ld\t%s\t%1.5f\n",matmap[j].bp_pos,matmap[j].rs,v_he[j]);}
		else{fprintf(o_eh,"%ld\t%s\tNA\n",matmap[j].bp_pos,matmap[j].rs);}
	}
	fclose(o_eh);
	printf("[done]\n");
}
else if (type[0]=='0'){//print only summary
	o_eh=fopen(iname,"w");
	printf(" - Writing expected het. summary to file\t");
	fprintf(o_eh,"%s_mean(H_EXP)\t%s_sd(H_EXP)\n",namepop,namepop);
	totv=totv/(col-ntmiss);
	totsd=0;
	for (j=0;j<col;j++){
		if (v_he[j]!=-9){totsd=totsd+(pow((v_he[j]-totv),2));}
	}
	fprintf(o_eh,"%1.5f\t%1.5f\n",totv,sqrt(totsd/(col-ntmiss-1)));
	fclose(o_eh);
	printf("[done]\n");
}
return 0;
}

int afs(struct ind_gen *matgen, struct ind_map *matmap, long int p1, long int p2, long int p3, long int col, long int row, int unfolded, int header, int nt)
/*This function compute the joint allele frequency spectrum of two or three populations, or the simple allele frequency spectrum if only
one population is provided. It computes both the folded or unfolded afs, based on the settings in the sumstat.par file.
p1, p2, p3 are the positions of the first individuals of maximum three populations, set -1 to exclude populations. */
{
long int n1, n2, n3;
long int i,k,j;
long int x,y,z;
long int tot[3];
char iname[100];
FILE *o_afs;
int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/

if ((p1!=-1)&(p2!=-1)&(p3!=-1)){
	if (unfolded==1){printf("-> Start computing unfolded joint AFS (pop1: [%s] - pop2: [%s] - pop3: [%s]) using %ld SNPs\t",matgen[p1].pop,matgen[p2].pop,matgen[p3].pop,col);}
	else{printf("-> Start computing folded joint AFS (pop1: [%s] - pop2: [%s] - pop3: [%s]) using %ld SNPs\t",matgen[p1].pop,matgen[p2].pop,matgen[p3].pop,col);}
}
else if ((p1!=-1)&(p2!=-1)&(p3==-1)){
	if (unfolded==1){printf("-> Start computing unfolded joint AFS (pop1: [%s] - pop2: [%s]) using %ld SNPs\t",matgen[p1].pop,matgen[p2].pop,col);}
	else{printf("-> Start computing folded joint AFS (pop1: [%s] - pop2: [%s]) using %ld SNPs\t",matgen[p1].pop,matgen[p2].pop,col);}
}
else if ((p1!=-1)&(p2==-1)&(p3==-1)){
	if (unfolded==1){printf("-> Start computing unfolded AFS (pop1: [%s]) using %ld SNPs\t",matgen[p1].pop,col);}
	else{printf("-> Start computing folded AFS (pop1: [%s]) using %ld SNPs\t",matgen[p1].pop,col);}
}

if((p1==-1)&(p2==-1)&(p3==-1)){printf("\nERROR: afs: at least one population required!");return 1;}
n1=-1;n2=-1;n3=-1;
if (p1!=-1){n1=0;}
if (p2!=-1){n2=0;}
if (p3!=-1){n3=0;}
//computing number of individuals for each population
for (i=0;i<row;i++){
	if (p1!=-1){if (strcmp(matgen[p1].pop,matgen[i].pop)==0){n1++;}}
	if (p2!=-1){if (strcmp(matgen[p2].pop,matgen[i].pop)==0){n2++;}}
	if (p3!=-1){if (strcmp(matgen[p3].pop,matgen[i].pop)==0){n3++;}}
	}
//set 3D array dimensions equals to the number of individuals plus one
if ((p1!=-1)&(p2!=-1)&(p3!=-1)){x=n1*2+1;y=n2*2+1;z=n3*2+1;}
else if((p1!=-1)&(p2!=-1)&(p3==-1)){x=n1*2+1;y=n2*2+1;z=1;}
else if((p1!=-1)&(p2==-1)&(p3==-1)){x=n1*2+1;y=1;z=1;}
//printf("Parameters: x=%ld, y=%ld, z=%ld, pop1=%s, pop2=%s\n",x,y,z,matgen[p1].pop, matgen[p2].pop);
//  Allocate 3D Array
long int *allElements = malloc(x * y * z * sizeof(long int));
if (allElements==NULL){printf("\nERROR: afs: error allocating memory (1)");return 1;}
long int ***af = malloc(x * sizeof(long int **));
if (af==NULL){printf("\nERROR: afs: error allocating memory (2)");return 1;}
for(i = 0; i < x; i++){
	af[i] = malloc(y * sizeof(long int *));
	if (af[i]==NULL){printf("\nERROR: afs: error allocating memory (3)");return 1;}
	for(j = 0; j < y; j++){
		af[i][j] = allElements + (i * y * z) + (j * z);
	}
}
// initialize 3D array
for(i = 0; i < x; i++){
	for(j = 0; j < y; j++){
		for(k = 0; k < z; k++){
			af[i][j][k]=0;
		}
	}
}
//compute afs
#pragma omp parallel for private(k,tot,i) shared(af) schedule(dynamic)
for (k=0;k<col;k++){
	tot[0]=0;tot[1]=0;tot[2]=0;
	//printf("SNP numero %ld\n", k);
	for (i=0;i<row;i++){
		if (p1!=-1){
			if (strcmp(matgen[p1].pop,matgen[i].pop)==0){
				if (unfolded==0){
					//printf("Pop1 - gen1 - obs:%c - ref:%c - ind:%ld\n",matgen[i].gen1[k], matmap[k].ref,i);
					//printf("Pop1 - gen2 - obs:%c - ref:%c - ind:%ld\n",matgen[i].gen2[k], matmap[k].ref,i);
					if (matgen[i].gen1[k]!=matmap[k].ref){tot[0]++;}
					if (matgen[i].gen2[k]!=matmap[k].ref){tot[0]++;}
				}
				if (unfolded==1){
					if (matgen[i].gen1[k]!=matmap[k].anc){tot[0]++;}
					if (matgen[i].gen2[k]!=matmap[k].anc){tot[0]++;}
				}
			}
		}
		if (p2!=-1){
			if (strcmp(matgen[p2].pop,matgen[i].pop)==0){
				if (unfolded==0){
					//printf("Pop2 - gen1 - obs:%c - ref:%c - ind:%ld\n",matgen[i].gen1[k], matmap[k].ref,i);
					//printf("Pop2 - gen2 - obs:%c - ref:%c - ind:%ld\n",matgen[i].gen2[k], matmap[k].ref,i);
					if (matgen[i].gen1[k]!=matmap[k].ref){tot[1]++;}
					if (matgen[i].gen2[k]!=matmap[k].ref){tot[1]++;}
				}
				if (unfolded==1){
					if (matgen[i].gen1[k]!=matmap[k].anc){tot[1]++;}
					if (matgen[i].gen2[k]!=matmap[k].anc){tot[1]++;}
				}
			}
		}
		if (p3!=-1){
			if (strcmp(matgen[p3].pop,matgen[i].pop)==0){
				if (unfolded==0){
					if (matgen[i].gen1[k]!=matmap[k].ref){tot[2]++;}
					if (matgen[i].gen2[k]!=matmap[k].ref){tot[2]++;}
				}
				if (unfolded==1){
					if (matgen[i].gen1[k]!=matmap[k].anc){tot[2]++;}
					if (matgen[i].gen2[k]!=matmap[k].anc){tot[2]++;}
				}
			}
		}
	}
	//printf("Tot1: %ld, tot2: %ld, tot3: %ld\n", tot[0],tot[1],tot[2]);
	//printf("x:%ld, y:%ld, z:%ld\n",x,y,z);
	//printf("a:%ld, b:%ld, c:%ld\n",a,b,c);
	#pragma omp critical
	if ((p1!=-1)&(p2!=-1)&(p3!=-1)){(af[tot[0]][tot[1]][tot[2]])++;}
	else if ((p1!=-1)&(p2!=-1)&(p3==-1)){(af[tot[0]][tot[1]][0])++;}
	else if ((p1!=-1)&(p2==-1)&(p3==-1)){(af[tot[0]][0][0])++;}
}
printf("[done]");
printf(" - Writing AFS to file\t");
/*printing allele frequency spectrum header*/
if (unfolded==1){strcpy(iname,"AFS-U_");}else{strcpy(iname,"AFS-F_");}
if (p1!=-1){strcat(iname,matgen[p1].pop);}
if (p2!=-1){strcat(iname,"_");strcat(iname,matgen[p2].pop);}
if (p3!=-1){strcat(iname,"_");strcat(iname,matgen[p3].pop);}
strcat(iname,".txt");
o_afs=fopen(iname,"w");
if (header==1){
	if (p1!=-1){
		for (i=0;i<x;i++){
			if (p2!=-1){
				for (k=0;k<y;k++){
					if (p3!=-1){
						for (j=0;j<z;j++){
							if (unfolded==1){fprintf(o_afs,"u_%s%ld_%s%ld_%s%ld\t",matgen[p1].pop,i,matgen[p2].pop,k,matgen[p3].pop,j);}
							else if (unfolded==0){fprintf(o_afs,"f_%s%ld_%s%ld_%s%ld\t",matgen[p1].pop,i,matgen[p2].pop,k,matgen[p3].pop,j);}
						}
					}
					else{
						if (unfolded==1){fprintf(o_afs,"u_%s%ld_%s%ld\t",matgen[p1].pop,i,matgen[p2].pop,k);}
						else if (unfolded==0){fprintf(o_afs,"f_%s%ld_%s%ld\t",matgen[p1].pop,i,matgen[p2].pop,k);}
					}
				}
			}
			else{
				if (unfolded==1){fprintf(o_afs,"u_%s%ld\t",matgen[p1].pop,i);}
				else if (unfolded==0){fprintf(o_afs,"f_%s%ld\t",matgen[p1].pop,i);}
			}
		}
	}
	fprintf(o_afs,"\n");
}
/*printing allele frequency spectrum*/
if ((p1!=-1)&(p2!=-1)&(p3!=-1)){
	for (i=0;i<x;i++){
		for (k=0;k<y;k++){
			for (j=0;j<z;j++){
				fprintf(o_afs,"%ld\t",af[i][k][j]);
			}
		}
	}
	fprintf(o_afs,"\n");
}
else if ((p1!=-1)&(p2!=-1)&(p3==-1)){
	for (i=0;i<x;i++){
		for (k=0;k<y;k++){
			fprintf(o_afs,"%ld\t",af[i][k][0]);
		}
	}
	fprintf(o_afs,"\n");
}
else if ((p1!=-1)&(p2==-1)&(p3==-1)){
	for (i=0;i<x;i++){
			fprintf(o_afs,"%ld\t",af[i][0][0]);
	}
	fprintf(o_afs,"\n");
}
//printf("x:%ld, y:%ld, z:%ld\n",x,y,z);
//printf("Tot1: %ld, tot2: %ld, tot3: %ld\n",n1,n2,n3);
//  Deallocate 3D array
free(allElements);
for(i = 0; i < x; i++){
	free(af[i]);
}
free(af);
fclose(o_afs);
printf("[done]\n");
return 0;
}

int dist(struct ind_gen *matgen, struct ind_map *matmap, char *namepop1, char *namepop2, long int col, long int row, double *outv, char *type)
/*This function compute the Nei Gst (Nei, PNAS 1973), Modified Nei Gst (Nei & Chesser 1983 AnnHumGen), G'st (Hedrick 2005 Evolution),
Jost D (Jost 2008 MolEcol), Fst (Weir and Cockerham 1984) between two populations and export the 5 statistics in the outv vector (the order is the same as explained here).
- matgen is the structure with genotypes
- namepop1 and namepop2 are the name of the population where perform the computation 
- col is the number or SNPs
- row is the number or individuals 
- outv is the results vector */
{
FILE *o_dist;
char iname[100];
long int j,m,miss;
double HT,HS,HT_tot,HS_tot,HE1,HE2,Nar;
double HTest,HSest,HTest_tot,HSest_tot;
double gstH,gstH_tot;
double dj, dj_tot;
double ht1[4],ht2[4],hbar[4],pbar[4],spbar[4],comb[4];
double nc,nbar,h1,h2,a,b,c,a_tot,b_tot,c_tot;
long int n1,n2,p1[4],p2[4];
double dist[5]={0,0,0,0,0};
strcpy(iname,"PAIR_DIST_");
strcat(iname,namepop1);
strcat(iname,"_");
strcat(iname,namepop2);
strcat(iname,".txt");
o_dist=fopen(iname,"w");
HT_tot=0.0;HS_tot=0.0;
HTest_tot=0.0;HSest_tot=0.0;
gstH_tot=0.0;dj_tot=0.0;
a_tot=0.0;b_tot=0.0;c_tot=0.0;
miss=0;
if (type[0]=='1'){//print on a file
	printf("-> Start computing pairwise distance between populations (Pop1: [%s] - Pop2: [%s]) using %ld SNPs\t",namepop1,namepop2,col);
	printf(" - Writing pairwise distances to file\t");
	fprintf(o_dist,"POS\tRS\tGSTNEI73\tGSTNEI83\tGSTHED05\tDJOST\tFSTWC84\n");
}
else if (type[0]=='0'){
	printf("-> Start computing pairwise distance between populations (Pop1: [%s] - Pop2: [%s]) using %ld SNPs\t",namepop1,namepop2,col);
}
for (j=0;j<col;j++){
	n1=0;n2=0;
	Nar=0.0;HT=0.0;HS=0.0;HE1=0.0;HE2=0.0;HTest=0.0;HSest=0.0;
	gstH=0.0;dj=0.0;
	//initialise allelic frequency vectors
	for (m=0;m<4;m++){p1[m]=0.0;p2[m]=0.0;ht1[m]=0.0;ht2[m]=0.0;hbar[m]=0.0;pbar[m]=0.0;spbar[m]=0.0;comb[4];}
	h1=0.0;h2=0.0;
	nc=0.0;nbar=0.0;
	/*allelic frequencies computation*/
	for (m=0;m<5;m++){dist[m]=0;}
	for (m=0;m<row;m++) {
		if (strcmp(matgen[m].pop,namepop1)==0) {
			if ((matgen[m].gen1[j]!='N')&&(matgen[m].gen1[j]!='0')) {
				n1++;
				if (matgen[m].gen1[j]=='A'){p1[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[0]++;}}
				else if (matgen[m].gen1[j]=='T'){p1[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[1]++;}}
				else if (matgen[m].gen1[j]=='C'){p1[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[2]++;}}
				else if (matgen[m].gen1[j]=='G'){p1[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[3]++;}}
			}
			if ((matgen[m].gen2[j]!='N')&&(matgen[m].gen2[j]!='0')) {
				n1++;
				if (matgen[m].gen2[j]=='A'){p1[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[0]++;}}
				else if (matgen[m].gen2[j]=='T'){p1[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[1]++;}}
				else if (matgen[m].gen2[j]=='C'){p1[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[2]++;}}
				else if (matgen[m].gen2[j]=='G'){p1[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[3]++;}}
			}
		}
		if (strcmp(matgen[m].pop,namepop2)==0) {
			if ((matgen[m].gen1[j]!='N')&&(matgen[m].gen1[j]!='0')) {
				n2++;
				if (matgen[m].gen1[j]=='A'){p2[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[0]++;}}
				else if (matgen[m].gen1[j]=='T'){p2[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[1]++;}}
				else if (matgen[m].gen1[j]=='C'){p2[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[2]++;}}
				else if (matgen[m].gen1[j]=='G'){p2[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[3]++;}}
			}
			if ((matgen[m].gen2[j]!='N')&&(matgen[m].gen2[j]!='0')) {
				n2++;
				if (matgen[m].gen2[j]=='A'){p2[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[0]++;}}
				else if (matgen[m].gen2[j]=='T'){p2[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[1]++;}}
				else if (matgen[m].gen2[j]=='C'){p2[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[2]++;}}
				else if (matgen[m].gen2[j]=='G'){p2[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[3]++;}}
			}
		}
	}
	if ((n1==0) || (n2==0)){
		miss++;
		printf("\n  WARNING: DIST: all genotypes are missing in a population at locus rs:\t[%s]\n",matmap[j].rs);
		if (type[0]=='1'){fprintf(o_dist,"%ld\t%s\tNA\tNA\tNA\tNA\tNA\n",matmap[j].bp_pos,matmap[j].rs);}
	}else{
		/*nei gst computation*/
		HE1=1.0-(pow(((double)p1[0]/(double)n1),2)+pow(((double)p1[1]/(double)n1),2)+pow(((double)p1[2]/(double)n1),2)+pow(((double)p1[3]/(double)n1),2));
		HE2=1.0-(pow(((double)p2[0]/(double)n2),2)+pow(((double)p2[1]/(double)n2),2)+pow(((double)p2[2]/(double)n2),2)+pow(((double)p2[3]/(double)n2),2));
		HS=((HE1+HE2)/2.0);
		HT=(1.0-(pow(((double)(p1[0]+p2[0])/(double)(n1+n2)),2)+pow(((double)(p1[1]+p2[1])/(double)(n1+n2)),2)+pow(((double)(p1[2]+p2[2])/(double)(n1+n2)),2)+pow(((double)(p1[3]+p2[3])/(double)(n1+n2)),2)));
		HS_tot=HS_tot+HS;
		HT_tot=HT_tot+HT;
		Nar=2.0/((1.0/((double)n1/2.0))+(1.0/((double)n2/2.0)));
		HSest=HS*((2.0*Nar)/(2.0*Nar-1));
		HTest=HT+(HSest/(4.0*Nar));
		HSest_tot=HSest_tot+HSest;
		HTest_tot=HTest_tot+HTest;
		if((HSest==0)&(HTest==0)){gstH=0.0;}
		else{gstH=(1.0-(HSest/HTest))/((1-HSest)/(1+HSest));}//originale
		gstH_tot=gstH_tot+gstH;
		dj=(((HTest-HSest)/(1-HSest))*2.0);
		if (dj>0){dj_tot=dj_tot+dj;}
		/*W&C84 computation*/
		nc = (((double)n1+(double)n2)/2.0)-((pow(((double)n1),2)/((double)(2.0*(n1+n2))))+(pow(((double)n2),2)/((double)(2.0*(n1+n2)))));
		nbar= ((double)(n1+n2))/4.0;
		a=0.0;b=0.0;c=0.0;
		for (m=0;m<4;m++){
			hbar[m]=(ht1[m]+ht2[m])*2.0/((double)(n1+n2));
			pbar[m]=((double)(p1[m]+p2[m]))/((double)(n1+n2));
			spbar[m]=((pow((((p1[m]/(double)n1))-pbar[m]),2)*n1) + (pow((((p2[m]/(double)n2))-pbar[m]),2)*n2)) / (((double)(n1+n2))/2.0);		
			comb[m]= pbar[m]*(1-pbar[m])-(spbar[m]/2.0);
			a = a + ((nbar/nc)*(spbar[m]-(1/(nbar-1))*(comb[m]-(hbar[m]/4.0))));
			b = b + ((nbar/(nbar-1))*(comb[m]-((2.0*nbar-1)/(4.0*nbar)*hbar[m])));
			c = c + (hbar[m]/2.0);
		}
		a_tot=a_tot+a;
		b_tot=b_tot+b;
		c_tot=c_tot+c;
		if (type[0]=='1'){
			if (HT==0){dist[0]=0.0;}else{dist[0]=1.0-(HS/HT);}//aggiunto
			if (HTest==0){dist[1]=0.0;}else{dist[1]=1-(HSest/HTest);}
			dist[2]=gstH;
			dist[3]=dj;
			if ((a==0)&(b==0)&(c==0)){dist[4]=0.0;}else{dist[4]=a/(a+b+c);}
			fprintf(o_dist,"%ld\t%s\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\n",matmap[j].bp_pos,matmap[j].rs,dist[0],dist[1],dist[2],dist[3],dist[4]);
		}
	}
}
if (type[0]=='0'){//print only summary
	printf("[done] - Writing pairwise distances summary to file\t");
	fprintf(o_dist,"GSTNEI73\tGSTNEI83\tGSTHED05\tDJOST\tFSTWC84\n");
	outv[0] = 1.0-(HS_tot/HT_tot);
	outv[1] = 1.0-(HSest_tot/HTest_tot);
	outv[2] = gstH_tot/(col-miss);
	outv[3] = dj_tot/(col-miss);
	outv[4] = a_tot / (a_tot+b_tot+c_tot);
	printf("\natot: %f - btot: %f - ctot: %f\n",a_tot,b_tot,c_tot);
	fprintf(o_dist,"%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\n",outv[0],outv[1],outv[2],outv[3],outv[4]);
}
fclose(o_dist);
printf("[done]\n");
//printf("\nGst(Nei73): %f - Gst(Nei83): %f - G'st(Hed05): %f - D(Jost08): %f - Fst(WC84): %f", outv[0], outv[1], outv[2], outv[3], outv[4]);
return 0;	
}

int dist_mp(struct ind_gen *matgen, struct ind_map *matmap, char *namepop1, char *namepop2, long int col, long int row, double *outv, char *type, int nt)
/*Parallel (OpenMP) version of the dist function*/
{
FILE *o_dist;
char iname[100];
long int i,j,m,miss;
double HT,HS,HT_tot,HS_tot,HE1,HE2,Nar;
double HTest,HSest,HTest_tot,HSest_tot;
double gstH,gstH_tot;
double dj, dj_tot;
double ht1[4],ht2[4],hbar[4],pbar[4],spbar[4],comb[4];
double nc,nbar,h1,h2,a,b,c,a_tot,b_tot,c_tot;
long int n1,n2,p1[4],p2[4];
/*variables for printing*/
//first five rows contain distanc measures;the sixth an 0/1 variable indicating if the distance can be computed or not
double **locDist;
int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/
HT_tot=0.0;HS_tot=0.0;
HTest_tot=0.0;HSest_tot=0.0;
gstH_tot=0.0;dj_tot=0.0;
a_tot=0.0;b_tot=0.0;c_tot=0.0;
miss=0;
//distances array allocation and init
locDist=(double**) malloc(6*sizeof(double*));
if (locDist==NULL){perror("ERROR: DIST: main memory allocation failure.");}
for (i=0;i<6;i++){
	locDist[i]=(double *) malloc(col*sizeof(double));
	if (locDist[i]==NULL){perror("ERROR: DIST: locus specific memory allocation failure.");}
}
for (j=0;j<col;j++){
	for (i=0;i<6;i++){locDist[i][j]=0.0;}
}
printf("-> Start computing pairwise distance between populations (Pop1: [%s] - Pop2: [%s]) using %ld SNPs\t",namepop1,namepop2,col);
#pragma omp parallel for private(j,m,n1,n2,HT,HS,HE1,HE2,Nar,HTest,HSest,gstH,dj,ht1,ht2,hbar,pbar,spbar,comb,nc,nbar,h1,h2,a,b,c,p1,p2) reduction(+:HS_tot,HT_tot,HSest_tot,HTest_tot,gstH_tot,dj_tot,a_tot,b_tot,c_tot,miss) schedule(dynamic)
for (j=0;j<col;j++){	
	n1=0;n2=0;
	Nar=0.0;HT=0.0;HS=0.0;HE1=0.0;HE2=0.0;HTest=0.0;HSest=0.0;
	gstH=0.0;dj=0.0;
	//initialise allelic frequency vectors
	for (m=0;m<4;m++){p1[m]=0.0;p2[m]=0.0;ht1[m]=0.0;ht2[m]=0.0;hbar[m]=0.0;pbar[m]=0.0;spbar[m]=0.0;comb[4];}
	h1=0.0;h2=0.0;
	nc=0.0;nbar=0.0;
	/*allelic frequencies computation*/
	for (m=0;m<row;m++) {
		if (strcmp(matgen[m].pop,namepop1)==0) {
			if ((matgen[m].gen1[j]!='N')&&(matgen[m].gen1[j]!='0')) {
				n1++;
				if (matgen[m].gen1[j]=='A'){p1[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[0]++;}}
				else if (matgen[m].gen1[j]=='T'){p1[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[1]++;}}
				else if (matgen[m].gen1[j]=='C'){p1[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[2]++;}}
				else if (matgen[m].gen1[j]=='G'){p1[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[3]++;}}
			}
			if ((matgen[m].gen2[j]!='N')&&(matgen[m].gen2[j]!='0')) {
				n1++;
				if (matgen[m].gen2[j]=='A'){p1[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[0]++;}}
				else if (matgen[m].gen2[j]=='T'){p1[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[1]++;}}
				else if (matgen[m].gen2[j]=='C'){p1[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[2]++;}}
				else if (matgen[m].gen2[j]=='G'){p1[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht1[3]++;}}
			}
		}
		if (strcmp(matgen[m].pop,namepop2)==0) {
			if ((matgen[m].gen1[j]!='N')&&(matgen[m].gen1[j]!='0')) {
				n2++;
				if (matgen[m].gen1[j]=='A'){p2[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[0]++;}}
				else if (matgen[m].gen1[j]=='T'){p2[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[1]++;}}
				else if (matgen[m].gen1[j]=='C'){p2[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[2]++;}}
				else if (matgen[m].gen1[j]=='G'){p2[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[3]++;}}
			}
			if ((matgen[m].gen2[j]!='N')&&(matgen[m].gen2[j]!='0')) {
				n2++;
				if (matgen[m].gen2[j]=='A'){p2[0]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[0]++;}}
				else if (matgen[m].gen2[j]=='T'){p2[1]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[1]++;}}
				else if (matgen[m].gen2[j]=='C'){p2[2]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[2]++;}}
				else if (matgen[m].gen2[j]=='G'){p2[3]++;if (matgen[m].gen1[j]!=matgen[m].gen2[j]){ht2[3]++;}}
			}
		}
	}
	if ((n1==0) || (n2==0)){
		miss++;
		printf("\n  -> WARNING: DIST: all genotypes are missing in a population at locus rs:\t[%s]",matmap[j].rs);
		locDist[5][j]=1;
	}else{
		/*nei gst computation*/
		HE1=1.0-(pow(((double)p1[0]/(double)n1),2)+pow(((double)p1[1]/(double)n1),2)+pow(((double)p1[2]/(double)n1),2)+pow(((double)p1[3]/(double)n1),2));
		HE2=1.0-(pow(((double)p2[0]/(double)n2),2)+pow(((double)p2[1]/(double)n2),2)+pow(((double)p2[2]/(double)n2),2)+pow(((double)p2[3]/(double)n2),2));
		HS=((HE1+HE2)/2.0);
		HT=(1.0-(pow(((double)(p1[0]+p2[0])/(double)(n1+n2)),2)+pow(((double)(p1[1]+p2[1])/(double)(n1+n2)),2)+pow(((double)(p1[2]+p2[2])/(double)(n1+n2)),2)+pow(((double)(p1[3]+p2[3])/(double)(n1+n2)),2)));
		HS_tot=HS_tot+HS;
		HT_tot=HT_tot+HT;
		Nar=2.0/((1.0/((double)n1/2.0))+(1.0/((double)n2/2.0)));
		HSest=HS*((2.0*Nar)/(2.0*Nar-1));
		HTest=HT+(HSest/(4.0*Nar));
		HSest_tot=HSest_tot+HSest;
		HTest_tot=HTest_tot+HTest;
		if((HSest==0)&(HTest==0)){gstH=0.0;}
		else{gstH=(1.0-(HSest/HTest))/((1-HSest)/(1+HSest));}
		gstH_tot=gstH_tot+gstH;
		dj=(((HTest-HSest)/(1-HSest))*2.0);
		if (dj>0){dj_tot=dj_tot+dj;}
		/*W&C84 computation*/
		nc = (((double)n1+(double)n2)/2.0)-((pow(((double)n1),2)/((double)(2.0*(n1+n2))))+(pow(((double)n2),2)/((double)(2.0*(n1+n2)))));
		nbar= ((double)(n1+n2))/4.0;
		a=0.0;b=0.0;c=0.0;
		for (m=0;m<4;m++){
			hbar[m]=(ht1[m]+ht2[m])*2.0/((double)(n1+n2));
			pbar[m]=((double)(p1[m]+p2[m]))/((double)(n1+n2));
			spbar[m]=((pow((((p1[m]/(double)n1))-pbar[m]),2)*n1) + (pow((((p2[m]/(double)n2))-pbar[m]),2)*n2)) / (((double)(n1+n2))/2.0);		
			comb[m]= pbar[m]*(1-pbar[m])-(spbar[m]/2.0);
			a = a + ((nbar/nc)*(spbar[m]-(1/(nbar-1))*(comb[m]-(hbar[m]/4.0))));
			b = b + ((nbar/(nbar-1))*(comb[m]-((2.0*nbar-1)/(4.0*nbar)*hbar[m])));
			c = c + (hbar[m]/2.0);
		}
		a_tot=a_tot+a;
		b_tot=b_tot+b;
		c_tot=c_tot+c;
		/*save locus specific distances*/
		if (HT==0){locDist[0][j]=0.0;}else{locDist[0][j]=1.0-(HS/HT);}
		if (HTest==0){locDist[1][j]=0.0;}else{locDist[1][j]=1.0-(HSest/HTest);}
		locDist[2][j]=gstH;
		locDist[3][j]=dj;
		if ((a==0)&(b==0)&(c==0)){locDist[4][j]=0.0;}else{locDist[4][j]=a/(a+b+c);};
	}
}
printf("[done]");
strcpy(iname,"PAIR_DIST_");
strcat(iname,namepop1);
strcat(iname,"_");
strcat(iname,namepop2);
strcat(iname,".txt");
if (type[0]=='1'){//print on a file
	o_dist=fopen(iname,"w");
	printf(" - Writing pairwise distances to file\t");
	fprintf(o_dist,"POS\tRS\tGSTNEI73\tGSTNEI83\tGSTHED05\tDJOST\tFSTWC84\n");
	for (i=0;i<col;i++){
		if (locDist[5][i]==0){
			fprintf(o_dist,"%ld\t%s\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\n",matmap[i].bp_pos,matmap[i].rs,locDist[0][i],locDist[1][i],locDist[2][i],locDist[3][i],locDist[4][i]);
		}else{
			fprintf(o_dist,"%ld\t%s\tNA\tNA\tNA\tNA\tNA\n",matmap[i].bp_pos,matmap[i].rs);
		}
	}
	fclose(o_dist);
}
else if (type[0]=='0'){//print only summary
	o_dist=fopen(iname,"w");
	printf(" - Writing pairwise distances summary to file\t");
	fprintf(o_dist,"GSTNEI73\tGSTNEI83\tGSTHED05\tDJOST\tFSTWC84\n");
	outv[0] = 1.0-(HS_tot/HT_tot);//gst73
	outv[1] = 1.0-(HSest_tot/HTest_tot);//gst83
	outv[2] = gstH_tot/(col-miss);
	outv[3] = dj_tot/(col-miss);//jostD
	outv[4] = a_tot / (a_tot+b_tot+c_tot);//fstW&C
	fprintf(o_dist,"%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\n",outv[0],outv[1],outv[2],outv[3],outv[4]);
	fclose(o_dist);
}
for (i=0;i<6;i++){free(locDist[i]);}
free(locDist);
//printf("\nGst(Nei73): %f - Gst(Nei83): %f - G'st(Hed05): %f - D(Jost08): %f - Fst(WC84): %f", outv[0], outv[1], outv[2], outv[3], outv[4]);
printf("[done]\n");
return 0;
}

int shAll(struct ind_gen *matgen, long int i, long int j, long int col, double *outv, int similarity)
/*This function compute the fraction of allelic states shared/not-shared between two individuals. This ratio is computed excluding missing positions */
{
long int nconf, k, tot;
*outv=0.0;
nconf=0;
tot=0;
for (k=0;k<col;k++){
		if (matgen[i].gen1[k]==matgen[j].gen1[k]){//1*if
			if ((matgen[i].gen1[k]!='N')&&(matgen[i].gen1[k]!='0')){//2*if
				tot++;
				nconf++;
				if (matgen[i].gen2[k]==matgen[j].gen2[k]){//3*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')){//4*if
						tot++;
						nconf++;
					}//end4*if
				}//end3*if
				else{//else3*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')&&(matgen[j].gen2[k]!='N')&&(matgen[j].gen2[k]!='0')){nconf++;}
				}//end else3*if
			}//end 2*if
			else{//else2*if
				if (matgen[i].gen2[k]==matgen[j].gen2[k]){//5*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')){//6*if
						tot++;
						nconf++;
					}//end6*if
				}//end 5* if
				else{//else5*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')&&(matgen[j].gen2[k]!='N')&&(matgen[j].gen2[k]!='0')){nconf++;}
				}//end else5*if
			}//end else2*if
		}//end1*if
		else{//else1*if
			if ((matgen[i].gen1[k]!='N')&&(matgen[i].gen1[k]!='0')&&(matgen[j].gen1[k]!='N')&&(matgen[j].gen1[k]!='0')){//7*if
				nconf++;
				if (matgen[i].gen2[k]==matgen[j].gen2[k]){//8*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')){//9*if
						tot++;
						nconf++;
					}//end9*if
				}//end8if
				else{//else8*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')&&(matgen[j].gen2[k]!='N')&&(matgen[j].gen2[k]!='0')){nconf++;}
				}//end else 8*if
			}//end7*if
		}//end else 1*if
	}
if (similarity==1){*outv=(double)tot/(double)nconf;}
else{*outv= 1.0-((double)tot/(double)nconf);}
//printf("\ntot: %ld - nconf: %ld\n",tot,nconf);
return 0;
}

int shAll_mp(struct ind_gen *matgen, long int i, long int j, long int col, double *outv, int similarity, int nt)
/*This function compute the fraction of allelic states shared/not-shared between two individuals. This ratio is computed excluding missing positions. OpenMP accellerated*/
{
long int nconf, k, tot;
int chunk_size;
int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/
*outv=0.0;
nconf=0;
tot=0;
chunk_size=ceil(col/nt);
#pragma omp parallel for private(k) reduction(+:nconf,tot) schedule(static,chunk_size)
for (k=0;k<col;k++){
		if (matgen[i].gen1[k]==matgen[j].gen1[k]){//1*if
			if ((matgen[i].gen1[k]!='N')&&(matgen[i].gen1[k]!='0')){//2*if
				tot++;
				nconf++;
				if (matgen[i].gen2[k]==matgen[j].gen2[k]){//3*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')){//4*if
						tot++;
						nconf++;
					}//end4*if
				}//end3*if
				else{//else3*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')&&(matgen[j].gen2[k]!='N')&&(matgen[j].gen2[k]!='0')){nconf++;}
				}//end else3*if
			}//end 2*if
			else{//else2*if
				if (matgen[i].gen2[k]==matgen[j].gen2[k]){//5*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')){//6*if
						tot++;
						nconf++;
					}//end6*if
				}//end 5* if
				else{//else5*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')&&(matgen[j].gen2[k]!='N')&&(matgen[j].gen2[k]!='0')){nconf++;}
				}//end else5*if
			}//end else2*if
		}//end1*if
		else{//else1*if
			if ((matgen[i].gen1[k]!='N')&&(matgen[i].gen1[k]!='0')&&(matgen[j].gen1[k]!='N')&&(matgen[j].gen1[k]!='0')){//7*if
				nconf++;
				if (matgen[i].gen2[k]==matgen[j].gen2[k]){//8*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')){//9*if
						tot++;
						nconf++;
					}//end9*if
				}//end8if
				else{//else8*if
					if ((matgen[i].gen2[k]!='N')&&(matgen[i].gen2[k]!='0')&&(matgen[j].gen2[k]!='N')&&(matgen[j].gen2[k]!='0')){nconf++;}
				}//end else 8*if
			}//end7*if
		}//end else 1*if
	}
if (similarity==1){*outv=((double)tot/(double)nconf);}
else{*outv= 1.0-((double)tot/(double)nconf);}
return 0;
}


