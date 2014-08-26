/*  		
		4P v1.0 - Parallel Processing of Polymorphism Panels
    Copyright (C) 2013  Andrea Benazzo, Alex Panziera and Giorgio Bertorelle 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "struct.h"

int loadPedSnp(struct ind_gen *matgen, long int col, long int row, FILE *i_f)
{
/* this function load plink ped input files (*.ped) and store the data into
ind_gen objects. matgen is the storing object, col is the number of SNPs,
row is the number of individuals and i_f is the pointer to the input file.
Note: don't use plink binary format*/
/*FARE UNA FUNZIONE CONVERTER PER ELIMINARE IL DOPPIO SPAZIO NEL FILE DI INPUT--*/
char *buf, *buf_split;
long int i,j,gc1,gc2;
int ok;
float perc;
buf=malloc(((col+1)*5+100)*sizeof(char));
for (i=0;i<row;i++){
	if (fgets(buf,((col+1)*5+100),i_f) !=NULL ){//leggo tutta una riga alla volta
	choppy(buf);
	perc=(((i+1.0)/row)*100);
	printf("\rReading genotypes from input file     (%1.0f%%)",perc);
	//printf("Row read: %s\n",buf);//debug
	//printf("Row length %u!\n",(unsigned)strlen(buf));
	ok=0;
	j=0;
	gc1=0;gc2=0;//counter for filling genotypes
	buf_split = strtok (buf," \t");
	while(buf_split != NULL){
		if (j==0){strcpy(matgen[i].pop,buf_split);}
		if (j==1){strcpy(matgen[i].name_ind,buf_split);}
		if (j==2){strcpy(matgen[i].paternal_id,buf_split);}
		if (j==3){strcpy(matgen[i].maternal_id,buf_split);}
		if (j==4){matgen[i].sex=atoi(buf_split);}
		if (j==5){matgen[i].phen=atoi(buf_split);}
		if (j>5){
			if (ok==0){
				strncpy(&matgen[i].gen1[gc1],buf_split,1);
				/*matgen[i].gen1[gc1]=(char)buf_split[0];//puts(buf_split);//puts(&matgen[i].gen1[gc1]);//printf("COL %d\n",j);strcpy(matgen[i].gen1[j],buf_split[0]);*/
				ok=1;
				gc1++;
				}
			else	{
				strncpy(&matgen[i].gen2[gc2],buf_split,1);
				/*matgen[i].gen2[gc2]=(char)buf_split[0];//puts(buf_split);//puts(&matgen[i].gen2[gc2]);//printf("COL %d\n",j);strcpy(matgen[i].gen2[j],buf_split[0]);*/
				ok=0;
				gc2++;
				}
			}
		buf_split = strtok(NULL," \t");
		j++;
		}
		matgen[i].gen1[gc1]='\0';
		matgen[i].gen2[gc2]='\0';
		if (strlen(matgen[i].gen1)!=strlen(matgen[i].gen2)){printf("ERROR: different number of snps in the two cromosomes! (check input file row %ld)\n",i+1);return 1;}
		if (((j-6)/2)!=col){printf("ERROR: loadPedSnp: a different SNPs number is present in the input file [-s: %ld / ped: %ld]\n",col,((j-6)/2));return 1;}
}
else{
	printf("There are only %ld individuals(rows) in the input file!\n",i+1);
	perror("Error attempting to read input file rows");
	return 1;
	}
}
printf("\t[%ld individuals loaded]\n",i);
if (fgets(buf,((col+1)*5+100),i_f) !=NULL ){printf("WARNING: loadPedSnp: In the input file there are more individuals than specified in the command line -n argument.\n");}
//printf("SNP vectors length: gen1:%ld, gen2:%ld\n",gc1,gc2);
free(buf);
return 0;
}

int loadMapSnp(struct ind_map *matmap, long int col, FILE *i_m)
{
/* this function load plink map input files (*.map) and store the data into
an ind_map object. matmap is the storing object, col is the number of SNPs
and i_m is the pointer to the input file.*/
long int i;
float perc;
char *bufmap, *bufmap_split;

bufmap=malloc(100*sizeof(char));
for (i=0;i<col;i++){
	perc=(((i+1.0)/col)*100);
	printf("\rReading snp info from input file      (%1.0f%%)",perc);
	if ( fgets(bufmap,100,i_m) !=NULL ){
		choppy(bufmap);
		//printf("la riga del mapfile è: %s\n",bufmap);
		bufmap_split=strtok(bufmap," \t");
		matmap[i].crom=(short unsigned)atoi(bufmap_split);
		bufmap_split=strtok(NULL," \t");
		strncpy(matmap[i].rs,bufmap_split,sizeof(matmap[i].rs));
		bufmap_split=strtok(NULL," \t");
		matmap[i].gen_dist=(short)atoi(bufmap_split);
		bufmap_split=strtok(NULL," \t");
		matmap[i].bp_pos=(long)atoi(bufmap_split);
		if (strtok(NULL," \t")!=NULL){printf("ERROR: when reading mapfile entry row nr %ld!\n",i);perror("more than 4 entries in mapfile");return 0;}
	}
	else{
		printf("ERROR:There are only %ld SNPs in the input mapfile!\n",i);
		perror("Error attempting to read input mapfile rows");
		return 0;
	}
}
printf("\t[info loaded for %ld SNPs]\n",i);
free(bufmap);
return 0;
}


int loadArlSnp(struct ind_gen *matgen, long int col, long int row, FILE *i_f)
/* this function load arlequin input files (*.arl) and store the data into
ind_gen objects. matgen is the storing object, col is the number of SNPs
and i_f is the pointer to the input file. 
NB: no space between snp permitted
NB: simulations with fastsimcoal neds -k and -g parameter */
{
char *buf;
float perc;
char pop[20]="0";
long int nind=0;
long int k,i=0;
int okpop=0;
int oknind=0;
int oksd=0;
int oksdc=0;
int r=0;
char *buf_split;
long int len;

buf=malloc(((col+1)+500)*sizeof(char));
while (fgets(buf,((col+1)+500),i_f) !=NULL){
	if (strcmp(buf,"\n")==0){continue;}//controllo se la riga è vuota
	choppy(buf);
	if (strpbrk(buf,"#,")!=NULL){continue;}
	if (i==0){printf("\rReading genotypes from arlequin input file     (0%%)");}
	else{	
		perc=((float)i/row)*100.0;
		printf("\rReading genotypes from arlequin input file     (%1.0f%%)",perc);
		}
	//printf("Row read: %s\n",buf);//debug
	buf_split = strtok(buf,"=");
	while (buf_split != NULL){
		if (strcmp(buf_split,"\t\tSampleName")==0){
			buf_split = strtok (NULL,"=");
			strcpy(pop,buf_split);
			okpop=1;
			//printf("Pop read: %s\n",pop);
		}
		if (strcmp(buf_split,"\t\tSampleSize")==0){
			buf_split = strtok (NULL,"=");
			nind=atol(buf_split);
			oknind=1;
			//printf("Number of individuals: %ld\n",nind);
		}
		if (strcmp(buf_split,"\t\tSampleData")==0){
			oksd=1;
			oksdc=1;
		}
		buf_split = strtok (NULL,"=");
	}
	if (oksdc==1){oksdc=0;continue;}
	/*if we get pop, number of individuals and skip "SampleData" line then...*/
	if (okpop==1&oknind==1&nind>0&oksd==1){
		if (i>(row-1)){perror("ERROR:Reading more individuals than these specified in the commandline!");return 1;}
		if (r==0){
			buf_split = strtok (buf,"\t");//get the nomeind
			strcpy(matgen[i].name_ind,buf_split);
			strcpy(matgen[i].pop,pop);
			buf_split = strtok (NULL,"\t");//get the frequency - unuseful
			buf_split = strtok (NULL,"\t ");//get all genotypes
			len=strlen(buf_split);
			for (k=0;k<len;k++){
				if (strncmp(&buf_split[k],"0",1)){strncpy(&buf_split[k],"A",1);}
				else if (strncmp(&buf_split[k],"1",1)){strncpy(&buf_split[k],"G",1);}
				else if (strncmp(&buf_split[k],"?",1)){strncpy(&buf_split[k],"N",1);}
				}
			strcpy(matgen[i].gen1,buf_split);
			len=strlen(matgen[i].gen1);
			r=1;
			buf_split = strtok (NULL,"\t");//terminate reading
			}
		else if (r==1){
			buf_split = strtok (buf,"\t ");//delete tab
			len=strlen(buf_split);
			for (k=0;k<len;k++){
				if (strncmp(&buf_split[k],"0",1)){strncpy(&buf_split[k],"A",1);}
				else if (strncmp(&buf_split[k],"1",1)){strncpy(&buf_split[k],"G",1);}
				else if (strncmp(&buf_split[k],"?",1)){strncpy(&buf_split[k],"N",1);}
				}
			strcpy(matgen[i].gen2,buf_split);
			r=0;
			nind--;
			matgen[i].gen1[col]='\0';
			matgen[i].gen2[col]='\0';
			//printf("Individual %ld genotype1 is: %s\n",i,matgen[i].gen1);
			//printf("Individual %ld genotype2 is: %s\n",i,matgen[i].gen2);
			i++;
			}
	}
	if(nind==0&oknind==1){okpop=0;oknind=0;oksd=0;}
}
if (i<row){perror("\nERROR:Reading less individuals than these specified in the commandline!\n");return 1;}
printf("\n");
free(buf);
return 0;
}

int loadVcfSnp(struct ind_gen *matgen, struct ind_map *matmap, long int col, long int row, FILE *i_f)
/* this function load not compressed vcf files (*.vcf) according to the 4.1 file format
from 1000 genomes project (http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)
and store the data into ind_gen objects. matgen is the storing object, col is the number of SNPs, row is the number of individuals
and i_f is the pointer to the input file.
The function strictly follow the 4.1 vcf format thus it does not store non recognised fields. 
It stores only the first 100 characters of the info field. */
{
char *buf;
char *buf_split, *bsf;
char *format, *spf;
char *g, *spg;
char *l, *spl;
long int j,k,i,ii;
int pos[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};//vector with the format fields order. 4.1 vcf standard has 13 fields. 
float perc;

k=0;
buf=malloc(((row*30*2)+500)*sizeof(char));
printf("Start reading vcf file...");
while (fgets(buf,((row*30*2)+500),i_f) !=NULL){
	choppy(buf);
	//printf("\n%s\n",buf);///////////////////////////
	/*skip first lines with non used informations*/
	if ((buf[0]=='#')&&(buf[1]=='#')){continue;}
	else{
		buf_split=strtok_r(buf,"\t",&bsf);
		if (strncmp(buf_split,"#CHROM",6)==0){
			i=0;
			/*skip first 9 elements of the header*/
			for (j=0;j<9;j++){buf_split=strtok_r(NULL,"\t",&bsf);}
			/*reading individuals names*/
			while (buf_split != NULL){//printf("\nindivid %s\n",buf_split);///////////////////////////
				i++;
				if (i>row){perror("ERROR: vcfread: attempting to read more inidividuals than specified in the command line!\n");return 1;}
				strncpy(matgen[(i-1)].name_ind,buf_split,sizeof(matgen[i].name_ind)/sizeof(char));
				matgen[(i-1)].name_ind[(strlen(matgen[(i-1)].name_ind)+1)]='\0';
				//printf("\nsource %s and sink %s , composed by %d characters.\n",buf_split,matgen[(i-1)].name_ind,(int)strlen(matgen[(i-1)].name_ind));
				buf_split=strtok_r(NULL,"\t",&bsf);
			}
			printf("\nChecking individuals name: %ld individuals present in vcf file. (%s , %s , %s, ...)\n",i,matgen[0].name_ind,matgen[1].name_ind,matgen[2].name_ind);
			if (i<row){printf("ERROR: vcfread: too few individuals in the vcf file (Nr Individuals - present: %ld, expected: %ld).\n",i,row);return 1;}
		}
		else{
			/*start reading snps*/
			if (k>(col-1)){perror("ERROR: vcfread: attempting to read more markers than specified in the command line!\n");return 1;}
			if (k==0){printf("\rReading genotypes from vcf input file     (0%%)");}
			else{	
				perc=((float)k/col)*100.0;
				printf("\rReading genotypes from vcf input file     (%1.0f%%)",perc);
				//printf("\nsnp: %ld\n",k);
				}
			i=0;
			while (buf_split != NULL){
				i++;
				if (i==1){matmap[k].crom=(short unsigned)atoi(buf_split);/*printf("\n%d",matmap[k].crom);*/}
				else if (i==2){matmap[k].bp_pos=atol(buf_split);/*printf("\n%ld",matmap[k].bp_pos);*/}
				else if (i==3){strcpy(matmap[k].rs,buf_split);/*printf("\n%s",matmap[k].rs);*/}
				else if (i==4){matmap[k].ref=*buf_split;/*strncpy(&matmap[k].ref,buf_split,1);printf("\n%c",matmap[k].ref);*/}
				else if (i==5){strcpy(matmap[k].alt,buf_split);/*printf("\n%s",matmap[k].alt);*/}
				else if (i==6){matmap[k].qual = atoi(buf_split);/*printf("\n%d",matmap[k].qual);*/}
				else if (i==7){strcpy(matmap[k].filt,buf_split);/*printf("\n%s",matmap[k].filt);*/}
				else if (i==8){
						if (strlen(buf_split)>99){
							strncpy(matmap[k].info,buf_split,99);
							matmap[k].info[99]='\0';
							//printf("\n%s",matmap[k].info);
							}
						else{
							strcpy(matmap[k].info,buf_split);
							//printf("\n%s",matmap[k].info);
							}
						}
				else if (i==9){
					/*detecting which fields (and the order) are present in the format field*/
					format=strtok_r(buf_split,":",&spf);
					for (ii=0;ii<13;ii++){pos[ii]=0;}
					ii=0;
					while (format!=NULL){
						ii++;
						if (strcmp(format,"GT")==0){pos[0]=ii;}
						else if (strcmp(format,"DP")==0){pos[1]=ii;}
						else if (strcmp(format,"FT")==0){pos[2]=ii;}
						else if (strcmp(format,"GL")==0){pos[3]=ii;}
						else if (strcmp(format,"GLE")==0){pos[4]=ii;}
						else if (strcmp(format,"PL")==0){pos[5]=ii;}
						else if (strcmp(format,"GP")==0){pos[6]=ii;}
						else if (strcmp(format,"GQ")==0){pos[7]=ii;}
						else if (strcmp(format,"HQ")==0){pos[8]=ii;}
						else if (strcmp(format,"PS")==0){pos[9]=ii;}
						else if (strcmp(format,"PQ")==0){pos[10]=ii;}
						else if (strcmp(format,"EC")==0){pos[11]=ii;}
						else if (strcmp(format,"MQ")==0){pos[12]=ii;}
						format=strtok_r(NULL,":",&spf);
					}
				}
				else if (i>9){
					g=strtok_r(buf_split,":",&spg);
					ii=0;
					while(g!=NULL){
						ii++;
						if (ii==pos[0]){
							if (strcmp(g,".")==0){/*The case of missing genotype indicated as "." (it is not in the vcf 4.1 standard but it happen) */
								matgen[(i-10)].gen1[k] = 'N';
								matgen[(i-10)].gen2[k] = 'N';
							}
							else{
								if (strchr(g,'|')!=NULL){//the phase is known
									if (i==10){matgen[(i-10)].phase=1;}
									l=strtok_r(g,"|",&spl);
									if (strcmp(l,".")==0){matgen[(i-10)].gen1[k] = 'N';}
									else{
										matgen[(i-10)].gen1[k] = getAll(matmap[k].ref,matmap[k].alt,atoi(l));/*matgen[(i-10)].gen1[k] = *l;*/
										if (matgen[(i-10)].gen1[k]=='n'){printf("ERROR: getAll: impossible to extract the right allele from ref/alt! (Check ind %ld marker %ld)\n",(i-9),(k+1));return 1;}
									}
									l=strtok_r(NULL,"|",&spl);
									if (strcmp(l,".")==0){matgen[(i-10)].gen2[k] = 'N';}
									else{
										matgen[(i-10)].gen2[k] = getAll(matmap[k].ref,matmap[k].alt,atoi(l));/*matgen[(i-10)].gen2[k] = *l;*/
										if (matgen[(i-10)].gen2[k]=='n'){printf("ERROR: getAll: impossible to extract the right allele from ref/alt! (Check ind %ld marker %ld)\n",(i-9),(k+1));return 1;}
									}
								}
								else{//the phase is unknown
									if (i==10){matgen[(i-10)].phase=0;}
									l=strtok_r(g,"/",&spl);
									if (strcmp(l,".")==0){matgen[(i-10)].gen1[k] = 'N';}
									else{
										matgen[(i-10)].gen1[k] = getAll(matmap[k].ref,matmap[k].alt,atoi(l));/*matgen[(i-10)].gen1[k] = *l;*/
										if (matgen[(i-10)].gen1[k]=='n'){printf("ERROR: getAll: impossible to extract the right allele from ref/alt! (Check ind %ld marker %ld)\n",(i-9),(k+1));return 1;}
									}
									l=strtok_r(NULL,"/",&spl);
									if (strcmp(l,".")==0){matgen[(i-10)].gen2[k] = 'N';}
									else{
										matgen[(i-10)].gen2[k] = getAll(matmap[k].ref,matmap[k].alt,atoi(l));/*matgen[(i-10)].gen2[k] = *l;*/
										if (matgen[(i-10)].gen2[k]=='n'){printf("ERROR: getAll: impossible to extract the right allele from ref/alt! (Check ind %ld marker %ld)\n",(i-9),(k+1));return 1;}
									}
								}
							}
						}//end pos1
						if (ii==pos[1]){matgen[(i-10)].dp=atoi(g);}
						if (ii==pos[2]){strcpy(matgen[(i-10)].ft,g);}
						if (ii==pos[3]){strcpy(matgen[(i-10)].gl,g);}
						if (ii==pos[4]){strcpy(matgen[(i-10)].gle,g);}
						if (ii==pos[5]){strcpy(matgen[(i-10)].pl,g);}
						if (ii==pos[6]){strcpy(matgen[(i-10)].gp,g);}
						if (ii==pos[7]){strcpy(matgen[(i-10)].gq,g);}
						if (ii==pos[8]){strcpy(matgen[(i-10)].hq,g);}
						if (ii==pos[9]){strcpy(matgen[(i-10)].ps,g);}
						if (ii==pos[10]){strcpy(matgen[(i-10)].pq,g);}
						if (ii==pos[11]){strcpy(matgen[(i-10)].ec,g);}
						if (ii==pos[12]){strcpy(matgen[(i-10)].mq,g);}
						g=strtok_r(NULL,":",&spg);
					}//end while g
					matgen[(i-10)].gen1[col]='\0';
					matgen[(i-10)].gen2[col]='\0';
				}//end else i>9
				buf_split=strtok_r(NULL,"\t",&bsf);
			}//end while buf_split
			k++;//nr snp counter
		}//end reading snps
	}//end else lettura
	printf(" - %ld SNPs loaded.",k);
}//end while
if (k<col){printf("ERROR: vcfread: too few markers in the vcf file (Nr markers - present: %ld, expected: %ld).\n",k,col);return 1;}
free(buf);
return 0;
}//end function

//vcf format example
//#CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
//20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
//20     17330   .         T      A,G       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3


int loadPop(struct ind_gen *matgen, long int row, FILE *i_p)
{
/*This function assign the population ID to each individual present in the matgen object reading the informations from a text file.
The file needs to have only two columns: the first with the individual ID and the second with the population ID. The column must be
spaced with TAB or SPACE, other character are not admitted.
Note: this function is called only with vcf files, in the other cases the pop info is included in the genotype input file 
Known Bugs: it does not work if there are multiple tabs or spaces between fields*/
char *buf;
char *buf_split, *spf;
long int line=0;
long int ok=0;
long int chk=0;
long int i;
long int size_buf=200;//size of the reading line buffer

buf=malloc(size_buf*sizeof(char));
printf("\nStart reading populations file...\n");
while (fgets(buf,size_buf,i_p) !=NULL){
	line++;
	if (strchr(buf,' ')==NULL & strchr(buf,'\t')==NULL ){printf("ERROR: loadPop: wrong separator character in line %ld. Please use only TAB or SPACE.\n",line);return 1;}
	choppy(buf);
	buf_split=strtok_r(buf," \t",&spf);
	chk=0;
	for (i=0;i<row;i++){
		if (strcmp(matgen[i].name_ind,buf_split)==0){
			buf_split=strtok_r(NULL," \t",&spf);
			//printf("\nAAA: %s\n",buf_split);
			if (buf_split!=NULL){strncpy(matgen[i].pop,buf_split,20);ok++;chk=1;break;printf("\nSI\n");}
			else{printf("ERROR: loadPop: only 1 field present in popfile. (Check line %ld)\n",line);return 1;}
		}
	}
	if (chk==0){printf("WARNING: loadPop: no match in the genotypes file for the individual called: %s\n",buf_split);}
}
if (ok!=row){printf("ERROR: loadPop: some individual have no population ID. Please check the population file or the individual ID in the genotype input file.\n");return 1;}
return 0;
printf("\nStart reading populations file...done!");
free(buf);
}

int loadAncAll(struct ind_map *matmap, long int col, FILE *i_a)
{
/*This function assign the ancestral allele information to each SNP present in the matmap object reading the informations from a text file.
The file needs to have only two columns: the first with the ID code of a SNP (the "rs" or a single number) and the second with the ancestral allele (only one allele allowed).
Note: if the first column contains only numbers, then the "rs" tag will be added for checking the correspondace.
Note:The columns must be spaced with TAB or SPACE, other character are not permitted.
Note: this function is called only with plink and vcf files.
Known Bugs: it does not work if there are multiple tabs or spaces between fields*/
char *buf;
char *buf_split, *spf;
long int line=0,tot=0;
long int i;
long int size_buf=200;//size of the reading line buffer
char *tmp;
short int *ind;

buf=malloc(size_buf*sizeof(char));
tmp=malloc(sizeof(buf)+(sizeof(char)*100));
ind=malloc(col*sizeof(short int));
for (i=0;i<col;i++){ind[i]=0;}
printf("Reading ancestral alleles file");
while (fgets(buf,size_buf,i_a) !=NULL){
	line++;
	//printf("\nAnc all line %ld",line);
	if (strchr(buf,' ')==NULL & strchr(buf,'\t')==NULL ){
		printf("ERROR: loadAncAll: wrong separator character in line %ld. Please use only TAB or SPACE.\n",line);
		free(buf);
		free(tmp);
		free(ind);
		return 1;
		}
	choppy(buf);
	buf_split=strtok_r(buf," \t",&spf);
	for (i=0;i<col;i++){
		//printf("\nRow : %ld",i);
		if (ind[i]!=0){continue;}
		if (strstr(matmap[i].rs,"rs")==NULL){
			//printf("\nAAMatmap: %s - Field: %s\n",matmap[i].rs,buf_split);
			if (strcmp(matmap[i].rs,buf_split)==0){
				buf_split=strtok_r(NULL," \t",&spf);
				if (buf_split!=NULL){
					matmap[i].anc=*buf_split;
					tot++;
					ind[i]=1;
					break;
				}
				else{
					printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);
					free(buf);
					free(tmp);
					free(ind);
					return 1;
				}
			}
			else{
				continue;
			}
		}
		else{
			strncpy(tmp,"rs",3);
			strcat(tmp, buf_split);
			//printf("\nFieldCat: %s\n",tmp);
			if (strcmp(matmap[i].rs,tmp)==0){
				//printf("\nMatmap: %s - Field: %s\n",matmap[i].rs,buf_split);
				buf_split=strtok_r(NULL," \t",&spf);
				if (buf_split!=NULL){
					matmap[i].anc=*buf_split;
					tot++;
					ind[i]=1;
					//printf("\nLine %ld - Found match! Tot match1:%ld",line,tot);
					break;
				}
				else{
					printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);
					free(buf);
					free(tmp);
					free(ind);
					return 1;
				}
			}
			//else{
			//	continue;
			//}
			
		}

	}
	if (tot==col){break;}
}
if (tot!=col){printf("WARNING: loadAncAll: some SNPs have no ancestral allele information. (%ld SNPs with allele information - %ld SNPs in total).\n",tot,col);}
else{printf("\t\t\t[ancestral state loaded for %ld/%ld SNPs]\n",tot,col);}
free(buf);
free(tmp);
free(ind);
return 0;
}

int loadAncAll1(struct ind_map *matmap, long int col, FILE *i_a)
{
/*This function assign the ancestral allele information to each SNP present in the matmap object reading the informations from a text file.
The file needs to have three columns: chromosome, position and the ancestral allele (only one allele allowed).
Note:The columns must be spaced with TAB or SPACE, other character are not permitted.
Note: this function is called only with plink and vcf files.
Known Bugs: it does not work if there are multiple tabs or spaces between fields*/
char *buf;
char *buf_split, *spf;
long int line=0,tot=0;
long int i;
long int size_buf=200;//size of the reading line buffer
char *tmp;
short int *ind;
short int crom;
long int bp_pos;
char all;

buf=malloc(size_buf*sizeof(char));
tmp=malloc(sizeof(buf)+(sizeof(char)*100));
ind=malloc(col*sizeof(short int));
for (i=0;i<col;i++){ind[i]=0;}
printf("Reading ancestral alleles file");
while (fgets(buf,size_buf,i_a) !=NULL){
	line++;
	//printf("\nAnc all line %s",buf);
	if (strchr(buf,' ')==NULL & strchr(buf,'\t')==NULL ){
		printf("ERROR: loadAncAll: wrong separator character in line %ld. Please use only TAB or SPACE.\n",line);
		free(buf);
		free(tmp);
		free(ind);
		return 1;
		}
	choppy(buf);
	all='N';
	for (i=0;i<3;i++){
		if (i==0){buf_split=strtok_r(buf," \t",&spf);}else{buf_split=strtok_r(NULL," \t",&spf);}
		if (buf_split!=NULL){
			if(i==0){ crom = (short unsigned)atoi(buf_split);}
			else if (i==1){ bp_pos = (long int)atol(buf_split);}
			else if (i==2){ all = *buf_split;}
		}
		else{
			printf("ERROR: loadAncAll: bad AncAllFile format. (Check line %ld)\n",line);
			free(buf);
			free(tmp);
			free(ind);
			return 1;
		}
	}
	//printf("\nCrom: %d, pos: %ld, ancall: %c",crom,bp_pos,all);
	for (i=0;i<col;i++){
		//printf("\nRow : %ld",i);
		if (ind[i]!=0){continue;}
		if (matmap[i].crom == crom ){
			//printf("\nAAMatmap: %s - Field: %s\n",matmap[i].rs,buf_split);
			if (matmap[i].bp_pos == bp_pos){
					matmap[i].anc=all;
					tot++;
					ind[i]=1;
					break;
	
			}
		}

	}
	if (tot==col){break;}
}
if (tot!=col){printf("WARNING: loadAncAll: some SNPs have no ancestral allele information. (%ld SNPs with allele information - %ld SNPs in total).\n",tot,col);}
else{printf("\t\t\t[ancestral state loaded for %ld/%ld SNPs]\n",tot,col);}
free(buf);
free(tmp);
free(ind);
return 0;
}




int loadAncAll_fast(struct ind_map *matmap, long int col, FILE *i_a)
{
/*This function assign the ancestral allele information to each SNP present in the matmap object reading the informations from a text file.
The file needs to have only two columns: the first with the ID code of a SNP (the "rs" or a single number) and the second with the ancestral allele (only one allele allowed). 
The columns must be spaced with TAB or SPACE, other character are not permitted.
Note: this function is called only with plink and vcf files, in the other cases the pop info is included in the genotype input file.
Known Bugs: it does not work if there are multiple tabs or spaces between fields*/
char *buf;
char *buf_split, *spf;
long int line=0,tot=0;
long int i,j,p;
long int size_buf=200;//size of the reading line buffer
char *tmp;
long int npos;
long int *pos; 

buf=malloc(size_buf*sizeof(char));
tmp=malloc(sizeof(buf)+(sizeof(char)*100));
pos=malloc(col*sizeof(long int));
for (i=0;i<col;i++){pos[i]=i;}
npos=col;
printf("\nStart reading ancestral allele state file...\n");
while (fgets(buf,size_buf,i_a) !=NULL){
	line++;
	printf("\nAnc all line %ld",line);
	if (strchr(buf,' ')==NULL & strchr(buf,'\t')==NULL ){printf("ERROR: loadAncAll: wrong separator character in line %ld. Please use only TAB or SPACE.\n",line);return 1;}
	choppy(buf);
	buf_split=strtok_r(buf," \t",&spf);
	for (i=0;i<npos;i++){
		p=pos[i];
		//printf("\necco la posizione p : %ld",p);
		if (strstr(matmap[p].rs,"rs")==NULL){
			//printf("\nAAMatmap: %s - Campo: %s\n",matmap[i].rs,buf_split);
			if (strcmp(matmap[p].rs,buf_split)==0){
				buf_split=strtok_r(NULL," \t",&spf);
				if (buf_split!=NULL){
					matmap[p].anc=*buf_split;
					tot++;
					npos--;
					for (j=p;j<npos;j++){
						pos[j]=pos[j+1];
					}
					break;
				}
				else{printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);return 1;}
			}
			else{
				continue;
			}
		}
		else{
			//printf("\nMatmap: %s - Campo: %s\n",matmap[i].rs,buf_split);
			strncpy(tmp,"rs",3);
			strcat(tmp, buf_split);
			//printf("\nCampoCat: %s\n",tmp);
			if (strcmp(matmap[p].rs,tmp)==0){
				printf("\nMatmap: %s - Campo: %s\n",matmap[p].rs,buf_split);
				buf_split=strtok_r(NULL," \t",&spf);
				if (buf_split!=NULL){
					matmap[p].anc=*buf_split;
					tot++;
					printf("\nLine %ld - Found match! Tot match:%ld",line,tot);
					npos--;
					for (j=p;j<npos;j++){
						pos[j]=pos[j+1];
					}
					break;
				}
				else{printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);return 1;}
			}
			else{
				continue;
			}
			
		}

	}
	if (tot==col){break;}
}
if (tot!=col){printf("WARNING: loadAncAll: some SNPs have no ancestral allele information. (%ld SNPs with allele information - %ld SNPs in total).\n",tot,col);}
else{printf("\nAncestral state loaded for %ld/%ld SNPs",tot,col);}
return 0;
printf("\nStart reading ancestral allele state file...done!");
free(buf);
free(tmp);
free(pos);
}

int loadAncAll_mp(struct ind_map *matmap, long int col, FILE *i_a, int nt)
{
/*This function assign the ancestral allele information to each SNP present in the matmap object reading the informations from a text file.
The file needs to have only two columns: the first with the ID code of a SNP (the "rs" or a single number) and the second with the ancestral allele (only one allele allowed). 
The columns must be spaced with TAB or SPACE, other character are not permitted.
Note: this function is called only with plink and vcf files, in the other cases the pop info is included in the genotype input file.
Known Bugs: it does not work if there are multiple tabs or spaces between fields*/
/*NB: this does not work, check mp pragma*/
char *buf;
char *buf_split, *spf;
long int line=0,tot=0;
long int i;
long int size_buf=200;//size of the reading line buffer
char tmp[200];
char al;
char rs[200];

int iCPU = omp_get_num_procs();/*maximum number of threads in the system*/
if (nt>iCPU){nt=iCPU;}
omp_set_num_threads(nt);/* Now set the number of threads*/

buf=malloc(size_buf*sizeof(char));
printf("\nStart reading ancestral allele state file...\n");
while (fgets(buf,size_buf,i_a) !=NULL){
	line++;
	printf("\nAnc all line %ld",line);
	if (strchr(buf,' ')==NULL & strchr(buf,'\t')==NULL ){printf("ERROR: loadAncAll: wrong separator character in line %ld. Please use only TAB or SPACE.\n",line);return 1;}
	choppy(buf);
	buf_split=strtok_r(buf," \t",&spf);
	strcpy(rs,buf_split);
	buf_split=strtok_r(NULL," \t",&spf);
	al=*buf_split;
	#pragma omp parallel for private(i,tmp,rs,al,line) reduction(+:tot)
	for (i=0;i<col;i++){
		//printf("\necco la linea i : %ld",i);
		strncpy(tmp,"rs",3);
		strcat(tmp, rs);
		if (strcmp(matmap[i].rs,tmp)==0){
			//printf("\nMatmap: %s - Campo: %s\n",matmap[i].rs,rs);
			matmap[i].anc=al;
			tot++;
			printf("\nLine %ld - Found match! rs:%s Tot match1",line,rs);
			//else{printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);return 1;}
		}
		//printf("\nciao\n");
		//else{
		//	continue;
		//}
	}
	//if (tot==col){break;}
if (line==10000){
printf("\ntot %ld",tot);
break;}
}
if (tot!=col){printf("WARNING: loadAncAll: some SNPs have no ancestral allele information. (%ld SNPs with allele information - %ld SNPs in total).\n",tot,col);}
else{printf("\nAncestral state loaded for %ld/%ld SNPs",tot,col);}
return 0;
printf("\nStart reading ancestral allele state file...done!");
free(buf);
}

int loadAncAll_alex(struct ind_map *matmap, long int col, FILE *i_a)
{
/*This function assign the ancestral allele information to each SNP present in the matmap object reading the informations from a text file.
The file needs to have only two columns: the first with the ID code of a SNP (the "rs" or a single number) and the second with the ancestral allele (only one allele allowed). 
The columns must be spaced with TAB or SPACE, other character are not permitted.
Note: this function is called only with plink and vcf files, in the other cases the pop info is included in the genotype input file.
Known Bugs: it does not work if there are multiple tabs or spaces between fields*/
/*fare una funzione ordinando preventivamente gli rs di matmap e ottimizzando la ricerca*/
char *buf;
char *buf_split, *spf;
long int line=0,tot=0;
long int i;
long int size_buf=200;//size of the reading line buffer
char *tmp;
short int *ind;

buf=malloc(size_buf*sizeof(char));
tmp=malloc(sizeof(buf)+(sizeof(char)*100));
//opos=malloc(col*sizeof());


for (i=0;i<col;i++){ind[i]=0;}

printf("\nStart reading ancestral allele state file...\n");
while (fgets(buf,size_buf,i_a) !=NULL){
	line++;
	printf("\nAnc all line %ld",line);
	if (strchr(buf,' ')==NULL & strchr(buf,'\t')==NULL ){printf("ERROR: loadAncAll: wrong separator character in line %ld. Please use only TAB or SPACE.\n",line);return 1;}
	choppy(buf);
	buf_split=strtok_r(buf," \t",&spf);
	for (i=0;i<col;i++){
		//printf("\necco la linea i : %ld",i);//
		if (ind[i]!=0){continue;}
		if (strstr(matmap[i].rs,"rs")==NULL){
			//printf("\nAAMatmap: %s - Campo: %s\n",matmap[i].rs,buf_split);//
			if (strcmp(matmap[i].rs,buf_split)==0){
				buf_split=strtok_r(NULL," \t",&spf);
				if (buf_split!=NULL){
					matmap[i].anc=*buf_split;
					tot++;
					ind[i]=1;
					break;
				}
				else{printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);return 1;}
			}
			else{
				continue;
			}
		}
		else{
			strncpy(tmp,"rs",3);
			strcat(tmp, buf_split);
			//printf("\nCampoCat: %s\n",tmp);//
			if (strcmp(matmap[i].rs,tmp)==0){
				printf("\nMatmap: %s - Campo: %s\n",matmap[i].rs,buf_split);
				buf_split=strtok_r(NULL," \t",&spf);
				if (buf_split!=NULL){
					matmap[i].anc=*buf_split;
					tot++;
					ind[i]=1;
					printf("\nLine %ld - Found match! Tot match1:%ld",line,tot);
					break;
				}
				else{printf("ERROR: loadAncAll: only 1 field present in AncAllFile. (Check line %ld)\n",line);return 1;}
			}
			//else{
			//	continue;
			//}
			
		}

	}
	if (tot==col){break;}
}
if (tot!=col){printf("WARNING: loadAncAll: some SNPs have no ancestral allele information. (%ld SNPs with allele information - %ld SNPs in total).\n",tot,col);}
else{printf("\nAncestral state loaded for %ld/%ld SNPs",tot,col);}
return 0;
printf("\nStart reading ancestral allele state file...done!");
free(buf);
free(tmp);
free(ind);
}

