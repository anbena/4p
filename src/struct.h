/*prototypes of structures used*/

struct ind_gen//PEDFILE-LIKE STRUCTURE: we allocate here informations about each individual present in the ped file
{
char pop[20];
char name_ind[50];
char paternal_id[20];
char maternal_id[20];
char sex;
int phen;
unsigned int phase;
int dp;
char ft[100];
char gl[100];// float but for now string 100
char gle[100];
char pl[100];
char gp[100];
char gq[100];
char hq[100];
char ps[100];
char pq[100];
char ec[100];
char mq[100];
char *gen1;//1°crom genotype
char *gen2;//2°crom genotype
};

struct ind_map//MAPFILE-LIKE STRUCTURE: we allocate here informations about each SNPs present in the map file
{
short unsigned int crom;
char rs[20];
short int gen_dist;
long int bp_pos;
char anc;
char ref;
char alt[20];//A,T,C,G,N
int qual;
char filt[10];
char info[100];
};

