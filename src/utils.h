/*prototypes of functions defined in utils.c */
struct ind_gen;
struct ind_map;

void logo(void);

void help(void);

void choppy(char *s);

int uniquePop(struct ind_gen *matgen, int n, int *pos, int *upop);

char flipBase(char base);

int outGen(struct ind_gen *matgen,long int row, long int col);

int outMap(struct ind_map *matmap, long int col);

int outDistMat(struct ind_gen *matgen, int upop, int *pos,double *outvt);

char getAll(char ref, char *alt, int x);

int updateAllInfo(struct ind_gen *matgen, struct ind_map *matmap, long int col, long int row, int inputType);

int updateAllInfo_mp(struct ind_gen *matgen, struct ind_map *matmap, long int col, long int row, int inputType, int nt);

void structDeall(struct ind_gen *matgen, struct ind_map *matmap, long int row);
