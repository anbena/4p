/*prototypes of functions defined in genStat.c*/
struct ind_gen;
struct ind_map;

int computeBaseFreq(struct ind_gen *matgen, struct ind_map *matmap, long int row, long int col, char *namepop, char *type, int nt);

int hetOss(struct ind_gen *matgen, struct ind_map *matmap, char *namepop, long int col, long int row, char *type, int nt);

int hetExp(struct ind_gen *matgen, struct ind_map *matmap, char *namepop, long int col, long int row, char *type, int nt);

int afs(struct ind_gen *matgen, struct ind_map *matmap, long int p1, long int p2, long int p3, long int col, long int row, int unfolded, int header, int nt);

int fst_mp(struct ind_gen *matgen, char *namepop1, char *namepop2, long int col, long int row, float *v_fst, int nt);

int dist(struct ind_gen *matgen, struct ind_map *matmap, char *namepop1, char *namepop2, long int col, long int row, double *outv, char *type);

int dist_mp(struct ind_gen *matgen, struct ind_map *matmap, char *namepop1, char *namepop2, long int col, long int row, double *outv, char *type, int nt);

int shAll(struct ind_gen *matgen, long int i, long int j, long int col, double *outv, int similarity);

int shAll_mp(struct ind_gen *matgen, long int i, long int j, long int col, double *outv, int similarity, int nt);
