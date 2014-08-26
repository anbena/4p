/*prototypes of functions defined in parsFiles.c */
struct ind_gen;
struct ind_map;

int loadPedSnp(struct ind_gen *matgen, long int col, long int row, FILE *i_f);

int loadMapSnp(struct ind_map *matmap, long int col, FILE *i_m);

int loadArlSnp(struct ind_gen *matgen, long int col, long int row, FILE *i_f);

int loadVcfSnp(struct ind_gen *matgen, struct ind_map *matmap, long int col, long int row, FILE *i_f);

int loadPop(struct ind_gen *matgen, long int row, FILE *i_p);

int loadAncAll(struct ind_map *matmap, long int col, FILE *i_a);

int loadAncAll1(struct ind_map *matmap, long int col, FILE *i_a);

int loadAncAll_fast(struct ind_map *matmap, long int col, FILE *i_a);

int loadAncAll_fast_mp(struct ind_map *matmap, long int col, FILE *i_a, int nt);

int loadAncAll_alex(struct ind_map *matmap, long int col, FILE *i_a);


