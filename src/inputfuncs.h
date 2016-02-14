#ifndef INP_F
#define INP_F

void set_input_params(int argc, char **argv, int *num_iters, int *print_freq, char *filename, double *stop_delta);
void print_input_help();
void print_input_params(int rank, int num_iters, int print_freq, char *filename, double stop_delta);

#endif
