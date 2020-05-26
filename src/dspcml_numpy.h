
#include "dspcml.h"

unsigned int print_matrix(unsigned long ulptr, char *cptr);
unsigned long zeros(size_t rows, size_t cols);
unsigned long rows(unsigned long ulptr);
unsigned long cols(unsigned long ulptr);
unsigned long numel(unsigned long ulptr);
void to_numpy(unsigned long ulptr, cml_real_t *out_data);

