typedef struct {
        matrix_TYP* element;
        int* product;
        int nproduct;
        int left, right;
        } derived_TYP;
typedef struct{
        derived_TYP** list;
        int sizemult;
        int firstfree;} derivedsg_TYP;
#define NIL (-1)
