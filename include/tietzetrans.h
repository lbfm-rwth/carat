typedef struct
        {
        int* lhsproduct;
        int lhsnproduct;
        int* rhsproduct;
        int rhsnproduct;
        } relator_TYP;

typedef struct 
        {
        derivedsg_TYP* generators;
        relator_TYP* relators;
        int norelators;
        int ext_factor;
        } presentation_TYP;

