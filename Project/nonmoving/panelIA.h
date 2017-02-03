
typedef struct {
	double *v0;
	double *v1;
        double *ny;
        double len;
        int idxv0;
        int idxv1;
} panel;

void panelIA(int d, int p, double h, panel *p1, panel *p2, double xLeg[], double wLeg[], double solution[]);
