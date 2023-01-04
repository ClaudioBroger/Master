#include "mex.h"

const int NRSTATES = 3;
const int NRPARAMETERS = 13;
const int NRVARIABLES = 3;
const int NRREACTIONS = 6;
const int NREVENTS = 1;

const int hasOnlyNumericICs = 1;
double defaultICs_num[3] = {
	0,0,0};
char *defaultICs_nonnum[1];

double defaultParam[13] = {
	50,50,0.5,0.0077,0.01,1,0.5,2.5,250,5000,0.0077,1,0.0173};
char *stateNames[3] = {
	"TetR","Citimmature","Citrine"};
char *parameterNames[13] = {
	"k7tetTetR","k7tetCit","dTetR","dCit","kL7tet","thetaTetR","KdTetR","nTetR","atcAdded","indTime","mu","nMperUnit","kmaturation"};
char *variableNames[3] = {
	"atc","KdTetR_InUnit","TetRfree"};
char *variableFormulas[3] = {
	"piecewiseIQM(3,atcAdded/nMperUnit,ge(time,indTime),0.0);","KdTetR/nMperUnit;","piecewiseIQM(1,TetR/2.0-KdTetR_InUnit/2.0-atc/2.0+pow(pow(KdTetR_InUnit,2.0)+2.0*KdTetR_InUnit*TetR+2.0*KdTetR_InUnit*atc+pow(TetR,2.0)-2.0*TetR*atc+pow(atc,2.0),1.0/2.0)/2.0),ge(TetR,1e-4),0.0);"};
char *reactionNames[6] = {
	"R1","R2","R3","D1","D2","D3"};
char *eventNames[1] = {
	"piecewise_event_1"};

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
