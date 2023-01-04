#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "MEX_EvoParam.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double TetR,Citimmature,Citrine;
    double k7tetTetR,k7tetCit,dTetR,dCit,kL7tet,thetaTetR,KdTetR,nTetR,atcAdded,indTime,mu,nMperUnit,kmaturation;
    double atc,KdTetR_InUnit,TetRfree;
    double R1,R2,R3,D1,D2,D3;
    double eventassign_1_1;

    time = time_local;

    TetR = stateVector[0];
    Citimmature = stateVector[1];
    Citrine = stateVector[2];
    k7tetTetR = paramdataPtr->parametervector[0]; /* 50 */
    k7tetCit = paramdataPtr->parametervector[1]; /* 50 */
    dTetR = paramdataPtr->parametervector[2]; /* 0.5 */
    dCit = paramdataPtr->parametervector[3]; /* 0.0077 */
    kL7tet = paramdataPtr->parametervector[4]; /* 0.01 */
    thetaTetR = paramdataPtr->parametervector[5]; /* 1 */
    KdTetR = paramdataPtr->parametervector[6]; /* 0.5 */
    nTetR = paramdataPtr->parametervector[7]; /* 2.5 */
    atcAdded = paramdataPtr->parametervector[8]; /* 250 */
    indTime = paramdataPtr->parametervector[9]; /* 5000 */
    mu = paramdataPtr->parametervector[10]; /* 0.0077 */
    nMperUnit = paramdataPtr->parametervector[11]; /* 1 */
    kmaturation = paramdataPtr->parametervector[12]; /* 0.0173 */
    atc = piecewiseIQM(3,atcAdded/nMperUnit,ge(time,indTime),0.0);;
    KdTetR_InUnit = KdTetR/nMperUnit;;
    TetRfree = piecewiseIQM(1,TetR/2.0-KdTetR_InUnit/2.0-atc/2.0+pow(pow(KdTetR_InUnit,2.0)+2.0*KdTetR_InUnit*TetR+2.0*KdTetR_InUnit*atc+pow(TetR,2.0)-2.0*TetR*atc+pow(atc,2.0),1.0/2.0)/2.0),ge(TetR,1e-4),0.0);;
    R1 = k7tetTetR*(kL7tet+(1.0-kL7tet)/(1.0+pow(TetRfree/thetaTetR,nTetR)+pow(Tup1free/thetaTup1,nTetR)));
    R2 = k7tetCit*(kL7tet+(1.0-kL7tet)/(1.0+pow(TetRfree/thetaTetR,nTetR)+pow(Tup1free/thetaTup1,nTetR)));
    R3 = kmaturation*Citimmature;
    D1 = (dTetR+dTetRadh1+mu)*TetR;
    D2 = (dCit+mu)*Citimmature;
    D3 = (dCit+mu)*Citrine;
    eventassign_1_1 = TetR;
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = R1-D1;
    	DDTvector[1] = R2-D2;
    	DDTvector[2] = R3-D3;
    } else if (DOflag == DOFLAG_VARREAC) {
        variableVector[0] = atc;
        variableVector[1] = KdTetR_InUnit;
        variableVector[2] = TetRfree;
        reactionVector[0] = R1;
        reactionVector[1] = R2;
        reactionVector[2] = R3;
        reactionVector[3] = D1;
        reactionVector[4] = D2;
        reactionVector[5] = D3;
    } else if (DOflag == DOFLAG_EVENTS) {
        gout[0] = ge(time,indTime)-0.5;
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
        if (eventVector[0] == 1 && gout[0] < 0) {
            DDTvector[0] = 1;
            stateVector[0] = eventassign_1_1;
        }
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double TetR,Citimmature,Citrine;
    double k7tetTetR,k7tetCit,dTetR,dCit,kL7tet,thetaTetR,KdTetR,nTetR,atcAdded,indTime,mu,nMperUnit,kmaturation;
    double atc,KdTetR_InUnit,TetRfree;
    k7tetTetR = paramdataPtr->parametervector[0]; /* 50 */
    k7tetCit = paramdataPtr->parametervector[1]; /* 50 */
    dTetR = paramdataPtr->parametervector[2]; /* 0.5 */
    dCit = paramdataPtr->parametervector[3]; /* 0.0077 */
    kL7tet = paramdataPtr->parametervector[4]; /* 0.01 */
    thetaTetR = paramdataPtr->parametervector[5]; /* 1 */
    KdTetR = paramdataPtr->parametervector[6]; /* 0.5 */
    nTetR = paramdataPtr->parametervector[7]; /* 2.5 */
    atcAdded = paramdataPtr->parametervector[8]; /* 250 */
    indTime = paramdataPtr->parametervector[9]; /* 5000 */
    mu = paramdataPtr->parametervector[10]; /* 0.0077 */
    nMperUnit = paramdataPtr->parametervector[11]; /* 1 */
    kmaturation = paramdataPtr->parametervector[12]; /* 0.0173 */
    atc = piecewiseIQM(3,atcAdded/nMperUnit,ge(time,indTime),0.0);;
    KdTetR_InUnit = KdTetR/nMperUnit;;
    TetRfree = piecewiseIQM(1,TetR/2.0-KdTetR_InUnit/2.0-atc/2.0+pow(pow(KdTetR_InUnit,2.0)+2.0*KdTetR_InUnit*TetR+2.0*KdTetR_InUnit*atc+pow(TetR,2.0)-2.0*TetR*atc+pow(atc,2.0),1.0/2.0)/2.0),ge(TetR,1e-4),0.0);;
    TetR = 0.0;
    Citimmature = 0.0;
    Citrine = 0.0;
    icVector[0] = TetR;
    icVector[1] = Citimmature;
    icVector[2] = Citrine;
}

