/*===========================================================================*
 * This file is part of the Branch, Constrain and Price Software (BiCePS)    *
 *                                                                           *
 * BiCePS is distributed under the Common Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 * Copyright (C) 2001-2007, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


//#############################################################################

#include "AlpsKnowledgeBroker.h"

#include "BcpsTreeNode.h"
#include "BcpsNodeDesc.h"

//#############################################################################

BcpsNodeDesc:: ~BcpsNodeDesc() { 
    int k;

    //------------------------------------------------------
    // Free variable related.
    //------------------------------------------------------

    delete [] vars_->posRemove; 
    vars_->posRemove = NULL;
    
    for (k = 0; k < vars_->numAdd; ++k) {
        delete vars_->objects[k];
    }
    delete [] vars_->objects; 
    vars_->objects = NULL;
    
    delete [] vars_->lbHard.posModify; 
    vars_->lbHard.posModify = NULL;
    delete [] vars_->lbHard.entries; 
    vars_->lbHard.entries = NULL;
    
    delete [] vars_->ubHard.posModify; 
    vars_->ubHard.posModify = NULL;
    delete [] vars_->ubHard.entries; 
    vars_->ubHard.entries = NULL;
    
    delete [] vars_->lbSoft.posModify; 
    vars_->lbSoft.posModify = NULL;
    delete [] vars_->lbSoft.entries; 
    vars_->lbSoft.entries = NULL;
    
    delete [] vars_->ubSoft.posModify; 
    vars_->ubSoft.posModify = NULL;
    delete [] vars_->ubSoft.entries; 
    vars_->ubSoft.entries = NULL;
    
    delete vars_;
    vars_ = NULL;
    
    //------------------------------------------------------
    // Free constraint related.
    //------------------------------------------------------
    
    delete [] cons_->posRemove; 
    cons_->posRemove = NULL;
    
    for (k = 0; k < cons_->numAdd; ++k) {
        delete cons_->objects[k];
    }
    delete [] cons_->objects;
    cons_->objects = NULL;
    
    delete [] cons_->lbHard.posModify; 
    cons_->lbHard.posModify = NULL;
    delete [] cons_->lbHard.entries; 
    cons_->lbHard.entries = NULL;
  
    delete [] cons_->ubHard.posModify; 
    cons_->ubHard.posModify = NULL;
    delete [] cons_->ubHard.entries; 
    cons_->ubHard.entries = NULL;
  
    delete [] cons_->lbSoft.posModify; 
    cons_->lbSoft.posModify = NULL;
    delete [] cons_->lbSoft.entries; 
    cons_->lbSoft.entries = NULL;
    
    delete [] cons_->ubSoft.posModify; 
    cons_->ubSoft.posModify = NULL;
    delete [] cons_->ubSoft.entries;
    cons_->ubSoft.entries = NULL;
    
    delete cons_;
    cons_ = NULL;
}
    

//#############################################################################

void BcpsNodeDesc::initToNull() {

    //------------------------------------------------------
    // Initialize the list of modified variable
    //------------------------------------------------------

    vars_ = new BcpsObjectListMod;
    cons_ = new BcpsObjectListMod;

    vars_->numRemove = 0;
    vars_->posRemove = NULL;
  
    vars_->numAdd = 0;
    vars_->objects = NULL;
  
    vars_->lbHard.relative = false;
    vars_->lbHard.numModify = 0;
    vars_->lbHard.posModify = NULL;
    vars_->lbHard.entries = NULL;
  
    vars_->ubHard.relative = false;
    vars_->ubHard.numModify = 0;
    vars_->ubHard.posModify = NULL;
    vars_->ubHard.entries = NULL;

    vars_->lbSoft.relative = false;
    vars_->lbSoft.numModify = 0;
    vars_->lbSoft.posModify = NULL;
    vars_->lbSoft.entries = NULL;
    
    vars_->ubSoft.relative = false;
    vars_->ubSoft.numModify = 0;
    vars_->ubSoft.posModify = NULL;
    vars_->ubSoft.entries = NULL;
    
    //------------------------------------------------------
    // Initialize the list of modified constraints.
    //------------------------------------------------------

    cons_->numRemove = 0;
    cons_->posRemove = NULL;

    cons_->numAdd = 0;
    cons_->objects = NULL;

    cons_->lbHard.relative = false;
    cons_->lbHard.numModify = 0;
    cons_->lbHard.posModify = NULL;
    cons_->lbHard.entries = NULL;
    
    cons_->ubHard.relative = false;
    cons_->ubHard.numModify = 0;
    cons_->ubHard.posModify = NULL;
    cons_->ubHard.entries = NULL;
    
    cons_->lbSoft.relative = false;
    cons_->lbSoft.numModify = 0;
    cons_->lbSoft.posModify = NULL;
    cons_->lbSoft.entries = NULL;
    
    cons_->ubSoft.relative = false;
    cons_->ubSoft.numModify = 0;
    cons_->ubSoft.posModify = NULL;
    cons_->ubSoft.entries = NULL;
}

//#############################################################################

void BcpsNodeDesc::setVars(int numRem,
			   const int    *posRem,
			   int numAdd,
			   const BcpsObject **objects,
			   bool relvlh,
			   int numvlh, 
			   const int    *vlhp,
			   const double *vlhe,
			   bool relvuh,
			   int numvuh,
			   const int    *vuhp,
			   const double *vuhe,
			   bool relvls,
			   int numvls,
			   const int    *vlsp,
			   const double *vlse,
			   bool relvus,
			   int numvus,
			   const int    *vusp,
			   const double *vuse)
{

    //------------------------------------------------------
    // Removed var.
    //------------------------------------------------------

    vars_->numRemove = numRem;
    if (numRem > 0) {
	int *posRemL = new int [numRem]; 
	memcpy(posRemL, posRem, sizeof(int) * numRem);
	vars_->posRemove = posRemL; 
	posRemL = NULL;
    }
    else {
	vars_->posRemove = NULL;
    }
    
    //------------------------------------------------------
    // Added var.
    //------------------------------------------------------

    vars_->numAdd = numAdd;
    if (numAdd > 0) {
	BcpsObject **objL = new BcpsObject* [numAdd]; 
	memcpy(objL, objects, sizeof(BcpsObject*) * numAdd);
	vars_->objects = objL; 
	objL = NULL;
    }
    else {
	vars_->objects = NULL;
    } 
    
    //------------------------------------------------------
    // Modified var hard lb.
    //------------------------------------------------------
    
    vars_->lbHard.relative = relvlh;
    vars_->lbHard.numModify = numvlh;
    if (numvlh > 0) {
	int    *vlhpL = new int [numvlh];
	double *vlheL = new double [numvlh];
	memcpy(vlhpL, vlhp, numvlh * sizeof(int));
	memcpy(vlheL, vlhe, numvlh * sizeof(double));
	vars_->lbHard.posModify = vlhpL;  vlhpL = NULL;
	vars_->lbHard.entries = vlheL;    vlheL = NULL;
    }
    else {
	vars_->lbHard.posModify = NULL;
	vars_->lbHard.entries = NULL;
    }  

    //------------------------------------------------------
    // Modified var hard ub.
    //------------------------------------------------------
    
    vars_->ubHard.relative = relvuh;
    vars_->ubHard.numModify = numvuh;
    if (numvuh > 0) {
	int    *vuhpL = new int [numvuh];
	double *vuheL = new double [numvuh];
	memcpy(vuhpL, vuhp, numvuh * sizeof(int));
	memcpy(vuheL, vuhe, numvuh * sizeof(double));       
    }
    else {
	vars_->ubHard.posModify = NULL;
	vars_->ubHard.entries = NULL;
    }
    
    //------------------------------------------------------
    // Modified var soft lb.
    //------------------------------------------------------

    vars_->lbSoft.relative = relvls;
    vars_->lbSoft.numModify = numvls;
    if (numvls > 0) {
	int    *vlspL = new int [numvls];
	double *vlseL = new double [numvls];
	memcpy(vlspL, vlsp, numvls * sizeof(int));
	memcpy(vlseL, vlse, numvls * sizeof(double));
	vars_->lbSoft.posModify = vlspL;  vlspL = NULL;
	vars_->lbSoft.entries = vlseL;    vlseL = NULL;
    }
    else {
	vars_->lbSoft.posModify = NULL;
	vars_->lbSoft.entries = NULL;
    }
      
    //------------------------------------------------------
    // Modified var soft ub.
    //------------------------------------------------------

    vars_->ubSoft.relative = relvus;
    vars_->ubSoft.numModify = numvus;
    if (numvus > 0) {
	int    *vuspL = new int [numvus];
	double *vuseL = new double [numvus];
	memcpy(vuspL, vusp, numvus * sizeof(int));
	memcpy(vuseL, vuse, numvus * sizeof(double));
	vars_->ubSoft.posModify = vuspL;  vuspL = NULL;
	vars_->ubSoft.entries = vuseL;    vuseL = NULL;
    }
    else {
	vars_->ubSoft.posModify = NULL;
	vars_->ubSoft.entries = NULL;
    }
}

//#############################################################################

void BcpsNodeDesc::assignVars(int numRem,
			      int    *&posRem,
			      int numAdd,
			      BcpsObject **&objects,
			      bool relvlh, int numvlh, int *&vlhp, double *&vlhe,
			      bool relvuh, int numvuh, int *&vuhp, double *&vuhe,
			      bool relvls, int numvls, int *&vlsp, double *&vlse,
			      bool relvus, int numvus, int *&vusp, double *&vuse)
{

    //------------------------------------------
    // Removed vars.
    //------------------------------------------

    vars_->numRemove = numRem;
    vars_->posRemove = posRem; posRem = NULL;

    //------------------------------------------
    // Added vars.
    //------------------------------------------
    
    vars_->numAdd = numAdd;
    vars_->objects = objects; objects = NULL; 

    //------------------------------------------
    // Modified var hard lb.
    //------------------------------------------

    vars_->lbHard.relative = relvlh;
    vars_->lbHard.numModify = numvlh;
    vars_->lbHard.posModify = vlhp;  vlhp = NULL;
    vars_->lbHard.entries = vlhe;    vlhe = NULL;

    //------------------------------------------  
    // Modified var hard u.
    //------------------------------------------

    vars_->ubHard.relative = relvuh;
    vars_->ubHard.numModify = numvuh;
    vars_->ubHard.posModify = vuhp;  vuhp = NULL;
    vars_->ubHard.entries = vuhe;    vuhe = NULL;
    
    //------------------------------------------
    // Modified var soft lb.
    //------------------------------------------

    vars_->lbSoft.relative = relvls;
    vars_->lbSoft.numModify = numvls;
    vars_->lbSoft.posModify = vlsp;  vlsp = NULL;
    vars_->lbSoft.entries = vlse;    vlse = NULL;
    
    //------------------------------------------      
    // Modified var soft ub.
    //------------------------------------------

    vars_->ubSoft.relative = relvus;
    vars_->ubSoft.numModify = numvus;
    vars_->ubSoft.posModify = vusp;  vusp = NULL;
    vars_->ubSoft.entries = vuse;    vuse = NULL;
}

//#############################################################################

void BcpsNodeDesc::setCons(int numRem,
			   const int    *posRem,
			   int numAdd,
			   const BcpsObject **objects,
			   bool relclh,
			   int numclh,
			   const int    *clhp,
			   const double *clhe,
			   bool relcuh,
			   int numcuh,
			   const int    *cuhp,
			   const double *cuhe,
			   bool relcls,
			   int numcls,
			   const int    *clsp,
			   const double *clse,
			   bool relcus,
			   int numcus,
			   const int    *cusp,
			   const double *cuse)
{
    //------------------------------------------------------
    // Removed con.
    //------------------------------------------------------
    
    cons_->numRemove = numRem;
    if (numRem > 0) {
	int *posRemL = new int [numRem]; 
	memcpy(posRemL, posRem, sizeof(int) * numRem);
	cons_->posRemove = posRemL; 
	posRemL = NULL;
    }
    else {
	cons_->posRemove = NULL;
    }

    //------------------------------------------------------
    // Added con.
    //------------------------------------------------------
    
    cons_->numAdd = numAdd;
    if (numAdd > 0) {
	BcpsObject **objL = new BcpsObject* [numAdd]; 
	memcpy(objL, objects, sizeof(BcpsObject*) * numAdd);
	cons_->objects = objL; 
	objL = NULL;
    }
    else {
	cons_->objects = NULL;
    } 
  
    //------------------------------------------------------
    // Modified col hard lb.
    //------------------------------------------------------

    cons_->lbHard.relative = relclh;
    cons_->lbHard.numModify = numclh;
    if(numclh > 0) {
	int    *clhpL = new int [numclh];
	double *clheL = new double [numclh];
	memcpy(clhpL, clhp, numclh * sizeof(int));
	memcpy(clheL, clhe, numclh * sizeof(double));
	cons_->lbHard.posModify = clhpL;  clhpL = NULL;
	cons_->lbHard.entries = clheL;    clheL = NULL;
    }
    else {
	cons_->lbHard.posModify = NULL;
	cons_->lbHard.entries = NULL;
    }
    
    //------------------------------------------------------
    // Modified col hard ub.
    //------------------------------------------------------

    cons_->ubHard.relative = relcuh;
    cons_->ubHard.numModify = numcuh;
    if (numcuh > 0) {
	int    *cuhpL = new int [numcuh];
	double *cuheL = new double [numcuh];
	memcpy(cuhpL, cuhp, numcuh * sizeof(int));
	memcpy(cuheL, cuhe, numcuh * sizeof(double));
	cons_->ubHard.posModify = cuhpL;  cuhpL = NULL;
	cons_->ubHard.entries = cuheL;    cuheL = NULL;
    }
    else {
	cons_->ubHard.posModify = NULL;
	cons_->ubHard.entries = NULL;
    }
    
    //------------------------------------------------------
    // Modified col soft lb.
    //------------------------------------------------------
    
    cons_->lbSoft.relative = relcls;
    cons_->lbSoft.numModify = numcls;
    if (numcls > 0) {
	int    *clspL = new int [numcls];
	double *clseL = new double [numcls];
	memcpy(clspL, clsp, numcls * sizeof(int));
	memcpy(clseL, clse, numcls * sizeof(double));
	cons_->lbSoft.posModify = clspL;  clspL = NULL;
	cons_->lbSoft.entries = clseL;    clseL = NULL;
    }
    else {
	cons_->lbSoft.posModify = NULL;
	cons_->lbSoft.entries = NULL;
    }
    
    //------------------------------------------------------
    // Modifed col soft ub.
    //------------------------------------------------------

    cons_->ubSoft.relative = relcus;
    cons_->ubSoft.numModify = numcus;
    if (numcus > 0) {
	int    *cuspL = new int [numcus];
	double *cuseL = new double [numcus];
	memcpy(cuspL, cusp, numcus * sizeof(int));
	memcpy(cuseL, cuse, numcus * sizeof(double));
	cons_->ubSoft.posModify = cuspL;  cuspL = NULL;
	cons_->ubSoft.entries = cuseL;    cuseL = NULL;
    }
    else {
	cons_->ubSoft.posModify = NULL;
	cons_->ubSoft.entries = NULL;
    }
}

//#############################################################################

void BcpsNodeDesc::assignCons(int numRem,
			      int    *&posRem,
			      int numAdd,
			      BcpsObject **&objects,
			      bool relclh,
			      int numclh,
			      int    *&clhp,
			      double *&clhe,
			      bool relcuh,
			      int numcuh,
			      int    *&cuhp,
			      double *&cuhe,
			      bool relcls,
			      int numcls,
			      int    *&clsp,
			      double *&clse,
			      bool relcus,
			      int numcus,
			      int    *&cusp,
			      double *&cuse)
{
    //------------------------------------------------------
    // Removed constraints.
    //------------------------------------------------------

    cons_->numRemove = numRem;
    cons_->posRemove = posRem; posRem = NULL;
    
    //------------------------------------------------------
    // Added constraints.
    //------------------------------------------------------

    cons_->numAdd = numAdd;
    cons_->objects = objects; objects = NULL;   

    //------------------------------------------------------
    // Modified col hard lb.
    //------------------------------------------------------

    cons_->lbHard.relative = relclh;
    cons_->lbHard.numModify = numclh;
    cons_->lbHard.posModify = clhp;  clhp = NULL;
    cons_->lbHard.entries = clhe;    clhe = NULL;

    //------------------------------------------------------
    // Modified col hard ub.
    //------------------------------------------------------

    cons_->ubHard.relative = relcuh;
    cons_->ubHard.numModify = numcuh;
    cons_->ubHard.posModify = cuhp;  cuhp = NULL;
    cons_->ubHard.entries = cuhe;    cuhe = NULL;

    //------------------------------------------------------
    // Modified col soft lb.
    //------------------------------------------------------

    cons_->lbSoft.relative = relcls;
    cons_->lbSoft.numModify = numcls;
    cons_->lbSoft.posModify = clsp;  clsp = NULL;
    cons_->lbSoft.entries = clse;    clse = NULL;

    //------------------------------------------------------
    // Modifed col soft ub.
    //------------------------------------------------------

    cons_->ubSoft.relative = relcus;
    cons_->ubSoft.numModify = numcus;
    cons_->ubSoft.posModify = cusp;  cusp = NULL;
    cons_->ubSoft.entries = cuse;    cuse = NULL;
}
  

//#############################################################################

void BcpsNodeDesc::setVarSoftBound(int numModSoftVarLB, 
				   const int *varLBi,
				   const double *varLBv,
				   int numModSoftVarUB,
				   const int *varUBi,
				   const double *varUBv)
{
    //------------------------------------------------------
    // Modified var soft lb.
    //------------------------------------------------------
    
    vars_->lbSoft.relative = true;
    vars_->lbSoft.numModify = numModSoftVarLB;
    if (vars_->lbSoft.posModify) delete [] vars_->lbSoft.posModify;
    if (vars_->lbSoft.entries) delete [] vars_->lbSoft.entries;
    
    if (numModSoftVarLB > 0) {
	int    *vlspL = new int [numModSoftVarLB];
	double *vlseL = new double [numModSoftVarLB];
	memcpy(vlspL, varLBi, numModSoftVarLB * sizeof(int));
	memcpy(vlseL, varLBv, numModSoftVarLB * sizeof(double));
	vars_->lbSoft.posModify = vlspL;  vlspL = NULL;
	vars_->lbSoft.entries = vlseL;    vlseL = NULL;
    }
    else {
	vars_->lbSoft.posModify = NULL;
	vars_->lbSoft.entries = NULL;
    }
    
    //------------------------------------------------------
    // Modified var soft ub.
    //------------------------------------------------------
    
    vars_->ubSoft.relative = true;
    vars_->ubSoft.numModify = numModSoftVarUB;
    if (vars_->ubSoft.posModify) delete [] vars_->ubSoft.posModify;
    if (vars_->ubSoft.entries) delete [] vars_->ubSoft.entries;
    
    if (numModSoftVarUB > 0) {
	int    *vuspL = new int [numModSoftVarUB];
	double *vuseL = new double [numModSoftVarUB];
	memcpy(vuspL, varUBi, numModSoftVarUB * sizeof(int));
	memcpy(vuseL, varUBv, numModSoftVarUB * sizeof(double));
	vars_->ubSoft.posModify = vuspL;  vuspL = NULL;
	vars_->ubSoft.entries = vuseL;    vuseL = NULL;
    }
    else {
	vars_->ubSoft.posModify = NULL;
	vars_->ubSoft.entries = NULL;
    }
}

//#############################################################################

void BcpsNodeDesc::assignVarSoftBound(int numModSoftVarLB, 
				      int *&varLBi,
				      double *&varLBv,
				      int numModSoftVarUB,
				      int *&varUBi,
				      double *&varUBv)
{
    if (vars_->lbSoft.posModify) delete [] vars_->lbSoft.posModify;
    if (vars_->lbSoft.entries) delete [] vars_->lbSoft.entries;

    vars_->lbSoft.relative = true;
    vars_->lbSoft.numModify = numModSoftVarLB;
    vars_->lbSoft.posModify = varLBi; varLBi = NULL;    
    vars_->lbSoft.entries = varLBv; varLBv = NULL;
    
    if (vars_->ubSoft.posModify) delete [] vars_->ubSoft.posModify;
    if (vars_->ubSoft.entries) delete [] vars_->ubSoft.entries;
    
    vars_->ubSoft.relative = true;
    vars_->ubSoft.numModify = numModSoftVarUB;
    vars_->ubSoft.posModify = varUBi; varUBi = NULL;
    vars_->ubSoft.entries = varUBv; varUBv = NULL;
}

//#############################################################################

void BcpsNodeDesc::setConSoftBound(int numModSoftConLB, 
				   const int *conLBi,
				   const double *conLBv,
				   int numModSoftConUB,
				   const int *conUBi,
				   const double *conUBv)
{
    //------------------------------------------------------
    // Modified con soft lb.
    //------------------------------------------------------

    cons_->lbSoft.relative = true;
    cons_->lbSoft.numModify = numModSoftConLB;
    if (cons_->lbSoft.posModify) delete [] cons_->lbSoft.posModify;
    if (cons_->lbSoft.entries) delete [] cons_->lbSoft.entries;
    
    if (numModSoftConLB > 0) {
	int    *vlspL = new int [numModSoftConLB];
	double *vlseL = new double [numModSoftConLB];
	memcpy(vlspL, conLBi, numModSoftConLB * sizeof(int));
	memcpy(vlseL, conUBv, numModSoftConLB * sizeof(double));
	cons_->lbSoft.posModify = vlspL;  vlspL = NULL;
	cons_->lbSoft.entries = vlseL;    vlseL = NULL;
    }
    else {
	cons_->lbSoft.posModify = NULL;
	cons_->lbSoft.entries = NULL;
    }
     
    //------------------------------------------------------ 
    // Modified con soft ub.
    //------------------------------------------------------
    
    cons_->ubSoft.relative = true;
    cons_->ubSoft.numModify = numModSoftConUB;
    if (cons_->ubSoft.posModify) delete [] cons_->ubSoft.posModify;
    if (cons_->ubSoft.entries) delete [] cons_->ubSoft.entries;

    if (numModSoftConUB > 0) {
	int    *vuspL = new int [numModSoftConUB];
	double *vuseL = new double [numModSoftConUB];
	memcpy(vuspL, conUBi, numModSoftConUB * sizeof(int));
	memcpy(vuseL, conUBv, numModSoftConUB * sizeof(double));
	cons_->ubSoft.posModify = vuspL;  vuspL = NULL;
	cons_->ubSoft.entries = vuseL;    vuseL = NULL;
    }
    else {
	cons_->ubSoft.posModify = NULL;
	cons_->ubSoft.entries = NULL;
    }
}

//#############################################################################

void BcpsNodeDesc::assignVarHardBound(int numModHardVarLB, 
				      int *&varLBi,
				      double *&varLBv,
				      int numModHardVarUB,
				      int *&varUBi,
				      double *&varUBv)
{
    if (vars_->lbHard.posModify) delete [] vars_->lbHard.posModify;
    if (vars_->lbHard.entries) delete [] vars_->lbHard.entries;

    vars_->lbHard.relative = true;
    vars_->lbHard.numModify = numModHardVarLB;
    vars_->lbHard.posModify = varLBi; varLBi = NULL;
    vars_->lbHard.entries = varLBv; varLBv = NULL;

    if (vars_->ubHard.posModify) delete [] vars_->ubHard.posModify;
    if (vars_->ubHard.entries) delete [] vars_->ubHard.entries;

    vars_->ubHard.relative = true;
    vars_->ubHard.numModify = numModHardVarUB;
    vars_->ubHard.posModify = varUBi; varUBi = NULL;
    vars_->ubHard.entries = varUBv; varUBv = NULL;
}

//#############################################################################

void BcpsNodeDesc::setVarHardBound(int numModHardVarLB, 
				   const int *varLBi,
				   const double *varLBv,
				   int numModHardVarUB,
				   const int *varUBi,
				   const double *varUBv)
{
    //------------------------------------------------------
    // Modified var hard lb.
    //------------------------------------------------------

    vars_->lbHard.relative = true;
    vars_->lbHard.numModify = numModHardVarLB;
    if (vars_->lbHard.posModify) delete [] vars_->lbHard.posModify;
    if (vars_->lbHard.entries) delete [] vars_->lbHard.entries;
    
    if (numModHardVarLB > 0) {
	int    *vlspL = new int [numModHardVarLB];
	double *vlseL = new double [numModHardVarLB];
	memcpy(vlspL, varLBi, numModHardVarLB * sizeof(int));
	memcpy(vlseL, varLBv, numModHardVarLB * sizeof(double));
	vars_->lbHard.posModify = vlspL;  vlspL = NULL;
	vars_->lbHard.entries = vlseL;    vlseL = NULL;
    }
    else {
	vars_->lbHard.posModify = NULL;
	vars_->lbHard.entries = NULL;
    }
  
    //------------------------------------------------------
    // Modified var hard ub.
    //------------------------------------------------------

    vars_->ubHard.relative = true;
    vars_->ubHard.numModify = numModHardVarUB;
    if (vars_->ubHard.posModify) delete [] vars_->ubHard.posModify;
    if (vars_->ubHard.entries) delete [] vars_->ubHard.entries;
    
    if (numModHardVarUB > 0) {
	int    *vuspL = new int [numModHardVarUB];
	double *vuseL = new double [numModHardVarUB];
	memcpy(vuspL, varUBi, numModHardVarUB * sizeof(int));
	memcpy(vuseL, varUBv, numModHardVarUB * sizeof(double));
	vars_->ubHard.posModify = vuspL;  vuspL = NULL;
	vars_->ubHard.entries = vuseL;    vuseL = NULL;
    }
    else {
	vars_->ubHard.posModify = NULL;
	vars_->ubHard.entries = NULL;
    }
}

//#############################################################################

void BcpsNodeDesc::setConHardBound(int numModHardConLB, 
				   const int *conLBi,
				   const double *conLBv,
				   int numModHardConUB,
				   const int *conUBi,
				   const double *conUBv)
{
    //------------------------------------------------------
    // Modified constraint hard lb.
    //------------------------------------------------------

    cons_->lbHard.relative = true;
    cons_->lbHard.numModify = numModHardConLB;
    if (cons_->lbHard.posModify) delete [] cons_->lbHard.posModify;
    if (cons_->lbHard.entries) delete [] cons_->lbHard.entries;
    
    if (numModHardConLB > 0) {
	int    *vlspL = new int [numModHardConLB];
	double *vlseL = new double [numModHardConLB];
	memcpy(vlspL, conLBi, numModHardConLB * sizeof(int));
	memcpy(vlseL, conUBv, numModHardConLB * sizeof(double));
	cons_->lbHard.posModify = vlspL;  vlspL = NULL;
	cons_->lbHard.entries = vlseL;    vlseL = NULL;
    }
    else {
	cons_->lbHard.posModify = NULL;
	cons_->lbHard.entries = NULL;
    }
    
    //------------------------------------------------------
    // Modified constraint hard ub.
    //------------------------------------------------------

    cons_->ubHard.relative = true;
    cons_->ubHard.numModify = numModHardConUB;
    if (cons_->ubHard.posModify) delete [] cons_->ubHard.posModify;
    if (cons_->ubHard.entries) delete [] cons_->ubHard.entries;
    
    if (numModHardConUB > 0) {
	int    *vuspL = new int [numModHardConUB];
	double *vuseL = new double [numModHardConUB];
	memcpy(vuspL, conUBi, numModHardConUB * sizeof(int));
	memcpy(vuseL, conUBv, numModHardConUB * sizeof(double));
	cons_->ubHard.posModify = vuspL;  vuspL = NULL;
	cons_->ubHard.entries = vuseL;    vuseL = NULL;
    }
    else {
	cons_->ubHard.posModify = NULL;
	cons_->ubHard.entries = NULL;
    }
}

//#############################################################################

AlpsReturnCode 
BcpsNodeDesc::encodeDblFieldMods(AlpsEncoded *encoded,
				 BcpsFieldListMod<double> *field) const
{
    AlpsReturnCode status = ALPS_OK;
    assert(encoded);
    
    encoded->writeRep(field->relative);
    encoded->writeRep(field->posModify, field->numModify);
    encoded->writeRep(field->entries, field->numModify);

#ifdef BCPS_DEBUG_MORE
    int k;
    for (k = 0; k < field->numModify; ++k) {
	std::cout << "BCPS encode() dbl: pos = " << field->posModify[k]
		  << ", value = " << field->entries[k] << std::endl;
    }
#endif
    
    return status;
}

//#############################################################################

AlpsReturnCode 
BcpsNodeDesc::encodeIntFieldMods(AlpsEncoded *encoded,
				 BcpsFieldListMod<int> *field) const
{
    AlpsReturnCode status = ALPS_OK;
    assert(encoded);
    
    encoded->writeRep(field->relative);
    encoded->writeRep(field->posModify, field->numModify);
    encoded->writeRep(field->entries, field->numModify);

#ifdef BCPS_DEBUG_MORE
    int k;
    for (k = 0; k < field->numModify; ++k) {
	std::cout << "BCPS encode() int: pos = " << field->posModify[k]
		  << ", value = " << field->entries[k] << std::endl;
    }
#endif

    return status;
}

//#############################################################################

AlpsReturnCode 
BcpsNodeDesc::encodeObjectMods(AlpsEncoded *encoded,
			       BcpsObjectListMod *objMod) const
{
    int k;
    AlpsReturnCode status = ALPS_OK;
    assert(encoded);
    
    // Pack removed object positions.
    encoded->writeRep(objMod->posRemove, objMod->numRemove);

    // Pack added objects.
    encoded->writeRep(objMod->numAdd);
    for (k = 0; k < objMod->numAdd; ++k) {
	//Pack a object to encoded.
	(objMod->objects)[k]->encode(encoded);
    }
    
    //std::cout << "---- BCPS encode lb hard:" << std::endl;
    status = encodeDblFieldMods(encoded, &(objMod->lbHard));
    //std::cout << "---- BCPS encode ub hard:" << std::endl;
    status = encodeDblFieldMods(encoded, &(objMod->ubHard));
    //std::cout << "---- BCPS encode lb soft:" << std::endl;
    status = encodeDblFieldMods(encoded, &(objMod->lbSoft));
    //std::cout << "---- BCPS encode ub soft:" << std::endl;
    status = encodeDblFieldMods(encoded, &(objMod->ubSoft));

    // Do know what's the use of status.
    //status = encodeIntFieldMods(encoded, &(objMod->status));


    return status;
}

//#############################################################################

/** Pack bcps node description into an encoded. */
AlpsReturnCode BcpsNodeDesc::encodeBcps(AlpsEncoded *encoded) const
{
    AlpsReturnCode status = ALPS_OK;    
    
    //std::cout << "---- BCPS encoded vars" << std::endl;
    status = encodeObjectMods(encoded, vars_);
    //std::cout << "---- BCPS encoded cons:" << std::endl;
    status = encodeObjectMods(encoded, cons_);
    
    return status;
}

//#############################################################################

AlpsReturnCode 
BcpsNodeDesc::decodeDblFieldMods(AlpsEncoded &encoded,
				 BcpsFieldListMod<double> *field)
{
    AlpsReturnCode status = ALPS_OK;
    
    encoded.readRep(field->relative);
    encoded.readRep(field->posModify, field->numModify);
    encoded.readRep(field->entries, field->numModify);

#ifdef BCPS_DEBUG_MORE
    int k;
    for (k = 0; k < field->numModify; ++k) {
	std::cout << "BCPS decode() dbl: pos = " << field->posModify[k]
		  << ", value = " << field->entries[k] << std::endl;
    }
#endif

    return status;
}

//#############################################################################

AlpsReturnCode 
BcpsNodeDesc::decodeIntFieldMods(AlpsEncoded &encoded,
				 BcpsFieldListMod<int> *field)
{
    AlpsReturnCode status = ALPS_OK;
    
    encoded.readRep(field->relative);
    encoded.readRep(field->posModify, field->numModify);
    encoded.readRep(field->entries, field->numModify);

    return status;
}

//#############################################################################

AlpsReturnCode 
BcpsNodeDesc::decodeObjectMods(AlpsEncoded &encoded,
			       BcpsObjectListMod *objMod)
{
    int k;
    AlpsReturnCode status = ALPS_OK;
    
    AlpsKnowledgeBroker *broker = model_->getKnowledgeBroker();
    assert(broker);

    // Pack removed object positions.
    encoded.readRep(objMod->posRemove, objMod->numRemove);

    // Pack added objects.
    encoded.readRep(objMod->numAdd);

    if (objMod->numAdd > 0) {
	objMod->objects = new BcpsObject* [objMod->numAdd];
	for (k = 0; k < objMod->numAdd; ++k) {
	    objMod->objects[k] = static_cast<BcpsObject *>
		( broker->decoderObject(BCPS_CONSTRAINT)->decode(encoded) );
	    
	    // Unpack a object from an encoded.
	    // (objMod->objects)[k]->encode(encoded);
	}
    }
    
    //std::cout << "---- BCPS decode lb hard:" << std::endl;
    status = decodeDblFieldMods(encoded, &(objMod->lbHard));
    //std::cout << "---- BCPS decode ub hard:" << std::endl;
    status = decodeDblFieldMods(encoded, &(objMod->ubHard));
    //std::cout << "---- BCPS decode lb soft:" << std::endl;
    status = decodeDblFieldMods(encoded, &(objMod->lbSoft));
    //std::cout << "---- BCPS decode ub soft:" << std::endl;
    status = decodeDblFieldMods(encoded, &(objMod->ubSoft));

    // Do know what's the use of status.
    //status = encodeIntFieldMods(encoded, &(objMod->status));

    return status;
}

//#############################################################################

/** Unpack bcps node description from an encoded. */
AlpsReturnCode BcpsNodeDesc::decodeBcps(AlpsEncoded &encoded)
{
    AlpsReturnCode status = ALPS_OK;    
    
    status = decodeObjectMods(encoded, vars_);
    status = decodeObjectMods(encoded, cons_);
    
    return status;
}

//#############################################################################
