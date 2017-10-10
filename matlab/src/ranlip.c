#include <time.h>
#include <string.h>

#include "mex.h"
#include "ranlipproc.h"

typedef struct _sParams
{
    char   m_pName[100];
    int    m_Dim;
    double *m_pLeft;
    double *m_pRight;
    int    m_Num;
    int    m_NumFine;
    double m_LC;
    double m_MinLC;
    char   m_pFileName[100];
    int    m_Seed;
    int     m_NumOfVecs;
//distribution callback
    mxArray *m_prhs_cb[3];
    mxArray *m_plhs_cb[1];
    double  *m_pDistParamArr;
// array dimensions
	int M,N;

}sParams;

static sParams gParams;

void debugParams(sParams *pParams)
{
    if (pParams->m_pName)
        mexPrintf("name=%s \n", pParams->m_pName);
    mexPrintf("dim=%d \n", pParams->m_Dim);
    mexPrintf("n=%d \n", pParams->m_Num);
}

bool parseName(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    int buflen = 0;
    bool bRes = false;
    if (nrhs)
    {
        bRes = mxIsChar(prhs[paramNo]);
        if (bRes)
        {
            buflen = (mxGetM(prhs[paramNo]) * mxGetN(prhs[paramNo])) + 1;
            memset(pParams->m_pName, 0, buflen);
            if (mxGetString(prhs[paramNo], pParams->m_pName, buflen))
            {
                mexPrintf("error extracting function name\n", paramNo);
                bRes = false;
            }
            else
                _strlwr(pParams->m_pName);
        }
        else 
            mexPrintf("%d input must be function name\n", paramNo);
    }
    else 
        mexPrintf("error: wrong number of parameters (%d). must be 1 or more \n", nrhs);
    return bRes;
}

bool parseDim(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) !=1 )
    {
        mexPrintf("%d input must be an integer (dimension)\n", paramNo);
        bRes = false;
    }
    else
        pParams->m_Dim = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseLeft(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_Dim ||mxGetN(prhs[paramNo])==1)  || !( mxGetM(prhs[paramNo])==1 || mxGetM(prhs[paramNo])==pParams->m_Dim  ))
    {
        mexPrintf("%d input must be a vector of doubles dimension size 1 x %d (left boundary)\n", paramNo, pParams->m_Dim);
        bRes = false;
    }
    else {
        pParams->m_pLeft = mxGetPr(prhs[paramNo]);
	  pParams->N=mxGetN(prhs[paramNo]);
	  pParams->M=mxGetM(prhs[paramNo]);
    }
    return bRes;
}


bool parseRight(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]) || !(mxGetN(prhs[paramNo])== pParams->m_Dim ||mxGetN(prhs[paramNo])==1)  || !( mxGetM(prhs[paramNo])==1 || mxGetM(prhs[paramNo])==pParams->m_Dim  ))
//    if (!mxIsDouble(prhs[paramNo]) || mxGetN(prhs[paramNo])!= pParams->m_Dim ||  mxGetM(prhs[paramNo])!=1)
    {
        mexPrintf("%d input must be a vector of doubles dimension size 1 x %d (right boundary)\n", paramNo, pParams->m_Dim);
        bRes = false;
    }
    else
        pParams->m_pRight = mxGetPr(prhs[paramNo]);
    return bRes;
}

bool parseNum(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) != 1)
    {
        mexPrintf("%d dst input must be an integer (number of cells the domain is subdivided in each direction)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_Num = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}


bool parseNumFine(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) != 1)
    {
        mexPrintf("%d dst input must be an integer (number of cells in the fine partition. should be a power of 2to speed up calculations)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_NumFine = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseLC(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]))
    {
        mexPrintf("%d input must be a double (Lipschitz constant)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_LC = mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseMinLC(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if (!mxIsDouble(prhs[paramNo]))
    {
        mexPrintf("%d input must be a double (min Lipschitz constant)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_MinLC = mxGetScalar(prhs[paramNo]);
    return bRes;
}


bool parseFileName(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    int buflen = 0;
    bool bRes = false;
    bRes = mxIsChar(prhs[paramNo]);
    if (bRes)
    {
        buflen = (mxGetM(prhs[paramNo]) * mxGetN(prhs[paramNo])) + 1;
        /* Allocate memory for input and output strings. */
        memset(pParams->m_pFileName, 0, buflen);
        if (mxGetString(prhs[paramNo], pParams->m_pFileName, buflen))
        {
            mexPrintf("error extracting file name\n", paramNo);
            bRes = false;
        }
    }
    else 
        mexPrintf("%d input must be file name\n", paramNo);
    return bRes;
}

bool parseSeed(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) !=1 )
    {
        mexPrintf("%d input must be an integer (seed)\n", paramNo);
        bRes = false;
    }
    else
        pParams->m_Seed = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

bool parseNumOfVecs(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    bool bRes = true;
    if ((!mxIsDouble(prhs[paramNo])) || mxGetN(prhs[paramNo])!= 1 ||  mxGetM(prhs[paramNo]) != 1)
    {
        mexPrintf("%d dst input must be an integer (number of random vectors)\n", paramNo);
        bRes = false;
    }
    else 
        pParams->m_NumOfVecs = (int)mxGetScalar(prhs[paramNo]);
    return bRes;
}

void freeDistMem(void)
{

    if (gParams.m_prhs_cb[0])
        mxDestroyArray(gParams.m_prhs_cb[0]);

    if (gParams.m_prhs_cb[1])
        mxDestroyArray(gParams.m_prhs_cb[1]);

    if (gParams.m_prhs_cb[2])
        mxDestroyArray(gParams.m_prhs_cb[2]);

    gParams.m_prhs_cb[0] = gParams.m_prhs_cb[1] = gParams.m_prhs_cb[2] = NULL;
}    

bool parseDistFuncRef(int nrhs, const mxArray *prhs[], sParams *pParams, int paramNo)
{
    int buflen = 0;
    double *pTmp;
    bool bRes = true;
    freeDistMem();

    gParams.m_prhs_cb[0] = mxDuplicateArray(prhs[paramNo]);
    mexMakeArrayPersistent(gParams.m_prhs_cb[0]);

    gParams.m_prhs_cb[1] = mxCreateDoubleMatrix(gParams.M, gParams.N, mxREAL);
    gParams.m_pDistParamArr = mxGetPr(gParams.m_prhs_cb[1]);
    mexMakeArrayPersistent(gParams.m_prhs_cb[1]);
    
    gParams.m_prhs_cb[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    pTmp = mxGetPr(gParams.m_prhs_cb[2]);
    pTmp[0] = gParams.m_Dim;
    mexMakeArrayPersistent(gParams.m_prhs_cb[2]);

    return bRes;
}

double distributionFunctionRef(double* p, int dim)
{ // example: multivariate normal distribution
    
    double res;
    memcpy(gParams.m_pDistParamArr, p, dim * sizeof(double));
    if (mexCallMATLAB(1 , gParams.m_plhs_cb, 3, gParams.m_prhs_cb, "feval")==0)
    {
        res = mxGetScalar(gParams.m_plhs_cb[0]);
    }
    else 
        mexPrintf("callback error \n");
    mxDestroyArray(gParams.m_plhs_cb[0]);
    return res;
}



bool parseParams(int nrhs, const mxArray *prhs[], sParams *pParams)
{
    bool bRes = true;
    int paramNo = 0;
    bRes = parseName(nrhs, prhs, pParams, paramNo);
    paramNo++;
    if (bRes)
    {
        if (strcmp(pParams->m_pName, "preparehatfunction") == 0)
        {       
            if (nrhs == 8)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLeft(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRight(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseDistFuncRef(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNum(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNumFine(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 8.\n", nrhs);
                bRes = false;
            }
        }
        else if (strcmp(pParams->m_pName, "init") == 0)
        {       
            if (nrhs == 4)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLeft(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRight(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 4.\n", nrhs);
                bRes = false;
            }
        }
        else if (strcmp(pParams->m_pName, "preparehatfunction0") == 0)
        {       
            if (nrhs == 4)
            {
                bRes = parseNum(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNumFine(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLC(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 4.\n", nrhs);
                bRes = false;
            }
        }
        else if (strcmp(pParams->m_pName, "preparehatfunctionauto") == 0)
        {       
            if (nrhs == 8)
            {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLeft(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRight(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseDistFuncRef(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNum(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNumFine(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseMinLC(nrhs, prhs, pParams, paramNo++) && bRes;
			} else if(nrhs==7) {
                bRes = parseDim(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseLeft(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseRight(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseDistFuncRef(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNum(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseNumFine(nrhs, prhs, pParams, paramNo++) && bRes;
				pParams->m_MinLC=0;
               
			}
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 7 or 8.\n", nrhs);
                bRes = false;
            }
        } 
        else if (strcmp(pParams->m_pName, "randomvec") == 0)
        {       
            if (nrhs == 1 || nrhs == 2)
            {
                pParams->m_NumOfVecs = 1;
                if (nrhs == 2)
                    bRes = parseNumOfVecs(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 1 or 2.\n", nrhs);
                bRes = false;
            }
        } 

        else if (strcmp(pParams->m_pName, "freemem") == 0 ||
                 strcmp(pParams->m_pName, "getseed") == 0 || strcmp(pParams->m_pName, "lipschitz") == 0 ||
                 strcmp(pParams->m_pName, "count_total") == 0 || strcmp(pParams->m_pName, "count_error") == 0)
        {       
            if (nrhs == 1)
            {
                bRes = true;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 1.\n", nrhs);
                bRes = false;
            }
        }
        else if (strcmp(pParams->m_pName, "savepartition") == 0 || strcmp(pParams->m_pName, "loadpartition") == 0)
        {       
            if (nrhs == 2)
            {
                bRes = bRes = parseFileName(nrhs, prhs, pParams, paramNo++) && bRes;;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 2.\n", nrhs);
                bRes = false;
            }
        }
        else if (strcmp(pParams->m_pName, "seed") == 0)
        {       
            if (nrhs == 2)
            {
                bRes = parseSeed(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 2.\n", nrhs);
                bRes = false;
            }
        }
        else if (strcmp(pParams->m_pName, "setdistfunction") == 0)
        {       
            if (nrhs == 2)
            {
                //bRes = parseDistFuncName(nrhs, prhs, pParams, paramNo++) && bRes;
                bRes = parseDistFuncRef(nrhs, prhs, pParams, paramNo++) && bRes;
            }
            else 
            {
                mexPrintf("error: wrong number of parameters (%d). must be 2.\n", nrhs);
                bRes = false;
            }
        }
        else
        {
            mexPrintf("error: unknown function name '%s'. \n", pParams->m_pName);
            bRes = false;
        }
    }
    return bRes;
}

void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] )
{
    int i,j;
    bool bRes = true;
    double *pRes = NULL;
    double* temp;   
    mexAtExit(freeDistMem);

    //mexLock();
    bRes = parseParams(nrhs, prhs, &gParams);
    if (bRes)
    {
        if (strcmp(gParams.m_pName, "preparehatfunction")==0)
        {
            InitRanLip(gParams.m_Dim, gParams.m_pLeft, gParams.m_pRight);
            SetDistFunctionRanLip(distributionFunctionRef);
            PrepareHatFunctionRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_LC);
			SeedRanLip(1);
        }
        else if (strcmp(gParams.m_pName, "init")==0)
        {
            InitRanLip(gParams.m_Dim, gParams.m_pLeft, gParams.m_pRight);
        }
        else if (strcmp(gParams.m_pName, "preparehatfunction0")==0)
        {
            PrepareHatFunctionRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_LC);
			SeedRanLip(1);
        }
        else if (strcmp(gParams.m_pName, "preparehatfunctionauto")==0)
        {
            InitRanLip(gParams.m_Dim, gParams.m_pLeft, gParams.m_pRight);
            SetDistFunctionRanLip(distributionFunctionRef);
            PrepareHatFunctionAutoRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_MinLC);
 			SeedRanLip(1);
       } 
        else if (strcmp(gParams.m_pName, "randomvec")==0)
        {
		if(gParams.M==1) {
            plhs[0] = mxCreateDoubleMatrix(gParams.m_NumOfVecs, gParams.m_Dim, mxREAL);
            pRes = mxGetPr(plhs[0]);
			if(gParams.m_NumOfVecs==1) RandomVecRanLip(pRes ); else
			{
				temp=(double*) malloc(sizeof(double)* gParams.m_Dim);
				for ( i = 0; i < gParams.m_NumOfVecs; i++)  //returns matrix in C style order
				{
					RandomVecRanLip(temp);
					for(j=0;j<gParams.m_Dim;j++)
						pRes[i + j*  gParams.m_NumOfVecs] = temp[j];  // transposed
				}
				free(temp);			
			}
		} else {
            	plhs[0] = mxCreateDoubleMatrix(gParams.m_Dim, gParams.m_NumOfVecs,  mxREAL);
            	pRes = mxGetPr(plhs[0]);
			for ( i = 0; i < gParams.m_NumOfVecs; i++)  //returns matrix in F style order
						RandomVecRanLip(pRes + i*gParams.m_Dim);	
			
		      }
        } 
        else if (strcmp(gParams.m_pName, "savepartition")==0)
        {
            switch (SavePartitionRanlip(gParams.m_pFileName))
            {
            case 0: 
                mexPrintf("hat function successfully saved");
                break;
            case 1: 
                mexPrintf("could not save hat function");
                break;
            case 2: 
                mexPrintf("hat function has not been created yet");
                break;
            default:
                mexPrintf("unknown error");
                break;
            }
        } 
        else if (strcmp(gParams.m_pName, "loadpartition")==0)
        {
            switch (LoadPartitionRanLip(gParams.m_pFileName))
            {
            case 0: 
                mexPrintf("hat function  loaded");
				gParams.m_Dim=GetDim();
				SeedRanLip(1);
                break;
            case 1: 
                mexPrintf("could not load hat function");
                break;
            case 2: 
                mexPrintf("not enough memory");
                break;
            default:
                mexPrintf("unknown error");
                break;
            }
        } 
        else if (strcmp(gParams.m_pName, "freemem")==0)
        {
            freeDistMem();
            FreeMemRanLip();
        } 
        else if (strcmp(gParams.m_pName, "seed")==0)
        {
            SeedRanLip(gParams.m_Seed);
        } 
        else if (strcmp(gParams.m_pName, "getseed")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = GetSeedRanLip();
        } 
        else if (strcmp(gParams.m_pName, "lipschitz")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = LipschitzRanLip();
        } 
        else if (strcmp(gParams.m_pName, "count_total")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = Count_totalRanLip();
        } 
        else if (strcmp(gParams.m_pName, "count_error")==0)
        {
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            pRes = mxGetPr(plhs[0]);
            *pRes = Count_errorRanLip();
        } 
        else if (strcmp(gParams.m_pName, "setdistfunction")==0)
        {
            SetDistFunctionRanLip(/*MyDist*/ distributionFunctionRef);
            //PrepareHatFunctionAutoRanLip(50,16,0);
        } 

    }

    if (!pRes && nlhs)
    {
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        pRes = mxGetPr(plhs[0]);
        *pRes = 0;
    }

}
