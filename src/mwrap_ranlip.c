#include "maplec.h"
#include "ranlipproc.h"

#ifndef bool
#define bool int
#endif
#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef NULL
#define NULL 0
#endif

#define DllExport   __declspec( dllexport )

typedef struct _sParams
{
    int    m_Dim;
    double *m_pLeft;
    double *m_pRight;
    int    m_Num;
    int    m_NumFine;
    double m_LC;
    double m_MinLC;
    char   *m_pFileName;
    int    m_Seed;
    double *m_pOutDoubleArr;
    int     m_NumOfVecs;
//distribution callback
    MKernelVector m_MapleKernelVec;
    ALGEB   m_CallBackFunc;
}sParams;

static sParams gParams;

//            MaplePrintf(kv, "QQQQQQQQQQQQ1111\n)");

bool parseDim(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_Dim = MapleToInteger32(kv,(ALGEB) args[paramNo]);
    return bRes;
}

bool parseLeft(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pLeft = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseRight(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pRight = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseNum(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_Num = MapleToInteger32(kv,(ALGEB) args[paramNo]);
    return bRes;
}


bool parseNumFine(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_NumFine = MapleToInteger32(kv,(ALGEB) args[paramNo]);
    return bRes;
}

bool parseLC(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_LC = MapleToFloat64( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseMinLC(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_MinLC = MapleToFloat64( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseFileName(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pFileName = MapleToString( kv, (ALGEB) args[paramNo]);
    return bRes;
}

bool parseSeed(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_Seed = MapleToInteger32(kv,(ALGEB) args[paramNo]);
    return bRes;
}

bool parseOutDoubleArr(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_pOutDoubleArr = (FLOAT64 *) RTableDataBlock( kv, (ALGEB) args[paramNo] );
    return bRes;
}

bool parseNumOfVecs(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{
    bool bRes = true;
    pParams->m_NumOfVecs = MapleToInteger32(kv,(ALGEB) args[paramNo]);
    return bRes;
}

bool parseDistFuncRef(MKernelVector kv,  ALGEB args, sParams *pParams, int paramNo)
{ 
    bool bRes = true;
    pParams->m_MapleKernelVec = kv;
	if( IsMapleProcedure(kv,(ALGEB)args[paramNo]) )
        pParams->m_CallBackFunc = (ALGEB)args[paramNo];
    else 
    {
        MapleRaiseError(kv,"procedure expected");
        bRes = false;
    }
    return bRes;
}

double distributionFunctionRef(double* p, int dim)
{ 
    return EvalhfMapleProc(gParams.m_MapleKernelVec,(ALGEB)gParams.m_CallBackFunc, dim, (p-1));
}

DllExport ALGEB MWRAP_PrepareHatFunction(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 7)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLeft(kv, args, pParams, paramNo++) && bRes;
        bRes = parseRight(kv, args, pParams, paramNo++) && bRes;
        bRes = parseDistFuncRef(kv, args, pParams, paramNo++) && bRes;
        bRes = parseNum(kv, args, pParams, paramNo++) && bRes;
        bRes = parseNumFine(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
        if (bRes)
        {
            InitRanLip(gParams.m_Dim, gParams.m_pLeft, gParams.m_pRight);
            SetDistFunctionRanLip(distributionFunctionRef);
            PrepareHatFunctionRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_LC);
			SeedRanLip(1);
        }
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 7. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}


DllExport ALGEB MWRAP_InitRanLip(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 3)
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLeft(kv, args, pParams, paramNo++) && bRes;
        bRes = parseRight(kv, args, pParams, paramNo++) && bRes;
        if (bRes)
            InitRanLip(gParams.m_Dim, gParams.m_pLeft, gParams.m_pRight);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 3. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}


DllExport ALGEB MWRAP_PrepareHatFunctionRanLip(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 3)
    {
        bRes = parseNum(kv, args, pParams, paramNo++) && bRes;
        bRes = parseNumFine(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLC(kv, args, pParams, paramNo++) && bRes;
		if (bRes) {
            PrepareHatFunctionRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_LC);
			SeedRanLip(1);
		}
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 3. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_PrepareHatFunctionAutoRanLip(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 3)
    {
        bRes = parseNum(kv, args, pParams, paramNo++) && bRes;
        bRes = parseNumFine(kv, args, pParams, paramNo++) && bRes;
        bRes = parseMinLC(kv, args, pParams, paramNo++) && bRes;
		if (bRes) {
            PrepareHatFunctionAutoRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_MinLC);
			SeedRanLip(1);
		}
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 3. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}


DllExport ALGEB MWRAP_PrepareHatFunctionAuto(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 7 || nargs == 6 )
    {
        bRes = parseDim(kv, args, pParams, paramNo++) && bRes;
        bRes = parseLeft(kv, args, pParams, paramNo++) && bRes;
        bRes = parseRight(kv, args, pParams, paramNo++) && bRes;
        bRes = parseDistFuncRef(kv, args, pParams, paramNo++) && bRes;
        bRes = parseNum(kv, args, pParams, paramNo++) && bRes;
        bRes = parseNumFine(kv, args, pParams, paramNo++) && bRes;
		if(nargs == 7 )
			bRes = parseMinLC(kv, args, pParams, paramNo++) && bRes;	else		
			gParams.m_MinLC=0;
        if (bRes)
        {
            InitRanLip(gParams.m_Dim, gParams.m_pLeft, gParams.m_pRight);
            SetDistFunctionRanLip(distributionFunctionRef);
            PrepareHatFunctionAutoRanLip(gParams.m_Num, gParams.m_NumFine, gParams.m_MinLC);
			SeedRanLip(1);
        }
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 6 or 7. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_RandomVecRanLip(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    int i = 0;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 1 || nargs == 2)
    {
        bRes = parseOutDoubleArr(kv, args, pParams, paramNo++) && bRes;
        pParams->m_NumOfVecs = 1;
        if (nargs == 2)
            bRes = parseNumOfVecs(kv, args, pParams, paramNo++) && bRes;
        if (bRes)
        {
            for (i = 0; i < pParams->m_NumOfVecs; i++)
                RandomVecRanLip(pParams->m_pOutDoubleArr + i * pParams->m_Dim);
        }
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 1 or 2. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_FreeMemRanLip(MKernelVector kv,  ALGEB args)
{
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        FreeMemRanLip();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 0. \n");
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_GetSeedRanLip(MKernelVector kv,  ALGEB args)
{
    int res = 0;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        res = GetSeedRanLip();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 0. \n");
    return ToMapleInteger(kv, res);
}


DllExport ALGEB MWRAP_LipschitzRanLip(MKernelVector kv,  ALGEB args)
{
    double res = 0;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        res = LipschitzRanLip();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 0. \n");
    return ToMapleFloat(kv, res);
}

DllExport ALGEB MWRAP_Count_totalRanLip(MKernelVector kv,  ALGEB args)
{
    int res = 0;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        res = Count_totalRanLip();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 0. \n");
    return ToMapleInteger(kv, res);
}

DllExport ALGEB MWRAP_Count_errorRanLip(MKernelVector kv,  ALGEB args)
{
    int res = 0;
    if (MapleNumArgs(kv,(ALGEB)args) == 0)
        res = Count_errorRanLip();
    else 
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 0. \n");
    return ToMapleInteger(kv, res);
}


bool parseSaveLoadParams(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 1)
    {
        bRes = parseFileName(kv, args, pParams, paramNo++) && bRes;;
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 1. \n");
        bRes = false;
    }
    return bRes;
}

DllExport ALGEB MWRAP_SavePartitionRanlip(MKernelVector kv,  ALGEB args)
{
    sParams *pParams = &gParams;
    if (parseSaveLoadParams(kv, args))
    {
        switch (SavePartitionRanlip(gParams.m_pFileName))
        {
            case 0: 
                MaplePrintf(kv, "Hat function successfully saved");
                break;
            case 1: 
                MaplePrintf(kv, "Could not save hat function");
                break;
            case 2: 
                MaplePrintf(kv, "Hat function has not been created yet");
                break;
            default:
                MaplePrintf(kv, "Unknown error");
                break;
        }
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_LoadPartitionRanLip(MKernelVector kv,  ALGEB args)
{
    sParams *pParams = &gParams;
    if (parseSaveLoadParams(kv, args))
    {
        switch (LoadPartitionRanLip(gParams.m_pFileName))
        {
            case 0: 
                MaplePrintf(kv, "Hat function  loaded");
				pParams->m_Dim=GetDim();
				SeedRanLip(1);
                break;
            case 1: 
                MaplePrintf(kv, "Could not load hat function");
                break;
            case 2: 
                MaplePrintf(kv, "Not enough memory");
                break;
            default:
                MaplePrintf(kv, "Unknown error");
                break;

        };
    }
    return ToMapleNULL(kv);
}


DllExport ALGEB MWRAP_SeedRanLip(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 1)
    {
        bRes = parseSeed(kv, args, pParams, paramNo++) && bRes;
        if (bRes)
            SeedRanLip(gParams.m_Seed);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 1. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}

DllExport ALGEB MWRAP_SetDistFunctionRanLip(MKernelVector kv,  ALGEB args)
{
    bool bRes = true;
    double res = 0;
    int paramNo = 1;
    sParams *pParams = &gParams;
    int nargs =  MapleNumArgs(kv,(ALGEB)args);
    if (nargs == 1)
    {
        bRes = parseDistFuncRef(kv, args, pParams, paramNo++) && bRes;
        if (bRes)
            SetDistFunctionRanLip(/*MyDist*/ distributionFunctionRef);
    }
    else 
    {
        MapleRaiseError(kv,"Error: wrong number of parameters . Must be 1. \n");
        bRes = false;
    }
    return ToMapleNULL(kv);
}
