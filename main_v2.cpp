/******************************************************************************
 * Project:  gisconf 2015 demo
 * Purpose:  gdal warp using demo
 * Author:   Dmitry Baryshnikov, dmitry.baryshnikov@nextigs.com
 ******************************************************************************
*   Copyright (C) 2015 NextGIS, info@nextgis.com
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include "cpl_string.h"
#include "gdalwarper.h"

#include <vector>

#define SEGMENT_STEPS 8
#define MULTI
#define APRROX_MAXERROR 0.125
#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
    do { if (iArg + nExtraArg >= nArgc) \
        Usage(CPLSPrintf("%s option requires %d argument(s)", papszArgv[iArg], nExtraArg)); } while(0)

static void Usage(const char* pszErrorMsg = NULL)

{
    printf( "Usage: warper [-c lat long] [-nw lat long] [-ne lat long]\n"
            "              [-sw lat long] [-se lat long] filename");

    if( pszErrorMsg != NULL )
        fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);

    exit( 1 );
}

void SaveGeometry(const CPLString &path, const OGRPolygon &polygon, const OGRSpatialReference &spaRef)
{
    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS = poDriver->Create( path, 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL )
    {
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    const char* pszLayerName = CPLGetBasename(path);

    OGRLayer *poLayer = poDS->CreateLayer( pszLayerName, spaRef.Clone(), wkbPolygon, NULL );
    if( poLayer == NULL )
    {
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    OGRFieldDefn oField( "Name", OFTString );
    oField.SetWidth(32);
    if( poLayer->CreateField( &oField ) != OGRERR_NONE )
    {
        printf( "Creating Name field failed.\n" );
        exit( 1 );
    }

    OGRFeature *poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
    //poFeature->SetField( "Name", szName );

    poFeature->SetGeometry( polygon.clone() );
    if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
    {
        printf( "Failed to create feature in shapefile.\n" );
        exit( 1 );
    }
    OGRFeature::DestroyFeature( poFeature );
    GDALClose( poDS );
}

int main( int nArgc, char ** papszArgv )
{
    // register drivers
    GDALAllRegister();

    if( nArgc < 2 )
        return EXIT_FAILURE;

    double dfaCornersX[5] = {0};
    double dfaCornersY[5] = {0};
    CPLString sFileName;

    // parse input values
    for( int iArg = 1; iArg < nArgc; iArg++ )
    {
        if( EQUAL(papszArgv[iArg],"-nw"))
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
            const char* pszCoord = papszArgv[++iArg];
            dfaCornersY[1] = CPLAtofM(pszCoord);
            pszCoord = papszArgv[++iArg];
            dfaCornersX[1] = CPLAtofM(pszCoord);
        }
        else if( EQUAL(papszArgv[iArg],"-ne"))
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
            const char* pszCoord = papszArgv[++iArg];
            dfaCornersY[2] = CPLAtofM(pszCoord);
            pszCoord = papszArgv[++iArg];
            dfaCornersX[2] = CPLAtofM(pszCoord);
        }
        else if( EQUAL(papszArgv[iArg],"-se"))
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
            const char* pszCoord = papszArgv[++iArg];
            dfaCornersY[3] = CPLAtofM(pszCoord);
            pszCoord = papszArgv[++iArg];
            dfaCornersX[3] = CPLAtofM(pszCoord);
        }
        else if( EQUAL(papszArgv[iArg],"-sw"))
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
            const char* pszCoord = papszArgv[++iArg];
            dfaCornersY[4] = CPLAtofM(pszCoord);
            pszCoord = papszArgv[++iArg];
            dfaCornersX[4] = CPLAtofM(pszCoord);
        }
        else if( EQUAL(papszArgv[iArg],"-c"))
        {
            CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(2);
            const char* pszCoord = papszArgv[++iArg];
            dfaCornersY[0] = CPLAtofM(pszCoord);
            pszCoord = papszArgv[++iArg];
            dfaCornersX[0] = CPLAtofM(pszCoord);
        }
        else if(sFileName.empty())
            sFileName = papszArgv[iArg];
    }



    // set default values in cm
    // KH-4(A,B) http://fas.org/irp/imint/docs/kh-4_camera_system.htm
    double dfFocus = 60.9602;
    double dfFilmHeight = 5.45338; // 2.147 in
    double dfMeanHeight = 156;

    // to meters
    double dfFocusM = dfFocus / 100;
    double dfFilmHalfWidth = dfFilmHeight / 200;

    OGRSpatialReference oOGRSpatialReference(SRS_WKT_WGS84);
    int nZoneNo = ceil( (180.0 + dfaCornersX[0]) / 6.0 );
    OGRSpatialReference oDstOGRSpatialReference(SRS_WKT_WGS84);
    oDstOGRSpatialReference.SetUTM(nZoneNo, dfaCornersY[0] > 0);

    // transform coordinates from WGS84 to UTM
    OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation( &oOGRSpatialReference, &oDstOGRSpatialReference);
    if(!poCT)
    {
        Usage("get coordinate transformation failed");
        return EXIT_FAILURE;
    }
        
    int nResult = poCT->Transform(5, dfaCornersX, dfaCornersY, NULL);
    if(!nResult)
    {
        Usage("transformation failed");
        return EXIT_FAILURE;
    }
    
    
    OGRPoint ptCenter(dfaCornersX[0], dfaCornersY[0]);
    OGRPoint pt1(dfaCornersX[1], dfaCornersY[1]); // NW Cormer
    OGRPoint pt2(dfaCornersX[2], dfaCornersY[2]); // NE Corner
    OGRPoint pt3(dfaCornersX[3], dfaCornersY[3]); // SE Corner
    OGRPoint pt4(dfaCornersX[4], dfaCornersY[4]); // SW Corner

    // check width
    double dfDist1 = pt1.Distance(&pt2);
    double dfDist2 = pt2.Distance(&pt3);
    double dfDist3 = pt3.Distance(&pt4);
    double dfDist4 = pt4.Distance(&pt1);

    // open dataset
    GDALDataset *poSrcDataset = (GDALDataset *) GDALOpen( sFileName, GA_ReadOnly ); // GA_Update
    char* pszSpaRefDef = NULL;
    if( oDstOGRSpatialReference.exportToWkt(&pszSpaRefDef) != OGRERR_NONE)
    {
        CPLFree( pszSpaRefDef );
        GDALClose( (GDALDatasetH) poSrcDataset );
        return EXIT_FAILURE;
    }

    double dfWidth = (dfDist2 + dfDist4) / 4;
    double dfLen = (dfDist1 + dfDist3) / 4;
    OGRPoint ptBeg, ptEnd;

    OGRLineString SubCenterLine1;
    SubCenterLine1.addPoint(&pt1);
    SubCenterLine1.addPoint(&pt4);
    SubCenterLine1.Centroid(&ptBeg);

    OGRLineString SubCenterLine2;
    SubCenterLine2.addPoint(&pt2);
    SubCenterLine2.addPoint(&pt3);
    SubCenterLine2.Centroid(&ptEnd);

    double dfHeightHyp = (dfWidth * dfFocusM) / dfFilmHalfWidth;
    //satellite height
    double dfHeight = sqrt( dfHeightHyp * dfHeightHyp - dfLen * dfLen);

    int nStepCount = SEGMENT_STEPS / 2;

    std::vector<OGRPoint> aPt1, aPt2, aPt3, aPt4;

    GDAL_GCP *paGSPs = (GDAL_GCP *) CPLMalloc ((SEGMENT_STEPS * 2 + 6) * sizeof(GDAL_GCP));
    GDALInitGCPs(SEGMENT_STEPS * 2 + 6, paGSPs);

    // add image corners
    int nGCPPos = 0;
    paGSPs[nGCPPos].dfGCPLine = 0;
    paGSPs[nGCPPos].dfGCPPixel = 0;
    paGSPs[nGCPPos].dfGCPX = pt1.getX();
    paGSPs[nGCPPos].dfGCPY = pt1.getY();
    paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
    paGSPs[nGCPPos].pszId = "nw";
    nGCPPos++;

    paGSPs[nGCPPos].dfGCPLine = 0;
    paGSPs[nGCPPos].dfGCPPixel = poSrcDataset->GetRasterXSize();
    paGSPs[nGCPPos].dfGCPX = pt2.getX();
    paGSPs[nGCPPos].dfGCPY = pt2.getY();
    paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
    paGSPs[nGCPPos].pszId = "ne";
    nGCPPos++;

    paGSPs[nGCPPos].dfGCPLine = poSrcDataset->GetRasterYSize();
    paGSPs[nGCPPos].dfGCPPixel = poSrcDataset->GetRasterXSize();
    paGSPs[nGCPPos].dfGCPX = pt3.getX();
    paGSPs[nGCPPos].dfGCPY = pt3.getY();
    paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
    paGSPs[nGCPPos].pszId = "se";
    nGCPPos++;

    paGSPs[nGCPPos].dfGCPLine = poSrcDataset->GetRasterYSize();
    paGSPs[nGCPPos].dfGCPPixel = 0;
    paGSPs[nGCPPos].dfGCPX = pt4.getX();
    paGSPs[nGCPPos].dfGCPY = pt4.getY();
    paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
    paGSPs[nGCPPos].pszId = "sw";
    nGCPPos++;

    //proceed left side of frame
    OGRLineString CenterLineLeft;
    CenterLineLeft.addPoint(&ptCenter);
    CenterLineLeft.addPoint(&ptBeg);
    double dfStepLen = CenterLineLeft.get_Length() / nStepCount;
    double dfImageStepLen = double(poSrcDataset->GetRasterXSize()) / SEGMENT_STEPS;
    double dfImageCenterX = double(poSrcDataset->GetRasterXSize()) / 2;

    //add center pt
    double dfOffsetC = dfHeight * dfFilmHalfWidth / dfFocusM;
        
    double dxC = ptCenter.getX() - ptBeg.getX();
    double dyC = ptCenter.getY() - ptBeg.getY();

    double uxC = -dfOffsetC * dxC / CenterLineLeft.get_Length();
    double uyC = -dfOffsetC * dyC / CenterLineLeft.get_Length();

    OGRPoint ptUp(ptCenter.getX() - uyC, ptCenter.getY() + uxC);
    aPt2.push_back(ptUp);

    paGSPs[nGCPPos].dfGCPLine = poSrcDataset->GetRasterYSize();
    paGSPs[nGCPPos].dfGCPPixel = dfImageCenterX;
    paGSPs[nGCPPos].dfGCPX = ptUp.getX();
    paGSPs[nGCPPos].dfGCPY = ptUp.getY();
    paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
    paGSPs[nGCPPos].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPPos));
    nGCPPos++;

    uxC = dfOffsetC * dxC / CenterLineLeft.get_Length();
    uyC = dfOffsetC * dyC / CenterLineLeft.get_Length();

    OGRPoint ptDown(ptCenter.getX() - uyC, ptCenter.getY() + uxC);
    aPt1.push_back(ptDown);

    paGSPs[nGCPPos].dfGCPLine = 0;
    paGSPs[nGCPPos].dfGCPPixel = dfImageCenterX;
    paGSPs[nGCPPos].dfGCPX = ptDown.getX();
    paGSPs[nGCPPos].dfGCPY = ptDown.getY();
    paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
    paGSPs[nGCPPos].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPPos));
    nGCPPos++;

    // generate points using defined step
    for(double i = dfStepLen; i < CenterLineLeft.get_Length(); i += dfStepLen)
    {
        OGRPoint pTmpPt;
        CenterLineLeft.Value(i, &pTmpPt);

        double dfHeightHypTmp = sqrt(i * i + dfHeight * dfHeight);
        double dfOffset = dfHeightHypTmp * dfFilmHalfWidth / dfFocusM;
        
        double dx = pTmpPt.getX() - ptCenter.getX();
        double dy = pTmpPt.getY() - ptCenter.getY();

        double ux = dfOffset * dx / i;
        double uy = dfOffset * dy / i;

        dfImageCenterX -= dfImageStepLen;

        OGRPoint ptUp(pTmpPt.getX() - uy, pTmpPt.getY() + ux);
        aPt2.push_back(ptUp);

        paGSPs[nGCPPos].dfGCPLine = poSrcDataset->GetRasterYSize();
        paGSPs[nGCPPos].dfGCPPixel = dfImageCenterX;
        paGSPs[nGCPPos].dfGCPX = ptUp.getX();
        paGSPs[nGCPPos].dfGCPY = ptUp.getY();
        paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
        paGSPs[nGCPPos].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPPos));
        nGCPPos++;

        ux = -dfOffset * dx / i;
        uy = -dfOffset * dy / i;

        OGRPoint ptDown(pTmpPt.getX() - uy, pTmpPt.getY() + ux);
        aPt1.push_back(ptDown);

        paGSPs[nGCPPos].dfGCPLine = 0;
        paGSPs[nGCPPos].dfGCPPixel = dfImageCenterX;
        paGSPs[nGCPPos].dfGCPX = ptDown.getX();
        paGSPs[nGCPPos].dfGCPY = ptDown.getY();
        paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
        paGSPs[nGCPPos].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPPos));
        nGCPPos++;
    }

    //proceed right side of frame
    OGRLineString CenterLineRight;
    CenterLineRight.addPoint(&ptCenter);
    CenterLineRight.addPoint(&ptEnd);

    dfStepLen = CenterLineRight.get_Length() / nStepCount;
    dfImageCenterX = double(poSrcDataset->GetRasterXSize()) / 2;

    for(double i = dfStepLen; i < CenterLineRight.get_Length(); i += dfStepLen)
    {
        OGRPoint pTmpPt;
        CenterLineRight.Value(i, &pTmpPt);

        double dfHeightHypTmp = sqrt(i * i + dfHeight * dfHeight);
        double dfOffset = dfHeightHypTmp * dfFilmHalfWidth / dfFocusM;
        
        double dx = pTmpPt.getX() - ptCenter.getX();
        double dy = pTmpPt.getY() - ptCenter.getY();

        double ux = dfOffset * dx / i;
        double uy = dfOffset * dy / i;

        dfImageCenterX += dfImageStepLen;

        OGRPoint ptUp(pTmpPt.getX() - uy, pTmpPt.getY() + ux);
        aPt3.push_back(ptUp);

        paGSPs[nGCPPos].dfGCPLine = 0;
        paGSPs[nGCPPos].dfGCPPixel = dfImageCenterX;
        paGSPs[nGCPPos].dfGCPX = ptUp.getX();
        paGSPs[nGCPPos].dfGCPY = ptUp.getY();
        paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
        paGSPs[nGCPPos].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPPos));
        nGCPPos++;

        ux = -dfOffset * dx / i;
        uy = -dfOffset * dy / i;

        OGRPoint ptDown(pTmpPt.getX() - uy, pTmpPt.getY() + ux);
        aPt4.push_back(ptDown);

        paGSPs[nGCPPos].dfGCPLine = poSrcDataset->GetRasterYSize();
        paGSPs[nGCPPos].dfGCPPixel = dfImageCenterX;
        paGSPs[nGCPPos].dfGCPX = ptDown.getX();
        paGSPs[nGCPPos].dfGCPY = ptDown.getY();
        paGSPs[nGCPPos].dfGCPZ = dfMeanHeight;
        paGSPs[nGCPPos].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPPos));
        nGCPPos++;
    }

    // add points to polygon
    OGRLinearRing Ring;

    Ring.addPoint(pt1.getX(), pt1.getY());
    for(int i = aPt1.size() - 2; i >= 0; --i)
        Ring.addPoint(aPt1[i].getX(), aPt1[i].getY());
    for(size_t i = 0; i < aPt3.size() - 1; ++i)
        Ring.addPoint(aPt3[i].getX(), aPt3[i].getY());
    Ring.addPoint(pt2.getX(), pt2.getY());
    Ring.addPoint(pt3.getX(), pt3.getY());
    for(int i = aPt4.size() - 2; i >= 0; --i)
        Ring.addPoint(aPt4[i].getX(), aPt4[i].getY());
    for(size_t i = 0; i < aPt2.size() - 1; ++i)
        Ring.addPoint(aPt2[i].getX(), aPt2[i].getY());

    Ring.addPoint(pt4.getX(), pt4.getY());

    Ring.closeRings();

    OGRPolygon Rgn;
    Rgn.addRingDirectly((OGRCurve*)Ring.clone());
    Rgn.assignSpatialReference(&oDstOGRSpatialReference);
    Rgn.flattenTo2D();

    OGREnvelope DstEnv;
    Rgn.getEnvelope(&DstEnv);

    //SaveGeometry(CPLResetExtension(sFileName, "shp"), Rgn, oDstOGRSpatialReference);

    // search point along image
    // add GCP to opened raster
    if(poSrcDataset->SetGCPs(nGCPPos, paGSPs, pszSpaRefDef) != CE_None)
    {
        printf( "Set GCPs failed\n" );
        exit( 1 );
    }

    // create warper
    char **papszTO = NULL;
    papszTO = CSLSetNameValue( papszTO, "METHOD", "GCP_TPS" );
    papszTO = CSLSetNameValue( papszTO, "NUM_THREADS", "4" );
    papszTO = CSLSetNameValue( papszTO, "DST_SRS", pszSpaRefDef );
    papszTO = CSLSetNameValue( papszTO, "SRC_SRS", pszSpaRefDef );
    papszTO = CSLSetNameValue( papszTO, "INSERT_CENTER_LONG", "FALSE" );

    GDALDriver *poOutputDriver = (GDALDriver *) GDALGetDriverByName( "GTiff" );
    CPLSetConfigOption( "CHECK_WITH_INVERT_PROJ", "TRUE" );
    void* hTransformArg = GDALCreateGenImgProjTransformer2( poSrcDataset, NULL, papszTO );
    GDALTransformerInfo* psInfo = (GDALTransformerInfo*)hTransformArg;

    double adfThisGeoTransform[6];
    double adfExtent[4];
    int nThisPixels, nThisLines;

    //
    if( GDALSuggestedWarpOutput2( poSrcDataset, psInfo->pfnTransform, hTransformArg, adfThisGeoTransform, &nThisPixels, &nThisLines, adfExtent, 0 ) != CE_None )
    {
        printf( "Suggest Output failed\n" );
        exit( 1 );
    }

    adfThisGeoTransform[0] = DstEnv.MinX;
    adfThisGeoTransform[3] = DstEnv.MaxY;

    int nPixels = (int) ((DstEnv.MaxX - DstEnv.MinX) / adfThisGeoTransform[1] + 0.5);
    int nLines = (int) ((DstEnv.MaxY - DstEnv.MinY) / -adfThisGeoTransform[5] + 0.5);

    GDALSetGenImgProjTransformerDstGeoTransform( hTransformArg, adfThisGeoTransform);

    CPLString sOutputRasterPath = CPLResetExtension(sFileName, "tif");
    GDALDataset  *poDstDataset = poOutputDriver->Create(sOutputRasterPath, nPixels, nLines, poSrcDataset->GetRasterCount(), GDT_Byte, NULL );
    if( NULL == poDstDataset )
    {
        printf( "Create Output failed\n" );
        exit( 1 );
    }
    poDstDataset->SetProjection( pszSpaRefDef );
    poDstDataset->SetGeoTransform( adfThisGeoTransform );

    //GDALDestroyGenImgProjTransformer( hTransformArg );
    //hTransformArg = GDALCreateGenImgProjTransformer2( poSrcDataset, poDstDataset, papszTO );

#ifdef APRROX_MAXERROR
    hTransformArg = GDALCreateApproxTransformer( GDALGenImgProjTransform,  hTransformArg, APRROX_MAXERROR);
    GDALTransformerFunc pfnTransformer = GDALApproxTransform;
    GDALApproxTransformerOwnsSubtransformer(hTransformArg, TRUE);
#else
    GDALTransformerFunc pfnTransformer = GDALGenImgProjTransform;
#endif // APRROX_MAXERROR

    // warp
    GDALWarpOptions *psWO = GDALCreateWarpOptions();

    psWO->eWorkingDataType = GDT_Byte;
    psWO->eResampleAlg = GRA_NearestNeighbour;

    psWO->hSrcDS = poSrcDataset;
    psWO->hDstDS = poDstDataset;

    psWO->pfnTransformer = pfnTransformer;
    psWO->pTransformerArg = hTransformArg;

    psWO->pfnProgress = GDALTermProgress;
    psWO->nBandCount = poSrcDataset->GetRasterCount();

    psWO->panSrcBands = (int *) CPLMalloc(psWO->nBandCount*sizeof(int));
    psWO->panDstBands = (int *) CPLMalloc(psWO->nBandCount*sizeof(int));

    for(int i = 0; i < psWO->nBandCount; ++i )
    {
        psWO->panSrcBands[i] = i+1;
        psWO->panDstBands[i] = i+1;
    }

    GDALWarpOperation oWO;
    if( oWO.Initialize( psWO ) == CE_None )
    {
#ifdef MULTI
        if( oWO.ChunkAndWarpMulti( 0, 0, poDstDataset->GetRasterXSize(), poDstDataset->GetRasterYSize() ) != CE_None)
#else //MULTI
        if( oWO.ChunkAndWarpImage( 0, 0, poDstDataset->GetRasterXSize(), poDstDataset->GetRasterYSize() ) != CE_None)
#endif //MULTI
        {
            const char* err = CPLGetLastErrorMsg();
            printf( "Warp failed.%s\n", err );
            exit( 1 );
        }
    }

    GDALDestroyWarpOptions( psWO );
    CSLDestroy( papszTO );

    CPLFree( pszSpaRefDef );
    GDALClose( (GDALDatasetH) poSrcDataset );
    GDALClose( (GDALDatasetH) poDstDataset );

    GDALDestroyDriverManager();

    return EXIT_SUCCESS;
}
    

