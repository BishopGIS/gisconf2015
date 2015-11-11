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

#define SEGMENT_STEPS 7
#define MULTI
#define APRROX_MAXERROR 0.125
#define CHECK_HAS_ENOUGH_ADDITIONAL_ARGS(nExtraArg) \
    do { if (iArg + nExtraArg >= nArgc) \
        Usage(CPLSPrintf("%s option requires %d argument(s)", papszArgv[iArg], nExtraArg)); } while(0)

void Usage(const char* pszErrorMsg = NULL)
{
    printf( "Usage: warper [-c lat long] [-nw lat long] [-ne lat long]\n"
            "              [-sw lat long] [-se lat long] filename");

    if( pszErrorMsg != NULL )
        fprintf(stderr, "\nFAILURE: %s\n", pszErrorMsg);

    exit( 1 );
}

// set default values in cm
// KH-4(A,B) http://fas.org/irp/imint/docs/kh-4_camera_system.htm
const double dfFocus = 60.9602;
const double dfFilmHeight = 5.45338; // 2.147 in
const double dfMeanHeight = 156;
const double dfMountAngle = 30.46 / 2; // camera mount angle
const double dfMountAngleRad = dfMountAngle * M_PI / 180;

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

OGRPoint GetCenterOfLine(OGRPoint* pt1, OGRPoint* pt2)
{
    OGRLineString SubCenterLine1;
    SubCenterLine1.addPoint(pt1);
    SubCenterLine1.addPoint(pt2);

    OGRPoint retPoint;
    SubCenterLine1.Centroid(&retPoint);

    return retPoint;
}

OGRPoint GetTangetPoint(OGRLineString &line, double dfDist)
{
    double dfLength = line.get_Length();
    OGRPoint ptBeg, ptEnd;
    line.StartPoint(&ptBeg);
    line.EndPoint(&ptEnd);

    double dx = ptEnd.getX() - ptBeg.getX();
    double dy = ptEnd.getY() - ptBeg.getY();

    double ux = dfDist * dx / dfLength;
    double uy = dfDist * dy / dfLength;

    OGRPoint resultPoint(ptEnd.getX() - uy, ptEnd.getY() + ux);

    return resultPoint;
}

double GetWidthForHeight(double dfHeight, double dfDist, double dfRatio)
{
    return sqrt(dfDist * dfDist + dfHeight * dfHeight) * dfRatio;
}

GDAL_GCP * PrepareGCP(const CPLString& sFileName, OGRPoint *pt1, OGRPoint *pt2, OGRPoint *pt3, OGRPoint *pt4, OGRPoint *ptCenter, const OGRSpatialReference &oDstOGRSpatialReference, const int nRasterSizeX, const int nRasterSizeY, int &nGCPCount, OGREnvelope &DstEnv)
{
    // to meters
    double dfFocusM = dfFocus / 100;
    double dfFilmHalfHeightM = dfFilmHeight / 200;
    double dfRatio = dfFilmHalfHeightM / dfFocusM;

    //create center point and line of scene
    OGRPoint ptShortSideBeg = GetCenterOfLine(pt1, pt4);
    OGRPoint ptShortSideEnd = GetCenterOfLine(pt2, pt3);

    OGRLineString lnTmp;
    lnTmp.addPoint(&ptShortSideBeg);
    lnTmp.addPoint(&ptShortSideEnd);

    lnTmp.Value(lnTmp.Project(pt1), &ptShortSideBeg);
    lnTmp.Value(lnTmp.Project(pt2), &ptShortSideEnd);

    double dfDist1 = pt1->Distance(pt2);
    double dfDist2 = pt2->Distance(pt3);
    double dfDist3 = pt3->Distance(pt4);
    double dfDist4 = pt4->Distance(pt1);
    double dfHalfWidth = (dfDist2 + dfDist4) / 4;
    double dfHalfHeight = (dfDist1 + dfDist3) / 4;

    double dfAltitudeAtSide = (dfHalfWidth * dfFocusM) / dfFilmHalfHeightM;
    double dfAltitudeAtSceneCenter = sqrt( dfAltitudeAtSide * dfAltitudeAtSide - dfHalfHeight * dfHalfHeight);
    double dfAltitude = dfAltitudeAtSceneCenter * cos(dfMountAngleRad); // 145 - 220 km
    double dfDistCenter = dfAltitudeAtSceneCenter * sin(dfMountAngleRad);

    OGRLineString lnCenterLeft;
    lnCenterLeft.addPoint(&ptShortSideBeg);
    lnCenterLeft.addPoint(ptCenter);
    OGRPoint ptSatCenter = GetTangetPoint(lnCenterLeft, dfDistCenter);

    std::vector<OGRPoint> aPt1, aPt2;
    int nTotalGCPCount = ((SEGMENT_STEPS + 1) * 2);
    GDAL_GCP *paGSPs = (GDAL_GCP *) CPLMalloc (nTotalGCPCount * sizeof(GDAL_GCP));
    GDALInitGCPs(nTotalGCPCount, paGSPs);

    double dfImageStepLen = double(nRasterSizeX) / SEGMENT_STEPS;
    double dfImageCenterX = 0;

    OGRLineString lnCenter;
    lnCenter.addPoint(&ptShortSideBeg);
    lnCenter.addPoint(ptCenter);
    lnCenter.addPoint(&ptShortSideEnd);

    double dfCenterLineLen = lnCenter.get_Length();
    double dfStepLen = dfCenterLineLen / SEGMENT_STEPS;
    double dfCenterLineHalfLen = dfCenterLineLen / 2;

    for(double i = 0; i <= dfCenterLineLen; i += dfStepLen)
    {
        OGRPoint ptTmp;
        lnCenter.Value(i, &ptTmp);

        double dfDist = fabs(dfCenterLineHalfLen - i);

        double dfWidthTmp = GetWidthForHeight(dfAltitudeAtSceneCenter, dfDist, dfRatio);

        OGRLineString lnTmpLine;
        lnTmpLine.addPoint(ptCenter);
        lnTmpLine.addPoint(&ptTmp);

        int direction = 1;
        if(dfCenterLineHalfLen < i)
            direction = -1;

        OGRPoint ptUp = GetTangetPoint(lnTmpLine, dfWidthTmp * direction);
        OGRPoint ptDown = GetTangetPoint(lnTmpLine, -dfWidthTmp * direction);

        //OGRPoint ptUp = GetValue(dfWidthTmp, lnTmpLine);
        //OGRPoint ptDown = GetValue(-dfWidthTmp, lnTmpLine);

        aPt1.push_back(ptUp);
        aPt2.push_back(ptDown);

        paGSPs[nGCPCount].dfGCPLine = 0;
        paGSPs[nGCPCount].dfGCPPixel = dfImageCenterX;
        paGSPs[nGCPCount].dfGCPX = ptDown.getX();
        paGSPs[nGCPCount].dfGCPY = ptDown.getY();
        paGSPs[nGCPCount].dfGCPZ = dfMeanHeight;
        paGSPs[nGCPCount].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPCount));
        nGCPCount++;

        paGSPs[nGCPCount].dfGCPLine = nRasterSizeY;
        paGSPs[nGCPCount].dfGCPPixel = dfImageCenterX;
        paGSPs[nGCPCount].dfGCPX = ptUp.getX();
        paGSPs[nGCPCount].dfGCPY = ptUp.getY();
        paGSPs[nGCPCount].dfGCPZ = dfMeanHeight;
        paGSPs[nGCPCount].pszId = CPLStrdup(CPLSPrintf("pt%d", nGCPCount));
        nGCPCount++;

        dfImageCenterX += dfImageStepLen;
    }

    // add points to polygon
    OGRLinearRing Ring;
    for(int i = 0; i < aPt1.size(); ++i)
        Ring.addPoint(aPt1[i].getX(), aPt1[i].getY());

    for(int i = aPt2.size() - 1; i >= 0; --i)
        Ring.addPoint(aPt2[i].getX(), aPt2[i].getY());

    Ring.closeRings();

    OGRPolygon Rgn;
    Rgn.addRingDirectly((OGRCurve*)Ring.clone());
    Rgn.assignSpatialReference(oDstOGRSpatialReference.Clone());
    Rgn.flattenTo2D();

    Rgn.getEnvelope(&DstEnv);

    SaveGeometry(CPLResetExtension(sFileName, "shp"), Rgn, oDstOGRSpatialReference);

    return paGSPs;
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

    OGRSpatialReference oOGRSpatialReference(SRS_WKT_WGS84);
    int nZoneNo = ceil( (180.0 + dfaCornersX[0]) / 6.0 );
    OGRSpatialReference oDstSpatialReference(SRS_WKT_WGS84);
    oDstSpatialReference.SetUTM(nZoneNo, dfaCornersY[0] > 0);

    // transform coordinates from WGS84 to UTM
    OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation( &oOGRSpatialReference, &oDstSpatialReference);
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

    // open input dataset
    GDALDataset *poSrcDataset = (GDALDataset *) GDALOpen( sFileName, GA_ReadOnly ); // GA_Update
    char* pszSpaRefDef = NULL;
    if( oDstSpatialReference.exportToWkt(&pszSpaRefDef) != OGRERR_NONE)
    {
        CPLFree( pszSpaRefDef );
        GDALClose( (GDALDatasetH) poSrcDataset );
        return EXIT_FAILURE;
    }

    // search point along image
    // add GCP to opened raster
    OGRPoint ptCenter(dfaCornersX[0], dfaCornersY[0]);
    OGRPoint pt1(dfaCornersX[1], dfaCornersY[1]); // NW Cormer
    OGRPoint pt2(dfaCornersX[2], dfaCornersY[2]); // NE Corner
    OGRPoint pt3(dfaCornersX[3], dfaCornersY[3]); // SE Corner
    OGRPoint pt4(dfaCornersX[4], dfaCornersY[4]); // SW Corner
    int nGCPCount = 0;
    OGREnvelope DstEnv;
    GDAL_GCP *paGSPs = PrepareGCP(sFileName, &pt1, &pt2, &pt3, &pt4, &ptCenter, oDstSpatialReference, poSrcDataset->GetRasterXSize(), poSrcDataset->GetRasterYSize(), nGCPCount, DstEnv);

    if(poSrcDataset->SetGCPs(nGCPCount, paGSPs, pszSpaRefDef) != CE_None)
    {
        Usage( "Set GCPs failed" );
        return EXIT_FAILURE;
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

    // suggest the raster output size
    if( GDALSuggestedWarpOutput2( poSrcDataset, psInfo->pfnTransform, hTransformArg, adfThisGeoTransform, &nThisPixels, &nThisLines, adfExtent, 0 ) != CE_None )
    {
        Usage( "Suggest Output failed" );
        return EXIT_FAILURE;
    }

    adfThisGeoTransform[0] = DstEnv.MinX;
    adfThisGeoTransform[3] = DstEnv.MaxY;

    int nPixels = (int) ((DstEnv.MaxX - DstEnv.MinX) / adfThisGeoTransform[1] + 0.5);
    int nLines = (int) ((DstEnv.MaxY - DstEnv.MinY) / -adfThisGeoTransform[5] + 0.5);

    GDALSetGenImgProjTransformerDstGeoTransform( hTransformArg, adfThisGeoTransform);

    // create new raster
    CPLString sOutputRasterPath = CPLResetExtension(sFileName, "tif");
    GDALDataset  *poDstDataset = poOutputDriver->Create(sOutputRasterPath, nPixels, nLines, poSrcDataset->GetRasterCount(), GDT_Byte, NULL );
    if( NULL == poDstDataset )
    {
        Usage( "Create Output failed" );
        return EXIT_FAILURE;
    }
    poDstDataset->SetProjection( pszSpaRefDef );
    poDstDataset->SetGeoTransform( adfThisGeoTransform );

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
            Usage( CPLSPrintf("Warp failed.%s", err) );
            return EXIT_FAILURE;
        }
    }

    // cleanup
    GDALDestroyWarpOptions( psWO );
    CSLDestroy( papszTO );

    CPLFree( pszSpaRefDef );
    GDALClose( (GDALDatasetH) poSrcDataset );
    GDALClose( (GDALDatasetH) poDstDataset );

    GDALDestroyDriverManager();

    return EXIT_SUCCESS;
}
    

