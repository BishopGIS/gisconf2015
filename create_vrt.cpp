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
 
 
#define OVERLAP 2900

typedef struct inFile{
    CPLString path;
    int width;
    int height;
} INPUTFILE;

int main( int nArgc, char ** papszArgv )
{
    // register drivers
    GDALAllRegister();

    // get path from config
    if( nArgc < 2 )
        exit( -1 );

    CPLString sFileName = papszArgv[1];
    CPLString sVrtFileName = CPLResetExtension(sFileName, "vrt");
    
    // read files from archive and fill struct

    std::vector<struct inFile> aFiles;

    sFileName = "/vsitar/" + sFileName;

    char **papszItems = CPLReadDir(sFileName);
    if(papszItems == NULL)
        exit( -1 );

    int nMaxWidth = 0, nMaxHeight = 0;
    for(int i = 0; i < CSLCount(papszItems); ++i )
    {
        if( EQUAL(papszItems[i],".") || EQUAL(papszItems[i],"..") )
            continue;

        INPUTFILE stInFile;
        stInFile.path = sFileName;
        stInFile.path += "/";
        stInFile.path += papszItems[i];

        // open image and read it size
        GDALDataset *poDataset = (GDALDataset *) GDALOpen( stInFile.path, GA_ReadOnly );
        if( poDataset != NULL ) {
            stInFile.width = poDataset->GetRasterXSize();
            stInFile.height = poDataset->GetRasterYSize();
            GDALClose( (GDALDatasetH) poDataset );

            nMaxWidth += stInFile.width - OVERLAP;

            if(nMaxHeight < stInFile.height)
                nMaxHeight = stInFile.height;

            aFiles.push_back(stInFile);
        }
    }
    CSLDestroy( papszItems );

    // create vrt file
    GDALDriver *poDriver = (GDALDriver *) GDALGetDriverByName( "VRT" );
    GDALDataset *poVRTDS = poDriver->Create( sVrtFileName, nMaxWidth, nMaxHeight, 1, GDT_Byte, NULL );
    char **papszMD = NULL;
    int nOff = 0;
    for(int i = aFiles.size() - 1; i > -1; --i){
        const char *pszXML = CPLSPrintf("source_%d=<SimpleSource>"
                "<SourceFilename relativeToVRT=\"0\">%s</SourceFilename>"
                "<SourceBand>1</SourceBand>"
                "<SrcRect xOff=\"0\" yOff=\"0\" xSize=\"%d\" ySize=\"%d\"/>"
                "<DstRect xOff=\"%d\" yOff=\"%d\" xSize=\"%d\" ySize=\"%d\"/>"
                "</SimpleSource>", i, aFiles[i].path.c_str(), aFiles[i].width, aFiles[i].height,
                                  nOff, 0, aFiles[i].width, aFiles[i].height);
        papszMD = CSLAddString(papszMD, pszXML);
        nOff += aFiles[i].width - OVERLAP;
    }
    poVRTDS->GetRasterBand(1)->SetMetadata( papszMD, "new_vrt_sources");

    CSLDestroy( papszMD );

    GDALClose( (GDALDatasetH) poVRTDS );
    
}
