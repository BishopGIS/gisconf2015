#*************************************************************************** 
# Project:  demo
# Purpose:  master class for gisconf 2015
# Author:   Dmitry Baryshnikov, dmitry.baryshnikov@nextgis.ru
#***************************************************************************
#   Copyright (C) 2015 NextGIS, info@nextgis.ru
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***************************************************************************

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

set(PRJ_NAME demo-gdal)
set(APP_NAME warper)

project(${PRJ_NAME})

#find needed packages
find_package(GDAL REQUIRED)
if(GDAL_FOUND)
    include_directories(${GDAL_INCLUDE_DIR})

    FOREACH(file_path ${GDAL_INCLUDE_DIR})       
        include_directories(${file_path}/ogr)
        include_directories(${file_path}/ogr/ogrsf_frmts)
        include_directories(${file_path}/port)
        include_directories(${file_path}/gcore)
        include_directories(${file_path}/alg)
        include_directories(${file_path}/frmts/vrt)
    ENDFOREACH()   
endif(GDAL_FOUND)

set(PROJECT_HEADERS
)

set(PROJECT_SOURCES
    main.cpp
)

add_executable(${APP_NAME} ${PROJECT_HEADERS} ${PROJECT_SOURCES})

if(GDAL_FOUND)
    target_link_libraries(${APP_NAME} ${GDAL_LIBRARIES})
endif() 


