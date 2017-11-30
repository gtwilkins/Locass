/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the Locass software package <https://github.com/gtwilkins/Locass>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EXPORTFILE_H
#define EXPORTFILE_H

#include <iostream>
#include <fstream>
#include "node.h"
#include "locus.h"

using namespace std;

class ExportFile {
public:
    ExportFile( string filePrefix );
    
    void createFile();
    virtual void closeFile();
   
    
protected:
    virtual void setFileName() = 0;
    ofstream outFile_;
    string filename_;
    string filePrefix_;
    int fileCount_;
    bool overWrite_;
};

class ExportFasta : public ExportFile
{
public:
    ExportFasta( string filePrefix );
    void exportContig( string header, string &seq );
    void exportExtends( string header, string &origin, string &left, string &right );
    void exportLoci( vector<Locus*> loci );
    void exportLocus( Locus* locus );
    
    void addNode( string header, string seq );
    
private:
    void setFileName();    
};

class ExportState : public ExportFile
{
public:
    ExportState( string filePrefix );
    void exportLoci( vector<Locus*> loci );
    void exportLocus( Locus* locus );
    void setInFileName( string filename );
    void updateLocus( Locus* locus, bool del=false );
    
private:
    void getFileName();
    void setFileName();  
};


#endif /* EXPORTFILE_H */

