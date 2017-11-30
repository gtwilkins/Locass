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

#ifndef IMPORTFILE_H
#define IMPORTFILE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include "locus.h"


using namespace std;

class ImportFile
{
public:
    ImportFile( string filename );
    virtual ~ImportFile(){};
    
protected:
    vector<string> splitString( string s, char delim );
    
    string filename_;
};

class ImportAlign : public ImportFile
{
public:
    ImportAlign( string filename );
    
    vector<string> getSeqs();
};

class ImportState : public ImportFile
{
public:
    ImportState( string filename );
    
    vector<Locus*> getLoci();
    
private:
    typedef tuple<string, string, int> InEdge;
    typedef vector<InEdge> Edges;
    typedef unordered_map<string, Node*> NodeMap;
    
    void addEdgesToNodes( NodeMap &nodes, Edges &edges );
};

#endif /* IMPORTFILE_H */

