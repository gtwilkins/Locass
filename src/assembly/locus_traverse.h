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

#ifndef LOCUS_TRAVERSE_H
#define LOCUS_TRAVERSE_H

#include "node.h"
#include "query.h"

struct Vertex
{
    Vertex( string seq, int coord=0, bool drxn=true );
    string seq_;
    vector< pair<ReadId, pair<int, int> > > reads_[2];
    int coord_[2];
};

struct TFork
{
    TFork( string seq, bool drxn );
    TFork( string& base, string seq, ReadId id, int ext, int ol, bool drxn );
    TFork( TFork* tf, string seq, ReadId id, int ext, int ol, bool drxn );
    ~TFork();
    bool add( string& seq, ReadId id, bool drxn );
    bool add( string& seq, ReadId id, int ol, bool drxn );
    bool update();
    string base_, ext_, seq_;
    vector<TFork*> exts_;
    vector< pair<ReadId, int> > reads_;
    int ol_;
    bool exted_;
};

class Traverse
{
    Traverse( Node* node, bool drxn );
    
    void print( TFork* fork, bool drxn );
    void setFlayed( Querier& bwt, bool drxn );
    void setPairs( Querier& bwt, bool drxn );
    void setTraces( bool drxn );
    
    vector<Node*> path_;
    unordered_map<Node*, int32_t> coords_;
    unordered_map<ReadId, int32_t> reads_;
    unordered_map<ReadId, string> pairs_;
    vector<TFork*> forks_;
    int32_t len_, ends_[2];
public:
    static bool leap( Querier& bwt, Node* node, bool drxn );
    static bool trace( Querier& bwt, Node* node, bool drxn );
};

#endif /* LOCUS_TRAVERSE_H */

