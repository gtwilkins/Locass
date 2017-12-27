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

#ifndef PATH_SEQUENCE_H
#define PATH_SEQUENCE_H

#include "types.h"
#include "node.h"

struct ReadEndMap
{
    ReadEndMap( string seq, ReadMark &mark, Node* hit, int32_t* hitCoords, int32_t offset, bool mapDrxn );
    
    string getSeq( unordered_set<ReadId> &usedIds, bool drxn );
    void move( bool mapDrxn );
    void setEdge( vector<ReadEndMap*> &reads, bool drxn );
    
    string seq;
    Node* node;
    ReadEndMap* edge;
    ReadId id;
    int32_t coords[3], off;
    int ol, edgeOl;
    bool drxn, edged, doMap;
};

struct SeqPath
{
    SeqPath( NodeList &path, Node* anchor, int32_t offset );
    string seq;
    NodeList nodes;
    vector<int32_t> nodeCoords[2];
    int32_t ends[2];
};

struct SeqPathReassemble : public SeqPath
{
    SeqPathReassemble( NodeList &path, Node* anchor, int32_t offset );
    ~SeqPathReassemble();
    
    bool doMap( PathVars &pv );
    int32_t getCoord( ReadEndMap* read );
    string getExtraSeq( ReadEndMap* read, int &Ol, bool drxn );
    void map( string &q, ReadMark &mark, int minOls[2] );
    void remap( string &q, ReadMark &mark );
    void setBridges( vector<SeqPathReassemble*> seqs, MapNode* mn, unordered_set<ReadId> &usedIds, float &bestScore, bool drxn, int fromDrxn );
    void setEdges();
    void setHalves( MapNode* mn, unordered_set<ReadId> &usedIds, int &bestScore, bool drxn );
    void sortReads();
    vector<ReadEndMap*> reads[2];
    vector< vector<ReadEndMap*> > paths[2];
    unordered_set<ReadId> usedIds[2];
    int32_t limits[2];
    bool complete;
};

struct SeqPathMerge : public SeqPath
{
    SeqPathMerge( NodeList &path, Node* anchor, int32_t offset );
    
    bool doMerge( PathVars &pv, NodeSet &delSet, bool drxn );
    Node* getSplitCoord( int32_t &coord, bool drxn );
    void merge( vector<SeqPathMerge*> seqs, bool endOnly, bool drxn );
    
    int32_t selfCoord, hitCoord, nodeCoord;
    int hitScore, ol;
    SeqPathMerge* hitSeq;
};


#endif /* PATH_SEQUENCE_H */

