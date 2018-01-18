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

#ifndef NODESTRUCTS_H
#define NODESTRUCTS_H

#include "types.h"
#include "parameters.h"
#include "query.h"
#include "node_types.h"

class Node;

struct Coords
{
    Coords( int32_t start, int32_t end, bool isRedundant=false );
    void offset( int32_t offset );
    int32_t& operator[]( bool i );
    int32_t coords[2];
    bool redundant;
};
typedef std::unordered_map<SeqNum, Coords> ReadCoords;

struct PairCoords
{
    int& operator[](int i){ return coords[i]; };
    int32_t coords[4];
};

struct Score
{
    Score();
    void clear();
    void operator += ( const Score &rhs );
    void operator -= ( const Score &rhs );
    bool operator > ( const Score &rhs );
    float& operator[](int i);
    float ratio( float hits=0, float misses=0 );
    
    float hits, misses;
};
typedef std::unordered_map<Node*, Score> ScoreMap;
typedef std::vector< std::pair<Node*, Score> > ScoreList;

struct PathVars
{
    PathVars( Querier &inBwt, NodeList* inNodes, unordered_set<ReadId> &remappedReads, bool inFinal, bool pathDrxn )
    : bwt( inBwt ), nds( inNodes ), usedIds( remappedReads )
    , finalise( inFinal ), drxn( pathDrxn ), furthest( 0 ), weak( NULL )
    , misassembled( NULL ), unspanned( NULL ), isInvalid( false ), rerun( false ) {};
    void addTarget( Node* t, NodeSet &fwdSet );
    void addTarget( NodeList &path );
    void resetPathReliable( NodeList &path );
    void setTarget( NodeSet &fwdSet, NodeList &islands );
    Querier &bwt;
//    NodeList &nodes, &islands;
    NodeList* nds;
    unordered_set<ReadId> &usedIds;
    NodeIntMap hits, reliable, islandHits;
    NodeFloatMap adjusted;
    NodeSet tSet, newSet;
    Node* weak,* misassembled,* unspanned;
    int32_t furthest, misassMark[2], misassEst[2];
    int32_t* reliLimits;
    bool drxn, isInvalid, rerun, finalise;
};

struct CloneTargetVars
{
    NodeSet tSetThis, tSetAlt;
    vector<int32_t> offsets;
};

struct CloneScore
{
    CloneScore(): selfPref( 0 ), altPref( 0 ), selfOnly( 0 ), altOnly( 0 ), selfHits( 0 ), selfScore( 0 ), altScore( 0 ), cumulScore( 0 ), cumulMisses( 0 ) {}
//    void init( Node* node, bool drxn );
    void addHit( Coords* coords, bool drxn );
    void mergeFurthest( const CloneScore &rhs, bool drxn );
    void operator += ( const CloneScore &rhs );
    bool operator > ( const CloneScore &rhs );
    void setScores();
    void setCumul( Node* node, int prevScore, int prevMisses );
    int selfPref, altPref, selfOnly, altOnly, selfHits, selfScore, altScore, cumulScore, cumulMisses, furthest[2];
};
typedef std::unordered_map<Node*, CloneScore> CloneScoreMap;

struct ReadMark
{
//    ReadMark();
//    ReadMark( const SeqNum &readId ): readId( readId ){}
    ReadMark( SeqNum &readId, Coords &coords, Lib* lib, bool drxn );
//    ReadMark( SeqNum num, Coords &inCoords, Lib &lib, bool markerRev, bool drxn );
//    void getPairCoords( UnpairedMap &pairs );
    void getMissScore( float &missScore, int32_t* limits );
    bool isValid( Coords &coords );
    bool isValid( Coords &coords, pair<int32_t, int32_t> &markOffset, pair<int32_t, int32_t> &hitOffset, bool drxn );
    void offset( int32_t offset );

    int32_t mark, estimate, coords[2];
//    PairCoords coords;
    SeqNum id;
//    uint8_t mode; // 0 = unknown direction and orientation, 1 = known orientation. 2 = both known
//    bool extDrxn;
};

struct PairingVars
{
    NodeIntMap pairs;
    unordered_set<SeqNum> hitIds;
    vector<ReadMark>* marks;
    NodeList tNodes;
};

struct ExtVars
{
    ExtVars( NodeList &nodes, NodeList &island, int32_t* limits, Querier &inBwt, bool doLoop=true, bool doOffset=true ) 
            : nodes( nodes ), island( island ), limits( limits ), bwt( inBwt ), doLoop( doLoop ), doOffset( doOffset ) {};
    NodeList &nodes, &island;
    int32_t* limits;
    NodeSet ante, dual, del, cloneSet, bypass, offset, rebranch;
    bool doLoop, doOffset;
    Querier &bwt;
};

struct LoopVars
{
    bool addLoopStart( Node* node, bool drxn );
    NodeList startList;
    int cutoff, target;
    NodeSet branches, clones, cloned, open, extended, rebranch;
    NodeSetList cloneFwdSets;
    CloneScoreMap scores;
    bool exitFound;
};

struct IslandVars
{
    IslandVars( ExtVars &ev, bool drxn ) : ev( ev ), drxn( drxn ){}
    
    Node* pathEnd;
    ExtVars ev;
    int round;
    int32_t limit;
    NodeList origin;
    NodeSet merged[2], ante;
    unordered_set<SeqNum> mpReads, peReadsReliable, peReads;
    bool drxn;
};

struct DrxnVars
{
    DrxnVars( Querier &inBwt, NodeList &inNodes, NodeList &inIslands, bool inDrxn )
    : bwt( inBwt ), nodes( inNodes ), islands( inIslands ), drxn( inDrxn ) {};
    Querier &bwt;
    NodeList &nodes, &islands;
    bool drxn;
};

struct MergeHit
{
    MergeHit() : node( NULL ), coords( NULL ), read( 0 ), overlap( 0 ) {};
    Node* node;
    SeqNum read;
    Coords* coords;
    uint16_t overlap;
};

struct PairHit
{
    PairHit(): node( NULL ), coords( NULL ){};
    PairHit( Node* node, Coords* coords ): node( node ), coords( coords ){};
    PairHit( Node* node, Coords* coords, int32_t diff ): node( node ), coords( coords ), diff( diff ){};
    Node* node;
    Coords* coords;
    int32_t diff;
};

struct Edge
{
    Edge( Node* node, int overlap, bool isLeap ): node( node ), overlap( overlap ), isLeap( isLeap ) {}
    float getOverlapWeakness();
    Node* node;
    int overlap;
    bool isLeap;
};

struct MapStruct
{
    MapStruct( string &inSeq, int32_t inMinLen, int32_t inCoord, int32_t target, bool drxn );
    string seq;
    NodeList nodes;
    int32_t coord;
    int32_t minLen;
    int32_t cutLen;
    int32_t len;
    bool isEnd;
};

struct MapResult
{
    MapResult(): len( 0 ), score( 0 ){ node[0] = node[1] = NULL; };
//    Node* l, * r;
    Node* node[2];
    int32_t coords[2];
    int32_t len;
    int score;
};

struct IslandRead
{
    IslandRead( string seq, ReadMark &mark );
    void setOverlap( IslandRead &read );
    string seq;
    IslandRead* edges[2];
    int overlaps[2];
    int32_t coord;
    ReadId id;
};

struct MapNode
{
    MapNode(): node( NULL ), estimate( 0 ){};
    void addEdge( MapNode* mapNode, int overlap, bool drxn );
    void addEdge( Node* edgeNode, int overlap, bool drxn );
    void checkLoop();
    static void collapse( vector<MapNode*> &mapNodes, vector<MapNode*> &mapEdges );
    static void fold( MapNode* mapNode, vector<MapNode*> &mapNodes, bool drxn );
    bool getSecondSeq( string &secondSeq, bool drxn );
    bool recoil();
    void removeEdge( MapNode* mn, bool drxn );
    void removeEdges( bool drxn );
    void setEdgeSeq();
    void setSecondSeq( bool drxn );
    void splitEdge( Node* node, int i, bool drxn );
    string seq;
    Node* node;
    vector<ReadId> ids;
    vector<int32_t> coords[2];
    vector<MapNode*> edges[2];
    vector<int> edgeOverlaps[2], bridgeOverlaps[2], bridgeCoords[2];
    vector<Node*> bridges[2];
    int32_t estimate;
};

struct NodeMapRead
{
    NodeMapRead( ReadId inId ): id( inId ) { lens[0] = lens[1] = 0; nodes[0] = nodes[1] = NULL; };
    bool checkDoubleHit();
    bool checkHit( ReadMark &mark, bool markDrxn, bool drxn );
    Node* nodes[2];
    ReadId id;
    string seq;
    int32_t coords[2][2];
    int lens[2];
};


struct NodeMapReadHits
{
    void add( NodeSet hitNodes, int32_t* hitCoords, bool doAdd, bool drxn );
    NodeSetList nodes;
    vector<int32_t> coords[2];
    vector<int> hits[2];
};

struct MappedReadEnd
{
    MappedReadEnd(){};
    MappedReadEnd( string inSeq, ReadId inId, int32_t* inCoords, int32_t offset, bool endDrxn, bool estDrxn );
    bool checkMap( int score, bool overDrxn );
    string getExtSeq( bool drxn );
    int getMinOverlap( int minlen, bool overDrxn );
    int getOverlap( bool overDrxn );
    void set( MappedReadEnd &read );
    bool withinLimits( int32_t* limits );
    string seq;
    ReadId id;
    int32_t coords[3], offset;
    int counts[2];
    int ol;
    bool drxn, doMap;
};

struct PathRead
{
    PathRead(){ edges[0] = edges[1] = NULL; ols[0] = ols[1] = 0; };
    static void offset( vector<PathRead*> reads, int offset );
    static void setOverlaps( vector<PathRead*> &reads );
    string seq;
    int32_t coord;
    vector<ReadId> ids;
    PathRead* edges[2];
    int ols[2];
};

#endif /* NODESTRUCTS_H */

