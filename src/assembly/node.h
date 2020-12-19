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

#ifndef NODE_H
#define NODE_H

#include "types.h"
#include "node_structs.h"
#include "node_pairs.h"
#include "node_groups.h"
#include "node_types.h"
#include "node_limits.h"
#include "calibrate_structs.h"
#include "shared_functions.h"
#include <cassert>
#include <iostream>

extern struct Parameters params;

using namespace std;

class Node
{
#define PAUSED -1
#define DEAD_END 1
#define MANY_END 2
#define SHORT_END 3
#define BACK_END 4
#define BLUNT_END 5
    
public:
    Node();
    Node( string seq );
    Node( string seq, int32_t coord, int subGraph );
    Node( vector<Overlap> &reads );
    Node( string seq, Extension &ext, int32_t prevEnd, bool drxn );
    Node( string seq, QueryNode* ext, int32_t prevEnd, bool drxn, int graph );
    Node( QueryNode* ext, int32_t prevEnd, bool drxn );
    Node( string seq, int32_t beginCoord, int* stop, bool drxn );
    Node( Node* toClone, ExtVars &ev, bool drxn );
    Node( Node* toClone );
    Node( Node* toClone, NodeRoll& nodes, int graph, bool bad );
    Node( string& seq, NodeMark& mark, bool graph );
    Node( string seq, ReadMark &mark, Extension &ext, bool extDrxn, bool islandDrxn );
    Node( string seq, Extension &ext, int32_t prevEnd, bool extDrxn, bool islandDrxn );
    Node( int32_t anchor, int32_t overlap, int subGraph, bool drxn );
    Node( string seq, vector<Overlap> &overlaps, int drxn );
    Node( Querier& bwt, string& seq, int graph );
    Node( string& seq, Node* flip, int graph );
    Node( MapNode* mapNode, int drxn );
    Node( ReadStruct &read );
    Node( NodeMapRead &mapRead, bool drxn );
    Node( MapNode* mn, int i, int j, int drxn );
    Node( string seq, ReadId id, int32_t estimate, int drxn );
    Node( string seq, ReadId id, int32_t estimate, int drxn, bool bad );
    ~Node();

    void edgeTest( bool drxn );
    void nodeTest();
    void readTest();
    void pairTest();
    void addCloned( Node* clone );
    void addEdge( Node* node, bool drxn );
    void addEdge( Node* node, int overlap, bool drxn, bool doOffset=true, bool isLeap=false );
    void addEdge( Edge edge, bool drxn, bool reciprocate );
    void blankEnd( int32_t len, bool drxn );
    void cleanEdges( NodeRoll& nodes, bool drxn );
    void clearEdges( bool drxn, bool stop=true );
    vector<Node*> clones( bool inclSelf=true );
//    int countEdges( bool drxn, bool inclClones=false );
    bool deleteTest( bool drxn );
    bool deleteTest( NodeList &tNodes, bool drxn );
    NodeRoll dismantle();
    vector<Edge> edges( bool drxn );
//    vector<Edge> getAltEdges( bool drxn );
//    vector<Node*> getAltForks( bool drxn );
    int32_t getBiggestOffset( bool drxn );
    int getBestOverlap( bool drxn, bool inclClones=false );
    NodeRoll getDependent();
    Edge getEdge( Node* node, bool drxn );
    static int32_t getFurthest( int32_t q, int32_t t, bool drxn );
    string getHeader( string header );
    static int32_t getLen( vector<Node*>& path );
    int32_t getLength();
    int getOverlap( Node* node, bool drxn );
    vector<Node*> getPath();
    void getPath( vector<Node*>& path, bool drxn );
    static string getSeq( vector<Node*>& path );
    string getSeqEnd( int length, bool drxn );
    NodeList getTargetNodes( bool drxn, bool inclSelf=false, bool inclClones=false );
    bool inheritEdges( bool drxn );
    bool isBeyond( int32_t bgn, int32_t nd, bool drxn );
    bool isClone( Node* node );
    bool isClone( vector<Node*>& path, int i, bool drxn );
    bool isContinue( bool drxn );
    bool isContinue( int readLimit, bool drxn );
    bool isContinue( int readLimit, int ol, int minOl, bool drxn );
    bool isDeadEnd( bool drxn );
    bool isEdge( Node* node, bool drxn, bool inclClone=false );
    bool isEnded( bool drxn );
    bool isFurther( int32_t coord, bool endDrxn, bool drxn );
    bool isNearer( int32_t coord, bool endDrxn, bool drxn );
    bool isInRange( Node* node );
    bool pause( bool drxn );
    void propagateOffset( bool drxn );
    void propagateOffset( NodeSet &propagated, bool drxn );
    bool removeEdge( Node* node, bool drxn, bool reciprocate=false );
    bool removeEdges( NodeSet &removeSet, bool drxn );
    void setOrigin();
    int size();
    void sortEdges( bool drxn );
    void stop( int stopCode, bool drxn );
    void trimEnd( bool drxn );
    void trimEnd( int32_t coord, NodeList &nodes, bool drxn );
    bool unpause( bool drxn );
    bool unpauseable( bool drxn );
    
    static void deleteNodes( NodeSet &delSet, NodeList &nodes );
    void dismantleNode();
    void dismantleNode( NodeSet &delSet, bool drxn );

private:
    void trimSeq( int32_t keepCoord, bool trimDrxn, bool stop=true );
    
// NodeBridging
public:
    static bool bridgeIsland( IslandVars &iv, NodeSetList &islandSets );
    static void bridgeIslandDump( IslandVars &iv, NodeSetList &islandSets, Querier &bwt );
    void bridgeIslandDump( Querier &bwt, ofstream &fh, bool drxn );
    static bool mapBridge( Node* target, PathVars &pv, MapNode* mn );
//    static bool mapBridge( PathVars &pv, Node* target, MapNode* mn, int32_t* coords, bool drxn );
    static void setBridge( PathVars &pv, NodeSet &newSet, NodeList &nodes, vector<int> &overlaps, MapNode* mn, bool drxn );
private:
    static void bridgeIslandGetEnds( IslandVars &iv, NodeIntMap &limitMap, NodeList &endList, NodeSet &endSet, bool isIsland, bool drxn );
    static NodeList bridgeIslandGetForks( IslandVars &iv, NodeIntMap &limitMap, NodeSet &endSet, bool drxn );
    static void bridgeIslandGetPairs( IslandVars &iv, NodeSetList &islandSets, NodeSet &peSet, NodeIntMap &qLimits, NodeIntMap &tLimits );
    void bridgeIslandOffset( IslandVars &iv, NodeSet &islandSet, bool drxn );
    static void bridgeIslandConsolidateEnds( NodeIntMap &limitMap, NodeList &candidates, NodeSet &bckSet, bool drxn );
    static void bridgeIslandGetIslandEnds( IslandVars &iv, NodeSet &peSet, NodeList &endList, NodeIntMap &limitMap, NodeSet &islandSet );
    static NodeList bridgeIslandGetMainEnds( IslandVars &iv, NodeSet &peSet, NodeIntMap &limitMap );
    static NodeList bridgeIslandGetMainEndsSort( IslandVars &iv, NodeList &endList, bool drxn );
//    static bool bridgeIslandPerfect( IslandVars &iv, NodeList &mainEnds, NodeList &islandEnds );
//    bool bridgeIslandPerfect( IslandVars &iv, Node* node, int32_t coord, int overlap, bool drxn );
    static bool bridgeIslandSet( IslandVars &iv, NodeList &mainEnds, NodeList &islandEnds, vector<bool> &iAttempted, NodeList &tNodes );
    static void bridgeIslandSetEndSets( IslandVars &iv, NodeList &mainEnds, NodeList &islandEnds, NodeSetList &mainSets, NodeSet &islandShared );
    static bool bridgeIslandSetOffsets( IslandVars &iv, NodeList &islandEnds, NodeList &tNodes );
    void bridgeIslandTrimEnd( ExtVars &ev, int32_t coord, bool doTruncate, bool drxn );
    void cloneBridge( PathVars &pv, NodeList &hitNodes );
    static bool originBridge( PathVars &pv, NodeList edgeNodes[2], NodeList &newNodes );
    
// NodeBypassing
public:
    static void resolveBypass( ExtVars &ev, bool drxn );
    static void resolveOffsets( ExtVars &ev, bool drxn );
    static void resolveRebranch( ExtVars &ev, bool drxn );
private:
    bool resolveBypass( ExtVars &ev, bool doExtend, bool drxn );
    bool resolveOffset( ExtVars &ev, bool doExtend, bool drxn );
    bool reviewOffset( int32_t diffOffset, NodeSet* offsetSets, int& hitsCount, bool drxn );
    bool reviewOffset( NodeSet* offsetSets, bool drxn );
    
// NodeCalibration
public:
    bool calibrateSeed( ExtVars &ev );
    static int calibrateCount( NodeList &nodes, SeedLibraryCount &libs );
    static void calibrate( NodeList &nodes, LocusLibraryCount &lib );
    
// NodeCloning
public:
    void declone( Node* cull, NodeRoll& nodes );
    Node* declone( NodeRoll& nodes );
    CloneScore getCloneComparison( Node* clone, bool drxn );
    NodeSet getCloneSet( bool inclSelf=false );
private:
    void addClone( Node* node, bool isFirst=true );
    bool anyCloneInSet( NodeSet &nodeSet );
    void cloneNode( NodeList &nodes, bool drxn );
//    void extendLoopNode( LoopVars &lv, ExtVars &ev, bool drxn );
    NodeSetList getCloneEdgeSets( bool drxn );
    vector<int32_t> getCloneOffsets();
    vector<int32_t> getCloneOffsets( NodeList &tNodes );
    CloneScore getCloneScore( bool drxn );
    CloneTargetVars getCloneTargetVars( bool drxn );
    CloneTargetVars getCloneTargetVars( NodeList &tNodes, bool drxn );
    CloneScore getLoopBranchScore( LoopVars &lv, bool drxn );
    void removeClone( Node* node );
    
// NodeCompletion
public:
    void complete();
private:

// NodeCoverage
public:
    float getAdjustedReadCount( int32_t* limits );
    static float getAlleleCoverage( Node* forks[2], NodeList paths[2], bool drxn );
    int getCoverage( bool drxn );
    static float getCoverage( vector<Node*>& nodes, int minLen=0 );
    static float getCoverageMedian( vector<Node*>& nodes );
    float getCoverageCoeff();
    float getCoverageDrxn( int32_t dist, bool drxn, bool inclSelf );
    float getMultiplicity();
    bool getReliability();
    bool isMultiple();
    bool isReliable( bool doSet=false );
    static void resetReliability( NodeSet &nodes, bool drxn );
    void setCoverage();
    void setReliable( bool doForce=false );
    bool setUnreliable();
private:
    vector< pair<int32_t, uint16_t> > getCoverageMarks( bool drxn );
    int32_t getLengthForCoverage();
    void resetReliability();
    void resetReliabilityPropagate( NodeSet &propagated, bool drxn );
//    void setCoverage( ExtVars &ev, bool subGraph, bool drxn, bool isIsland=false );
    void setReliability();
    
// NodeExport
public:
    void exportNodeAlign( ofstream &fh, int32_t &coords );
    void exportNodeDump( ofstream &fh );
    static NodeList importNodes( ifstream &fh );
private:
    static NodeList importAlign( ifstream &fh, bool drxn );
    static void importNodes( NodeList &nodes, unordered_map<int, Node*> &nodeMap, ifstream &fh, bool drxn );
    static vector<string> importNodesGetSubStrings( string line );
    
// NodeExtension
public:
    bool extendable( Querier& bwt, bool drxn );
    void extendEdge( Querier& bwt, NodeRoll &nodes, int minOl, bool drxn );
    bool extendFork( Querier &bwt, NodeRoll &nodes, int32_t dist, int branchLimit, bool drxn );
    void extendNode( ExtVars &ev, bool drxn );
    void extendNode( Querier &bwt, NodeRoll &nodes, bool drxn );
    void extendSeed( ExtVars &ev, bool drxn );
    vector< pair<Node*, int32_t> > getExtendable( int32_t dist, bool drxn );
    void reEnd( ExtVars &ev, bool drxn );
    void rebranchNode( ExtVars &ev, bool drxn );
    Node* splitNode( int32_t splitBegin, NodeList &nodes, bool subGraph, bool drxn, bool isOriginal=true );
    vector<Node*> splitNode( int32_t splitBegin, NodeRoll& nodes, bool drxn );
private:
    void addAlts( vector<QueryNode*> &alts, NodeRoll &nodes, bool drxn, int graph );
    void addExtension( Extension &ext, ExtVars &ev, vector<MergeHit> &selfMerges, bool doesBranch, bool drxn );
    void addExtension( Extension &ext, ExtVars &ev, NodeSet &acceptable, NodeSet &ignore, bool doesBranch, bool drxn );
    void addExtensionMerge( MergeHit &merge, Extension &ext, ExtVars &ev, bool doesBranch, bool drxn );
    void addExtensions( vector<Extension> &exts, ExtVars &ev, bool doesBranch, bool drxn );
    void addExtensions( vector<QueryNode*> &exts, NodeRoll &nodes, bool drxn, int graph );
    void addSelfLoop( vector<MergeHit> &selfMerges, ExtVars &ev, bool drxn );
    void appendNode( Extension &ext, bool drxn );
    void appendNode( QueryNode* ext, bool drxn );
    void appendSeq( string seq, bool drxn );
    bool checkExtensionMerge( Extension &ext, MergeHit &merge );
    bool checkMerge( QueryNode* ext, MergeHit &merge, bool drxn );
    bool checkMergeRev( QueryNode* ext, MergeHit &merge, bool drxn );
    static bool checkMergeRev( QueryNode* ext, NodeRoll& nodes, bool drxn );
    void debriefExtension( ExtVars &ev, bool drxn );
    void extendComplete( ExtVars &ev, bool drxn );
    bool extendForward( ExtVars &ev, int base, int &best, int target, bool drxn );
    void extendForward( ExtVars &ev, NodeSet &acceptable, bool drxn );
    int32_t getBestMergeClone( MergeHit &merge, int32_t fromEnd, bool drxn );
    bool isMergeBypass( ExtVars &ev, MergeHit &merge, bool drxn );
    bool isMergeCis( ExtVars &ev, MergeHit &merge, bool drxn );
    void isMergeDualLocus( ExtVars &ev, MergeHit &merge, bool drxn );
    bool isMergeFwd( ExtVars &ev, MergeHit &merge, bool drxn );
    bool isMergeSelf( ExtVars &ev, MergeHit &merge );
    Node* mergeNode( NodeList &nodes, int32_t coord, bool subGraph );
    Node* mergeNode( NodeList &nodes, Coords* coords, bool subGraph, bool drxn );
    void pauseBad( bool graph );
    bool reEnd( Querier &bwt, NodeRoll &nodes, bool drxn );
    void removeReExtensions( vector<Extension> &exts, bool drxn, bool inclSelf=false );
    void setExtVars( ExtVars &ev, bool drxn );
    Node* splitNode2( int32_t splitBegin, NodeRoll& nodes, bool drxn );
    Node* splitNodeDual( int32_t* coords, NodeList &nodes, int subGraph );
    
// NodeFilling
public:
    void checkClones();
    void clearDupes();
    void clearReads();
    void fillReads( Querier &bwt, NodeSet &delSet );
    static void fillEdges( Querier& bwt, NodeRoll& nodes );
    void fillEdge( Querier& bwt, NodeRoll& nodes, Edge e );
    static void graphCover( string filename, NodeList &nodes );
    static void graphPairs( string filename, NodeList &nodes );
    void mapMates( Querier &bwt, int &count );
    static void mergeAll( NodeList* nodes, NodeSet &delSet );
//    Node* merge( bool drxn );
//    void merge( Nodes &nodes, bool drxn );
    void mergeDrxn( NodeSet &delSet, bool drxn );
    void recoil();
    void recoil( int32_t diff, bool drxn );
    static void recoordinate( NodeRoll& nodes );
    void recoordinate( Nodes& tested, bool drxn );
    static void refill( Querier& bwt, NodeRoll& nodes );
    static void remap( Querier& bwt, NodeRoll& nodes );
    bool remap( Querier &bwt );
    bool remap( NodeRoll& tar );
    
// NodeFolding
public:
    static bool foldAlleles( NodeList &nodes, Node* forks[2], NodeList paths[2], NodeSet &delSet, bool drxn );
    static void foldBranch( NodeList &nodes, Node* merges[2], int32_t coords[2], int ol, bool drxn );
    bool getNextReadCoord( int32_t &coord, bool coordDrxn, bool readDrxn );
private:
    void bluntEnd( bool drxn );
    NodeSet foldEdge( ExtVars &ev, Node* targetNode, bool drxn );
    static NodeSet foldEdgeHit( ExtVars &ev, MapResult &result, bool drxn );
    NodeSet foldEdgeMiss( ExtVars &ev, Node* targetNode, int32_t limit, bool drxn );
    NodeListList foldEdgeGetTargets( Node* targetNode, int32_t* limits, bool drxn );
    Node* foldEnd( ExtVars &ev, Node* altNode, bool drxn );
    Node* foldEndHit( ExtVars &ev, Node* hitNode, MapResult &result, bool drxn );
    Node* foldEndMiss( ExtVars &ev, NodeList &altPath, int32_t qLimit, int32_t tLimit, bool drxn );
    NodeListList foldEndGetAltPath( Node* foldEnd, Node* tFork, NodeList &altPath, int32_t* limits, bool drxn );
    NodeList foldEndGetFoldPath( NodeSet &altBckSet, bool drxn );
    Node* foldEndGetPairs( NodeIntMap &limitMap, NodeList &foldPath, NodeSet &tSet, bool drxn );
    void foldNodesIn( NodeSet &folded, NodeSet &foldable, NodeSet &acceptors, bool drxn );
    bool foldPrep( ExtVars &ev, int32_t coord, int &overlap, bool isEnd, bool drxn );
    void getFoldPath( NodeListList &targetPaths, NodeList &thisPath, NodeSet &usedSet, NodeSet &allowSet, int32_t* limits, bool allowFork, bool drxn );
    MapStruct getMapStructEnd( int32_t target, int32_t* limits, bool drxn );
    vector<MapStruct> getMapStructQuery( int32_t target, bool drxn );
    void getMapStructQuery( vector<MapStruct> &qStructs, MapStruct &qStruct, int32_t* limits, bool drxn );
    vector<MapStruct> getMapStructTarget( Node* targetNode, int32_t target, bool drxn );
    vector<MapStruct> getMapStructTarget( NodeListList &targetPaths, int32_t target, int32_t* limits, bool drxn );
    MapStruct getMapStructTarget( NodeList &targetPath, int32_t target, int32_t* limits, bool drxn );
    void getMapStructTargetCheckPaths( NodeListList &targetPaths, int32_t* limits, bool drxn );
    void interEdge( ExtVars &ev, Node* node, bool drxn );
//    static void joinEnds( IslandVars &iv, Node** nodes, int overlap );
//    static bool joinEnds( IslandVars &iv, Node** nodes, vector<Overlap> &overlaps, int overlap );
    void mapFold( MapResult &result, MapStruct &q, MapStruct &t, int32_t distWiggle, int minLen, bool fixedDist, bool drxn );
    void mapSeqs( MapStruct &l, MapStruct &r, int &bestScore, int &iBest, int &jBest, int &bestLen, int minLen, int32_t distWiggle, bool fixedDist, bool &isSet );
    void mapSetResult( MapResult &result, MapStruct &l, MapStruct &r, int i, int j, int len );
    void truncateNode( int32_t coord, bool drxn );

// NodeFurthest
public:
    void setFurthest( NodeSet &tSet, bool drxn );
private:
    void resetFurthest( Node* node );
    void resetFurthest( NodeSet &resetSet );
    
// NodeGroupings
public:
    static void getLimits( int32_t* limits, NodeList &nodes, bool doReset=true );
    static NodeList getNodeListIntersection( NodeList &a, NodeList &b );
    static void sortNodeListByFurthest( NodeList &nodes, bool endDrxn, bool drxn );
    
// NodeIslands
public:
    void extendIsland( IslandVars &iv, bool drxn );
    static bool islandCompare( IslandVars &iv, NodeSet &islandSet, NodeSet &otherSet );
    static void islandDelete( IslandVars &iv, Node* node );
    bool islandPairToBase( IslandVars &iv, NodeList &tNodes );
    void islandUpdateValid( SeqNum &readId );
    static bool islandReviewUnbridged( IslandVars &iv, NodeSet &islandSet );
    static bool islandReview( IslandVars &iv, NodeSet &islandSet );
    static void islandReview( IslandVars &iv, NodeSet &islandSet, NodeSetList &peIslands, NodeSetList &mpIslands );
    static void islandSetExtend( IslandVars &iv, NodeSet &islandSet, NodeSet* &extSets );
    static NodeList islandSetExtendGetAnchors( IslandVars &iv, NodeIntMap &pairMap, NodeSet &fwdSet, int &totalHits );
    static NodeIntMap islandSetExtendGetBridges( IslandVars &iv, NodeSet &islandSet, NodeList &peNodes, NodeList &mpNodes, NodeSet &fwdSet );
    bool islandSetExtendMulti( IslandVars &iv, NodeSet* &extSets, NodeSet &usedSet, int32_t* cutoffs, bool isBestAnchor, int &bestHits, int drxn );
    bool islandSetExtendSingle( IslandVars &iv, NodeSet* &extSets, int32_t* cutoffs, int drxns );
    static void islandSetValidate( IslandVars &iv, NodeList &testNodes, NodeSet &tested );
    void mergeIsland( ExtVars &ev, bool drxn, bool trimBack=false );
    bool overlapExtend( NodeList &nodes, int32_t* coords, NodeList &hitNodes, vector<int32_t>* hitCoords, bool subGraph, bool drxn );
    static void reviewMerged( ExtVars &ev, NodeSet &mergeSet, bool drxn );
    static bool seedIslandBridge( IslandVars &iv, vector<IslandRead*> &path, Node* hitNode, int32_t* coords, bool drxn );
    static void seedIslandsClump( IslandVars &iv, vector<ReadMark> &marks, unordered_set<SeqNum> &seeds, bool drxn );
//    static void seedIslandsClumps( IslandVars &iv, vector<ReadMark> &clumps );
    static bool seedIslandsClumpsCheck( IslandVars &iv, vector<MapNode*> &mns, MapNode* mapNode );
    static bool seedIslandsSingle( IslandVars &iv, ReadMark &mark, unordered_set<SeqNum> &seeds, bool drxn );
    static void trimIsland( IslandVars &iv, NodeSet &bgnSet );
    static void trimIsland( IslandVars &iv, NodeSet &fwdSet, NodeSet &hitBckSet, NodeIntMap &hitMap, NodeSet &goodSet );
    bool validate( IslandVars &iv );
private:
    void addExtension( Extension &ext, IslandVars &iv, bool doesBranch, bool drxn );
    void addExtensionMerge( MergeHit &merge, Extension &ext, IslandVars &iv, bool doesBranch, bool drxn );
    void addExtensions( vector<Extension> &exts, IslandVars &iv, bool drxn );
    bool canValidate( IslandVars &iv );
    static void checkExtension( Extension &ext, IslandVars &iv, MergeHit &merge );
    void pushValidLimits( IslandVars &iv, Node* markNode, Node* hitNode, int32_t &markCoord, Coords* coords );
    static bool seedIslandsCheckRead( IslandVars &iv, ReadMark &mark );
    static bool seedIslandsCheckSeq( IslandVars &iv, string &seq, ReadMark &mark );
    bool seedIslandsConfirm( IslandVars &iv, unordered_set<SeqNum> &seeds, bool drxn );
    bool setBlank( IslandVars &iv, NodeSet &foldable, bool drxn );
    
// NodeLeaping
//public:
//    static bool leap( Querier& bwt, NodeRoll& nodes, NodeList& path, bool drxn );
//    bool leap( Querier& bwt, NodeRoll& nodes, bool drxn );
//    bool leapSeed( Querier& bwt, NodeRoll& nodes );
//private:
//    void leapExtend( Querier& bwt, Node* seed, NodeRoll& nodes, bool drxn );
//    NodeIntList leapScore( NodeScores& scores, Nodesx& allowed, bool drxn );
//    Nodesx leapTarget( bool resetGood, bool drxn );
//    void leapTarget( Nodesx& tar, int extReadMax, bool resetGood, bool drxn );
    
// NodeLooping
public:
    static void resolveLoops( ExtVars &ev, bool drxn );
private:
    static void loop( LoopVars &lv, ExtVars &ev, bool drxn );
    void loopDelete( LoopVars &lv, ExtVars &ev, bool drxn );
//    bool loopExtend( LoopVars &lv, ExtVars &ev, int &best, int base, bool drxn );
//    void loopExtend( LoopVars &lv, ExtVars &ev, NodeSet &acceptable, bool drxn );
    static void loopReview( LoopVars &lv, ExtVars &ev, bool drxn );
    void loopWalk( LoopVars &lv, ExtVars &ev, int &best, bool drxn );
    bool loopWalkDoWalk( LoopVars &lv, bool drxn );
    void loopWalkEdges( LoopVars &lv, ExtVars &ev, bool drxn );
    void loopWalkEdgesCatalog( LoopVars &lv, NodeSet &usedEdges, vector<Edge> &edges, bool drxn );
    void loopWalkEdgesSet( LoopVars &lv, ExtVars &ev, NodeSet &usedEdges, vector<Edge> &edges, bool drxn );
    static Node* loopWalkGetNext( LoopVars &lv, bool drxn );
    bool loopWalkIsValidClone( LoopVars &lv, ExtVars &ev, int &cutoff, bool drxn );
    
// NodeNavigation
public:
    NodeSet getBetweenNodes( Node* node, bool drxn );
    NodeSet getConnectedNodes( bool sameGraph );
    void getConnectedNodes( NodeSet &nodes, bool sameGraph );
    NodeSet getDrxnNodes( bool drxn, bool sameGraph=false, bool inclSelf=false );
    void getDrxnNodes( NodeSet &nodes, bool drxn, bool sameGraph=false );
    void getDrxnNodes( NodeSet &nodes, bool drxn, int32_t limit );
    NodeSet getDrxnNodesInSet( NodeSet &inSet, bool drxn, bool inclSelf=false );
    void getDrxnNodesInSet( NodeSet &nodes, NodeSet &inSet, bool drxn );
    NodeSet getDrxnNodesNotInSet( NodeSet &notSet, bool drxn, bool inclSelf=false );
    void getDrxnNodesNotInSet( NodeSet &nodes, NodeSet &notSet, bool drxn );
    NodeOffsetMap getDrxnNodesOffset( bool drxn, int32_t limitDist=0, bool inclSelf=false );
    void getDrxnNodesOffset( NodeOffsetMap &nodes, bool drxn, int32_t &limit );
    NodeSet getEndNodes( bool drxn, bool inclStopped=true );
    static NodeSetList getNodeSetsExclusive( Node* a, Node* b, bool drxn );
    NodeSet getInvalidNodes( bool drxn );
    void getInvalidNodes( NodeSet &nodes, bool drxn );
    void getNextNodes( NodeSet &nextNodes, bool drxn );
    NodeSet getNextNodes( bool drxn);
    void getNextNodesInSet( NodeSet &nextSet, NodeSet &inSet, bool drxn );
    NodeSet getNextNodesInSet( NodeSet &inSet, bool drxn );
    static NodeSet getNotForwardSet( NodeSet &tmpCurrSet, bool drxn );
    int32_t getOffset( Node* node, bool drxn );
    void getRedundantNodes( NodeSet &nodes, int32_t* coords, bool drxn );
    bool offset( int32_t off );
    bool offset( Node* prv, bool drxn );
    void offsetEdge( Edge &e, bool drxn );
    void offsetForward( bool drxn, bool sameGraph=false, bool notSelf=false );
    bool setBad( bool drxn );
    bool setState();
    void setOffsetMap( NodeIntMap &offsets, NodeSet useSet, int32_t limit, bool drxn );
private:
    NodeIntList getOffsetEdges( bool drxn );
    bool isOffset( bool drxn );
    static void offsetForward( NodeSet &tmpCurrSet, bool drxn, bool notSelf=false );
    static void offsetForward( NodeSet &currSet, NodeSet &alreadySet, NodeSet &fwdSet, bool drxn );
    void offsetIsland( NodeSet &propagated, bool drxn );
    void setNotBad( bool drxn );
    
// NodePairing
public:
    void clearPaired( bool readd );
    Score getBranchScore( int32_t* limits, bool drxn );
    int getBridgeCount();
    int getFurthestAndReliable( NodeList &tNodes, int32_t &furthest, bool drxn );
    int32_t getFurthestPair( bool drxn );
    void getFurthestReliable( int32_t &furthest, int32_t &furthestReliable, NodeList &tNodes, bool drxn );
    int32_t getLengthForScoring( float hits );
    Score getPairScore( NodeList &tNodes, int32_t* limits, bool drxn );
    void getPairScore( Score &score, NodeList &tNodes, int32_t* limits, bool drxn );
    int getPairHits( Node* node );
    int getPairHitsTmp( bool drxn );
    int getPairHitsTotal();
    static void getReliablePairNodes( NodeList &nodes, NodeSet &reliableSet, bool drxn );
    static ScoreMap getScoreMap( NodeSet &qSet, int32_t* limits, bool drxn );
    void getUnpairedMarks( NodeList &nodes, vector<ReadMark> &peMarks, vector<ReadMark> &mpMarks, unordered_set<SeqNum> &usedIds, unordered_set<SeqNum> &readIds, int32_t* limits, bool inclMp, bool drxn );
    void limitScore( Score &score, bool inclNotValidSelf );
    static void resetPairing( NodeSet &nodes );
    void updatePairs();
//    void resetPairing( bool drxn );
private:
    void addPairs( NodeIntMap &pairs, unordered_set<SeqNum> &hitIds, bool drxn );
    void clearPairs();
    void clearPairsForward( bool drxn );
    void getFurthestPair( int32_t &furthest, bool drxn );
    void getFurthestPair( int32_t &furthest, int32_t &furthestReliable, NodeList &tNodes, bool drxn, bool getReliable=false );
    void getFurthestPairs( int32_t* furthest, NodeSet &qSet, CloneScoreMap &score, bool drxn );
    PairHit findReadPair( ReadMark &mark, bool &valid, bool drxn );
    PairHit findReadPairBestClone( NodeList &tNodes, ReadMark &mark, bool &valid, bool drxn );
    float getMissScore( int32_t* limits, bool drxn );
    int getPairsHits( NodeList &tNodes );
    vector<ReadId> getPairsIdsBase();
    bool isPaired( Node* node );
    CloneScore setCloneScore( PairingVars &pv, CloneTargetVars &ctv, bool drxn );
    CloneScore setCloneScore( NodeIntMap &pairs, NodeList &tNodes, vector<int32_t> offsets, vector<ReadMark> &marks, unordered_set<SeqNum> &hitIds, bool drxn, bool checkValid=true );
    void setPaired( Node* node );
    int setPairs( NodeList &tNodes, bool drxn );
    void setPairs( PairingVars &pv, bool drxn );
    void setPairs( NodeOffsetMap &fwdMap, NodeOffsetMap &revMap, bool drxn );
    void setPairs( Node* focus, Node* t, pair<int32_t, int32_t> &qOffset, pair<int32_t, int32_t> &tOffset, bool drxn );
    
// NodePruning
public:
    static void prune( Querier& bwt, NodeRoll& nodes );
    bool claimBranch( Querier& bwt, NodeRoll& nodes, bool drxn );
    bool isBlunt( int readCount, int readLimit, bool drxn );
    bool isBranchComplete( int readCount, int readLimit, bool drxn );
    bool isForkComplete( int32_t dist, int readLimit, bool drxn );
    bool pruneBranch( Querier& bwt, NodeRoll& nodes, int readsLeft, bool drxn );
    bool pruneFork( Querier& bwt, NodeRoll& nodes, bool drxn );
    bool pruneFwd( Querier& bwt, NodeRoll& nodes, Nodes* good, Nodes* tested, bool drxn );
    static void pruneLoops( NodeRoll& nodes );
    static void prunePaths( Querier& bwt, NodeRoll& nodes );
    static bool merge( NodeRoll &nodes );
    static void updateStates( NodeRoll& nodes );
private:
//    bool isFoldable( Node* from, bool drxn );
    bool isMergeable( Node* alt, vector<Node*> merges[2], bool drxn );
    bool isSubstantial( int readCount, int readLimit, bool drxn );
    void merge( NodeRoll &nodes, Nodes& tested, int& merged );
    bool merge( NodeRoll& nodes, bool drxn, bool verify=false );
    void updateState( Nodes& tested, Nodes& oppose, bool drxn );
    static void pruneBad( NodeRoll& nodes );
    void pruneBad( Nodes keep[3], Nodes tested[2], int32_t rev, int32_t dblrev, bool drxn );
    void pruneBad( Nodes keep[3], bool drxn );
    bool prepFork( Querier &bwt, NodeRoll& nodes, int32_t dist, bool drxn );
    static void pruneBlunt( Querier &bwt, NodeRoll& nodes );
    static void pruneIslands( Querier &bwt, NodeRoll& nodes );
//    bool pruneClone( bool drxn );
//    static void pruneClones( Querier& bwt, NodeRoll& nodes );
//    bool pruneClones( NodeRoll& nodes, int& i, bool drxn );
    bool prunePaths( Querier& bwt, NodeRoll& nodes, Nodes& tested, Nodes tried[2], int branched, bool drxn );
    void setStrong( Nodes& strong, int readCount, int readLimit, bool drxn );
    bool updateBad();
    bool updateBad( Nodes& tested, bool drxn );
//    int setStrong( int readLimit, bool drxn );
    
// NodeReads
public:
    void add( ReadId id, int32_t i, int32_t j, bool redundant, bool ignore=false, int unaligned=0 );
    bool add( ReadId id, string& seq );
    void addMatchedReads( vector<Overlap> &reads );
    bool addRead( SeqNum readId, int32_t bgn, int32_t nd, bool isRedundant, bool notClones=false );
    void addRead( NodeMapRead &mapRead, bool drxn );
    int countReads( bool notRedundant=false );
    static bool findOverlap( Node* &hitNode, int32_t* coords, string &seq, NodeList &nodes, int minOl, bool drxn );
    int getEndMarks( bool drxn );
    void getMarksCount( int counts[2] );
    Coords* getRead( ReadId id );
    ReadId getTerminalRead( bool drxn );
    bool isRedundant( int i, int j );
    bool isRedundant( Coords* coords );
    bool offsetNode( bool drxn );
    void resetMarks();
    void resetUnmarked( bool drxn );
    bool rmvMark( ReadId id, bool drxn );
private:
//    void addMark( ReadId id, int i, int j );
    void addMark( SeqNum readId, Coords &coords );
    void addMarks( vector<ReadId> &ids );
    void addRead( Overlap &read, int32_t anchor, bool olDrxn );
    bool anyReadBeyondCoord( int32_t coord, bool coordDrxn, bool drxn );
    bool anyReadInNode( unordered_set<SeqNum> &readIds );
    void cullMarks();
    void cullMarks( vector<NodeMark> &marks, bool drxn );
    int32_t findNextRead( int32_t mark, bool drxn );
    bool findRead( ReadId id, Coords*& hit, bool inclRedundant=true );
    vector<NodeMark> getMarks( bool drxn, bool pe );
    vector<ReadMark> getMarksBase( int drxn );
    bool getSplitCoords( int32_t coords[2], int split, bool drxn );
    void reAddMark( SeqNum readId, Coords &coords );
    void reAddMarks( vector<SeqNum> &readIds );
    void requery( Querier& bwt );
    void remark();
    void removeMark( SeqNum &readId );
    void removeMarks( unordered_set<SeqNum> &readIds, bool pushLimits, bool isPair, bool drxn );
    void resort();
    void sortMarks( vector<ReadMark> &marks, bool drxn );
    int32_t split( Node* node, int32_t cut, bool drxn );
    int32_t splitReads( Node* node, int32_t splitBegin, bool drxn );
    void trimReads( int32_t endCoord, bool drxn );

// Node Reliability
public:
//    bool checkMisassembly( PathVars &pv );
    bool checkMisassembly( PathVars &pv, NodeList* paths );
    bool getMisassembled( PathVars &pv, bool drxn );
    bool isMisassembled();
//    bool isMisassembled( int32_t reliLimit, bool drxn );
//    bool isMisassembled( NodeList &nodes, int32_t* validLimits, int32_t reliLimit, bool drxn );
    bool isMisassembled( PathVars &pv, bool drxn );
    bool isMisassembled( PathVars &pv, Node* farFork, NodeSet &forkSet, bool drxn );
    bool isMisassembledRev( int32_t markCoords[2], bool drxn );
//    bool reassemble( PathVars &pv, bool drxn );
private:
//    vector<PathSeq> getCombinedSeqs( int dist, bool drxn );
    bool isMisassembled( vector<int32_t> &marks, vector<int32_t> &ests, int32_t markLimits[2], int32_t estLimits[2], bool drxn );
//    bool isMisassembled( int &iBest, int &jBest, int32_t reliLimit, bool drxn );
    void mapRead( NodeMapRead &nmr, NodeSet &tSet );
//    bool reassemble( PathVars &pv, Node* targ, int32_t estimates[2], int32_t dist, bool remap, bool drxn );
    
// NodeSeed
public:
    static void addSeed( NodeRoll& nodes, ReadStruct& read );
//    void printSeed( Querier &bwt, NodeSet nodes[2] );
//    bool fixSeed( IslandVars &iv, bool drxn );
//    bool fixSeedGetNodes( IslandVars &iv, vector<PairPath> &paths, bool drxn );
//    void fixSeedGetPairs( Querier &bwt, vector<GoodPair> &gps, vector<BadPair> &bps, bool drxn );
//    vector<PairPath> fixSeedGraphPairs( vector<BadPair> &bps, bool drxn );
//    bool fixSeedSetEdge( IslandVars &iv, vector<GoodPair> &gp, bool drxn );
//    bool isSeed( int32_t seedLen );
//    bool plotSeed( IslandVars &iv, NodeSet &delSet, bool drxn, bool finished );
//    static void plotSeed( NodeList &ends, bool drxn );
//    void seedAdd( ReadStruct &read );
//    bool seedCongruent( ReadStruct &read, int32_t &coord );
//    bool seedDiminutive();
//    static void seedGetExtend( NodeList* extendNodes, NodeSet &seedSet, NodeSet &delSet, int32_t* limits );
//    static bool seedLeap( Querier &bwt, NodeList &nodes, int minCoord, int maxCoord );
//    void seedLoop( bool drxn );
//    bool seedLoop( NodeSet &loopSet, NodeSet &usedSet, int len, int limit, bool drxn );
//    void seedSetDrxnNodes( Node* fork, NodeList &nodes, bool drxn );
//    static Node* seedSetOrigin( NodeList &forkList );
//    void seedSplit( NodeList &nodes, int32_t coord );
//    static void seedValidate( NodeList &base, NodeSet &delSet, bool drxn );
//    void seedValidate( NodeSet &delSet, bool drxn );
//    void setIds( int &id, bool drxn );
    static void seedNode( Querier& bwt, NodeRoll& nodes, vector<Node*>& added, vector<Node*>& seeded, string& seq, ReadId id, int32_t coord, bool drxn );
//    static void seedNode( Querier& bwt, NodeRoll& nodes, ReadId id, int32_t coord, int drxn );
    static vector<Node*> seedNode( Querier& bwt, NodeRoll& nodes, string s, int lOl, int rOl, int drxn, int32_t coord=0 );
    static void trimSeed( Querier &bwt, NodeRoll &nodes );
private:
    bool addSeed( ReadStruct &read );
    void addSeed( ReadStruct &read, int ol );
    void addSeed( ReadStruct& read, NodeRoll& nodes, NodeIntList& ols, NodeIntIntList& splits );
//    bool seedJoin( Node* node, int32_t coord, bool drxn );
//    void seedJoinLoci( Node** nodes );
//    bool seedValidate( bool drxn );
//    bool seedValidate( NodeList &tNodes, NodeOffsetMap &fwdMap, NodeOffsetMap &revMap, bool drxn );
    
// NodeSlicing
public:
    void slice( PathVars &pv, NodeSet &delSet, bool drxn );
    bool slice( PathVars &pv, bool misassembled, bool drxn );
    bool sliceOrBridge( PathVars &pv, NodeList &hitNodes, vector<int32_t> hitCoords[2], NodeSet &delSet );
private:
    bool sliceExtend( IslandVars &iv, int32_t limit, bool drxn );
    bool splitExtend( Querier &bwt, Node* target, NodeList &nodes, bool drxn );
    
// NodeValidation
public:
    void forceValidate( int32_t* limits, NodeSetList &delSets, NodeSetList &goodSets );
    bool isValidated();
    bool isValidated( bool drxn );
    void propagateValidation( int32_t* limits, bool drxn );
    static void propagateValidation( NodeSet &qSet, int32_t* limits, bool drxn );
    void resetValid();
    void reviewEnd( ScoreMap &scores, NodeSet &validSet, NodeSet &goodSet, NodeSet &delSet, bool drxn );
    void reviewFork( ScoreMap &scores, NodeSet &goodSet, NodeSet &badSet, bool drxn );
//    bool trimValid( NodeSet &delSet, bool drxn );
    void setValid();
    bool validate( int32_t* limits );
    bool validate( bool drxn );
    static void validate( NodeSet &seedSet, NodeSet &delSet, int32_t* validLimits, int32_t* ends, bool doDel=true );
    void validateMerge( Node* mergeNode, bool drxn );
private:
    bool anyValid( bool drxn );
    int32_t getValidLimit( bool drxn );
    bool isValid( int32_t coord, bool drxn );
    bool isValidCoords( Coords &coords, int drxn );
    bool isValidHit( Coords* coords, bool drxn );
    void propagateValidation( int32_t* limits, NodeSet &testedSet, NodeSet &backedSet, bool drxn );
    void pullValidLimits();
    void pushValidLimits( int32_t mark, bool drxn );
    void pushValidLimits( Node* t, int32_t qMark, int32_t tMark, bool drxn );
    void pushValidLimts( NodeSet &bckSet, int hits, bool drxn );
    void setValid( bool drxn );
    void setValid( int32_t* limits );
    bool testInvalid( int32_t* limits, bool drxn );
    
// NodeVerify
public:
    static void reverify( NodeRoll& nodes );
    bool reverify();
    void setAllPairs();
    bool setVerified();
    int testVerify();
    static void verify( NodeRoll& nodes );
    void verify( bool drxn );
    void verify( Nodes &tested, bool drxn );
    bool verify();
    void verifyFork( int32_t dist, bool force, bool drxn );
private:
    void add( Node* tar, NodeMark& mark, bool drxn );
    void add( Node* tar, NodeMark& mark, Coords* hit, bool drxn );
    bool add( NodeRoll& nodes, string& seq, NodeMark& mark, bool drxn );
    bool bridgeVerified( Nodes& verified, bool drxn );
    bool canSetVerified();
    bool canVerify();
    void confirmVerified();
    void getVerified( Nodes& verified, bool bridge, bool drxn );
    void setUnverified();
    void setVerifyLimits( int32_t limits[2] );
    int verify( Node* tar, vector<NodeMark> &marks, int fwdHits, bool drxn );
    int verifyClone( Node* tar, vector<NodeMark> &marks, int fwdHits, bool drxn );
    int verifyClone( NodeRoll clones[2], vector<NodeOffsets>& offs, vector<NodeMark>& marks, int fwdHits, bool drxn );
    void verifyFill();
    

public:
    Node* allele_;
    NodeRoll* cloned_;
    ReadCoords reads_;
    vector<ReadMark> marks_[2];
    vector<NodeMark> pe_[2], mp_[2];
    int drxn_, id2_;
    NodeLimits ends_;
    string id_;
    string seq_;
    vector<Edge> edges_[2];
    NodePairs hits_;
    NodeIntMap pairs_;
    NodeList* clones_;
    Node* farPairNodes_[2]; // 0 = overall, 1 = reliable
    int32_t farPairCoords_[2];
    float coverage_;
    int extendCount_;
    int edgeCount_[2];
    int stop_[2];
    int32_t validLimits_[4];
    bool verified_, unreliable_, dontExtend_, bad_, branch_;
    bool assembled_[2], misassembled_[2];
    bool validated_, reliable_, mapped_, culled_;
private:
    NodeSet* paired_;
    Nodes* pairedNodes_;

};



#endif /* NODE_H */

