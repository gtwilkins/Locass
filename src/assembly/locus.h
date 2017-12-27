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

#ifndef LOCUS_H
#define LOCUS_H

#include <iostream>
#include <cassert>
#include "query.h"
#include "node.h"
#include "locus_pathing_structs.h"
#include "path_review.h"

extern struct Parameters params;

using namespace std;

class Locus
{
public:
    Locus( Querier &bwt );
    Locus( Querier &bwt, Node* origin );
    Locus( Querier &bwt, NodeList* nodes );
    Locus( Querier &bwt, ifstream &fh );
    virtual ~Locus();
    
    void locusTest();
    void addNode( Node* node, int drxn );
    void deleteNodes( NodeSet &delSet, bool drxn );
    void deleteNodes( NodeSet &delSet, NodeSet &goodSet, bool drxn );
    NodeListList exportNodes();
    void rebootLocus();
    NodeList getAllNodes();
    NodeSet getEnds( bool drxn );
    uint getLen();
    ScoreMap getScoreMap( Node* bgn, NodeList &tNodes, bool drxn );
    NodeSet getValidSet();
    NodeList getValidNodes();
    int getWeakestEdge( Node* begin, Node* end, NodeSet &fwdSet, bool drxn );
    void setForkLimits();
    void setHeader( string &header );
    void setOriginEnds();
    
    
private:
    NodeListList getEndLists( NodeSet &fwdSet, bool drxn );
    
    void summarise( bool drxn );
    
// LocusCalibration
public:
    void calibrate( LocusLibraryCount &lib );
    
// LocusExport
public:
    void exportLocus( ofstream &align, ofstream &dump );
    void getContig( string &header, string &seq, uint &locusId );
    void getExtends( string &header, string &origin, string &lSeq, string &rSeq );
private:
    void blankMismatchedEnds( vector< pair<Node*, int> > &nodePath );
    vector< pair<Node*, int> > getExportPath();
    void getExportPathBridge( vector< pair<Node*, int> > &nodePath, Node* lft, Node* rght );
    
// LocusExtension
public:
    void extendLocus();
private:
    Node* advanceEnd( Node* node, ScoreMap &scores, NodeList &sideNodes, NodeSet &delSet, bool drxn, bool isEnd=false );
    void advanceEndScoreNext( Node* node, ScoreMap &scores, ScoreList &nxtScores, ScoreList &stopScores, bool drxn );
    Node* advanceEndSetNext( ScoreList &nxtScores, ScoreList &stopScores, NodeList &sideNodes, NodeSet &delSet, bool isEnd, bool drxn );
    bool canExtend( bool drxn );
    bool debriefExtend( ExtVars &ev, bool drxn, bool rePlot=false );
    bool doExtendNode( Node* &node, bool drxn );
    bool doReplot( bool drxn );
    bool getForwardLoops( Node* &loop, ExtVars &ev, bool drxn );
    bool extendNodes( bool drxn );
    void resolveLoop( ExtVars &ev, bool drxn );
    void reviewSet( NodeSet &qSet, bool drxn );
    bool setBestEnds( Node* bgn, ScoreMap &scores, NodeSet &extSet, NodeSet &delSet, bool drxn, bool isEnd );
    NodeSet setBestEndsGetForks( Node* bgn, ScoreMap &scores, bool drxn );
    bool setBestEndsSetEnd( NodeList &endList, float &forkScore, float &endScore, ScoreMap &scores, NodeFloatList &endScores, float &cutoff, NodeSet &delSet, bool drxn );
    bool setBestEndsSetExtend( NodeFloatList &endScores, NodeSet &extSet, bool isEnd );
    bool setBestEndsSetFork( Node* fork, float &forkScore, ScoreMap &scores, NodeFloatList &endScores, float &endCutoff, NodeSet &delSet, bool drxn );
    void setExtend( bool drxn );
    void setExtendEnds( NodeSet &extSet, NodeSet &delSet, bool drxn );
    void setExtendLimits();
    void setExtendSides( NodeSet &extSet, NodeSet &delSet, bool drxn );
    bool updateExtension( bool drxn );
    
// LocusFill
public:
    void fill();
    void fill2();
    void fill3();
    void fill4();
    
// LocusLeaping
public:
private:
//    void leapTest( bool drxn );
    void resetIslandVars( bool drxn );
    void setIslandVars( Node* pathEnd, bool drxn );
    void leap();
    bool leap( Path &path, bool drxn );
    void leapManual( IslandVars &iv );
    bool leapBridge( IslandVars &iv );
//    void leapTest();
    void leapCleanup( IslandVars &iv );
    void leapCleanup( IslandVars &iv, NodeSet &goodSet );
    bool leapExtend( IslandVars &iv );
    void leapGetNodes( Node* node, NodeList &reliableNodes, int32_t* limits, bool drxn );
    void leapGetPeClumps( vector<ReadMark> &peClumps, vector<ReadMark> &peMarks );
    void leapGetReads( NodeList &reliableNodes, vector<ReadMark> &peClumps, vector<ReadMark> &peMarks, vector<ReadMark> &mpMarks, unordered_set<SeqNum> &readIds, int32_t* limits, bool drxn );
    bool leapReview( IslandVars &iv );
    bool leapSetExtend( IslandVars &iv, NodeSet* extSets );
    void leapSetExtendBridge( IslandVars &iv );
    void leapSetExtendEnds( IslandVars &iv, NodeSet* &extSets );
    void leapSetExtendMerge( IslandVars &iv );
    void leapSetExtendValidate( IslandVars &iv, NodeList &validNodes );
    bool leapSetIslandClumps( IslandVars &iv, vector<ReadMark> &peClumps, unordered_set<SeqNum> &readIds );
    void leapSetIslandSingles( IslandVars &iv, vector<ReadMark> &peMarks, vector<ReadMark> &mpMarks, unordered_set<SeqNum> &readIds, int32_t* limits );
    
// LocusPathing
public:
    void finalise();
private:
    bool plot();
    bool plotPath( Path &path, bool drxn );
    bool plotPathAdd( PathVars &pv, PathBranch &best, BranchList &branches, Path &path );
    bool plotPathConverge( PathVars &pv, PathBranch &convFork, Path &path, bool drxn );
//    bool plotPathConverge( PathVars &pv, Node* convEnd, Path &path, bool drxn );
//    bool plotPathConverge( PathVars &pv, PathBranch &convBranch, Node* convEnd, BranchList &newDiv, NodeSet &convSet, Path &path, bool drxn );
//    bool plotPathConvergeClones( PathScore* scores, bool drxn );
//    bool plotPathConvergeClonesAnyBetter( Node* origBgn, NodeSet &origSet, bool drxn );
//    bool plotPathConvergeCompare( PathScore* scores, Path &path, bool drxn );
//    void plotPathConvergeCompareReliable( PathScore* scores, Path &path, bool* onlyReliable, bool* moreReliable, bool drxn );
//    bool plotPathConvergeConflict( PathScore &prime, PathScore &alt, Path &path, bool drxn );
//    bool plotPathConvergePath( PathVars pv, PathBranch &convBranch, Node* convEnd, BranchList &newDiv, Path &path, NodeList &convPath, NodeSet &convSet, bool drxn );
//    bool plotPathConvergeRedundant( PathScore &alt, Path &path, bool drxn );
//    void plotPathConvergeScore( PathVars &pv, int32_t bgnCoord, int32_t bestOls[2], Node* postFork, PathScore &score, bool drxn );
    PathBranch plotPathGetBestNext( BranchList &branches, bool isFirst );
    PathBranch plotPathGetBestNext( PathBranch &convBranch, Node* convEnd, BranchList &branches, BranchList &newDiv, NodeSet &convSet, Path &path, NodeSet &delSet, bool drxn );
    BranchList plotPathGetBranches( PathVars &pv, bool drxn );
    BranchList plotPathGetBranches( Node* target, bool drxn );
    void plotPathGetFirst( PathVars &pv, PathBranch &best, Path &path, BranchList &branches, bool drxn );
    void plotPathGetNext( PathVars &pv, PathBranch &best, PathBranch &last, Path &path, BranchList &branches, bool drxn );
    void plotPathSetEnds( PathVars &pv, PathBranch &best, PathBranch &last, Path &path, bool drxn );
//    bool plotPathSetForks( PathVars &pv, Path &path, NodeSet &goodSet, NodeList &notReliableNodes, bool drxn );
    void plotPathSetSpans( bool drxn );
    void plotPathTrimBranch( PathVars &pv, PathBranch &branch, NodeSet &goodSet, NodeSet &delSet, bool isEnd, bool drxn );
    void plotPathTrimDivergent( PathVars &pv, Path &path, NodeSet &goodSet, NodeSet &delSet, bool drxn );
    void plotPathTrimEnds( PathVars &pv, PathBranch &best, BranchList &branches, NodeSet &goodSet, NodeSet &delSet, bool drxn );
    void plotPathTrimPrep( PathVars &pv, Path &path, bool drxn );
    void plotPathUpdateFork( PathVars &pv, Path &path, bool drxn );
    void plotPathUpdateMultiplicity( PathBranch &best, Path &path, NodeList &notReliableNodes, bool drxn );
    void plotPathUpdateReliability( PathBranch &best, Path &path, NodeList &notReliableNodes, bool drxn );
    void plotPathUpdateSpans( PathBranch &best, Path &path, bool drxn );
    void plotPrep( bool drxn );
    void reviewSpans( bool drxn );
    void setReliable( NodeList &path, bool drxn );
    bool updateForks( bool drxn );
    
public:
    uint8_t stopCodes_[2]; // 1: Couldn't extend, 2: Bad node, 3: Bad branches, 4: Excess coverage, 5: Hit limit, 6: Problematic
    string header_;
    double duration_, leapTime_, revTime_;
    int32_t ends_[2];
    uint desperation_[2], multiplicity_, forkCount_[2], loopCount_[2];
    ofstream* debugLog_, *export_;
private:
    NodeList originEnds_[2], toExtend_[2], endNodes_[2], sideNodes_[2], forkNodes_[2];
    NodeList nodes_[5];
    IslandVars* ivs_[2];
    DrxnVars* dvs_[2];
    unordered_set<SeqNum> leptReads_[2], remappedReads_[2];
    vector<Path> paths_[2];
    int32_t limits_[2], validLimits_[2], reliable_[2], forkLimits_[2];
    bool completed_[2], harsh_[2], finished_[2], leapFar_[2], finalise_;
    Querier &bwt_;
};

#endif /* LOCUS_H */

