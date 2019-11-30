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

#ifndef LEAP_H
#define LEAP_H

#include "node.h"
#include "types.h"
#include "query.h"

struct LeapHit
{
    LeapHit( Node* node, ReadId id, int32_t coord, int mp );
    static bool add( Node* q, vector<LeapHit> hits[2], NodeRoll& nodes, ReadId ids[2], int32_t coords[2], int mp );
    static void add( vector<LeapHit>& hits, Node* node, ReadId id, int32_t coord, int mp );
    static void affirm( vector<LeapHit>& hits, NodeRoll& nodes );
    Node* node;
    vector< pair<ReadId, int32_t> > hit[2];
};

struct LeapMark
{
    LeapMark( Node* node, string seq, ReadId qId, ReadId tId, int32_t qCoord, int32_t tCoord, bool mp );
    bool add( NodeRoll& nodes );
    void join( Node* node, bool drxn );
    void overlap( LeapMark* lm );
    bool seed( Querier& bwt );
    Node* node;
    string seq;
    ReadId ids[2];
    int32_t coords[2];
    int ols[3], joined;
    bool mp, redundant;
};

struct LeapEnd
{
    LeapEnd( LeapHit& lh, vector<LeapHit>& hits, bool orient, bool drxn );
    bool add( ReadId id, int32_t coord, bool mp );
    void add( LeapEnd* le, int32_t off );
    int32_t estimate();
    int getBranchReads( Node* node, int readCount, int readLimit, bool drxn );
    NodeDists getForks( bool drxn );
    bool obsolete( LeapEnd* le );
    bool retract( LeapEnd* le, bool drxn );
    void setForks( Node* node, NodeDists& forks, int readsLeft, bool drxn );
    Node* fork;
    NodeDists dists;
    unordered_set<ReadId> base, reads;
    vector<Node*> branches;
    vector< pair<ReadId, int32_t> > estimates[2];
    int score;
};

struct LeapEnds
{
    LeapEnds( vector<LeapHit>& hits, NodeRoll& nodes, bool consolidate, bool orient, bool drxn );
    ~LeapEnds();
    void cull( unordered_set<LeapEnd*>& culled, unordered_map<LeapEnd*, unordered_set<LeapEnd*> > forks[2] );
    vector< pair<Node*, int> > get( bool orient, bool drxn );
    vector<LeapEnd*> ends;
};

class Leap
{
public:
    static bool leapBranch( Querier& bwt, NodeRoll& nodes, Node* node, bool drxn );
    static Node* leapEnd( Querier& bwt, NodeRoll& nodes, Node* node, unordered_map<Node*, int32_t>& offs, bool drxn );
private:
    Leap( Querier& bwt, NodeRoll& nodes, unordered_map<Node*, int32_t>& offs, int32_t coord, bool drxn );
    ~Leap();
    void add( NodeRoll& nodes );
    Node* bridge( NodeRoll& nodes, bool drxn );
    static void confirmBranch( Node* node, unordered_map<Node*, int32_t>& offs, Nodes& branch, bool drxn );
    void extend( Querier& bwt, NodeRoll& nodes, bool drxn );
    void extend( Querier& bwt, NodeRoll& nodes, Node* node, bool drxn );
    static void getBranch( Node* node, unordered_map<Node*, int32_t>& offs, int32_t dist, int32_t limit, bool neg, bool drxn );
    Node* join( Querier& bwt, NodeRoll& nodes, bool drxn );
    Node* leap( Querier& bwt, NodeRoll& nodes, bool drxn );
    void remap( Querier& bwt, NodeRoll& nodes, bool drxn );
    bool seed( Querier& bwt, NodeRoll& nodes, bool drxn );
    vector<LeapHit> hits_[2];
    vector<LeapMark*> marks_;
    vector<Node*> joins_;
    int32_t coord_;
    bool seeded_;
};

class KmerLeap
{
    KmerLeap( Node* node, bool drxn );
};

#endif /* LEAP_H */

