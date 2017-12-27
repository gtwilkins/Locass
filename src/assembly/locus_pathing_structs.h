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

#ifndef LOCUSSTRUCTS_H
#define LOCUSSTRUCTS_H

#include "node.h"

struct Span
{
    Span( Node* bgn );
    
    Node *node;
    bool spanned, ended, complete;
};

struct PathBranch;
typedef vector<PathBranch> BranchList;
struct PathBranch
{
    PathBranch();
    PathBranch( Node* node, bool drxn );
    PathBranch( Node* node, Node* prev, bool drxn );
    bool isContinue( bool drxn );
    void setAlt( PathBranch &alt );
    static void setScores( PathVars &pv, BranchList &branches, vector<Span> &spans );
    void setScores( PathVars &pv );
    void setSpans( vector<Span> &spans, bool drxn );

    Node* branch,* altBranch,* fork;
//    Node* far;
    NodeSet fwdSet;
    Score score, altScore;
    bool valid, invalid, doesContinue;
    int overlap, spanned;
    float adjusted, hits, reliable, reads, islands;
    Node* farNodes[2];
    int32_t farCoords[2];
    uint8_t state; // 0 = has converged, 1 = will converge, 2 = divergent, 3 = divergent from convergent
};

struct PathScore
{
    PathScore(): bgnFork( NULL ), endFork( NULL ), len( 0 ), furthest( 0 ), readCount( 0 ), reliable( 0 ) {};
    Node* bgnFork,* endFork;
    NodeList path, altPath;
    int32_t len, furthest;
    int readCount, reliable;
    bool isHaploid, notHaploid;
};

struct Allele
{
    Allele( PathScore* scores );
    bool anyInSet( NodeSet &nodes );
    void anyInSet( NodeSet &nodes, bool* found );
    NodeList paths[2];
    bool complete;
};

struct Path
{
    void completeSpan( Node* node, bool drxn );
    void getAllelesInSet( vector<Allele> &rAlleles, NodeSet &nodes );
    void reset( NodeList &nodes );
    
    Node* fork;
    NodeList path, convFork;
    vector<Allele> alleles;
    vector<PathScore> convergents;
    BranchList divergent;
    vector<Span> spans;
    bool completed, isReliable, isMulti;
};

#endif /* LOCUSSTRUCTS_H */

