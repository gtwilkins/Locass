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

#ifndef PATHCONVERGE_H
#define PATHCONVERGE_H

#include "node.h"
#include "node_structs.h"
#include "locus_pathing_structs.h"
#include "types.h"

struct AltPath
{
    AltPath( Node* node ): fork( node ), doesBridge( false ), doesContinue( false ) { hits = reads = score = 0; };
    bool doesConflict( AltPath &bck, bool drxn );
    Node* fork;
    NodeList path;
    float hits, reads, score;
    int32_t limits[2], far;
    bool doesBridge, doesContinue;
};

struct ConPath
{
    ConPath():diverged( false ){};
    bool doesConflict( AltPath &bck, bool drxn );
    Node* forks[2];
    NodeList paths[2];
    float score;
    bool diverged;
};

class PathReview {
public:
    PathReview( PathVars &pv, NodeList &path, int32_t* reliable, int32_t* forkLimits, bool anteFinished );
    bool review( Path &path, NodeList &sideNodes, NodeSet &delSet );
    
private:
    bool resolveAlleles( ConPath &con, NodeSet &delSet, bool &diverged );
    AltPath resolveBranch( Node* node, NodeSet &tSet, NodeSet &delSet );
    void resolveConverge( Node* forks[2], NodeSet &convSet, NodeSet &delSet );
    bool resolveConverge( Node* fork, NodeSet &delSet );
    AltPath resolveDiverge( Node* fork, NodeSet &delSet );
    bool resolveEnd( NodeSet &delSet );
    bool resolveForks( NodeSet &delSet );
    bool resolveMisassembly( NodeSet &delSet );
    bool resolveUnspanned( Path &path, NodeSet &delSet );
    void reviewSpans( Path &path );
    
    PathVars &pv_;
    Node* fork_,* truncate_;
    NodeList path_, forks_[2];
    NodeSet pathSet_, convSet_, sideSet_, bridgeSet_, usedSet_;
    NodeFloatMap scores_;
    vector<ConPath> cons_;
    vector<AltPath> divs_;
    int32_t reliLimits_[2], forkLimits_[2];
    bool drxn_, anteFinished_;
};

#endif /* PATHCONVERGE_H */

