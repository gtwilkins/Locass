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

#ifndef SEED_H
#define SEED_H

#include "node.h"
#include "locus.h"
#include "query.h"
#include "seed_fork.h"

class Seed 
{
public:
    Seed( string &header, string &seq, Querier &inBwt, int errorCount, bool bestOnly=true );
    void assemble();
    void assembleGraph();
    vector<Locus*> getLoci();
    bool warning();
    virtual ~Seed();
private:
    bool extend( NodeRoll &exts, bool drxn );
    bool restart( vector<SeedFork> forks[2] );
    void assembleHaploid();
    void checkDivergent( NodeList &path );
    void checkDivergentBack( NodeList &div, NodeSet &pathSet, bool drxn );
    void deleteNodes( NodeSet &delSet );
    NodeList getLociGetPath( bool doForce );
    void getLociGetPath( Node* curr, NodeList &path, int &score, bool drxn );
    int getLociGetPathCurrScore( Node* curr, NodeList &path, bool drxn );
    NodeList getLociGetPathGetStarts( bool doForce );
    void getLociResolveDivergent( NodeIntMap &scores, NodeList &path, Node** forks );
    bool getLociSetConvergent( NodeList &path, Node** forks );
    bool getLociSetDivergent( NodeList &path, Node** forks );
    void merge();
    void plot();
    bool resolveBackFork( Node** forks, NodeSet &delSet );
    bool resolveBackForkBypass( Node* fork, NodeSet &delSet );
    bool resolveBackForkDouble( Node* fork, NodeSet &delSet );
    void resolveBackForkDouble( Node** forks, NodeSet &delSet );
    void resolveBackForks();
    
    string header_, seq_;
    NodeRoll seed_;
    NodeList nodes_;
    int32_t validLimits_[2], ends_[2], tether_[2];
    int readCount_;
    Querier &bwt_;
};

#endif /* SEED_H */

