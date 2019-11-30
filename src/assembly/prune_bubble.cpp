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

#include "prune_bubble.h"
#include <algorithm>

//BubbleBranch::BubbleBranch( Node* fork, vector<Node*> branch, int32_t coord, int32_t reads )
//: fork( fork ), branch( branch ), coord( coord ), reads( reads )
//{}
//
//bool BubbleBranch::pop( Node* tar, Nodes& fwd, bool drxn )
//{
//    for ( Edge& e : fork->edges_[!drxn] )
//    {
//        vector<Node*> path;
//        set( e.node, path, -e.ol, tar, fwd, !drxn );
//    }
//    
//    bool best = true;
//    for ( BubbleBranch& bb : alts )
//    {
//        for ( Node* node : bb.branch ) for ( auto& read : node->reads_ ) if ( !read.second.redundant && !read.second.coords[2] ) bb.reads++;
//        if ( bb.reads >= reads ) best = false;
//    }
//    if ( alts.empty() ) return false;
//    
//    return true;
//}

Bubble::Bubble( Node* fork, Edge& e, Nodes& fwd, bool drxn )
: coord( -e.ol ), reads( 0 ), branched( false )
{
    forks[!drxn] = fork;
    forks[drxn] = e.node;
    
    while ( forks[drxn]->edges_[0].size() == 1 && forks[drxn]->edges_[1].size() == 1 )
    {
        for ( auto& read : forks[drxn]->reads_ ) if ( !read.second.redundant && !read.second.coords[2] ) reads++;
        if ( forks[drxn]->cloned_ ) branched = true;
        branch.insert( drxn ? branch.end() : branch.begin(), forks[drxn] );
        coord += forks[drxn]->size() - forks[drxn]->edges_[drxn][0].ol;
        forks[drxn] = forks[drxn]->edges_[drxn][0].node;
    }
    
    if ( branched || forks[drxn]->edges_[!drxn].size() < 2 || branch.empty() ) forks[drxn] = NULL;
    
    if ( forks[drxn] ) for ( Edge& e : forks[drxn]->edges_[!drxn] )
    {
        vector<Node*> path;
        set( e.node, path, -e.ol, fwd, !drxn );
    }
}

Bubble::Bubble( vector<Node*>& path, int32_t dist, bool drxn )
: branch( drxn ? path : vector<Node*>( path.rbegin(), path.rend() ) ), coord( dist ), reads( 0 ), branched( false )
{
    forks[0] = forks[1] = NULL;
    for ( Node* node : branch )
    {
        for ( auto& read : node->reads_ ) if ( !read.second.redundant && !read.second.coords[2] ) reads++;
        if ( node->edges_[0].size() > 1 || node->edges_[1].size() > 1 || node->cloned_ ) branched = true;
    }
}

bool Bubble::confirm()
{
    if ( coord >=0 ) return false;
    int x = 0;
    return false;
}

bool Bubble::pop( Nodes& prunes, bool drxn )
{
    if ( !forks[0] || !forks[1] || alts.empty() || branched ) return false;
    
    bool bad = false;
    for ( Bubble& alt : alts ) if ( alt.reads > reads && alt.coord != coord ) bad = true;
    if ( bad && reads > 5 && !confirm() ) return false;
    if ( !bad ) return false;
    assert( coord < 0 );
    
    cout << "Bubble from " << forks[0]->ends_[1] << " to " << forks[1]->ends_[0] << endl;
    cout << "    Main branch: " << reads << " reads, " << coord << " distance." << endl;
    for ( Bubble& alt : alts ) cout << "    Alt branch: " << alt.reads << " reads, " << alt.coord << " distance." << endl;
    
    for ( Node* node : branch )
    {
        node->dismantle();
        prunes += node;
    }
    
    return true;
}

bool Bubble::prune( NodeRoll& nodes )
{
    bool pruned = false;
    for ( int again = 1; again-- > 0; )
    {
        Nodes prunes;
        for ( Node* node : nodes.nodes ) if ( resolve( node, prunes, !node->drxn_ ) ) pruned = true;
        for ( Node* node : prunes.nodes ) nodes.erase( node );
        again = !prunes.empty();
    }
    return pruned;
}

bool Bubble::resolve( Node* fork, Nodes& prunes, bool drxn )
{
    int count = 0;
    if ( !fork->bad_ ) for ( Edge& e : fork->edges_[drxn] ) if ( !e.node->bad_ ) count++;
    if ( count < 2 ) return false;
    
    Nodes fwd( fork, params.readLen*2, drxn, false );
    if ( count > 1 ) for ( Edge& e : fork->edges_[drxn] ) if ( !e.node->bad_ && e.node->edges_[!drxn].size() < 2 )
    {
        Bubble bubble( fork, e, fwd, drxn );
        if ( bubble.pop( prunes, drxn ) )
        {
            resolve( fork, prunes, drxn );
            return true;
        }
    }
    return false;
}

void Bubble::set( Node* node, vector<Node*>& path, int32_t dist, Nodes& fwd, bool drxn )
{
    if ( find( branch.begin(), branch.end(), node ) != branch.end() ) return;
    if ( node == forks[drxn] ) alts.push_back( Bubble( path, dist, 0 ) );
    if ( !fwd.find( node ) ) return;
    
    path.push_back( node );
    for ( Edge& e : node->edges_[drxn] ) set( e.node, path, dist + node->size() - e.ol, fwd, drxn );
    path.pop_back();
}
