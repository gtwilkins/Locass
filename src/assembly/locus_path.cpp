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

#include "locus_path.h"
#include "leap.h"
#include <algorithm>

LocusPath::LocusPath( PathEdge* edge, unordered_set<NodePath*> forks[2], unordered_set<PathEdge*>& edged )
: path_{ edge->edge[0], edge->edge[1] }
{
    edged.insert( edge );
    for ( int d : { 0, 1 } ) extend( forks, edged, d );
}

bool LocusPath::extend( unordered_set<NodePath*> forks[2], unordered_set<PathEdge*>& edged, bool drxn )
{
    bool extended = false;
    branch_[drxn] = fork_[drxn] = blocked_[drxn] = false;
    PathEdge* pe = NULL;
    for ( NodePath* np = drxn ? path_.back() : path_[0]; np && !np->edges_[drxn].empty() && ( pe = np->edges_[drxn][0] ); )
    {
        if ( blocked_[drxn] = !isBest( pe ) ) break;
        if ( branch_[drxn] = ( forks[drxn].find( np ) != forks[drxn].end() ) ) break;
        if ( fork_[drxn] = ( forks[!drxn].find( pe->edge[drxn] ) != forks[!drxn].end() ) ) break;
        np = pe->edge[drxn];
        path_.insert( drxn ? path_.end() : path_.begin(), np );
        edged.insert( pe );
        extended = true;
    }
    return extended;
}

bool LocusPath::isBest( PathEdge* edge )
{
    for ( int d : { 0, 1 } ) if ( edge->edge[d]->edges_[!d][0]->edge[!d] != edge->edge[!d]) return false;
    return true;
}

bool LocusPath::isBridge( NodePath* np )
{
    if ( np->edges_[0].empty() || np->edges_[1].empty() ) return false;
    for ( int d : { 0, 1 } )
    {
        if ( np->edges_[d].empty() ) return false;
        if ( isBest( np->edges_[d][0] ) ) return false;
        int score[2]{ np->edges_[d][0]->score, np->edges_[d][0]->edge[d]->edges_[!d][0]->score };
        if ( score[0] >= score[1] / 3 ) return false;
    }
    return true;
}

bool LocusPath::isBridge( unordered_set<NodePath*> forks[2], unordered_set<NodePath*>& bridges )
{
    if ( !blocked_[0] || !blocked_[1] ) return false;
    PathEdge* edge[2]{ path_[0]->edges_[0][0], path_.back()->edges_[1][0] };
    int score[2]{ max( edge[0]->score, edge[1]->score ), edge[0]->edge[0]->edges_[1][0]->score + edge[1]->edge[1]->edges_[0][0]->score };
    for ( int i = 1; i < path_.size(); i++ ) score[0] = max( score[0], path_[i-1]->edges_[1][0]->score );
    if ( score[0] >= score[1] / 8 ) return false;
    bridges.insert( path_.begin(), path_.end() );
    for ( int d : { 0, 1 } ) if ( !isFork( edge[d]->edge[d], bridges, !d ) ) forks[!d].erase( edge[d]->edge[d] );
    return true;
}

bool LocusPath::isFork( NodePath* np, unordered_set<NodePath*>& bridges, bool drxn )
{
    for ( int i = 1; i < np->edges_[drxn].size(); i++ ) if ( bridges.find( np->edges_[drxn][i]->edge[drxn] ) == bridges.end() )
    {
        return np->edges_[drxn][i]->score + max( 0, np->edges_[drxn][i]->score-1 ) * 5 >= ( np->edges_[drxn][0]->score-2 ) / 2;
    }
    return false;
}

bool LocusPath::isRedundant( LocusPath* t )
{
    if ( t->path_.size() < path_.size() ) for ( int i = 0; i < path_.size(); i++ ) if ( path_[i] == t->path_[0] )
    {
        int j = 0;
        while ( j < t->path_.size() && path_[i+j] == t->path_[j] ) j++; 
        if ( j == t->path_.size() ) return true;
    }
    return false;
}

bool LocusPath::setBranches( vector<LocusPath*>& paths, unordered_set<NodePath*> forks[2], unordered_set<PathEdge*>& edged )
{
    vector< pair<NodePath*, int> > forked[2];
    for ( int d : { 0, 1 } ) for ( NodePath* np : forks[d] ) forked[d].push_back( make_pair( np, np->edges_[d][1]->score ) ); 
    for ( int d : { 0, 1 } ) sort( forked[d].begin(), forked[d].end(), []( pair<NodePath*, int>&a , pair<NodePath*, int>& b ){ return a.second > b.second; } );
    
    for ( int d : { 0, 1 } ) for ( pair<NodePath*, int>& fork : forked[d] )
    {
        
    }
}

bool LocusPath::setBridges( vector<LocusPath*>& paths, unordered_set<NodePath*> forks[2], unordered_set<NodePath*>& bridges, unordered_set<PathEdge*>& edged )
{
    bool bridged = false;
    for ( int i = 0; i < paths.size(); i++ ) if ( paths[i]->isBridge( forks, bridges ) )
    {
        bridges.insert( paths[i]->path_.begin(), paths[i]->path_.end() );
        delete paths[i];
        paths.erase( paths.begin() + i-- );
        bridged = true;
        for ( int j = 0; j < paths.size(); j++ ) if ( i != j ) for ( int d : { 0, 1 } ) if ( paths[j]->extend( forks, edged, d ) )
        {
            for ( int k = j+1; k < paths.size(); k++ ) if ( paths[j]->isRedundant( paths[k] ) )
            {
                delete paths[k];
                paths.erase( paths.begin() + k-- );
                if ( k <= i ) i--;
            }
        }
    }
    
    return bridged;
}

bool LocusPath::setRedundant( vector<LocusPath*>& paths )
{
    for ( int i = 0; i < paths.size(); i++ ) for ( int j = 0; j < paths.size(); j++ ) if ( i != j )
    {
        
    }
}

vector<LocusPath*> LocusPath::create( vector<NodePath*>& paths )
{
    vector<PathEdge*> edges;
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) PathEdge::sort( np->edges_[d] );
    for ( NodePath* np : paths ) for ( PathEdge* pe : np->edges_[1] ) edges.push_back( pe );
    PathEdge::sort( edges );
    
    vector<LocusPath*> loci;
    unordered_set<NodePath*> forks[2], bridges;
    unordered_set<PathEdge*> edged;
    for ( NodePath* np : paths ) if ( isBridge( np ) ) bridges.insert( np );
    for ( NodePath* np : paths ) for ( int d : { 0, 1 } ) if ( isFork( np, bridges, d ) ) forks[d].insert( np );
    
//    for ( int again = 1; again-- > 0; )
//    {
//        bool pathed = false;
//        for ( PathEdge* pe : edges ) if ( pe->score > 0 && isBest( pe ) 
//                && edged.find( pe ) == edged.end()
//                && forks[0].find( pe->edge[1] ) == forks[0].end() 
//                && forks[1].find( pe->edge[0] ) == forks[1].end() )
//        {
//            paths.push_back( new LocusPath( pe, forks, edged ) );
//            pathed = true;
//        }
//        if ( again = setBridges( paths, forks, bridges, edged ) ) continue;
//        if ( again = setBranches( paths, forks, edged ) ) continue;
//        if ( pathed ) again = true;
//    }
//    
//    int x = 0;
    return loci;
}

Diploid::Diploid( NodePath* a, NodePath* b )
: homo_( a == b )
{
    alleles_[0] = { a };
    alleles_[1] = { b };
}

void Diploid::easyExtend( bool drxn )
{
    vector<NodePath*> paths[2];
}

vector<NodePath*> Diploid::getEdges( bool i, bool drxn )
{
    vector<NodePath*> edges;
    if ( alleles_[i].empty() ) return edges;
    for ( PathEdge* pe : ( drxn ? alleles_[i].back() : alleles_[i][0] )->edges_[drxn] ) edges.push_back( pe->edge[drxn] );
    return edges;
}

void Diploid::create( NodePath* a, NodePath* b, vector<NodePath*> alleles[2] )
{
    for ( int i : { 0, 1 } ) alleles[i].clear();
    vector<Diploid> diploid{ Diploid( a, b ) };
    for ( int d : { 0, 1 } ) for ( int again = 1; again-- > 0; )
    {
        Diploid* dip = d ? &diploid.back() : &diploid[0];
        vector<NodePath*> edges[2]{ dip->getEdges( 0, d ), dip->getEdges( 1, d ) };
        assert( edges[0].empty() == edges[1].empty() );
        assert( edges[0].size() < 3 && edges[1].size() < 3 );
        for ( int i : { 0, 1 } ) if ( !again && ( again = edges[i].size() == 1 ) )
        {
            if ( edges[!i].size() == 1 ) diploid.insert( d ? diploid.end() : diploid.begin(), Diploid( edges[0][0], edges[1][0] ) );
            else if ( edges[!i].size() == 2 )
            {
                int j = edges[!i][0] == edges[i][0];
                assert( edges[!i][!j] == edges[i][0] );
                dip->alleles_[i].insert( d ? dip->alleles_[i].end() : dip->alleles_[i].begin(), edges[i][0] );
                dip->alleles_[!i].insert( d ? dip->alleles_[!i].end() : dip->alleles_[!i].begin(), edges[!i][j] );
            }
        }
        if ( !again && ( again = ( edges[0].size() == 2 && edges[1].size() == 2 ) ) )
        {
            for ( int i : { 0, 1 } ) assert( edges[0][i] == edges[1][0] || edges[0][i] == edges[1][1] );
            diploid.insert( d ? diploid.end() : diploid.begin(), Diploid( edges[0][0], edges[0][1] ) );
        }
    }
    
    struct DipPair{ Diploid* dip[2]; int net; };
    vector<DipPair> pairs, conflicts;
    for ( int i = 0; i < diploid.size(); i++ ) if ( !diploid[i].homo_ ) for ( int j = i+1; j < diploid.size(); j++ ) if ( !diploid[j].homo_ )
    {
        unordered_set<ReadId> ids[2];
        for ( int ii: { 0, 1 } ) for ( int jj: { 0, 1 } ) for ( NodePath* l : diploid[i].alleles_[ii] ) for ( NodePath* r : diploid[j].alleles_[jj] )
        {
            for ( PathPairs* pp : l->paired_[1] ) if ( pp->node_[1] == r )
            {
                ids[ii == jj].insert( pp->ids_.begin(), pp->ids_.end() );
                ids[ii == jj].insert( pp->mates_.begin(), pp->mates_.end() );
                break;
            }
        }
        DipPair paired;
        paired.dip[0] = &diploid[i];
        paired.dip[1] = &diploid[j];
//        for ( int k : { 0, 1 } )for ( ReadId id : ids[k] ) if ( ids[!k].find( id ) == ids[!k].end() ) paired.net += k ? 1 : -1;
        int scores[2]{0};
        for ( int k : { 0, 1 } )for ( ReadId id : ids[k] ) if ( ids[!k].find( id ) == ids[!k].end() ) scores[k]++;
        paired.net = scores[1] > scores[0] ? scores[1] - scores[0]*2 : scores[1]*2 - scores[0];
        if ( max( scores[0], scores[1] ) > min( scores[0], scores[1] )*2 ) pairs.push_back( paired );
    }
    sort( pairs.begin(), pairs.end(), []( DipPair& a, DipPair& b ){ return abs( a.net ) > abs( b.net ); } );
    unordered_map<Diploid*, bool> used;
    while ( !pairs.empty() )
    {
        used.insert( make_pair( pairs[0].dip[0], 0 ) );
        used.insert( make_pair( pairs[0].dip[1], pairs[0].net < 0 ) );
        pairs.erase( pairs.begin() );
        for ( int i = 0; i < pairs.size(); i++ )
        {
            int found[2]{ 2, 2 };
            for ( int j : { 0, 1 } )
            {
                auto it = used.find( pairs[i].dip[j] );
                if ( it != used.end() ) found[j] = it->second;
            }
            if ( found[0] == 2 && found[1] == 2 ) continue;
            bool added = true;
            if ( found[0] != 2 && found[1] != 2 )
            {
                if ( ( pairs[i].net > 0 ) != ( found[0] == found[1] ) )
                {
                    conflicts.push_back( pairs[i] );
                }
                added = false;
            }
            for ( int j : { 0, 1 } ) if ( found[j] == 2 ) used.insert( make_pair( pairs[i].dip[j], pairs[i].net > 0 ? found[!j] : !found[!j] ) );
            pairs.erase( pairs.begin() + i-- );
            if ( added ) i = -1;
        }
    }
    
    for ( int i = 0; i < diploid.size(); i++ )
    {
        bool flip = false;
        auto it = used.find( &diploid[i] );
        if ( it != used.end() ) flip = it->second;
        alleles[0].insert( alleles[0].end(), diploid[i].alleles_[flip].begin(), diploid[i].alleles_[flip].end() );
        alleles[1].insert( alleles[1].end(), diploid[i].alleles_[!flip].begin(), diploid[i].alleles_[!flip].end() );
    }
}

