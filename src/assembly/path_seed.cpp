/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "path_seed.h"
#include <algorithm>

//SeedBranch::SeedBranch( Edge &e, NodeList &tNodes, NodeSet &shared, bool drxn )
//: n( e.node ), ol( e.ol )
//{
//    score[0] = score[1] = reads = 0;
//    for ( Node* fwd : e.node->getDrxnNodes( drxn, false, true ) )
//    {
//        reads += fwd->reads_.size();
//        for ( ReadMark &mark : fwd->marks_[drxn] )
//        {
//            for ( Node* t : tNodes )
//            {
//                auto it = t->reads_.find( mark.id );
//                if ( it == t->reads_.end() ) continue;
//                score[ shared.find( t ) != shared.end() ]++;
//                if ( shared.find( t ) == shared.end() )
//                {
//                    int32_t tCoords[2] = { it->second[0], it->second[1] };
//                    int x = 0;
//                }
//            }
//        }
//    }
//}
//
//bool SeedBranch::operator >( SeedBranch &rhs )
//{
//    if ( score[0] > rhs.score[0] ) return true;
//    if ( score[0] < rhs.score[0] ) return false;
//    if ( score[1] > rhs.score[1] ) return true;
//    if ( score[1] < rhs.score[1] ) return false;
//    if ( ol > rhs.ol ) return true;
//    if ( ol < rhs.ol ) return false;
//    if ( reads > rhs.reads ) return true;
//    return false;
//}

//SeedPath::SeedPath( Node* node, NodeSet &usedSet, bool drxn )
//: ended( false ), doBranch( true )
//{
//    NodeList tNodes;
//    for ( Node* fwd : node->getDrxnNodes( drxn, false, true ) )
//    {
//        if ( abs( fwd->ends_[drxn] - node->ends_[drxn] ) < params.maxPeMax ) tNodes.push_back( fwd );
//    }
//    
//    Node* curr = node;
//    while ( curr )
//    {
//        path.push_back( curr );
//        SeedBranch* best = NULL;
//        for ( Edge &e : curr->edges_[!drxn] )
//        {
//            if ( usedSet.find( e.node ) != usedSet.end() ) continue;
//            SeedBranch* branch = new SeedBranch( e, tNodes, shared, !drxn );
//            if ( !best || *branch > *best )
//            {
//                if ( best ) delete best;
//                best = branch;
//            }
//            else delete branch;
//        }
//        curr = best ? best->n : NULL;
//        if ( best ) delete best;
//    }
//    reverse( path.begin(), path.end() );
//    
//    for ( Node* bck : node->getDrxnNodes( !drxn ) )
//    {
//        if ( find( path.begin(), path.end(), bck ) != path.end() ) continue;
//        if ( usedSet.find( bck ) != usedSet.end() ) continue;
//        starts.insert( bck );
//    }
//}
//
//NodeList SeedPath::getTargets( bool drxn )
//{
//    NodeList tNodes( path.begin(), path.end() );
//    for ( Node* node : starts )
//    {
//        if ( find( tNodes.begin(), tNodes.end(), node ) == tNodes.end() ) tNodes.push_back( node );
//    }
//    for ( NodeList &alt : alts )
//    {
//        for ( Node* node : alt )
//        {
//            if ( find( tNodes.begin(), tNodes.end(), node ) == tNodes.end() ) tNodes.push_back( node );
//        }
//    }
//    for ( Node* bck : path[0]->getDrxnNodes( !drxn ) )
//    {
//        if ( find( tNodes.begin(), tNodes.end(), bck ) == tNodes.end() ) tNodes.push_back( bck );
//    }
//    
//    return tNodes;
//}
//
//NodeSet SeedPath::getUsed( vector<SeedPath> &paths, bool drxn )
//{
//    NodeSet usedSet;
//    for ( SeedPath &path : paths )
//    {
//        path.path[0]->getDrxnNodes( usedSet, !drxn );
//        usedSet.insert( path.path.begin(), path.path.end() );
//        for ( NodeList &alt : path.alts )
//        {
//            usedSet.insert( alt.begin(), alt.end() );
//        }
//    }
//    return usedSet;
//}
//
//bool SeedPath::resolve( vector<SeedPath> &paths, NodeSet &altFwd, NodeSet &usedSet, int i, bool drxn )
//{
//    if ( paths[i].ended ) return false;
//    vector<SeedPath*> sameEnds;
//    Node* fork = paths[i].path.back();
//    NodeSet endSet;
//    for ( int j = 0; j < paths.size(); j++ )
//    {
//        if ( fork == paths[j].path.back() ) sameEnds.push_back( &paths[j] );
//        if ( i == j || paths[j].ended ) continue;
//        endSet.insert( paths[j].path.back() );
//        paths[j].path.back()->getDrxnNodes( altFwd, drxn );
//    }
//    
//    if ( sameEnds.size() > 1 && altFwd.find( fork ) != altFwd.end() ) return false;
//    altFwd.insert( endSet.begin(), endSet.end() );
//    
//    if ( sameEnds.size() == 1 )
//    {
//        usedSet = SeedPath::getUsed( paths, drxn );
//        return false;
//    }
//    
//    vector< vector<SeedBranch*> > pathEdges;
//    for ( int j = 0; j < sameEnds.size(); j++ )
//    {
//        sameEnds[j]->shared.insert( fork );
//        vector<SeedBranch*> edges;
//        NodeList tNodes = sameEnds[j]->getTargets( drxn );
//        for ( Edge &e : fork->edges_[drxn] )
//        {
//            edges.push_back( new SeedBranch( e, tNodes, sameEnds[j]->shared, drxn ) );
//        }
//        sort( edges.begin(), edges.end(), []( SeedBranch* a, SeedBranch* b ){ return *a > *b; } );
//        pathEdges.push_back( edges );
//    }
//    
//    for ( vector<SeedBranch*> &e1 : pathEdges )
//    {
//        for ( int j = 0; j < e1.size(); j++ )
//        {
//            bool doErase = false;
//            for ( vector<SeedBranch*> &e2 : pathEdges )
//            {
//                if ( !j && e1[j]->score[0] ) break;
//                if ( &e1 == &e2 ) continue;
//                for ( int k = 0; k < e2.size(); k++ )
//                {
//                    if ( e1[j]->n != e2[k]->n ) continue;
//                    if ( *e2[k] > *e1[j] ) doErase = true;
//                    if ( j && k && !e1[j]->score[0] && !e2[k]->score[0] )
//                    {
//                        doErase = true;
//                        while( k < e2.size() )
//                        {
//                            delete e2[k];
//                            e2.erase( e2.begin() + k );
//                        }
//                    }
//                }
//            }
//            if ( doErase )
//            {
//                while( j < e1.size() )
//                {
//                    delete e1[j];
//                    e1.erase( e1.begin() + j );
//                }
//            }
//        }
//    }
//    
//    int x = 0;
//    
//    for ( int j = 0; j < sameEnds.size(); j++ )
//    {
//        if ( pathEdges[j].empty() ) sameEnds[j]->ended = true;
//        else
//        {
//            for ( int k = 0; k < j; k++ )
//            {
//                if ( pathEdges[j][0]->score[0] ) break;
//                if ( pathEdges[k].empty() ) continue;
//                if ( pathEdges[j][0]->n == pathEdges[k][0]->n ) sameEnds[j]->ended = true;
//            }
//            if ( sameEnds[j]->ended ) continue;
//            sameEnds[j]->path.push_back( pathEdges[j][0]->n );
//        }
//    }
//    
//    return true;
//}
//
//void SeedPath::setBranches( Node* node, NodeSet &altFwd, NodeSet &usedSet, NodeList &tNodes, bool drxn )
//{
//    for ( int i = 0; i < branches.size(); )
//    {
//        if ( branches[i] == node || usedSet.find( branches[i] ) != usedSet.end() )
//        {
//            branches.erase( branches.begin() + i );
//            branchFwd.erase( branchFwd.begin() + i );
//        }
//        else if ( altFwd.find( node ) != altFwd.end() && branchFwd[i].find( node ) != branchFwd[i].end() )
//        {
//            NodeSet pathSet = node->getDrxnNodesInSet( branchFwd[i], !drxn, true );
//            NodeList allele;
//            Node* curr = branches[i];
//            while ( curr && usedSet.find( curr ) == usedSet.end() )
//            {
//                if ( find( tNodes.begin(), tNodes.end(), node ) == tNodes.end() ) tNodes.push_back( node );
//                allele.push_back( curr );
//                SeedBranch* best = NULL;
//                SeedBranch* bestAlt = NULL;
//                for ( Edge &e : curr->edges_[drxn] )
//                {
//                    if ( pathSet.find( e.node ) == pathSet.end() || altFwd.find( curr ) != altFwd.end() ) continue;
//                    SeedBranch* branch = new SeedBranch( e, tNodes, shared, drxn );
//                    if ( !best || *branch > *best )
//                    {
//                        if ( best ) delete best;
//                        if ( bestAlt ) delete bestAlt;
//                        best = branch;
//                        NodeList tAlt;
//                        NodeSet altShared;
//                        for ( Node* bck : best->n->getDrxnNodes( !drxn ) )
//                        {
//                            if ( find( tNodes.begin(), tNodes.end(), bck ) == tNodes.end() ) tAlt.push_back( bck );
//                        }
//                        bestAlt = new SeedBranch( e, tAlt, altShared, drxn );
//                    }
//                    else delete branch;
//                }
//                curr = NULL;
//                if ( best )
//                {
//                    if ( best > bestAlt ) curr = best->n;
//                    delete best;
//                    delete bestAlt;
//                }
//            }
//            
//            alts.push_back( allele );
//            usedSet.insert( allele.begin(), allele.end() );
//            branches.erase( branches.begin() + i );
//            branchFwd.erase( branchFwd.begin() + i );
//        }
//        else i++;
//    }
//}
//
//bool SeedPath::setLoop( Node* node, NodeList &tNodes, bool drxn )
//{
//    for ( int i = 0; i < path.size(); i++ )
//    {
//        if ( path[i] != node ) continue;
//        SeedBranch* best = NULL;
//        NodeList branchPath;
//        for ( int j = i; j < path.size(); j++ )
//        {
//            for ( Edge &e: path[j]->edges_[drxn] )
//            {
//                if ( find( tNodes.begin(), tNodes.end(), e.node ) != tNodes.end() ) continue;
//                SeedBranch* branch = new SeedBranch( e, tNodes, shared, drxn );
//                if ( !best || *branch > *best )
//                {
//                    if ( best ) delete best;
//                    best = branch;
//                    branchPath = { path.begin() + i, path.begin() + j + 1 };
//                    branchPath.push_back( best->n );
//                }
//                else delete branch;
//            }
//        }
//        if ( best )
//        {
//            for ( Node* node : branchPath )
//            {
//                if ( find( tNodes.begin(), tNodes.end(), node ) == tNodes.end() ) tNodes.push_back( node );
//                path.push_back( node );
//            }
//            delete best;
//        }
//        else ended = true;
//        
//        return true;
//    }
//    
//    return false;
//}
//
//bool SeedPath::setNext( NodeSet &altFwd, NodeSet &usedSet, NodeSet &bckSet, NodeList &tNodes, bool drxn )
//{
//    if ( altFwd.find( path.back() ) != altFwd.end() ) return false;
//    if ( ended ) return false;
//    
//    vector<SeedBranch*> edges;
//    SeedBranch* best = NULL;
//    for ( Edge &e : path.back()->edges_[drxn] )
//    {
//        SeedBranch* branch = new SeedBranch( e, tNodes, shared, drxn );
//        edges.push_back( branch );
//        if ( !best || *branch > *best ) best = branch;
//    }
//    if ( !best ) ended = true;
//    else if ( setLoop( best->n, tNodes, drxn ) )
//    {
//        branches.clear();
//        branchFwd.clear();
//    }
//    else
//    {
//        setBranches( best->n, altFwd, usedSet, tNodes, drxn );
//        path.push_back( best->n );
//        if ( find( tNodes.begin(), tNodes.end(), best->n ) == tNodes.end() ) tNodes.push_back( best->n );
//        usedSet.insert( best->n );
//        for ( SeedBranch* edge : edges )
//        {
//            if ( edge == best || usedSet.find( edge->n ) != usedSet.end() ) continue;
//            NodeList tAlt;
//            NodeSet altShared;
//            for ( Node* bck : best->n->getDrxnNodes( !drxn ) )
//            {
//                if ( find( tNodes.begin(), tNodes.end(), bck ) == tNodes.end() ) tAlt.push_back( bck );
//            }
//            Edge e( edge->n, edge->ol, false );
//            SeedBranch alt( e, tAlt, altShared, drxn );
//            if ( *edge > alt )
//            {
//                branches.push_back( edge->n );
//                branchFwd.push_back( edge->n->getDrxnNodes( drxn, false, true ) );
//            }
//        }
//        for ( SeedBranch* edge : edges ) delete edge;
//    }
//    
//    if ( bckSet.find( path.back() ) != bckSet.end() )
//    {
//        branches.clear();
//        branchFwd.clear();
//    }
//    
//    return true;
//}
//
//PathSeed::PathSeed( NodeList &nodes )
//: nodes_( nodes )
//{
//    setEnds();
//}
//
//void PathSeed::exportAlign( ofstream &fp )
//{
//    int32_t base = nodes_[0]->ends_[0];
//    for ( Node* node : nodes_ )
//    {
//        node->id2_ = -1;
//        base = min( base, node->ends_[0] );
//    }
//    
//    int i = 1, k = 1;
//    for ( SeedPath &path : paths_ )
//    {
//        for ( Node* node : path.starts ) if ( node->id2_ < 0 ) node->id2_ = i++;
//        if ( path.path[0]->id2_ < 0 ) path.path[0]->id2_ = i++;
//        path.path[0]->offsetNode( 1 );
//        string seq = path.path[0]->seq_;
//        for ( int j = 0; j < path.path.size()-1; j++ )
//        {
//            for ( Edge &e : path.path[j]->edges_[1] )
//            {
//                if ( e.node->id2_ < 0 )
//                {
//                    path.path[j]->offsetEdge( e, 1 );
//                    e.node->id2_ = i++;
//                }
//                if ( e.node == path.path[j + 1] )
//                {
//                    seq += path.path[j + 1]->seq_.substr( e.ol );
//                }
//                else
//                {
//                    for ( int k = 0; k < path.alts.size(); )
//                    {
//                        if ( e.node == path.alts[k][0] )
//                        {
//                            if ( path.alts[k][0]->id2_ < 0 )
//                            {
//                                path.alts[k][0]->id2_ = i++;
//                            }
//                            NodeSet pathSet( path.alts[k].begin(), path.alts[k].end() );
//                            for ( Node* node : path.alts[k] )
//                            {
//                                for ( Edge &f : node->edges_[1] )
//                                {
//                                    if ( pathSet.find( f.node ) == pathSet.end() ) continue;
//                                    if ( f.node->id2_ < 0 )
//                                    {
//                                        node->offsetEdge( f, 1 );
//                                        f.node->id2_ = i++;
//                                    }
//                                }
//                            }
//                            path.alts.erase( path.alts.begin() + k );
//                        }
//                        else k++;
//                    }
//                }
//            }
//        }
//        
//        NodeSet sets[2] = { path.path[0]->getDrxnNodes( 0 ), path.path.back()->getDrxnNodes( 1 ) };
//        bool dim[2] = { false, false };
//        for ( SeedPath &alt : paths_ )
//        {
//            if ( &path == &alt ) continue;
//            if ( sets[0].find( alt.path[0] ) != sets[0].end() ) dim[0] = true;
//            for ( Node* node : alt.starts ) if ( sets[0].find( node ) != sets[0].end() ) dim[0] = true;
//            if ( sets[1].find( alt.path.back() ) != sets[1].end() ) dim[1] = true;
//        }
//        if ( ( dim[0] || dim[1] ) && seq.length() < 200 ) continue;
//        if ( !dim[0] && !dim[1] ) fp << ">Complete ";
//        else if ( dim[0] && dim[1] ) fp << ">Middle ";
//        else if ( dim[1] ) fp << ">Left ";
//        else fp << ">Right ";
//        fp << to_string( k++ ) << endl;
//        fp << string( path.path[0]->ends_[0] - base, '-' ) + seq << endl;
//    }
//}
//
//void PathSeed::plot( bool drxn )
//{
//    NodeSet bckSet;
//    for ( Node* node : ends_[!drxn] ) paths_.push_back( SeedPath( node, usedSet_, drxn ) );
//    for ( SeedPath &path : paths_ ) path.path.back()->getDrxnNodes( bckSet, !drxn );
//    
//    bool complete = false;
//    while ( !complete )
//    {
//        complete = true;
//        for ( int i = 0; i < paths_.size(); i++ )
//        {
//            if ( paths_[i].ended ) continue;
//            complete = false;
//            NodeSet altFwd, usedSet;
//            while ( SeedPath::resolve( paths_, altFwd, usedSet, i, drxn ) );
//            if ( paths_[i].ended ) continue;
//            NodeList tNodes = paths_[i].getTargets( drxn );
//            while ( paths_[i].setNext( altFwd, usedSet, bckSet, tNodes, drxn ) )
//            {
//                if ( usedSet_.find( paths_[i].path.back() ) != usedSet_.end() )
//                {
//                    paths_[i].ended = true;
//                    paths_[i].path.erase( paths_[i].path.end() - 1 );
//                }
//            }
//        }
//    }
//    
//    NodeSet fwdSet;
//    for ( Node* node : SeedPath::getUsed( paths_, drxn ) ) usedSet_.insert( node );
//    setEnds();
//    for ( Node* node : ends_[!drxn] ) node->getDrxnNodes( fwdSet, drxn );
//    for ( int i = 0; i < ends_[!drxn].size(); )
//    {
//        if ( fwdSet.find( ends_[!drxn][i] ) != fwdSet.end() )
//        {
//            ends_[!drxn].erase( ends_[!drxn].begin() + i );
//        }
//        else i++;
//    }
//    if ( !ends_[!drxn].empty() ) plot( drxn );
//}
//
//void PathSeed::setEnds()
//{
//    for ( bool drxn : { 0, 1 } )
//    {
//        ends_[drxn].clear();
//        vector< pair<Node*, Node*> > forks;
//        for ( Node* node : nodes_ )
//        {
//            if ( usedSet_.find( node ) != usedSet_.end() ) continue;
//            bool isEnd = true;
//            for ( Node* nxt : node->getNextNodes( drxn ) ) if ( usedSet_.find( nxt ) == usedSet_.end() ) isEnd = false;
//            if ( !isEnd ) continue;
//            if ( !usedSet_.empty() )
//            {
//                ends_[drxn].push_back( node );
//                continue;
//            }
//            Node* curr = node;
//            while ( curr->edges_[!drxn].size() == 1 && curr->edges_[drxn].size() < 2 ) curr = node->edges_[!drxn][0].node;
//            forks.push_back( make_pair( node, curr ) );
//        }
//        for ( int i = 0; i < forks.size(); i++ )
//        {
//            bool useFork = false;
//            for ( int j = i + 1; j < forks.size(); )
//            {
//                if ( forks[i].second == forks[j].second )
//                {
//                    useFork = true;
//                    forks.erase( forks.begin() + j );
//                }
//                else j++;
//            }
//            ends_[drxn].push_back( ( useFork ? forks[i].second : forks[i].first ) );
//        }
//    }
//}
