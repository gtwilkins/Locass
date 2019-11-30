/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "query_graph.h"
#include "query_state.h"
#include "constants.h"
#include "shared_functions.h"
#include "node.h"
#include <algorithm>
#include <cassert>

extern Parameters params;

void QueryRead::offset( int off )
{
    ext -= off;
    ol += off;
}

QueryNode::QueryNode( string &ext, int maxOl )
: ext( ext ), len( 0 ), base( 0 ), maxOl( maxOl ), maxExt( 0 ), readCount( 0 ), merged( 0 ), fixed( false ), good( false ), bad( false )
{}

// Base branch
QueryNode::QueryNode( string &q, Overlap &ol, bool perfect )
: QueryNode( ol.seq, ol.overLen )
{
    seq = q.substr( 0, maxOl );
    reads.push_back( QueryRead( ol.readId, ol.overLen, ol.extLen ) );
    if ( perfect ) len = ol.overLen;
}

// Back fork
QueryNode::QueryNode( Overlap &ol, vector<QueryNode*> &hits, bool perfect )
: QueryNode( ol.seq, ol.overLen )
{
    for ( QueryNode* hit : hits ) hit->addEdge( this );
    reads.push_back( QueryRead( ol.readId, ol.overLen, ol.extLen ) );
    if ( perfect ) len = ol.overLen;
}

QueryNode::QueryNode( QueryNode* node, int j )
: QueryNode( node->ext, node->reads[j].ol )
{
    reads.insert( reads.end(), node->reads.begin() + j, node->reads.end() );
    node->reads.erase( node->reads.begin() + j, node->reads.end() );
    for ( QueryNode* e : node->edges[1] )
    {
        addEdge( e );
        e->removeEdge( node, 0 );
    }
    node->edges[1].clear();
    if ( node->len ) len = maxOl;
}

QueryNode::~QueryNode()
{
    for ( int d : { 0, 1 } ) for ( QueryNode* node : edges[d] ) if ( node ) node->removeEdge( this, !d );
    for ( QueryNode* node : edges[1] ) if ( node->edges[0].empty() ) delete node;
}

void QueryNode::addEdge( QueryNode* qn )
{
    edges[1].push_back( qn );
    qn->edges[0].push_back( this );
}

//void QueryNode::cull( vector<QueryNode*> &nodes, int cutoff )
//{
//    for ( int i = 0; i < nodes.size(); i++ )
//    {
//        if ( nodes[i]->confirm() ) continue;
//        nodes[i]->discard();
//        nodes.erase( nodes.begin() + i-- );
//    }
//    
//    int maxReads = 0;
//    for ( QueryNode* node : nodes ) maxReads = max( maxReads, node->readCount );
//    if ( maxReads < 10 ) return;
//    
//    for ( int i = 0; i < nodes.size(); i++ )
//    {
//        if ( nodes[i]->readCount > 1 ) continue;
//        nodes[i]->discard();
//        nodes.erase( nodes.begin() + i-- );
//    }
//    
//    for ( QueryNode* node : nodes ) node->cull( cutoff );
//}
//
//void QueryNode::cull( int cutoff )
//{
////    if ( edges[1].empty() )
////    {
////        uint16_t exts[3]{0}, ols[3]{0};
////        for ( QueryRead &qr : reads )
////        {
////            for ( int i = 0; i < 3; i++ )
////            {
////                if ( qr.ext + ( qr.ol > ols[i] ) <= exts[i] ) continue;
////                for ( int j = 3; --j > i; )
////                {
////                    exts[j] = exts[i];
////                    ols[j] = ols[i];
////                }
////                exts[i] = qr.ext;
////                ols[i] = qr.ol;
////                break;
////            }
////        }
////        if ( exts[2] == exts[0] || !exts[1] ) return;
////        // Cut 2
////        if ( ols[1] < maxOl && exts[2] < exts[1] && ols[1] >= ols[0] ) trim( ols[1], 1 );
////        // Cut 1
////        else if ( ols[0] < maxOl ) trim( ols[0], 1 );
////    }
////    
//    int maxReads = 0;
//    for ( QueryNode* node : edges[1] ) maxReads = max( maxReads, node->readCount );
//    for ( int i = 0; maxReads >= 10 && ext.size() < cutoff && i < edges[1].size(); i++ )
//    {
//        if ( edges[1][i]->readCount > 1 ) continue;
//        if ( !edges[1][i]->removeEdge( this, 0 ) ) edges[1][i]->discard();
//        edges[1].erase( edges[1].begin() + i-- );
//    }
//    for ( QueryNode* node : edges[1] ) node->cull( cutoff );
//    merge();
//    
//    bool blunt = false;
//    for ( QueryNode* node : edges[1] )
//    {
//        if ( node->edges[0].size() > 1 || node->readCount > 3 ) return;
//        if ( node->readCount < 2 ) blunt = true;
//    }
//    if ( !blunt ) return;
//    for ( QueryNode* node : edges[1] )
//    {
//        node->removeEdge( this, 0 );
//        node->discard();
//    }
//    edges[1].clear();
//}

bool QueryNode::confirm()
{
    if ( len ) return true;
    for ( QueryNode* edge : edges[1] ) if ( edge->confirm() ) return true;
    if ( edges[0].empty() ) delete this;
    return false;
}

//void QueryNode::discard()
//{
//    for ( QueryNode* node : edges[1] ) if ( !node->removeEdge( this, 0 ) ) node->discard();
//    delete this;
//}

void QueryNode::extend( Overlap &ol )
{
    assert( ol.seq.size() > ext.size() );
    ext += ol.seq.substr( ext.size() );
    for ( QueryRead &qr : reads ) assert( qr.id != ol.readId );
    reads.push_back( QueryRead( ol.readId, ol.overLen, ol.extLen ) );
}

string QueryNode::getSeq( bool drxn )
{
    assert( merged );
    uint16_t extLen = reads[0].ext;
    for ( int i = 1; i < merged; i++ ) extLen = max( extLen, reads[i].ext );
    return drxn ? ext.substr( base, extLen-base ) : ext.substr( ext.size() - extLen, extLen - base );
}

vector<QueryNode*> QueryNode::graph( QueryBinaries* qb, QState &qs, string seq, bool drxn )
{
    for ( uint8_t i : qs.q ) seq += drxn ? intToChar[i] : intToCharComp[i];
    
    vector<QueryNode*> nodes;
    for ( QState &edge : qs.edges )
    {
        vector<QueryNode*> branches = graph( qb, edge, seq, drxn );
        nodes.insert( nodes.end(), branches.begin(), branches.end() );
    }
    
    for ( int i = qs.ols.size(); --i >= 0; )
    {
        vector<Overlap> ols;
        qb->getOverlaps( ols, qs.ols[i].rank, qs.ols[i].count, qs.ols[i].ol, 1 );
        Overlap::sortByExt( ols );
        for ( Overlap &ol : ols )
        {
            bool added = false;
            vector<QueryNode*> hits;
            for ( QueryNode* qn : nodes ) qn->match( ol, hits, 0, added );
            if ( added ) continue;
            if ( hits.empty() ) nodes.push_back( new QueryNode( seq, ol, qs.perfect ) );
            else if ( hits.size() > 1 || !hits[0]->edges[1].empty() ) new QueryNode( ol, hits, qs.perfect );
            else if ( !hits[0]->len && qs.perfect ) new QueryNode( ol, hits, qs.perfect );
            else hits[0]->extend( ol );
        }
    }
    
    return nodes;
}

bool QueryNode::isBlunt()
{
    if ( edges[0].size() > 1 ) return false;
    for ( QueryNode* node : edges[1] ) if ( !node->isBlunt() ) return false;
    return true;
}

bool QueryNode::match( Overlap &ol, vector<QueryNode*> &hits, int i, bool &added )
{
    if ( find( hits.begin(), hits.end(), this ) != hits.end() ) return true;
    int limit = min( ext.size(), ol.seq.size() );
    while ( i < limit && ext[i] == ol.seq[i] ) i++;
    
    // Perfect fit within this node
    if ( i == ol.seq.size() )
    {
        if ( !reads.empty() && ol.readId == reads.back().id ) return true;
        QueryRead qr( ol.readId, ol.overLen, ol.extLen );
        qr.redundant = i < ext.size() || reads.back().ol > qr.ol;
        for ( int j = 0; j < reads.size() && !qr.redundant; j++ ) qr.redundant = reads[j].ol > qr.ol && reads[j].ext == qr.ext;
        for ( QueryRead &qr : reads ) assert( qr.id != ol.readId );
        reads.push_back( qr );
        added = true;
        return true;
    }
    
    // Extension goes beyond this node
    if ( i == ext.size() )
    {
        bool branched = false;
        for ( QueryNode* edge : edges[1] ) if ( edge->match( ol, hits, i, added ) ) branched = true;
        if ( !branched ) hits.push_back( this );
        return true;
    }
    
    // Mismatch in this node
    uint16_t j = 0, maxLen = 0;
    for ( ; j < reads.size() && reads[j].ext <= i; j++ ) maxLen = max( maxLen, reads[j].ext );
    if ( !j ) return false;
    assert( j < reads.size() );
    addEdge( new QueryNode( this, j ) );
    ext.erase( ext.begin() + maxLen, ext.end() );
    hits.push_back( this );
    return true;
}

bool QueryNode::merge()
{
    for ( int i = 0; i < edges[1].size(); i++ ) if ( edges[1][i]->bad && !edges[1][i]->good ) delete edges[1][i--];
    if ( edges[1].size() != 1 || edges[1][0]->edges[0].size() != 1 ) return false;
    QueryNode* node = edges[1][0];
    if ( bool( len ) != bool( node->len ) ) return false;
    assert( !fixed );
    edges[1].clear();
    ext = node->ext;
    reads.insert( reads.end(), node->reads.begin(), node->reads.end() );
    for ( QueryNode* e : node->edges[1] ) addEdge( e );
    delete node;
    return true;
}

void QueryNode::removeEdge( QueryNode* node, int i )
{
    edges[i].erase( remove( edges[i].begin(), edges[i].end(), node ), edges[i].end() );
}

void QueryNode::resetBase()
{
    base = 0;
    for ( QueryNode* node : edges[1] ) node->resetBase();
}

int QueryNode::setExt()
{
    maxExt = ext.size();
    for ( QueryNode* node : edges[1] ) maxExt = max( maxExt, node->setExt() );
    return maxExt;
}

int QueryNode::setReads( int cutoff )
{
    assert( !reads.empty() );
    if ( readCount ) return readCount;
    
    // Set read count
    for ( QueryRead &qr : reads ) if ( !qr.redundant ) readCount++;
    for ( QueryNode* node : edges[1] ) readCount += node->setReads( cutoff );
    if ( edges[1].empty() ) trim();
    
    int maxReads = 0;
    for ( QueryNode* node : edges[1] ) maxReads = max( maxReads, node->readCount );
    for ( QueryNode* node : edges[1] ) ( maxReads < 10 || node->readCount > 1 || ext.size() >= cutoff ? node->good : node->bad ) = true;
    if ( maxReads < 4 )
    {
        bool blunt = false;
        for ( QueryNode* node : edges[1] ) if ( node->readCount < 2 ) blunt = true;
        if ( blunt ) for ( QueryNode* node : edges[1] ) if ( !node->isBlunt() ) blunt = false;
        while ( blunt && !edges[1].empty() ) delete edges[1].back();
    }
    return readCount;
}

bool QueryNode::setSeq( bool drxn )
{
    if ( fixed ) return false;
    if ( !drxn ) revComp( ext );
    else reverse( seq.begin(), seq.end() );
    if ( !drxn ) for ( QueryRead &qr : reads ) qr.id = params.getRevId( qr.id );
    return fixed = true;
}

//void QueryNode::setNonAlt( vector<QueryNode*> &nodes )
//{
//    if ( !len ) for ( QueryNode* edge : edges[1] ) edge->setNonAlt( nodes );
//    else if ( len < maxOl ) assert( false );
//    else if ( find( nodes.begin(), nodes.end(), this ) == nodes.end() )
//    {
//        for ( QueryNode* node : edges[0] ) if ( node->len ) return;
//        nodes.push_back( this );
//        edges[0].push_back( NULL );
//    }
//}

void QueryNode::trim()
{
    uint16_t exts[3]{0}, ols[3]{0};
    for ( QueryRead &qr : reads )
    {
        for ( int i = 0; i < 3; i++ )
        {
            if ( qr.ext + ( qr.ol > ols[i] ) <= exts[i] ) continue;
            for ( int j = 3; --j > i; )
            {
                exts[j] = exts[i];
                ols[j] = ols[i];
            }
            exts[i] = qr.ext;
            ols[i] = qr.ol;
            break;
        }
    }
    if ( exts[2] == exts[0] || !exts[1] ) return;
    // Cut 2
    if ( ols[1] < maxOl && exts[2] < exts[1] && ols[1] >= ols[0] ) trim( ols[1] );
    // Cut 1
    else if ( ols[0] < maxOl ) trim( ols[0] );
}

void QueryNode::trim( uint16_t ol )
{
//    for ( int i = 0; i < edges[1].size(); i++ )
//    {
//        if ( !edges[1][i]->removeEdge( this, 0 ) ) edges[1][i]->discard();
//        edges[1].erase( edges[1].begin() + i-- );
//    }
//    
    uint16_t extLen = base;
    for ( int i = 1; i < reads.size(); i++ ) if ( reads[i].ol <= ol ) reads.erase( reads.begin()+i, reads.end() );
    for ( int i = 0; i < reads.size(); i++ ) extLen = max( extLen, reads[i].ext );
    assert( !reads.empty() && reads[0].ol > ol );
    ext = ext.substr( 0, extLen );
}

