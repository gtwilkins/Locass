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

#include "node.h"
#include <algorithm>

void Node::addPairs( NodeIntMap &pairs, unordered_set<SeqNum> &hitIds, bool drxn )
{
    if ( pairs.empty() ) return;
    
    NodeSet bckSet = getDrxnNodes( !drxn, false, params.getFurthestMpDist( ends_[!drxn], !drxn ) );
    
    // Update pairs and remove markers from pair if applicable
    for ( const pair<Node*, int> &np : pairs )
    {
        // Set or update hits in this
        auto result = pairs_.insert( np );
        if ( !result.second ) result.first->second += np.second;
        
        // Set hits in pair
        np.first->pairs_[this] = result.first->second;
        
        // Remove markers from pair
        np.first->removeMarks( hitIds, true, true, !drxn );
        
        // Set valid limits in between
        np.first->pushValidLimts( bckSet, np.second, drxn );
        if ( np.first != this )
        {
            for ( int i( 0 ); i < max( np.second, 2 ); i++ )
            {
                pushValidLimits( ends_[!drxn], !drxn );
                np.first->pushValidLimits( np.first->ends_[drxn], drxn );
            }
        }
        
        if ( np.first->drxn_ != drxn )
        {
            np.first->farPairNodes_[0] = np.first->farPairNodes_[1] = NULL;
        }
    }
    
    // Remove markers from this
    farPairNodes_[0] = farPairNodes_[1] = NULL;
    removeMarks( hitIds, true, false, drxn );
}

void Node::clearPairs()
{
    NodeList tNodes;
    for ( pair<Node*, int> np : pairs_ )
    {
        if ( np.first != this )
        {
            tNodes.push_back( np.first );
            np.first->pairs_.erase( this );
            np.first->resetFurthest( this );
        }
        
        if ( np.first->paired_ )
        {
            np.first->paired_->erase( this );
        }
        
        if ( np.first->farPairNodes_[0] == this )
        {
            np.first->farPairNodes_[0] = NULL;
        }
        
        if ( np.first->farPairNodes_[1] == this )
        {
            np.first->farPairNodes_[1] = NULL;
        }
    }
    
    // Reset pairs in this node
    pairs_.clear();
    if ( paired_ )
    {
        paired_->clear();
    }
    
    farPairNodes_[0] = farPairNodes_[1] = NULL;
    
    if ( !tNodes.empty() )
    {
        vector<SeqNum> pairIds = getPairsIdsBase();
        for ( Node* t : tNodes )
        {
            t->reAddMarks( pairIds );
        }
    }
}

void Node::clearPairsForward( bool drxn )
{
    NodeSet fwdSet = getDrxnNodes( drxn, false, true );
    
    for ( Node* fwd : fwdSet )
    {
        NodeList tNodes;
        for ( auto &np : pairs_ )
        {
            if ( fwdSet.find( np.first ) == fwdSet.end() )
            {
                tNodes.push_back( np.first );
            }
        }
        
        if ( !tNodes.empty() )
        {
            vector<SeqNum> fwdReads;
            vector<ReadMark> marks = getMarksBase( drxn );
            for ( Node* t : tNodes )
            {
                vector<SeqNum> tReads;
                fwd->pairs_.erase( t );
                fwd->resetFurthest( t );
                t->pairs_.erase( fwd );
                t->resetFurthest( fwd );
                for ( ReadMark &mark : marks )
                {
                    if ( t->reads_.find( mark.id ) != t->reads_.end() )
                    {
                        tReads.push_back( mark.id );
                        fwdReads.push_back( params.getPairId( mark.id ) );
                    }
                }
                t->reAddMarks( tReads );
            }
            fwd->reAddMarks( fwdReads );
        }
        
        if ( fwd->paired_ )
        {
            delete fwd->paired_;
        }
        fwd->paired_ = new NodeSet();
        fwd->validated_ = false;
        fwd->farPairNodes_[0] = fwd->farPairNodes_[1] = NULL;
    }
}

PairHit Node::findReadPair( ReadMark &mark, bool &valid, bool drxn )
{
    auto hit = reads_.find( mark.id );
    if ( hit != reads_.end() )
    {
        valid = mark.isValid( hit->second );
        return PairHit( this, &hit->second );
    }
    return PairHit( NULL, NULL );
}

PairHit Node::findReadPairBestClone( NodeList &tNodes, ReadMark &mark, bool &valid, bool drxn )
{
    auto hit = reads_.find( mark.id );
    if ( hit != reads_.end() )
    {
        int32_t diff = abs( mark.estimate - hit->second[!drxn] );
        valid = mark.isValid( hit->second );
        Node* bestNode = this;
        Coords* bestCoords = &hit->second;
        if ( clones_ )
        {
            for ( Node* clone : *clones_ )
            {
                hit = clone->reads_.find( mark.id );
                if ( find( tNodes.begin(), tNodes.end(), clone ) != tNodes.end() && hit != clone->reads_.end() )
                {
                    bool thisValid = mark.isValid( hit->second );
                    int32_t thisDiff =  abs( mark.estimate - hit->second[!drxn] );
                    if ( thisValid && ( !valid || thisDiff < diff ) )
                    {
                        bestNode = clone;
                        bestCoords = &hit->second;
                        diff = thisDiff;
                        valid = thisValid;
                    }
                }
            }
        }
        return PairHit( bestNode, bestCoords, diff );
    }
    return PairHit( NULL, NULL, 0 );
}

Score Node::getBranchScore( int32_t* limits, bool drxn )
{
    Score score;
    NodeSet fwdSet = { this }, tSet;
    NodeList tNodes;
    getDrxnNodes( fwdSet, drxn );
    getDrxnNodes( tSet, !drxn, params.getFurthestMpDist( ends_[!drxn], !drxn ) );
    for ( Node* t : tSet )
    {
        if ( t->validated_ )
        {
            tNodes.push_back( t );
        }
    }
    
    for ( Node* fwd : fwdSet )
    {
        Score fwdScore = fwd->getPairScore( tNodes, limits, drxn );
        fwd->limitScore( fwdScore, false );
        score += fwdScore;
    }
    return score;
}

int Node::getBridgeCount()
{
    int bridgeCount = 0;
    for ( auto &np : pairs_ )
    {
        if ( ( drxn_ > 2 && np.first->drxn_ <= 2 ) || ( drxn_ <= 2 && np.first->drxn_ > 2 ) )
        {
            bridgeCount += np.second;
        }
    }
    return bridgeCount;
}

int Node::getFurthestAndReliable( NodeList &tNodes, int32_t &furthest, bool drxn )
{
    int reliable = 0;
    NodeList hitNodes;
    for ( Node* t : tNodes )
    {
        auto hit = pairs_.find( t );
        if ( hit != pairs_.end() )
        {
            reliable += t->reliable_ ? hit->second : 0;
            if ( drxn ? furthest < t->ends_[1] : t->ends_[0] < furthest )
            {
                hitNodes.push_back( t );
            }
        }
    }
    int32_t dummy;
    getFurthestPair( furthest, dummy, hitNodes, drxn );
    return reliable;
}

int32_t Node::getFurthestPair( bool drxn )
{
    int32_t furthest = ends_[drxn];
    getFurthestPair( furthest, drxn );
    return furthest;
}

void Node::getFurthestPair( int32_t &furthest, bool drxn )
{
    // Get all paired nodes that are further than furthest
    NodeList tNodes;
    for ( pair<Node*, int> pn : pairs_ )
    {
        if ( drxn ? pn.first->ends_[1] > max( furthest, ends_[1] ) : pn.first->ends_[0] < min( furthest, ends_[0] ) )
        {
            tNodes.push_back( pn.first );
        }
    }
    
    int32_t dummy;
    // Update furthest read pair from paired nodes
    getFurthestPair( furthest, dummy, tNodes, drxn );
}

void Node::getFurthestPair( int32_t &furthest, int32_t &furthestReliable, NodeList &tNodes, bool drxn, bool getReliable )
{
    if ( !tNodes.empty() )
    {
        bool valid, cloneValid;
        for ( ReadMark &mark : getMarksBase( !drxn ) )
        {
            for ( Node* t : tNodes )
            {
                PairHit r = t->findReadPairBestClone( tNodes, mark, valid, !drxn );
                if ( clones_ && r.node && valid )
                {
                    for ( Node* clone : *clones_ )
                    {
                        ReadMark cloneMark = mark;
                        cloneMark.offset( clone->ends_[drxn] - ends_[drxn] );
                        PairHit cloneR = t->findReadPairBestClone( tNodes, cloneMark, cloneValid, !drxn );
                        if ( cloneR.diff < r.diff && cloneValid && cloneR.node )
                        {
                            r.node = NULL;
                        }
                    }
                }

                if ( r.node && valid )
                {
                    if ( getReliable && r.node->isReliable() )
                    {
                        furthestReliable = getFurthest( (*r.coords)[drxn], furthestReliable, drxn );
                    }
                    furthest = getFurthest( (*r.coords)[drxn], furthest, drxn );
                }
            }
        }
    }
}

void Node::getFurthestPairs( int32_t* furthest, NodeSet &qSet, CloneScoreMap &score, bool drxn )
{
    furthest[0] = ends_[drxn];
    furthest[1] = ends_[drxn];
    
    for ( Node* q : qSet )
    {
        auto it = score.find( q );
        if ( it != score.end() )
        {
            if ( drxn )
            {
                if ( it->second.furthest[0] < furthest[0] )
                {
                    furthest[1] = furthest[0];
                    furthest[0] = it->second.furthest[0];
                    furthest[1] = min( furthest[1], it->second.furthest[1] );
                }
                else
                {
                    furthest[1] = min( furthest[1], it->second.furthest[0] );
                }
            }
            else
            {
                if ( it->second.furthest[0] > furthest[0] )
                {
                    furthest[1] = furthest[0];
                    furthest[0] = it->second.furthest[0];
                    furthest[1] = max( furthest[1], it->second.furthest[1] );
                }
                else
                {
                    furthest[1] = max( furthest[1], it->second.furthest[0] );
                }
            }
        }
    }
}

void Node::getFurthestReliable( int32_t &furthest, int32_t &furthestReliable, NodeList &tNodes, bool drxn )
{
    NodeList tNodesTmp;
    for ( Node* t : tNodes )
    {
        if ( t->isFurther( furthestReliable, drxn, drxn ) && pairs_.find( t ) != pairs_.end() )
        {
            tNodesTmp.push_back( t );
        }
    }
    getFurthestPair( furthest, furthestReliable, tNodes, drxn, true );
}

float Node::getMissScore( int32_t* limits, bool drxn )
{
    float missScore = 0;
    for ( ReadMark &mark : marks_[drxn] )
    {
        mark.getMissScore( missScore, limits );
    }
    return missScore;
}

int Node::getPairHits( Node* node )
{
    auto hit = pairs_.find( node );
    if ( hit != pairs_.end() )
    {
        return hit->second;
    }
    return 0;
}

int Node::getPairsHits( NodeList &tNodes )
{
    int hits = 0;
    for ( Node* t : tNodes )
    {
        hits += getPairHits( t );
    }
    return hits;
}

int Node::getPairHitsTmp( bool drxn )
{
    NodeSet tSet = getNextNodes( !drxn );
    for ( Node* nxt : getNextNodes( !drxn ) )
    {
        nxt->getDrxnNodes( tSet, !drxn, params.getFurthestMpDist( nxt->ends_[drxn], !drxn ) );
    }
    NodeList tNodes( tSet.begin(), tSet.end() );
    
    int hits = 0;
    for ( ReadMark &mark : getMarksBase( drxn ) )
    {
        for ( Node* t : tNodes )
        {
            if ( t->reads_.find( mark.id ) != t->reads_.end() )
            {
                hits++;
                break;
            }
        }
    }
    
    return hits;
}

int Node::getPairHitsTotal()
{
    int hits = 0;
    for ( auto &pn : pairs_ )
    {
        hits += pn.second;
    }
    return hits;
}

vector<SeqNum> Node::getPairsIdsBase()
{
    vector<SeqNum> pairIds;
    for ( auto &read : reads_ )
    {
        SeqNum pairId = params.getPairId( read.first );
        if ( reads_.find( pairId ) == reads_.end() )
        {
            pairIds.push_back( pairId );
        }
    }
    return pairIds;
}

int32_t Node::getLengthForScoring( float hits )
{
    int32_t maxLen = getLengthForCoverage();
    int32_t minLen = max( params.readLen, params.maxPeMean - params.readLen );
    int32_t adjLen = maxLen;
    if ( maxLen > minLen && hits > 0 )
    {
        float coverage = max( min( params.cover, coverage_ ), params.cover / 2 );
        float minHits = float( minLen * params.cover / 4 ) / float(params.readLen);
        float maxHits = float( maxLen * coverage / 2 ) / float(params.readLen);
        adjLen = max( minLen * ( minHits / hits ), maxLen * ( hits / maxHits ));
    }
    return min( maxLen, adjLen );
}

Score Node::getPairScore( NodeList &tNodes, int32_t* limits, bool drxn )
{
    Score score;
    getPairScore( score, tNodes, limits, drxn );
    return score;
}

void Node::getPairScore( Score &score, NodeList &tNodes, int32_t* limits, bool drxn )
{
    score[1] += getPairsHits( tNodes );
    score[0] += getMissScore( limits, drxn );
}

void Node::getReliablePairNodes( NodeList &nodes, NodeSet &reliableSet, bool drxn )
{
    for ( Node* node : nodes )
    {
        for ( const pair<Node*, int> &np : node->pairs_ )
        {
            if ( np.first->reliable_ 
                    && ( drxn ? np.first->ends_[1] > node->ends_[1] : np.first->ends_[0] < node->ends_[0] )
                    && find( nodes.begin(), nodes.end(), np.first ) == nodes.end() )
            {
                reliableSet.insert( np.first );
            }
        }
    }
}

ScoreMap Node::getScoreMap( NodeSet &qSet, int32_t* limits, bool drxn )
{
    ScoreMap scores;
    for ( Node* q : qSet )
    {
        Score score;
        score[1] = q->getPairHitsTotal();
        score[0] = q->getMissScore( limits, drxn );
        q->limitScore( score, true );
        scores[q] = score;
    }
    return scores;
}

void Node::getUnpairedMarks( NodeList &nodes, vector<ReadMark> &peMarks, vector<ReadMark> &mpMarks, unordered_set<SeqNum> &usedIds, unordered_set<SeqNum> &readIds, int32_t* limits, bool inclMp, bool drxn )
{
    for ( ReadMark &mark : getMarksBase( !drxn ) )
    {
        bool isPe = params.isReadPe( mark.id );
        bool found = ( !inclMp && !isPe )
                || ( usedIds.find( mark.id ) != usedIds.end() ) 
                || ( mark.estimate < limits[0] )
                || ( limits[1] < mark.estimate );
        
        if ( !found )
        {
            for ( Node* node : nodes )
            {
                if ( node->reads_.find( mark.id ) != node->reads_.end() )
                {
                    found = true;
                    break;
                }
            }
        }
        if ( !found )
        {
            readIds.insert( mark.id );
            ( isPe ? peMarks : mpMarks ).push_back( mark );
        }
    }
}

void Node::limitScore( Score &score, bool inclNotValidSelf )
{
    int32_t len = getLengthForScoring( score.hits );
    float coverage = max( min( params.cover, coverage_ ), params.cover / float(2) );
    float expected[2];
    expected[0] = ( ( coverage * params.peRatio ) / (float)params.readLen ) * len * 0.5;
    expected[1] = ( coverage / (float)params.readLen ) * len * 0.5;
    
    if ( inclNotValidSelf && !validated_ )
    {
        auto it = pairs_.find( this );
        if ( it != pairs_.end() )
        {
            float expectedBranch = ( params.peCover / 8 ) * float(min( len, params.maxPeMean - params.readLen )) / (float)params.readLen;
            score.hits += ( score.hits / expectedBranch ) * it->second;
        }
    }
    
    if ( score.hits && score.hits < expected[0] )
    {
        expected[0] *= max( score.hits / expected[0], (expected[0] - score.hits) / expected[0] );
    }
    
    score.misses = max( float( 0 ), expected[0] - score.hits );
    score.hits = score.hits <= expected[1] ? score.hits : max( 1, int( expected[1] ) );
}

void Node::updatePairs()
{
    for ( bool drxn : { 0, 1 } )
    {
        if ( drxn && ( drxn_ == 0 || drxn_ == 4 ) ) continue;
        if ( !drxn && ( drxn_ == 1 || drxn_ == 3 ) ) continue;
        if ( paired_ ) paired_->clear();

        NodeList tNodes = getTargetNodes( drxn, true );
        setPairs( tNodes, drxn );

        for ( Node* fwd : getDrxnNodes( drxn ) )
        {
            NodeList tSelf = { this };
            fwd->setPairs( tSelf, drxn );
        }
    }
}

void Node::resetPairing( NodeSet &nodes )
{
    for ( Node* node : nodes )
    {
        if ( node->paired_ )
        {
            delete node->paired_;
        }
        node->paired_ = new NodeSet();
        node->pairs_.clear();
        node->pullValidLimits();
        node->resetMarks();
        node->validated_ = false;
    }
}

CloneScore Node::setCloneScore( PairingVars &pv, CloneTargetVars &ctv, bool drxn )
{
    CloneScore score;
    score.furthest[0] = ends_[!drxn];
    score.furthest[1] = ends_[!drxn];
    
    for ( ReadMark &mark : *pv.marks )
    {
        bool valid = false;
        for ( Node* t : pv.tNodes )
        {
            PairHit r = t->findReadPairBestClone( pv.tNodes, mark, valid, drxn );
            if ( r.node )
            {
                bool isBest = true;
                bool isGood = valid && ctv.tSetThis.find( r.node ) != ctv.tSetThis.end();
                bool doRecord = isGood && r.node->isValidHit( r.coords, drxn );
                bool altBest = !valid;
                bool altGood = false;
                
                for ( int32_t &off : ctv.offsets )
                {
                    bool cloneValid = false;
                    mark.offset( off );
                    PairHit r2 = t->findReadPairBestClone( pv.tNodes, mark, cloneValid, drxn );
                    mark.offset( -off );
                    altBest = altBest || ( cloneValid && r2.diff < r.diff );
                    isBest = isBest && ( !cloneValid || r2.diff > r.diff );
                    altGood = altGood || ( cloneValid && ctv.tSetAlt.find( r2.node ) != ctv.tSetAlt.end() );
                }
                
                score.selfPref += isBest && isGood;
                score.selfOnly += isGood && !altGood;
                score.altPref += altBest && altGood;
                score.altOnly += !valid && altGood;
                
                if ( isGood && !altBest )
                {
                    if ( doRecord )
                    {
                        auto result = pv.pairs.insert( make_pair( r.node, 1 ) );
                        if ( result.second == false )
                        {
                            result.first->second++;
                        }
                    }
                    score.selfHits++;
                    if ( isBest )
                    {
                        score.addHit( r.coords, drxn );
                    }
                }
                
                if ( ( isGood || altGood ) && doRecord )
                {
                    pv.hitIds.insert( mark.id );
                }
                
                break;
            }
        }
    }
    return score;
}

CloneScore Node::setCloneScore( NodeIntMap &pairs, NodeList &tNodes, vector<int32_t> offsets, vector<ReadMark> &marks, unordered_set<SeqNum> &hitIds, bool drxn, bool checkValid )
{
    CloneScore score;
    score.furthest[0] = ends_[!drxn];
    score.furthest[1] = ends_[!drxn];
    for ( ReadMark mark : marks )
    {
        bool valid = false;
        for ( Node* t : tNodes )
        {
            PairHit r = t->findReadPairBestClone( tNodes, mark, valid, drxn );
            if ( r.node )
            {
                bool isBest = true, isGood = valid && ( !checkValid || r.node->isValidHit( r.coords, drxn ) );
                bool altBest = !valid, altGood = false;
                
                for ( int32_t off : offsets )
                {
                    bool cloneValid = false;
                    mark.offset( off );
                    PairHit r2 = t->findReadPairBestClone( tNodes, mark, cloneValid, drxn );
                    mark.offset( -off );
                    altBest = altBest || ( cloneValid && r2.diff > r.diff );
                    isBest = isBest && ( !cloneValid || r2.diff < r.diff );
                    altGood = altGood || ( cloneValid );
                }
                
                score.selfPref += isBest && isGood;
                score.selfOnly += isGood && !altGood;
                score.altPref += altBest && altGood;
                score.altOnly += !valid && altGood;
                
                if ( isGood && !altBest )
                {
                    auto result = pairs.insert( make_pair( r.node, 1 ) );
                    if ( result.second == false )
                    {
                        result.first->second++;
                    }
                    score.selfHits++;
                    if ( isBest )
                    {
                        score.addHit( r.coords, drxn );
                    }
                }
                
                if ( isGood || altGood )
                {
                    hitIds.insert( mark.id );
                }
                
                break;
            }
        }
    }
    return score;
}

int Node::setPairs( NodeList &tNodes, bool drxn )
{
    PairingVars pv;
    pv.marks = &marks_[drxn];
    
    NodeSet pairedSet;
    if ( paired_ )
    {
        pairedSet = *paired_;
    }
    
    // Get untested target nodes
    for ( Node* t : tNodes )
    {
        if ( ( !paired_ || paired_->find( t ) == paired_->end() ) )
        {
            pv.tNodes.push_back( t );
            if ( paired_ && t->validated_ )
            {
                paired_->insert( t );
            }
        }
    }
    
    if ( pv.tNodes.empty() )
    {
        return 0;
    }
    
    // Set pairs
    if ( clones_ )
    {
        CloneTargetVars ctv = getCloneTargetVars( drxn );
        setCloneScore( pv, ctv, drxn );
    }
    else
    {
        setPairs( pv, drxn );
    }
    
    addPairs( pv.pairs, pv.hitIds, drxn );
    
    if ( pv.pairs.find( this ) != pv.pairs.end() )
    {
        NodeList tSelf = { this };
        int i = setPairs( ( drxn_ == 2 ? pv.tNodes : tSelf ), drxn );
    }
    
    return pv.hitIds.size();
}

void Node::setPairs( PairingVars &pv, bool drxn )
{
    bool valid;
    Node* bestHit = NULL;
    for ( ReadMark &mark : *pv.marks )
    {
        if ( drxn_ != 2 || ( drxn ? validLimits_[1] <= mark.mark - params.readLen
                                  : mark.mark + params.readLen <= validLimits_[2] ) )
        {
            bestHit = NULL;
            for ( Node* t : pv.tNodes )
            {
                // Find pair in target node
                if ( t->clones_ )
                {
                    PairHit r = t->findReadPairBestClone( pv.tNodes, mark, valid, drxn );
                    if ( r.node && valid && ( validated_ || drxn_ == !r.node->drxn_ || r.node->isValidHit( r.coords, drxn ) ) ) bestHit = r.node;
                }
                else
                {
                    PairHit r = t->findReadPair( mark, valid, drxn );
                    if ( r.node && valid && ( validated_ || drxn_ == !r.node->drxn_ || r.node->isValidHit( r.coords, drxn ) ) ) bestHit = r.node;
                }
            }

            // Record pair hit
            if ( bestHit )
            {
                auto result = pv.pairs.insert( make_pair( bestHit, 1 ) );
                if ( result.second == false )
                {
                    result.first->second++;
                }
                pv.hitIds.insert( mark.id );
            }
        }
    }
}

void Node::setPairs( NodeOffsetMap &fwdMap, NodeOffsetMap &revMap, bool drxn )
{
    for ( auto &fwd : fwdMap )
    {
        for ( auto &rev : revMap )
        {
            fwd.first->setPairs( this, rev.first, fwd.second, rev.second, drxn );
        }
    }
    
    pair<int32_t, int32_t> offset = make_pair( 0, 0 );
    int32_t validLimits[2] = { numeric_limits<int32_t>::max(), numeric_limits<int32_t>::max() };
    while( validLimits[0] != validLimits_[1] || 
           validLimits[1] != validLimits_[2] )
    {
        validLimits[0] = validLimits_[1];
        validLimits[1] = validLimits_[2];
        setPairs( this, this, offset, offset, drxn );
        
        for ( auto &fwd : fwdMap )
        {
            fwd.first->setPairs( this, this, fwd.second, offset, drxn );
        }
        for ( auto &rev : revMap )
        {
            setPairs( this, rev.first, offset, rev.second, drxn );
        }
    }
    
    for ( auto &fwd : fwdMap )
    {
        for ( auto &rev : revMap )
        {
            if ( fwd.first->paired_ )
            {
                fwd.first->paired_->insert( rev.first );
            }
            if ( rev.first->paired_ )
            {
                rev.first->paired_->insert( fwd.first );
            }
        }
    }
}

void Node::setPairs( Node* focus, Node* t, pair<int32_t, int32_t> &qOffset, pair<int32_t, int32_t> &tOffset, bool drxn )
{
    if ( !paired_ || paired_->find( t ) == paired_->end() )
    {
        unordered_set<SeqNum> usedIds;
        Coords* coords;
        
        for ( ReadMark &mark : marks_[drxn] )
        {
            if ( this != focus || ( drxn ? validLimits_[0] + params.readLen <= mark.mark 
                                         : mark.mark <= validLimits_[3] - params.readLen ) )
            {
                if ( t->findRead( mark.id, coords, true ) )
                {
                    if ( ( t != focus || t->isValidHit( coords, drxn ) )
                            && mark.isValid( *coords, qOffset, tOffset, drxn ) )
                    {
                        pushValidLimits( t, mark.mark, (*coords)[!drxn], drxn );

                        auto r1 = pairs_.insert( make_pair( t, 1 ) );
                        if ( !r1.second )
                        {
                            r1.first->second++;
                        }
                        auto r2 = t->pairs_.insert( make_pair( this, 1 ) );
                        if ( !r2.second )
                        {
                            r2.first->second++;
                        }
                        
                        usedIds.insert( mark.id );
                    }
                }
            }
        }
        
        removeMarks( usedIds, false, false, drxn );
        t->removeMarks( usedIds, false, true, !drxn );
    }
}

