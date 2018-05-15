/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   path_reassembly.h
 * Author: glen
 *
 * Created on 18 December 2017, 2:43 PM
 */

#ifndef PATH_REASSEMBLY_H
#define PATH_REASSEMBLY_H

#include "types.h"
#include "node.h"
#include "shared_functions.h"
#include "path_sequence.h"

class Reassemble
{
public:
    Reassemble( Node* node, PathVars &pv, bool invalid, bool calibrate );
    ~Reassemble();
    
    bool reassemble( PathVars &pv, NodeSet &delSet, bool isAlleleFork=false );
    
private:
    bool doesSpanOrigin( PathVars &pv );
    void map( PathVars &pv );
    void removeDubious( PathVars &pv );
    void removeRedundant( vector<ReadEndMap*> &reads );
    void setWeakspot();
    bool tryBridge( PathVars &pv );
    bool tryComplete( PathVars &pv );
    bool tryExact( PathVars &pv );
    bool tryGap( PathVars &pv );
    bool tryHalf( PathVars &pv, NodeSet &delSet, bool drxn );
    bool tryMap( PathVars &pv );
    bool trySlice( PathVars &pv, NodeSet &delSet, bool isAlleleFork );
    
    NodeIntMap offsets_;
    Node* fork_;
    vector<SeqPathReassemble*> seqs_;
    int32_t est_, minCover_;
    int32_t estLimits_[2], markLimits_[2], good_[2], high_[2], limits_[2];
    bool allHigh_, allGood_, invalid_, calibrate_;
};

#endif /* PATH_REASSEMBLY_H */

