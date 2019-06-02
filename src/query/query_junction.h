/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   query_junction.h
 * Author: glen
 *
 * Created on 27 December 2018, 6:59 PM
 */

#ifndef QUERY_JUNCTION_H
#define QUERY_JUNCTION_H

#include "types.h"
#include "constants.h"
#include "index_reader.h"
#include "query_structs.h"
#include "query_state.h"
#include "query_graph.h"

class QueryJunction
{
public:
    QueryJunction( IndexReader* ir, QueryBinaries* qb, string &seq, ReadId* counts, int* ols, int cutoff, bool drxn );
    ~QueryJunction();
//    void setAlts();
    vector<QueryNode*> nodes_, alts_;
    int failure_;
    
private:
    bool add( QueryNode* node, bool cull, bool first );
    uint8_t getChar( int ol );
    void query( IndexReader* ir, QState &qs, ReadId* counts, int* ols );
    string &seq_;
    bool drxn_, branch_;
};


#endif /* QUERY_JUNCTION_H */

