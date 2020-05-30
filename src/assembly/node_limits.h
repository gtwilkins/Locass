/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   node_limits.h
 * Author: glen
 *
 * Created on 31 December 2018, 3:47 AM
 */

#ifndef NODE_LIMITS_H
#define NODE_LIMITS_H

#include "types.h"

struct NodeLimits
{
    NodeLimits();
    int32_t& operator[]( int i );
    void clear();
    void inherit( int32_t (&alt)[2][3] );
    void init( int32_t limit );
    void init( int32_t limit, bool drxn );
    bool isValidMark( int32_t mark, int fwdHits, bool drxn );
    void merge( NodeLimits &alt );
    void offset( int32_t off );
    void push( int32_t limit, bool drxn );
    void recoil();
    void reset( int drxn );
    void splitOrigin( NodeLimits &split, int &curDrxn, int &splitDrxn );
    bool targetable( bool drxn );
    bool targetable( int32_t limit, bool drxn );
    bool verified();
    bool verified( bool drxn );
    int32_t ends[2], limits[2][3], origin[2];
};

#endif /* NODE_LIMITS_H */

