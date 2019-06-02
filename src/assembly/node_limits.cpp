/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "node_limits.h"
#include <string.h>
#include <limits>
#include <cassert>

NodeLimits::NodeLimits()
{
    clear();
}
    
int32_t& NodeLimits::operator[]( int i )
{
    return ends[i];
}

void NodeLimits::clear()
{
    for ( int i = 0; i < 3; i++ ) limits[0][i] = numeric_limits<int32_t>::max();
    for ( int i = 0; i < 3; i++ ) limits[1][i] = numeric_limits<int32_t>::min();
    origin[0] = origin[1] = 0;
}

void NodeLimits::inherit( int32_t (&alt)[2][3] )
{
    memcpy( limits, alt, sizeof( alt ) );
    recoil();
}

void NodeLimits::init( int32_t limit )
{
    init( limit, 0 );
    init( limit, 1 );
}

void NodeLimits::init( int32_t limit, bool drxn )
{
    for ( int i = 0; i < 3; i++ )
    {
        if ( drxn ) limits[1][i] = max( limits[1][i], limit );
        else limits[0][i] = min( limits[0][i], limit );
    }
}

bool NodeLimits::isValidMark( int32_t mark, int fwdHits, bool drxn )
{
    if ( fwdHits > 2 ) return true;
    return drxn ? mark <= limits[drxn][2-fwdHits] : limits[drxn][2-fwdHits] <= mark;
}

void NodeLimits::merge( NodeLimits &alt )
{
    ends[0] = min( ends[0], alt[0] );
    ends[1] = max( ends[1], alt[1] );
    for ( int i = 0; i < 3; i++ )
    {
        push( alt.limits[0][i], 0 );
        push( alt.limits[1][i], 1 );
    }
    recoil();
}

void NodeLimits::offset( int32_t off )
{
    ends[0] += off;
    ends[1] += off;
    for ( int i = 0; i < 2; i++ ) for ( int j = 0; j < 3; j++ ) limits[i][j] += off;
    if ( origin[0] != origin[1] )
    {
        origin[0] += off;
        origin[1] += off;
    }
}

void NodeLimits::push( int32_t limit, bool drxn )
{
    for ( int i = 3; --i >= 0 && ( drxn ? limits[1][i] < limit : limit < limits[0][i] ); )
    {
        if ( !i || ( drxn ? limit <= limits[1][i-1] : limits[0][i-1] <= limit ) )
        {
            limits[drxn][i] = limit;
            return;
        }
        limits[drxn][i] = limits[drxn][i-1];
    }
}

void NodeLimits::recoil()
{
    for ( int i = 0; i < 3; i++ )
    {
        limits[0][i] = max( limits[0][i], ends[0] );
        limits[1][i] = min( limits[1][i], ends[1] );
    }
}

void NodeLimits::reset( int drxn )
{
    assert( drxn < 2 || origin[0] != origin[1] );
    for ( int i = 0 ; i < 3; i++ )
    {
        if ( drxn < 2 ) limits[0][i] = limits[1][i] = ends[!drxn];
        else for ( int d : { 0, 1 } ) limits[d][i] = origin[d];
    }
}

void NodeLimits::splitOrigin( NodeLimits &split, int &curDrxn, int &splitDrxn )
{
    if ( origin[0] == origin[1] ) return;
    if ( origin[1] > ends[1] )
    {
        split.origin[1] = origin[1];
        split.origin[0] = max( origin[0], split[0] );
        if ( split[0] < origin[0] )
        {
            origin[0] = origin[1] = 0;
            curDrxn = 0;
        }
        else origin[1] = ends[1];
        return;
    }
    if ( origin[0] < ends[0] )
    {
        split.origin[0] = origin[0];
        split.origin[1] = min( origin[1], split[1] );
        if ( origin[1] < split[1] )
        {
            origin[0] = origin[1] = 0;
            curDrxn = 1;
        }
        else origin[0] = ends[0];
        return;
    }
    splitDrxn = split[1] > ends[1];
}

bool NodeLimits::targetable( bool drxn )
{
    return targetable( ends[!drxn], drxn );
}

bool NodeLimits::targetable( int32_t limit, bool drxn )
{
    return drxn ? limit < limits[1][2] : limits[0][2] < limit;
}

bool NodeLimits::verified()
{
    return limits[0][2] <= ends[0] && ends[1] <= limits[1][2];
}

bool NodeLimits::verified( bool drxn )
{
    return drxn ? ends[1] <= limits[1][2] : limits[0][2] <= ends[0];
}