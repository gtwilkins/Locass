/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "local_alignment.h"
#include <algorithm>
#include <cassert>
#include <iostream>

// Glocal = b can align anywhere along a

LocalAlignment::LocalAlignment( std::string &a, std::string &b, bool glocal, bool freePolymer )
: a_( a ), b_( b ), m_( a.size()+1, std::vector<int>( b.size()+1 ) ), freePolymer_( freePolymer )
{
    bool aSame = freePolymer && a_[0] == b_[0], bSame = freePolymer && a_[0] == b_[0];
    for ( int i = 0; i < a_.size()+1; i++ )
    {
        if ( i && a_[i] != b_[0] ) aSame = false;
        m_[i][0] = glocal || aSame ? 0 : -i;
    }
    for ( int j = 0; j < b_.size()+1; j++ )
    {
        if ( j && a_[0] != b_[j] ) bSame = false;
        m_[0][j] = bSame ? 0 : -j;
    }
}

int LocalAlignment::isGapPoly( std::string (&a)[2], int d, int i, int gap )
{
    assert( a[0].size() == a[1].size() && a[d][i] == '-' && i + gap <= a[0].size() );
    char c = a[!d][i];
    bool n = c == 'N';
    int poly = gap;
    for ( int j = i+1; j < i + gap; j++ )
    {
        if ( a[!d][j] == 'N' ) n = true;
        else if ( c == 'N' ) c = a[!d][j];
        else if ( a[!d][j] != c ) poly = std::min( poly, j - i );
    }
    
    if ( poly < gap && !n )
    {
        int len = 1;
        bool perf = false;
        while ( !perf && ++len <= gap && ( i + gap + len*2 ) < a[0].size() )
        {
            if ( gap % len ) continue;
            perf = true;
            for ( int j = 0; j < len; j++ ) for ( int k = len; k < gap; k += len ) if ( a[!d][i+j] != a[!d][i+j+k] ) perf = false;
        }
        for ( int j = 0; perf && j < len; j++ ) if ( a[!d][i+j] != a[!d][i+gap+j] || a[!d][i+j] != a[!d][i+gap+len+j] ) perf = false;
        if ( perf ) return gap;
    }
    
    for ( int j = i; poly && --j >= 0; )
    {
        if ( a[0][j] == c && a[1][j] == c ) break;
        if ( a[0][j] != 'N' && a[0][j] != c ) poly = 0;
        if ( a[1][j] != 'N' && a[1][j] != c ) poly = 0;
    }
    
    return poly;
}

void LocalAlignment::print( int iMax, int jMax )
{
    if ( a_.size() < iMax ) iMax = a_.size();
    if ( b_.size() < jMax ) jMax = b_.size();
    for ( int i = 0; i <= iMax; i++ )
    {
        for ( int j = 0; j <= jMax; j++ )
        {
            if ( !i && !j ) std::cout << "* ";
            else if ( !i ) std:: cout << b_[j-1] << " ";
            else if ( !j ) std:: cout << a_[i-1] << " ";
            else std:: cout << std::to_string( m_[i][j] ) + " ";
        }
        std::cout << std::endl;
    }
}

void LocalAlignment::realign( std::string &a, std::string &b, bool conform, bool bluntStart, bool bluntEnd, bool trimEnd )
{
    int iMax = 0, jMax = 0, coords[2];
    bool polyMax = true;
    char c;
    for ( int i = 0; i < a_.size(); i++ )
    {
        for ( int j = 0; j < b_.size(); j++ )
        {
            m_[i+1][j+1] = score( i, j, c, coords );
            if ( i < a_.size()-1 && j < b_.size()-1 ) continue;
            if ( !coords[0] || !coords[1] ) continue;
            bool poly = coords[0] != coords[1], closer = abs( i - j ) < abs( iMax - jMax );
            if ( iMax && jMax && m_[iMax+1][jMax+1] > m_[i+1][j+1] ) continue;
            if ( iMax && jMax && m_[iMax+1][jMax+1] == m_[i+1][j+1] && ( poly > polyMax ? : !closer ) ) continue;
            iMax = i;
            jMax = j;
            polyMax = poly;
        }
    }
    
    int i = a_.size()-1, j = b_.size()-1;
    a = std::string( j-jMax, '-' );
    b = std::string( i-iMax, '-' );
    while ( i > iMax ) a += a_[i--];
    while ( j > jMax ) b += b_[j--];
    while ( i >= 0 && j >= 0 )
    {
        score( i, j, c, coords );
        if ( coords[1] > coords[0] ) a += std::string( coords[1] - coords[0], '-' );
        if ( coords[0] > coords[1] ) b += std::string( coords[0] - coords[1], '-' );
        if ( conform && c != 'N' )
        {
            i -= coords[0];
            j -= coords[1];
            a += std::string( coords[0], c );
            b += std::string( coords[1], c );
            continue;
        }
        while ( coords[0]-- ) a += a_[i--];
        while ( coords[1]-- ) b += b_[j--];
    }
    if ( j >= 0 ) a += std::string( j+1, '-' );
    if ( i >= 0 ) b += std::string( i+1, '-' );
    for ( ; i >= 0; i-- ) a += a_[i];
    for ( ; j >= 0; j-- ) b += b_[j];
    std::reverse( a.begin(), a.end() );
    std::reverse( b.begin(), b.end() );
    
    if ( bluntStart && ( a[0] == '-' || b[0] == '-' ) )
    {
        std::string &x = ( a[0] == '-' ? a : b ), &y = ( a[0] == '-' ? b : a );
        int i = 0, len = 1;
        for ( ; len < x.size() && x[len] == '-'; len++ );
        for ( ; i+len < x.size(); i++ )
        {
            if ( x[i+len] != 'N' && y[i] != 'N' && x[i+len] != y[i] ) break;
            x[i] = x[i+len];
            if ( y[i+len] == '-' )
            {
                assert( false );
                x.erase( x.begin() + i+len );
                y.erase( y.begin() + i+len-- );
            }
            else x[i+len] = '-';
        }
    }
    if ( bluntEnd && ( a.back() == '-' || b.back() == '-' ) )
    {
        std::string &x = ( a.back() == '-' ? a : b ), &y = ( a.back() == '-' ? b : a );
        int i = 1, j = 0;
        while ( i < x.size() && x.end()[ -i-1 ] == '-' ) i++;
        while ( j < i && i+j < y.size() && y.end()[ -i-j-1 ] == '-' ) j++;
        if ( j )
        {
            x.erase( x.end()-j, x.end() );
            y.erase( y.end()-i-j, y.end()-i );
        }
    }
    if ( trimEnd && ( a.back() == '-' || b.back() == '-' ) )
    {
        std::string &x = ( a.back() == '-' ? a : b ), &y = ( a.back() == '-' ? b : a );
        int i = 1;
        while ( i < x.size() && x.end()[ -i-1 ] == '-' ) i++;
        x.erase( x.end()-i, x.end() );
        y.erase( y.end()-i, y.end() );
    }
//    if ( bluntEnd && ( a.back() == '-' || b.back() == '-' ) )
//    {
//        std::string &x = ( a.back() == '-' ? a : b ), &y = ( a.back() == '-' ? b : a );
//        for ( int i = x.size()-1; i > 0 && x[i] == '-'; i--) if ( x[i-1] != '-' ) x.erase( x.begin() + i, x.end() );
//        int gap = 0;
//        for ( int i = x.size()-1; i >= 0 && y[i] == '-'; i-- ) gap++;
//        gap = std::min( gap, (int)y.size() - (int)x.size() );
//        if ( gap > 0 ) y.erase( y.begin() + (int)x.size() - gap, y.begin() + (int)x.size() );
//        if ( y.size() > x.size() ) y.erase( y.begin() + x.size(), y.end() );
//    }
}

int LocalAlignment::score( int i, int j, char &c, int* coords )
{
    int best = m_[i][j] + ( a_[i] == 'N' || b_[j] == 'N' ? 0 : ( a_[i] == b_[j] ? 1 : -1 ) );
    coords[0] = coords[1] = 1;
    c = a_[i] == b_[j] ? a_[i] : ( a_[i] == 'N' ? b_[j] : ( b_[j] == 'N' ? a_[i] : 'N' ) );
    
    int s = m_[i][j+1] - 1;
    if ( s > best )
    {
        best = s;
        coords[0] = 1;
        coords[1] = 0;
    }
    
    s = m_[i+1][j] - 1;
    if ( s > best )
    {
        best = s;
        coords[0] = 0;
        coords[1] = 1;
    }
    
    if ( c != 'N' && ( i+1 == a_.size() || a_[i+1] != c ) && ( j+1 == b_.size() || b_[j+1] != c ) )
    {
        int lens[2] = { 1, 1 }, ns[2]{ a_[i] == 'N', b_[j] == 'N' };
        while ( true )
        {
            for ( int &k = lens[0]; k <= i && a_[i-k] == c; k++ );
            for ( int &k = lens[1]; k <= j && b_[j-k] == c; k++ );
            s = m_[ i+1-lens[0] ][ j+1-lens[1] ] + std::min( lens[0] - ns[0], lens[1] - ns[1] );
            if ( lens[0] != lens[1] && s > best )
            {
                best = s;
                coords[0] = lens[0];
                coords[1] = lens[1];
            }
            if ( lens[0] <= i && a_[ i-lens[0] ] == 'N' ) { lens[0]++; ns[0]++; }
            else if ( lens[1] <= j && b_[ j-lens[1] ] == 'N' ){ lens[1]++; ns[1]++; lens[0] = 1; ns[0] = a_[i] == 'N'; }
            else break;
        }
    }
    
    return best;
}

//void LocalAlignment::setAlign( std::string &a, std::string &b )
//{
//    int i = a_.size()-1, j = b_.size()-1, run[2];
//    while ( i > i_ )
//    {
//        a += a_[i--];
//        b += '-';
//    }
//    while ( j > j_ )
//    {
//        a += '-';
//        b += b_[j--];
//    }
//    
//    while ( i >= 0 && j >= 0 )
//    {
//        if ( m_[i+1][j+1] == m_[i][j] + ( a_[i] == b_[j] ? 1 : ( ( a_[i] == 'N' || b_[j] == 'N' ) ? 0 : -1 ) ) )
//        {
//            a += a_[i--];
//            b += b_[j--];
//        }
//        else
//        {
//            int tmpRuns[2]{0};
//            if ( freePolymer_ )
//            {
//                setRuns( run, i, j );
//                tmpRuns[0] = run[0];
//                tmpRuns[1] = run[1];
//                if ( run[0] && m_[i+1][j+1] != m_[ i+1-run[0] ][j+1] ) run[0] = 0;
//                if ( run[1] && m_[i+1][j+1] != m_[i+1][ j+1-run[1] ] ) run[1] = 0;
//            }
//            
//            if ( !run[0] && !run[1] && m_[i+1][j+1] == m_[i][j+1] - 1 ) run[0] = 1;
//            if ( !run[0] && !run[1] && m_[i+1][j+1] == m_[i+1][j] - 1 ) run[1] = 1;
//            if ( run[0] && run[1] ) run[1] = 0;
//            
//            
//            while( run[0]-- )
//            {
//                a += a_[i--];
//                b += '-';
//            }
//            while( run[1]-- )
//            {
//                a += '-';
//                b += b_[j--];
//            }
//        }
//    }
//    
//    while ( i >= 0 )
//    {
//        a += a_[i--];
//        b += '-';
//    }
//    while ( j >= 0 )
//    {
//        a += '-';
//        b += b_[j--];
//    }
//    
//    std::reverse( a.begin(), a.end() );
//    std::reverse( b.begin(), b.end() );
//}

//int LocalAlignment::setCoords( int* coords )
//{
//    int i = i_, j = b_.size()-1, run[2];
//    
//    while ( i >= 0 && j >= 0 )
//    {
//        if ( m_[i+1][j+1] == m_[i][j] + ( a_[i] == b_[j] ? 1 : ( ( a_[i] == 'N' || b_[j] == 'N' ) ? 0 : -1 ) ) )
//        {
//            i--;
//            j--;
//        }
//        else
//        {
//            if ( a_[i] != b_[j] ) setRuns( run, i, j );
//            else run[0] = run[1] = 0;
//            
//            if ( run[0] && m_[i+1][j+1] != m_[ i+1-run[0] ][j+1] ) run[0] = 0;
//            if ( run[1] && m_[i+1][j+1] != m_[i+1][ j+1-run[1] ] ) run[1] = 0;
//            if ( !run[0] && !run[1] && m_[i+1][j+1] == m_[i][j+1] - 1 ) run[0] = 1;
//            if ( !run[0] && !run[1] && m_[i+1][j+1] == m_[i+1][j] - 1 ) run[1] = 1;
//            
//            assert( ( run[0] || run[1] ) && !( run[0] && run[1] ) );
//            
//            while( run[0]-- ) i--;
//            while( run[1]-- ) j--;
//        }
//    }
//    
//    coords[0] = i+1;
//    coords[1] = i_+1;
//    
//    return m_[i_+1][j_+1];
//}

void LocalAlignment::setRuns( int* run, int i, int j )
{
    run[0] = run[1] = 0;
    char c = a_[i+1] == 'N' ? b_[j+1] : a_[i+1];
    if ( c == 'N' ) return;
    if ( a_[i+1] != b_[j+1] && a_[i+1] != 'N' && b_[j+1] != 'N' ) return;
    
    for ( int k = 1; k <= i && ( a_[i-k+1] == c || a_[i-k+1] == 'N' ); k++ )
    {
        if ( a_[i-k] != c && m_[i+1-k][j+1] >= m_[ i+1-run[0] ][j+1] ) run[0] = k;
    }
    for ( int k = 1; k <= j && ( b_[j-k+1] == c || b_[j-k+1] == 'N' ); k++ )
    {
        if ( b_[j-k] != c && m_[i+1][j+1-k] >= m_[i+1][ j+1-run[1] ] ) run[1] = k;
    }
}