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

#include "calibrate_structs.h"
#include <cassert>
#include <iostream>
#include <algorithm>

SeedLibraryCount::SeedLibraryCount( Parameters &inParams )
: params( inParams ), ratioScore( 0 )
{
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        ffSum.push_back( 0 );
        frSum.push_back( 0 );
        rfSum.push_back( 0 );
        ffCounts.push_back( 0 );
        frCounts.push_back( 0 );
        rfCounts.push_back( 0 );
        libCounts.push_back( 0 );
    }
}

void SeedLibraryCount::add( ReadId id, int dist, bool leftFwd, bool rightFwd )
{
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        if ( id < params.libs[i].endCount )
        {
            if ( leftFwd && rightFwd )
            {
                ffSum[i] += dist;
                ffCounts[i]++;
            }
            else if ( leftFwd && !rightFwd )
            {
                frSum[i] += dist;
                frCounts[i]++;
            }
            else if ( !leftFwd && rightFwd )
            {
                rfSum[i] += dist;
                rfCounts[i]++;
            }
            break;
        }
    }
}

vector<int> SeedLibraryCounts::set( Parameters &params )
{
    for ( SeedLibraryCount &lib : libs )
    {
        ReadId libsTotal = 0;
        for ( ReadId &libCount : lib.libCounts )
        {
            libsTotal += libCount;
        }
        for ( int i ( 0 ); i < lib.libCounts.size(); i++ )
        {
            float thisRatio = (float)lib.libCounts[i] / (float)libsTotal;
            lib.ratios.push_back( thisRatio );
            libRatios[i].push_back( thisRatio );
        }
    }
    
    for ( vector<float> &libRatio : libRatios )
    {
        sort( libRatio.begin(), libRatio.end() );
    }
    vector<float> medianRatios;
    for ( vector<float> &libRatio : libRatios )
    {
        medianRatios.push_back( libRatio[ libRatio.size() / 2 ] );
    }
    
    vector<float> ratioScoreList;
    for ( SeedLibraryCount &lib : libs )
    {
        for ( int i ( 0 ); i < lib.ratios.size(); i++ )
        {
            lib.ratioScore += sqrt( abs( lib.ratios[i] - medianRatios[i] ) );
            ratioScoreList.push_back( lib.ratioScore );
        }
    }
    sort( ratioScoreList.begin(), ratioScoreList.end() );
    float ratioCutoff = ratioScoreList[ 1 + ratioScoreList.size() / 2 ];
    
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        ffSum.push_back( 0 );
        frSum.push_back( 0 );
        rfSum.push_back( 0 );
        ffCounts.push_back( 0 );
        frCounts.push_back( 0 );
        rfCounts.push_back( 0 );
    }
    
    vector<int> goodIds;
    for ( SeedLibraryCount &lib : libs )
    {
        if ( lib.ratioScore <= ratioCutoff )
        {
            goodIds.push_back( lib.id );
            for ( int i ( 0 ); i < params.libs.size(); i++ )
            {
                ffSum[i] += lib.ffSum[i];
                frSum[i] += lib.frSum[i];
                rfSum[i] += lib.rfSum[i];
                ffCounts[i] += lib.ffCounts[i];
                frCounts[i] += lib.frCounts[i];
                rfCounts[i] += lib.rfCounts[i];
            }
        }
    }
    
    int iBest = -1;
    CharId bestCount = 0;
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        CharId libCount = max( ffCounts[i], max( frCounts[i], rfCounts[i] ) );
        if ( libCount > bestCount )
        {
            iBest = i;
            bestCount = libCount;
        }
    }
    
    bool anyGood = false, anyBad = false;
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        CharId totalCount = ffCounts[i] + frCounts[i] + rfCounts[i];
        CharId cutoff = totalCount * 0.8;
        if ( totalCount < 100 ) continue;
        
        if ( ffCounts[i] > cutoff )
        {
            params.libs[i].size = (float)ffSum[i] / (float)ffCounts[i];
            params.libs[i].orientation = 1;
        }
        else if ( frCounts[i] > cutoff )
        {
            params.libs[i].size = (float)frSum[i] / (float)frCounts[i];
            params.libs[i].orientation = 2;
        }
        else if ( rfCounts[i] > cutoff )
        {
            params.libs[i].size = (float)rfSum[i] / (float)rfCounts[i];
            params.libs[i].orientation = 3;
        }
        else anyBad = true;
        
        if ( params.libs[i].size )
        {
            params.libs[i].isPe = params.libs[i].size < 1000;
            params.libs[i].setMinMax();
            anyGood = true;
        }
    }
    
    if ( !anyGood )
    {
        if ( anyBad ) cerr << endl << "Error: Read pairs appear to be randomly oriented, which locass does not support." << endl;
        else cerr << endl << "Error: Could not find sufficient number of pairs." << endl;
        exit( EXIT_FAILURE );
    }
    
//    if ( bestCount > 100 )
//    {
//        CharId totalCount = ffCounts[iBest] + frCounts[iBest] + rfCounts[iBest];
//        CharId cutoff = totalCount * 0.8;
//        
//        if ( ffCounts[iBest] > cutoff )
//        {
//            params.libs[iBest].size = (float)ffSum[iBest] / (float)ffCounts[iBest];
//            params.libs[iBest].orientation = 1;
//        }
//        else if ( frCounts[iBest] > cutoff )
//        {
//            params.libs[iBest].size = (float)frSum[iBest] / (float)frCounts[iBest];
//            params.libs[iBest].orientation = 2;
//        }
//        else if ( rfCounts[iBest] > cutoff )
//        {
//            params.libs[iBest].size = (float)rfSum[iBest] / (float)rfCounts[iBest];
//            params.libs[iBest].orientation = 3;
//        }
//        else
//        {
//            cerr << endl << "Error: Read pairs appear to be randomly oriented, which locass does not support." << endl;
//            exit( EXIT_FAILURE );
//        }
//        
//        params.libs[iBest].setMinMax();
//        params.libs[iBest].isPe = true;
//    }
//    else
//    {
//        cerr << endl << "Error: Could not find sufficient number of pairs." << endl;
//        exit( EXIT_FAILURE );
//    }
    
    params.set();
    
    return goodIds;
}

int32_t LocusLibraryCount::getMedian( int i, bool isFf, bool isFr, bool isRf )
{
    vector<int32_t>* whichVec = ( isFf ? &ffPairs[i] : ( isFr ? &frPairs[i] : &rfPairs[i] ) );
    
    if ( !whichVec->empty() )
    {
        return (*whichVec)[ whichVec->size() / 2 ];
    }
    
    return 0;
}

void LocusLibraryCount::setMedians( Parameters &params )
{
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        sort( ffPairs[i].begin(), ffPairs[i].end() );
        sort( frPairs[i].begin(), frPairs[i].end() );
        sort( rfPairs[i].begin(), rfPairs[i].end() );
    }
}

void LocusLibraryCounts::set( Parameters &params )
{
    for ( int i ( 0 ); i < params.libs.size(); i++ )
    {
        CharId ffTotal = 0, frTotal = 0, rfTotal = 0;
        CharId ffCounts = 0, frCounts = 0, rfCounts = 0;
        for ( LocusLibraryCount &count : libs )
        {
            count.setMedians( params );
            ffCounts += count.ffPairs[i].size();
            for ( const int32_t &ffPair : count.ffPairs[i] )
            {
                ffTotal += ffPair - min( ffPair, params.readLen );
            }
            frCounts += count.frPairs[i].size();
            for ( const int32_t &frPair : count.frPairs[i] )
            {
                frTotal += frPair - min( frPair, params.readLen );
            }
            rfCounts += count.rfPairs[i].size();
            for ( const int32_t &rfPair : count.rfPairs[i] )
            {
                rfTotal += rfPair - min( rfPair, params.readLen );
            }
        }

        bool isFf = ffTotal > frTotal && ffTotal > rfTotal;
        bool isFr = frTotal > ffTotal && frTotal > rfTotal;
        bool isRf = rfTotal > ffTotal && rfTotal > frTotal;
        if ( !isFf && !isFr && !isRf ) continue;
        int32_t maxMedian = 0;
        int32_t maxLen = 0;
        for ( LocusLibraryCount &count : libs )
        {
            maxMedian = max( maxMedian, count.getMedian( i, isFf, isFr, isRf ) );
            maxLen = max( maxLen, count.lens[i] );
        }
        int32_t lenCutoff = min( maxLen, int32_t(maxMedian * 1.5) );

        double divisor = 0;
        double medianSum = 0;
        vector<int32_t> minPairs, maxPairs;
        int limitsCutoff = max( (CharId)1, ( isFf ? ffCounts : ( isFr ? frCounts : rfCounts ) ) / 40 );
        for ( LocusLibraryCount &count : libs )
        {
            vector<int32_t>* whichVec = ( isFf ? &count.ffPairs[i] : ( isFr ? &count.frPairs[i] : &count.rfPairs[i] ) );
            minPairs.insert( minPairs.end(), whichVec->begin(), whichVec->begin() + min( (int)whichVec->size(), limitsCutoff ) );
            maxPairs.insert( maxPairs.end(), whichVec->end() - min( (int)whichVec->size(), limitsCutoff ), whichVec->end() );
            if ( count.lens[i] >= lenCutoff )
            {
                divisor += whichVec->size();
                medianSum += whichVec->size() * (double)count.getMedian( i, isFf, isFr, isRf );
            }
        }
        sort( minPairs.begin(), minPairs.end() );
        sort( maxPairs.begin(), maxPairs.end() );
        params.libs[i].orientation = ( isFf ? 1 : ( isFr ? 2 : 3 ) );
        params.libs[i].size = medianSum / divisor;
        params.libs[i].minDist = min( params.libs[i].size, minPairs[ min( (int)maxPairs.size() - 1, limitsCutoff * 2 ) ] ) * 0.8;
        params.libs[i].maxDist = max( params.libs[i].size, maxPairs[ max( 0, (int)maxPairs.size() - limitsCutoff ) ] ) * 1.2;
    }
    
    vector< pair<float, int> > coverages;
    int totalLen = 0;
    for ( LocusLibraryCount &count : libs )
    {
        if ( count.len )
        {
            coverages.push_back( make_pair( count.coverage, count.len ) );
            totalLen += count.len;
        }
    }
    
    sort( coverages.begin(), coverages.end(), []( pair<float, int> &a, pair<float, int> &b ){
        return a.first < b.first;
    });
    
    int cutoffLen = totalLen * 0.67;
    for ( int i ( 0 ); i < coverages.size(); i++ )
    {
        if ( cutoffLen < coverages[i].second )
        {
            params.cover = coverages[i].first;
            break;
        }
        cutoffLen -= coverages[i].second;
    }
    
    params.set();
}

