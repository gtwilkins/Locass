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

#ifndef QUERY_H
#define QUERY_H

#include "types.h"
#include "constants.h"
#include "index_reader.h"
#include "query_binary.h"
#include "query_state.h"
#include "query_structs.h"

class Querier {
public:
    Querier( IndexReader* ir, QueryBinaries* qb );
    ~Querier();
    
    string getSequence( ReadId id );
    vector<Extension> mapExtensions( bool &noMatches, string &seq, bool drxn, uint16_t minOver=1 );
    vector<Extension> mapExtensions( string &seq, bool drxn, uint16_t minOver=1 );
    vector<Extension> mapExtensions( string &seq, bool drxn, unordered_set<SeqNum> &seeds, uint16_t minOver=1 );
    vector<Overlap> mapJoin( string seq1, string seq2, uint8_t overLen );
    MappedSeqs mapSeed( string &seq, int errorRate, bool bestMatch );
    void mapSequence( string &seq, vector<ReadId> &ids, vector<int32_t>* coords );
   
private:
    vector<Overlap> getOverlaps( string &seq, uint16_t minOver, bool drxn );
    vector<Overlap> getOverlaps( string &seq, uint16_t minOver, uint8_t &maxConsidered, bool drxn );
    vector<Extension> compileExtensions( vector<Overlap> &overlaps, bool drxn, bool doTrim );
    bool isVector( Extension &ext );
    void mapExtensionsCull( string &seq, vector<Extension> &exts );
    void mapExtensionsCull( string &seq, vector<Extension> &exts, vector<Overlap> &overlaps, unordered_set<SeqNum> &seeds );
    void mapExtensionsTrim( vector<Extension> &exts, vector<Overlap> &overlaps );
    void mapReads( QueryState &q, uint8_t it, CharId rank, CharId count );
    bool mapSeed( QuerySeedState &qs, int errors, bool drxn );
    void setQuery( string &seq, uint8_t *query, int length, bool drxn );
    
    IndexReader* reader_;
    QueryBinaries* bin_;
    
    uint16_t minOver_, maxSeqs_, constCutoff_;
    
    CharId baseCounts[4];
    CharId baseOverlapRanks_[4][4], baseOverlapCounts_[4][4];
    
//    LetterCountEachPile countsCumulative_;
//    LetterCountEachPile countsPerPile_;
//    LetterCountEachPile baseRanges_;
};

#endif /* QUERY_H */

