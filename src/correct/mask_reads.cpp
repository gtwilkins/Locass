/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <fstream>

#include "mask_reads.h"
#include "error.h"
#include "timer.h"
#include <cassert>
#include <iostream>

//MaskReads::MaskReads( string fn, bool input )
//: input( input )
//{
//    fns[0] = fn;
//    size_t it = fn.find_last_of( '.' );
//    if ( it != fn.npos ) fn = fn.substr( 0, it );
//    // TEMP
//    {
//        it = fn.find_last_of( '/' );
//        fn = fn.substr( it );
//        fn = "/media/glen/ssd/Hp" + fn;
//    }
//    fns[1] = fn + "_masked.fastq";
//    fns[2] = fn + "_masked.bin";
//    fns[3] = fn + "_tmp1.bin";
//    fns[4] = fn + "_tmp2.bin";
//    
//    bf = new BitFile( fns[4], 2 );
//    
//    fpBin = NULL;
//    ifstream ifsReads( fns[0] );
//    if ( !ifsReads.good() ) failure( "Could not open reads file: \"" + fns[0] + "\"" );
//    ifsReads.close();
//    ifstream ifsTest( fns[1] );
//    if ( ifsTest.good() ) failure( "Intended output filename \"" + fns[1] + "\" already exists." );
//    ifsTest.close();
//    
//    readCount = kmersTotal = readCounted[0] = posBin[0] = posBit[0] = 0;
//    filled[0] = filled[1] = false;
//}
//
//void MaskReads::create( bool overwrite )
//{
//    posBin[0] = 16;
//    if ( !bf->open( false, true, overwrite ) )
//    {
//        fpBin = fopen( fns[3].c_str(), "rb" );
//        if ( !fpBin ) failure( "Incomplete temporary files detected." );
//        fread( &readCount, 8, 1, fpBin );
//        fread( &kmersTotal, 8, 1, fpBin );
//        fclose( fpBin );
//        cout << "Found temporary files for file with " << to_string( readCount ) << " reads and " << to_string( kmersTotal ) << " k-mers." << endl;
//        return;
//    }
//    
//    double startTime = clock();
//    uint8_t charTable[4][256];
//    for ( int i = 0; i < 4; i++ )
//    {
//        for ( int j = 0; j < 256; j++ ) charTable[i][j] = 0;
//        charTable[i]['C'] = charTable[i]['c'] = 1 << ( ( 3 - i ) * 2 );
//        charTable[i]['G'] = charTable[i]['g'] = 2 << ( ( 3 - i ) * 2 );
//        charTable[i]['T'] = charTable[i]['t'] = 3 << ( ( 3 - i ) * 2 );
//    }
//    
//    fpBin = fopen( fns[3].c_str(), "wb" );
//    fwrite( &readCount, 8, 1, fpBin );
//    fwrite( &kmersTotal, 8, 1, fpBin );
//    string lines[4];
//    ifstream ifsReads( fns[0] );
//    while ( getline( ifsReads, lines[0] ) && !lines[0].empty() )
//    {
//        if ( lines[0][0] != '@' ) failure( "Unrecoginsed fastq header: " + lines[0][0] );
//        for ( int i = 1; i < 4; i++ ) getline( ifsReads, lines[i] );
//        if ( lines[1].size() != lines[3].size() ) failure( "Mismatch between sequence length and phred score length" );
//        int lastBad = -1;
//        uint8_t len = lines[1].size(), lastGood = lines[1].size(), byte = 0;
//        while ( lastGood > 0 && ( lines[1][lastGood-1] == 'N' || lines[3][lastGood-1] < 38 ) ) lastGood--;
//        fwrite( &len, 1, 1, fpBin );
//        fwrite( &lastGood, 1, 1, fpBin );
//        int i = 0;
//        for ( int j = 0; j < len; j++ )
//        {
//            byte += charTable[i][ lines[1][j] ];
//            if ( ++i == 4 )
//            {
//                fwrite( &byte, 1, 1, fpBin );
//                byte = 0;
//                i = 0;
//            }
//            if ( lines[1][j] == 'N' ) lastBad = j;
//            if ( j < 31 ) continue;
//            bf->write( j - lastBad > 31 ? 0 : 3 );
//            kmersTotal++;
//        }
//        if ( i ) fwrite( &byte, 1, 1, fpBin );
//        readCount++;
//    }
//    
//    bf->close();
//    fclose( fpBin );
//    fpBin = fopen( fns[3].c_str(), "rb+" );
//    fwrite( &readCount, 8, 1, fpBin );
//    fwrite( &kmersTotal, 8, 1, fpBin );
//    fclose( fpBin );
//    
//    cout << "Created temporary files for file with " << to_string( readCount ) << " reads and " << to_string( kmersTotal ) << " k-mers in " << getDuration( startTime ) << endl;
//}
//
//void MaskReads::complete()
//{
//    ifstream ifsReads( fns[0] );
//    ofstream ofsReads( fns[1] );
//    bf->open( true, false );
//    BitFile obf( fns[2], 1 );
//    obf.open( false, true );
//    string lines[4];
//    while ( getline( ifsReads, lines[0] ) && !lines[0].empty() )
//    {
//        for ( int i = 1; i < 4; i++ ) getline( ifsReads, lines[i] );
//        int lastGood = lines[1].size(), interval[2]{0};
//        while ( lastGood > 0 && ( lines[1][lastGood-1] == 'N' || lines[3][lastGood-1] < 38 ) ) lastGood--;
//        
//        for ( int i = 31; i < lines[1].size(); i++ )
//        {
//            uint8_t status = bf->read();
//            assert( status );
//            obf.write( status == 2 );
//            if ( status > 2 ) continue;
//            if ( interval[1] == i ) interval[1]++;
//            else if ( interval[1] > lastGood ) break;
//            else
//            {
//                interval[0] = i - 31;
//                interval[1] = i + 1;
//            }
//        }
//        
//        int unique = 1;
//        for ( int i = interval[0]+1; i < min( interval[1], lastGood ); i++ ) if ( lines[1][i] != lines[1][i-1] ) unique++;
//        if ( lastGood - interval[0] < 32 && unique < 10 ) interval[1] = lastGood;
//        for ( int i = lastGood; i < interval[1]; i++ ) lines[3][i] = 43;
//        for ( int i = max( lastGood, interval[1] ); i < lines[1].size(); i++ ) lines[3][i] = 0;
//        for ( int i = 0; i < 4; i++ ) ofsReads << lines[i] + "\n";
//    }
//}
//
//void MaskReads::mask( vector<string> &fnInputs, vector<string> &fnAudits, int minKmers )
//{
//    double startTime = clock();
//    
//    vector<MaskReads> mrs;
//    for ( string &fn : fnInputs ) mrs.push_back( MaskReads( fn, true ) );
//    for ( string &fn : fnAudits ) mrs.push_back( MaskReads( fn, false ) );
//    
//    cout << "Creating temporary files..." << endl;
//    uint64_t kmersLeft, kmersTotal = 0;
//    for ( MaskReads &mr : mrs ) mr.create( true );
//    for ( MaskReads &mr : mrs ) kmersTotal += mr.kmersTotal;
//    cout << endl << "Commencing sequence parsing cycles..." << endl;
//    
//    int filled[2]{0};
//    unordered_map<uint64_t, uint32_t>* kmers = new unordered_map<uint64_t, uint32_t>( 500 * 1000 * 1000 );
//    bool full = true;
//    while ( full )
//    {
//        full = false;
//        double iterTime = clock();
//        kmersLeft = 0;
//        for ( int i = filled[0]; i < mrs.size(); i++ )
//        {
//            if ( mrs[i].parse( kmers, minKmers, true, full ) ) filled[1] = i;
//            kmersLeft += mrs[i].kmersLeft;
//        }
//        assert( mrs[ filled[1] ].input );
//        
//        for ( int i = filled[0]; i < mrs.size() && mrs[i].input; i++ )
//        {
//            mrs[i].parse( kmers, minKmers, false, full );
//        }
//        kmers->clear();
//        filled[0] = filled[1];
//        cout << "Cycle complete in " << getDuration( iterTime ) << ". " << to_string( kmersLeft ) << " k-mers still unresolved of " << to_string( kmersTotal ) << " total." << ( full ? " Commencing next cycle..." : "" ) << endl;
//    }
//    
//    delete kmers;
//    
//    cout << endl << "Masking complete after " << getDuration( startTime ) << "." << endl;
//}
//
//bool MaskReads::parse( unordered_map<uint64_t, uint32_t>* kmers, int minKmer, bool add, bool &full )
//{
//    assert( input || add );
//    uint64_t kmerLimit = kmers->bucket_count() / 5, kmer, fKmer = 0, rKmer = 0;
//    bf->open( true, true );
//    fpBin = fopen( fns[3].c_str(), "rb+" );
//    if ( add && filled[1] )
//    {
//        posBin[0] = posBin[1];
//        posBit[0] = posBit[1];
//        readCounted[0] = readCounted[1];
//        iBit[0] = iBit[1];
//        filled[0] = true;
//        filled[1] = false;
//    }
//    if ( filled[0] ) bf->setPos( posBit[0], iBit[0] );
//    fseek( fpBin, posBin[0], SEEK_SET );
//    kmersLeft = 0;
//    
//    uint64_t fTable[4][256], rTable[4][256];
//    for ( int i = 0; i < 4; i++ )
//    {
//        for ( int j = 0; j < 256; j++ )
//        {
//            fTable[i][j] = ( j >> ( ( 3 - i ) * 2 ) ) & 0x3;
//            rTable[i][j] = ( 3 - fTable[i][j] ) << 62;
//        }
//    }
//    
//    uint8_t len, lastGood, byte, seq[1024];
//    for ( uint64_t i = readCounted[0]; i < readCount; i++ )
//    {
//        assert( fread( &len, 1, 1, fpBin ) == 1 );
//        assert( fread( &lastGood, 1, 1, fpBin ) == 1 );
//        
//        if ( len ) fread( seq, 1, ( 1 + ( len - 1 ) / 4 ), fpBin );
//        if ( len < 32 ) continue;
//        int limit = len - 31, last = -31;
//        for ( int j = 0; j < limit; j++ )
//        {
//            uint8_t status = bf->read();
//            if ( add && status == 1 ) assert( false );
//            if ( add ? status : status != 1 ) continue;
//            last = max( last, j );
//            for ( ; last < j + 32; last++ )
//            {
//                fKmer = ( fKmer << 2 ) + fTable[last%4][ seq[last/4] ];
//                rKmer = ( rKmer >> 2 ) + rTable[last%4][ seq[last/4] ];
//            }
//            kmer = min( fKmer, rKmer );
//            bool good = add && input && last <= lastGood;
//            auto it = kmers->find( kmer );
//            int edit = !full;
//            if ( it != kmers->end() )
//            {
//                if ( it->second >= minKmer ) edit = 2;
//                else if ( good ) it->second++;
//                else if ( !input && !full ) edit = 3;
//            }
//            else if ( input ? status : !full ) edit = 3;
//            else if ( input && add && !full ) kmers->insert( make_pair( kmer, good ) );
//            
//            if ( edit ) bf->edit( edit );
//            if ( !edit ) kmersLeft++;
//        }
//        
//        if ( !full && add && kmers->size() >= kmerLimit )
//        {
//            posBin[1] = ftell( fpBin );
//            bf->tellPos( posBit[1], iBit[1] );
//            readCounted[1] = i + 1;
//            full = filled[1] = true;
//        }
//    }
//    bf->close();
//    fclose( fpBin );
//    
//    return filled[1];
//}
//
