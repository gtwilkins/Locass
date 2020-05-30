/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "correct_kmers.h"
#include "timer.h"
#include <string.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

//MarksFile::MarksFile( string fn, string mode )
//{
//    out = mode[0] == 'w';
//    assert( !mode.empty() && ( mode[0] == 'w' || mode[0] == 'r' ) );
//    fp = fopen( fn.c_str(), mode.c_str() );
//    assert( fp );
//    i = p = 0;
//    for ( int j = 0; j < 8; j++ ) x[j] = 1 << ( 7-j );
//    if ( out ) memset( b, 0, 4096 );
//    else fread( b, 1, 4096, fp );
//}
//
//MarksFile::~MarksFile()
//{
//    if ( out && ( i || p ) ) fwrite( b, 1, i+1, fp );
//    fclose( fp );
//}
//
//void MarksFile::advance( size_t posAdv, int pAdv )
//{
//    pos = posAdv;
//    fseek( fp, pos, SEEK_SET );
//    fread( b, 1, 4096, fp );
//    i = 0;
//    p = pAdv;
//}
//
//void MarksFile::edit( bool mark )
//{
//    if ( mark ) assert( !( b[i] & x[p] ) );
//    if ( mark ) b[i] += x[p];
//    increment();
//    if ( i < 4096 ) return;
//    fseek( fp, pos, SEEK_SET );
//    fwrite( b, 1, 4096, fp );
//    pos = ftell( fp );
//    fread( b, 1, 4096, fp );
//    i = 0;
//}
//
//void MarksFile::editSkip( int count )
//{
//    p += count;
//    while ( p > 7 )
//    {
//        p -= 8;
//        if ( ++i < 4096 ) continue;
//        fseek( fp, pos, SEEK_SET );
//        fwrite( b, 1, 4096, fp );
//        pos = ftell( fp );
//        fread( b, 1, 4096, fp );
//        i = 0;
//    }
//}
//
//vector< pair<int,int> > MarksFile::getIntervals( int len )
//{
//    vector< pair<int,int> > ils;
//    int first = 0;
//    len -= 31;
//    for ( int i = 0; i <= len; i++ )
//    {
//        if ( i < len && read() ) continue;
//        if ( first < i ) ils.push_back( make_pair( first, i+31 ) );
//        first = i+1;
//    }
//    
//    return ils;
//}
//
//void MarksFile::increment()
//{
//    if ( p == 7 ){ p = 0; i++; }
//    else p++;
//}
//
//bool MarksFile::read()
//{
//    bool r = b[i] & x[p];
//    increment();
//    if ( i == 4096 )
//    {
//        fread( b, 1, 4096, fp );
//        i = 0;
//    }
//    return r;
//}
//
//void MarksFile::write( bool mark )
//{
//    if ( mark ) b[i] += x[p];
//    increment();
//    if ( i == 4096 )
//    {
//        fwrite( b, 1, 4096, fp );
//        i = 0;
//        memset( b, 0, 4096 );
//    }
//}
//
//KmerFileStruct::KmerFileStruct( string &fn, int worth, bool seed )
//: kmers( 2048, 0 ), worth( worth ), seed( seed )
//{
//    string fnBase = fn.substr( 0, fn.find_last_of( '.' ) );
//    fnBase = "/media/glen/ssd/Lv/" + fnBase.substr( fnBase.find_last_of( '/' ) + 1 );
//    fnSeq = fnBase + "_seqs.bin";
//    fnMark = fnBase + "_marks.bin";
//    ifstream checks[3] = { ifstream( fn ), ifstream( fnSeq ), ifstream( fnMark ) };
//    assert( checks[0].good() || ( checks[1].good() && checks[2].good() ) );
//    assert( ( checks[1].good() && checks[2].good() ) || ( !checks[1].good() && !checks[2].good() ) );
//    CharId kmerCount = 0;
//    charCount = 0;
//    seqCount = 0;
//    fp = NULL;
//    
//    for ( int i = 0; i < 4; i++ )
//    {
//        rChar[i] = CharId(3-i) << 62;
//        for ( int j = 0; j < 256; j++ )
//        {
//            rKmer[i][j] = ( j >> ( ( 3 - i ) * 2 ) ) & 0x3;
//        }
//    }
//    
//    if ( checks[1] && checks[2] ) return;
//    
//    ifstream ifs( fn );
//    fp = fopen( fnSeq.c_str(), "wb" );
//
//    string lines[2];
//    getline( ifs, lines[0] );
//    int mode = lines[0][0] == '@' ? 3 : ( lines[0][0] == '>' ? 1 : 0 );
//    ifs.seekg( 0, ifs.beg );
//    
//    while ( getline( ifs, lines[0] ) )
//    {
//        for ( int i = 0; i < mode; i++ ) getline( ifs, lines[bool(i)] );
//        int first = 0;
//        vector< pair<uint16_t,uint16_t> > ils;
//        for ( int i = 0; i <= lines[0].size(); i++ )
//        {
//            if ( i < lines[0].size() && lines[0][i] != 'N' ) continue;
//            if ( i - first > 31 ) ils.push_back( make_pair( first, i ) );
//            first = i+1;
//        }
//
//        uint16_t seqLen = lines[0].size();
//        kmerCount += seqLen < 32 ? 0 : seqLen - 31;
//        charCount += seqLen;
//        uint8_t ilCount = ils.size();
//        fwrite( &seqLen, 2, 1, fp );
//        fwrite( &ilCount, 1, 1, fp );
//
//        for ( pair<uint16_t, uint16_t> &il : ils )
//        {
//            int writeLen = 1 + ( ( il.second - il.first - 1 ) / 4 );
//            uint8_t seq[writeLen]{0};
//            int j = 0;
//            for ( int i = il.first; i < il.second; i++ )
//            {
//                char c = lines[0][i] == 'A' ? 0 : ( lines[0][i] == 'C' ? 1 : ( lines[0][i] == 'G' ? 2 : 3 ) );
//                seq[j/4] += c << ( ( 3 - ( j & 0x3 ) ) * 2 );
//                j++;
//            }
//            fwrite( &il.first, 2, 1, fp );
//            fwrite( &il.second, 2, 1, fp );
//            fwrite( seq, 1, writeLen, fp );
//        }
//
//        seqCount++;
//    }
//
//    fclose( fp );
//    fp = fopen( fnMark.c_str(), "wb" );
//
//    uint8_t buf[4096]{0};
//    CharId kmersLeft = kmerCount;
//    while ( kmersLeft )
//    {
//        CharId thisTotal = min( kmersLeft, (CharId)4096*8 );
//        CharId thisWrite = 1 + ( thisTotal-1 ) / 8;
//        fwrite( buf, 1, thisWrite, fp );
//        kmersLeft -= thisTotal;
//    }
//    
//    cout << "Read " << seqCount << " sequences and " << charCount << " nucleotides, including " << kmerCount << " kmers." << endl;
//
//    fclose( fp );
//}
//
//void KmerFileStruct::clean()
//{
//    struct stat st;
//    assert( !stat( fnMark.c_str(), &st ) );
//    size_t fileSize = st.st_size;
//    fileSize = 1166368800;
//    uint8_t buf[16384]{0};
//    fp = fopen( fnMark.c_str(), "wb" );
//    while ( fileSize )
//    {
//        size_t thisWrite = min( fileSize, (size_t)16384 );
//        fwrite( buf, 1, thisWrite, fp );
//        fileSize -= thisWrite;
//    }
//    fclose( fp );
//}
//
//void KmerFileStruct::close()
//{
//    if ( fp ) fclose( fp );
//    fp = NULL;
//}
//
//void KmerFileStruct::open( size_t posStart )
//{
//    fp = fopen( fnSeq.c_str(), "rb" );
//    if ( posStart ) fseek( fp, posStart, SEEK_CUR );
//    seqsCounted = 0;
//    lastEnd = 0;
//    ilCount = 0;
//}
//
//bool KmerFileStruct::read()
//{
//    if ( !ilCount )
//    {
//        if ( fread( &seqLen, 2, 1, fp ) != 1 ) return false;
//        fread( &ilCount, 1, 1, fp );
//        lastEnd = 31;
//        seqsCounted++;
//        if ( !ilCount )
//        {
//            il[0] = seqLen > 31 ? seqLen - 31 : 0; 
//            il[1] = seqLen;
//            kmerCount = 0;
//            return true;
//        }
//    }
//    else lastEnd = il[1];
//    ilCount--;
//    
//    CharId fk, rk;
//    fread( il, 2, 2, fp );
//    int len = il[1]-il[0];
//    kmerCount = len - 31;
//    int readLen = 1 + ( ( len-1 ) / 4 );
//    uint8_t buf[readLen];
//    fread( &buf, 1, readLen, fp );
//    int j = 0;
//    for ( int i = 0; i < len; i++ )
//    {
//        uint8_t c = rKmer[ i & 0x3 ][ buf[i/4] ];
//        fk <<= 2;
//        rk >>= 2;
//        fk |= c;
//        rk |= rChar[c];
//        if ( i >= 31 ) kmers[j++] = fk < rk ? fk : rk;
//    }
//    
//    return true;
//}
//
//KmerCount::KmerCount( vector<string> &filenames, string fnBase )
//: fnBase_( fnBase )
//{
//    assert( false );
//    for ( int i = 0; i < filenames.size(); i++ )
//    {
//        kfss_.push_back( KmerFileStruct( filenames[i], ( i < 2 ? 2 : 1 ), i < 4 ) );
//    }
//}
//
//void KmerCount::catalog()
//{
//    iKfs_[0] = iKfs_[1] = pMarks_[0] = pMarks_[1] = 0;
//    posKfs_[0] = posKfs_[1] = posMarks_[0] = posMarks_[1] = 0;
//    mf_[0] = mf_[1] = NULL;
//    
//    bool complete = false;
//    int attempt = 0;
//    string fns[2] = { "", "" };
//    
////    {
////        for ( int i = 0; i < kfss_.size(); i++ ) kfss_[i].clean();
//////        fns[1] = fnBase_ + "tmp" + to_string( ++attempt );
//////        mf_[1] = new MarksFile( fns[1], "rb" );
//////        posKfs_[0] = 82763679;
//////        posMarks_[0] = 21401160;
////////        write();
////////        posMarks_[0] = posMarks_[1];
////////        pMarks_[0] = pMarks_[1];
////////        posKfs_[0] = posKfs_[1];
////////        iKfs_[0] = iKfs_[1];
////    }
//    
//    while ( !complete )
//    {
//        fns[0] = fns[1];
//        if ( mf_[1] ) { delete mf_[1]; mf_[0] = new MarksFile( fns[0], "rb" ); }
//        fns[1] = fnBase_ + "tmp" + to_string( ++attempt );
//        mf_[1] = new MarksFile( fns[1], "wb" );
//        
//        complete = read();
//        
//        if ( mf_[0] ) { delete mf_[0]; mf_[0] = new MarksFile( fns[0], "rb" ); }
//        
//        write();
//        
//        if ( mf_[0] ) { delete mf_[0]; remove( fns[0].c_str() ); }
//        if ( complete ) { delete mf_[1]; remove( fns[1].c_str() ); }
//        
//        posMarks_[0] = posMarks_[1];
//        pMarks_[0] = pMarks_[1];
//        posKfs_[0] = posKfs_[1];
//        iKfs_[0] = iKfs_[1];
//    }
//    
//}
//
//bool KmerCount::read()
//{
//    bool complete = true;
//    unordered_map<uint64_t, uint8_t>* kmers = new unordered_map<uint64_t, uint8_t>( 500 * 1000 * 1000 );
//    size_t kmerLimit = kmers->bucket_count() / 4;
//    
//    CharId kmersCounted = 0;
//    double readStart = clock();
//    for ( int i = iKfs_[0]; i < kfss_.size(); i++ )
//    {
//        bool doAdd = complete && kfss_[i].seed;
//        kfss_[i].open( ( i == iKfs_[0] ? posKfs_[0] : 0 ) );
//        CharId seqCount = 0, kmerCount = 0;
//        
//        while ( kfss_[i].read() )
//        {
//            if ( complete ) kmerCount += kfss_[i].kmerCount;
//            for ( int j = 0; j < kfss_[i].kmerCount; j++ )
//            {
//                bool used = true;
//                if ( !mf_[0] || !mf_[0]->read() )
//                {
//                    auto it = kmers->find( kfss_[i].kmers[j] );
//                    if ( it != kmers->end() ) it->second = min( 4, it->second + kfss_[i].worth );
//                    else if ( doAdd ) kmers->insert( make_pair( kfss_[i].kmers[j], kfss_[i].worth ) );
//                    else used = false;
//                }
//
//                if ( !complete ) mf_[1]->write( used );
//            }
//            if ( kfss_[i].ilCount ) continue;
//            if ( complete && kmers->size() >= kmerLimit )
//            {
//                iKfs_[1] = i;
//                posKfs_[1] = ftell( kfss_[i].fp );
//                doAdd = false;
//                complete = false;
//            }
//            if ( complete ) seqCount++;
//        }
//        kfss_[i].close();
//        
//        kmersCounted += kmerCount;
//    }
//    
//    cout << "Read " << to_string( kmersCounted ) << " kmers counted in " << getDuration( readStart ) << endl;
//    
//    FILE* fpKmer = fopen( ( fnBase_ + "kmers" ).c_str(), "wb" );
//    for ( const pair<uint64_t, uint8_t> &kmer : *kmers ) if ( kmer.second > 3 ) fwrite( &kmer.first, 8, 1, fpKmer );
//    fclose( fpKmer );
//    delete kmers;
//    
//    return complete;
//}
//
//void KmerCount::write()
//{
//    unordered_set<uint64_t>* kmers = new unordered_set<uint64_t>( 500 * 1000 * 1000 );
//    FILE* fpKmer = fopen( ( fnBase_ + "kmers" ).c_str(), "rb" );
//    uint64_t inKmer;
//    while ( fread( &inKmer, 8, 1, fpKmer ) == 1 ) kmers->insert( inKmer );
//    fclose( fpKmer );
//    
//    bool ended = false;
//    CharId kmerCount = 0;
//    double writeStart = clock();
//    for ( int i = iKfs_[0]; i < kfss_.size(); i++ )
//    {
//        kfss_[i].open( ( i == iKfs_[0] ? posKfs_[0] : 0 ) );
//        MarksFile* mf = new MarksFile( kfss_[i].fnMark, "rb+" );
//        mf->advance( i == iKfs_[0] ? posMarks_[0] : 0, i == iKfs_[0] ? pMarks_[0] : 0 );
//        CharId seqCount = 0;
//        
//        while ( kfss_[i].read() )
//        {
//            mf->editSkip( kfss_[i].il[0] - ( kfss_[i].lastEnd - 31 ) );
//            for ( int j = 0; j < kfss_[i].kmerCount; j++ )
//            {
//                if ( !mf_[0] || !mf_[0]->read() )
//                {
//                    bool mark = kmers->find( kfss_[i].kmers[j] ) != kmers->end();
//                    mf->edit( mark );
//                    if ( mark ) kmerCount++;
//                }
//                else mf->editSkip( 1 );
//            }
//            if ( kfss_[i].ilCount ) continue;
//            if ( kfss_[i].il[1] < kfss_[i].seqLen ) mf->editSkip( kfss_[i].seqLen - kfss_[i].il[1] );
//            if ( !ended && i == iKfs_[1] && ftell( kfss_[i].fp ) == posKfs_[1] )
//            {
//                pMarks_[1] = mf->p;
//                posMarks_[1] = mf->pos + mf->i;
//                ended = true;
//            }
//            seqCount++;
//        }
//        
//        kfss_[i].close();
//        delete mf;
//    }
//    
//    cout << "Wrote " << kmerCount << " marks from " << to_string( kmers->size() ) << " kmers, counted in " << getDuration( writeStart ) << endl;
//    
//    delete kmers;
//}
