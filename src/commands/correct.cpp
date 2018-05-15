/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "correct.h"
#include "correct_structs.h"
#include "shared_functions.h"
#include "filenames.h"
#include <cassert>
#include <iostream>
#include <string.h>
#include "timer.h"
#include "index.h"

extern Parameters params;

Correct::Correct( int argc, char** argv )
{
    PreprocessFiles* ppf = NULL;
    Filenames* fns = NULL;
    ifstream infile;
    string oprefix;
    
    for ( int i ( 2 ); i < argc; )
    {
        if ( !strcmp( argv[i], "-i" ) )
        {
            if ( infile.is_open() )
            {
                cerr << "Error: multiple inputs provided." << endl;
                exit( EXIT_FAILURE );
            }
            infile.open( argv[i+1] );
            i += 2;
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns )
            {
                cerr << "Error: more than one output prefix provided." << endl;
                exit( EXIT_FAILURE );
            }
            fns = new Filenames( argv[i+1] );
            ppf = new PreprocessFiles( argv[i+1] );
            i += 2;
        }
        else if ( !strcmp( argv[i], "-o" ) )
        {
            oprefix = argv[i+1];
            i += 2;
        }
    }
    
    IndexWriter idx( ppf, 64, 1024 );
    delete ppf;
    
    IndexReader* ir = new IndexReader( fns );
    QueryBinaries* qb = new QueryBinaries( fns );
    Querier bwt( ir, qb );
    delete fns;
    
    correct( bwt, infile, oprefix );
}

void Correct::correct( Querier &bwt, ifstream &infile, string oprefix )
{
    assert( infile.is_open() && infile.good() );
    string line;
    int i = 1;
    int parsed = 0, corrected = 0, report = 100000;
    double startTime = clock();
    while ( getline( infile, line ) )
    {
        Fastq fq( line, 44 );
        ofstream lib( oprefix + "_lib-" + to_string( i++ ) + ".seq" );
        while ( fq.setNext() )
        {
            int nCoords[2], qCoords[2];
            string seq;
            if ( fq.getSeq( seq, nCoords, qCoords, 0 ) && correctRead( bwt, seq, fq.seqs[0], nCoords, qCoords ) ) corrected++;
            if ( fq.getSeq( seq, nCoords, qCoords, 1 ) && correctRead( bwt, seq, fq.seqs[1], nCoords, qCoords ) ) corrected++;
            lib << fq.seqs[0] + '\n' + fq.seqs[1] + '\n';
            parsed += 2;
            if ( parsed >= report )
            {
                cout << to_string( report ) + " | " + to_string( corrected ) + ": " + getDuration( startTime ) << endl;
                report += 100000;
            }
        }
        lib.close();
    }
}

bool Correct::correctRead( Querier &bwt, string &qSeq, string &finalSeq, int nCoords[2], int qCoords[2] )
{
    if ( nCoords[1] - nCoords[0] < 45 ) return false;
    int seqLen = finalSeq.length();
    uint8_t query[seqLen];
    int coords[2] = { nCoords[0], nCoords[1] };
    if ( coords[0] == 0 && coords[1] == seqLen )
    {
        coords[0] = qCoords[0];
        coords[1] = qCoords[1];
    }
    
    if ( coords[0] == 0 && coords[1] == seqLen ) return false;
    
    bool rev = coords[1] == seqLen;
    int start = rev ? coords[0] : coords[1] - 1;
    int len = rev ? seqLen - start : start + 1;
    
    for ( int i = 0; i < len; i++ )
    {
        int j = rev ? i + start : start - i;
        if ( qSeq[j] == 'A' ) query[i] = rev ? 3 : 0;
        else if ( qSeq[j] == 'a' ) query[i] = rev ? 8 : 5;
        else if ( qSeq[j] == 'C' ) query[i] = rev ? 2 : 1;
        else if ( qSeq[j] == 'c' ) query[i] = rev ? 7 : 6;
        else if ( qSeq[j] == 'G' ) query[i] = rev ? 1 : 2;
        else if ( qSeq[j] == 'g' ) query[i] = rev ? 6 : 7;
        else if ( qSeq[j] == 'T' ) query[i] = rev ? 0 : 3;
        else if ( qSeq[j] == 't' ) query[i] = rev ? 5 : 8;
        else query[i] = 4;
        
    }
    
    if ( query[0] > 4 ) query[0] -= 5;
    if ( query[1] > 4 ) query[1] -= 5;
    if ( query[0] > 3 || query[1] > 3 ) return false;
    vector<Overlap> ols = bwt.mapCorrection( query, len );
    
    if ( ols.empty() ) return false;
    
    string extSeq = rev ? qSeq.substr( 0, coords[0] ) : qSeq.substr( coords[1] );
    if ( rev ) revComp( extSeq );
    
    qSeq.clear();
    for ( int i = len; --i >= 0; )
    {
        if ( query[i] == 0 ) qSeq += 'A';
        else if ( query[i] == 1 ) qSeq += 'C';
        else if ( query[i] == 2 ) qSeq += 'G';
        else if ( query[i] == 3 ) qSeq += 'T';
        else break;
    }
    
    bool changed = false;
    for ( uint16_t i = 0; i < extSeq.length(); i++ )
    {
        bool lowQ = islower( extSeq[i] );
        if ( lowQ ) extSeq[i] = toupper( extSeq[i] );
        if ( ols.empty() ) continue;
        int baseCount = 0;
        bool diff = false;
        for ( int j = 0; j < ols.size(); j++ )
        {
            if ( ols[j].extLen <= i ) ols.erase( ols.begin() + j, ols.end() );
            else
            {
                if ( ols[j].seq[i] == extSeq[i] ) baseCount++;
                if ( ols[j].seq[i] != ols[0].seq[i] ) diff = true;
            }
        }
        
        if ( baseCount && diff )
        {
            for ( int j = 0; j < ols.size(); )
            {
                if ( ols[j].seq[i] != extSeq[i] ) ols.erase( ols.begin() + j );
                else j++;
            }
        }
        
        if ( !baseCount )
        {
            if ( diff || ols.size() < 4 )
            {
                if ( changed ) extSeq = extSeq.substr( 0, i );
                ols.clear();
            }
            else extSeq[i] = ols[0].seq[i];
            if ( lowQ ) changed = true;
        }
    }
    
    qSeq += extSeq;
    if ( qSeq.length() < 45 ) return false;
    finalSeq = qSeq;
    if ( rev ) revComp( finalSeq );
    
    return true;
}

//void Correct::correct()
//{
//    vector<TrimFile> files;
//    files.push_back( TrimFile( prefix_ + "mp.fastq", true ) );
//    files.push_back( TrimFile( prefix_ + "pe.fastq", false ) );
//    files.push_back( TrimFile( prefix_ + "single.fastq", false, true ) );
//    
//    ofstream pe( prefix_ + "pe.seq" );
//    ofstream mp( prefix_ + "mp.seq" );
//    ofstream se( prefix_ + "single.seq" );
//    
//    string lines[8];
//    ReadId report = 20000;
//    double startTime = clock();
//    int good = 0, bad = 0;
//    for ( TrimFile &file : files )
//    {
//        while ( file.getNext( &lines[0], 4, 43 ) )
//        {
//            if ( good + bad > report )
//            {
//                cout << getDuration( startTime ) << endl;
//                report += 20000;
//            }
//            ( bwt_.correct( lines[1] ) ? good : bad )++;
//            if ( !file.single && file.getNext( &lines[4], 4, 43 ) )
//            {
//                ( bwt_.correct( lines[5] ) ? good : bad )++;
//                if ( lines[1].find( lines[5] ) != lines[1].npos ) lines[5].clear();
//                if ( lines[5].find( lines[1] ) != lines[5].npos ) lines[1].clear();
//                if ( lines[5].empty() ) 
//                {
//                    if ( !lines[1].empty() ) writeFile( se, &lines[0], 1 );
//                }
//                else if ( lines[1].empty() )
//                {
//                    if ( !lines[1].empty() ) writeFile( se, &lines[4], 1 );
//                }
//                else
//                {
//                    writeFile( ( file.mp ? mp : pe ), &lines[0], 1 );
//                    writeFile( ( file.mp ? mp : pe ), &lines[4], 1 );
//                }
//            }
//            else
//            {
//                if ( !lines[1].empty() ) writeFile( se, &lines[0], 1 );
//            }
//        }
//    }
//    
//    pe.close();
//    mp.close();
//    se.close();
//    
//    assert( false );
//}

//void Correct::decontaminate()
//{
//    vector<TrimFile> files;
//    files.push_back( TrimFile( prefix_ + "mp.seq", true ) );
//    files.push_back( TrimFile( prefix_ + "pe.seq", false ) );
//    
//    vector<SamFile> sams;
//    sams.push_back( prefix_ + "mp.sam" );
//    sams.push_back( prefix_ + "pe.sam" );
//    
//    vector<ofstream> ofps;
//    ofps.push_back( ofstream( prefix_ + "mp.seq-tmp" ) );
//    ofps.push_back( ofstream( prefix_ + "pe.seq-tmp" ) );
//    ofstream se( prefix_ + "single.seq", ios::out | ios::app );
//    
//    int pairs = 0, badPairs = 0, bads = 0, singles = 0;
//    for ( int i : { 0, 1 } )
//    {
//        int j = 0;
//        int k = sams[i].getNext();
//        string seqs[2];
//        bool isBad[2];
//        while ( files[i].getNext( &seqs[0], 1, 0 ) && files[i].getNext( &seqs[1], 1, 0 ) )
//        {
//            isBad[0] = k == j++;
//            if ( isBad[0] ) k = sams[i].getNext();
//            isBad[1] = k == j++;
//            if ( isBad[1] ) k = sams[i].getNext();
//            
//            if ( isBad[0] && isBad[1] ) badPairs++;
//            else if ( isBad[0] ) bads++;
//            else if ( isBad[1] ) bads++;
//            
//            if ( !isBad[0] && !isBad[1] )
//            {
//                ofps[i] << seqs[0] << '\n';
//                ofps[i] << seqs[1] << '\n';
//                pairs += 2;
//            }
//            else if ( !isBad[0] ) se << seqs[0] << '\n';
//            else if ( !isBad[1] ) se << seqs[1] << '\n';
//        }
//        ofps[i].close();
//    }
//    
//    se.close();
//}

//void Correct::trim( vector<string> &pes, vector<string> &mps )
//{
//    bool isFastq = true;
//    
//    vector<TrimFile> files;
//    for ( string &mp : mps ) files.push_back( TrimFile( mp, true ) );
//    for ( string &pe : pes ) files.push_back( TrimFile( pe, false ) );
//    
//    string lines[8];
//    
//    ofstream pe( prefix_ + "pe.fastq" );
//    ofstream mp( prefix_ + "mp.fastq" );
//    ofstream se( prefix_ + "single.fastq" );
//    
//    for ( TrimFile &file : files )
//    {
//        while ( file.getNext( &lines[0], 4, 43 ) && file.getNext( &lines[4], 4, 43 ) )
//        {
//            string extra[4];
//            bool didJunc = false;
//            if ( trimPrimers( &lines[0], &extra[0], isFastq, file.mp ) ) didJunc = true;
//            if ( trimPrimers( &lines[4], &extra[2], isFastq, file.mp ) ) didJunc = true;
//
//            if ( !extra[0].empty() && !extra[2].empty() )
//            {
//                trimPair( lines[1], extra[2] );
//                trimPair( lines[5], extra[0] );
//            }
//            
//            bool didWrite = false;
//            for ( int i : { 1, 5 } )
//            {
//                int j = i == 5 ? 1 : 5;
//                if ( ( lines[i].length() >= 35 && lines[j].length() < 13 )
//                        || lines[i].find( lines[j] ) != lines[i].npos  )
//                {
//                    writeFile( se, &lines[i-1], isFastq ? 4 : 1 );
//                    didWrite = true;
//                    break;
//                }
//            }
//            
//            if ( didWrite || ( lines[1].length() < 35 && lines[5].length() < 35 ) ) continue;
//            
//            writeFile( ( file.mp || didJunc ? mp : pe ), &lines[0], isFastq ? 4 : 1 );
//            writeFile( ( file.mp || didJunc ? mp : pe ), &lines[4], isFastq ? 4 : 1 );
//        }
//        
//        file.fp.close();
//    }
//    
//    pe.close();
//    mp.close();
//    se.close();
//}

//void Correct::trimPair( string &seq, string &rev )
//{
//    revComp( rev );
//    for ( int i = 0; i < min( seq.length(), rev.length()  ); i++ )
//    {
//        if ( seq[seq.length()-i-1] == 'N' && rev[rev.length()-i-1] != 'N' )
//        {
//            seq[seq.length()-i-1] = rev[rev.length()-i-1];
//        }
//        else if ( seq[seq.length()-i-1] != rev[rev.length()-i-1] ) return;
//    }
//}
//
//bool Correct::trimPrimers( string lines[4], string extra[2], bool isFastq, bool isMp )
//{
//    size_t it, len;
//    it = lines[1].find( "AGATCGGAAGAGCA" );
//    if ( it == lines[1].npos ) it = lines[1].find( "AGATCGGAAGAGCG" );
//    if ( it == lines[1].npos )
//    {
//        len = getEndTrim( lines[1], "AGATCGGAAGAGC", 1 );
//        if ( len ) it = lines[1].length() - len;
//    }
//    if ( it == lines[1].npos && lines[1].length() < 235 && lines[1].back() == 'A' )
//    {
//        it = lines[1].length() - 1;
//    }
//    if ( it != lines[1].npos )
//    {
//        lines[1] = lines[1].substr( 0, it );
//        if ( isFastq ) lines[3] = lines[3].substr( 0, it );
//    }
//    
//    bool didJunc = false;
//    it = lines[1].find( "CTGTCTCTTATACACATCT" );
//    if ( it != lines[1].npos )
//    {
//        len = min( lines[1].length() - it, size_t(38) );
//    }
//    else
//    {
//        it = lines[1].find( "AGATGTGTATAAGAGACAG" );
//        if ( it != lines[1].npos )
//        {
//            if ( it > 19 ) { it -= 19; len = 38; }
//            else { len = it + 19; it = 0; }
//        }
//    }
//    
//    if ( it == lines[1].npos && isMp )
//    {
//        len = getEndTrim( lines[1], "CTGTCTCTTATACACATC", 1 );
//        if ( len ) it = lines[1].length() - len;
//    }
//    
//    if ( it != lines[1].npos )
//    {
//        extra[0] = lines[1].substr( it + len );
//        lines[1] = lines[1].substr( 0, it );
//        if ( isFastq )
//        {
//            extra[1] = lines[3].substr( it + len );
//            lines[3] = lines[3].substr( 0, it );
//        }
//        didJunc = true;
//    }
//    
//    return didJunc;
//}

//void Correct::writeFile( ofstream &fp, string lines[4], int lineCount )
//{
//    if ( lineCount == 1 ) fp << lines[1] << '\n';
//    else for ( int i = 0; i < lineCount; i++ ) fp << lines[i] << '\n';
//}
