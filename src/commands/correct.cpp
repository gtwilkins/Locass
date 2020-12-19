/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "correct.h"
#include "shared_functions.h"
#include "filenames.h"
#include <cassert>
#include <iostream>
#include <string.h>
#include "timer.h"
#include "index.h"
#include "local_alignment.h"
#include "error.h"
#include "deadapter.h"
#include "constants.h"
#include "correct_structs.h"
#include "correct_amplicon.h"

extern Parameters params;

Correct::Correct( int argc, char** argv )
: fileType_( 0 )
{
    ir_ = NULL;
    da_ = NULL;
    fns_ = NULL;
    df_ = NULL;
    used_ = NULL;
//    PreprocessFiles* ppf = NULL;
    string oprefix;
    demate_ = deadapter_ = deamplify_ = false;
    pairCount_ = trimCount_ = adpCount_ = discardCount_ = overlapCount_ = reoverlapCount_ = connectCount_ = connectPairCount_ = 0;
    cap_ = 0;
    
    for ( int i ( 2 ); i < argc; i++ )
    {
        if ( !strcmp( argv[i], "-i" ) )
        {
            if ( !ifn_.empty() ) failure( "Multiple inputs provided." );
            ifn_ = argv[++i];
            if ( i + 1 < argc && argv[i+1][0] != '-' ) assert( false );
        }
        else if ( !strcmp( argv[i], "-p" ) )
        {
            if ( fns_ ) failure( "More than one index prefix provided." );
            fns_ = new Filenames( argv[++i] );
//            ppf = new PreprocessFiles( argv[i], true );
        }
        else if ( !strcmp( argv[i], "-t" ) )
        {
            cap_ = atoi( argv[++i] );
        }
        else if ( !strcmp( argv[i], "-o" ) ) oprefix = argv[++i];
        else if ( !strcmp( argv[i], "--de-mate" ) ) demate_ = true;             // Looks for mated pairs that either contain an adapter or are chimeric
        else if ( !strcmp( argv[i], "--de-adapter" ) ) deadapter_ = true;       // Looks for reads that run through to a 3' adapter
        else if ( !strcmp( argv[i], "--de-amplify" ) ) deamplify_ = true;       // Looks for reads for reads that are redundant amplicons
        else failure( "Unrecoginsed command: " + string( argv[i] ) );
    }
//    IndexWriter idx( ppf, 64, 1024 );
//    IndexWriter idx( ppf, 1024, 20000 );
//    delete ppf;
//    {
//        double startTime = clock();
//        string fn = "/media/glen/ssd/Sp/Sp-mer.dat";
//        ir_ = new IndexReader( fns_ );
//        ir_->createSeeds( fn, 12 );
//        cout << getDuration( startTime );
//        exit( 1 );
//    }
    if ( !ifn_.empty() )
    {
        ifstream ifs( ifn_ );
        string line;
        getline( ifs, line );
        if ( line.empty() ) failure( "Couldn't read sequences from file." );
        if ( line[0] == '>' ) fileType_ = 2;
        else if ( line[0] == '@' ) fileType_ = 4;
        else fileType_ = 1;
        ifs.close();
    }
    else failure( "Input file not provided." );
    if ( !Filenames::exists( ifn_ ) ) failure( "Could not open input file." );
    if ( oprefix.empty() ) oprefix = ifn_.substr( 0, ifn_.find_last_of( '.' ) );
    ofns_[0] = oprefix + "_corrected." + ( fileType_ == 4 ? "fastq" : ( fileType_ == 2 ? "fasta" : "seq") );
    ofns_[1] = oprefix + "_tmp_seqs.bin";
    ofns_[2] = oprefix + "_tmp_stat.bin";
    ofns_[3] = oprefix + "_tmp_kmer.bin";
    ofns_[4] = oprefix + "_tmp_amps.bin";
//    for ( int i = 1; i < 5; i++ ) if ( Filenames::exists( ofns_ [i] ) ) failure( "Output file \"" + ofns_[i] + "\" already exists." );
    
    if ( !fns_ && demate_ ) failure( "Index files must be provided for --de-mate command." );
    
    
    double startTime = clock();
    if ( deadapter_ ) da_ = new Deadapter( ifn_ );
    ofs_.open( ofns_[0] );
    
    
    if ( deamplify_ ) deamplify();
    correct();
    
    ofs_.close();
    
    uint64_t discardCount = discardCount_ + Amp::ampCount - Pile::pileCount;
    bool corrected = adpCount_ || trimCount_ || CorrectAlign::nCorrectCount || CorrectAlign::polyCorrectCount || CorrectAlign::errorCorrectCount;
    
    cout << "Overall summary:" << endl;
    cout << pairCount_ << " read pairs were scrutinised in " << getDuration( startTime ) << endl;
    if ( Amp::ampCount ) cout << to_string( Amp::ampCount - Pile::pileCount ) << " pairs were discarded as redundant amplicons." << endl;
    if ( discardCount_ ) cout << discardCount_ << " pairs were discarded as self-complements." << endl;
    if ( corrected && discardCount ) cout << ( pairCount_ - discardCount ) << " pairs remained undiscarded." << endl;
    if ( adpCount_ ) cout << adpCount_ << " pairs were consolidated and trimmed as self-compliments." << endl;
    if ( connectCount_ ) cout << connectPairCount_ << " pairs contained connecting adapters. Consequently, " << connectCount_ << " reads were trimmed" << endl;
    if ( overlapCount_ ) cout << overlapCount_ << " paired-end reads overlapped, " << reoverlapCount_ << " of those reads were corrected " << endl;
    if ( trimCount_ ) cout << trimCount_ << " chimeric mated reads were trimmed." << endl;
    if ( CorrectAlign::nCorrectCount ) cout << CorrectAlign::nCorrectCount << " undetermined bases were corrected." << endl;
    if ( CorrectAlign::polyCorrectCount ) cout << CorrectAlign::polyCorrectCount << " homopolymer bases were corrected." << endl;
    if ( CorrectAlign::errorCorrectCount ) cout << CorrectAlign::errorCorrectCount << " incorrect bases were corrected." << endl;
    
    remove( ofns_[1].c_str() );
    remove( ofns_[2].c_str() );
    remove( ofns_[3].c_str() );
    remove( ofns_[4].c_str() );
}

Correct::~Correct()
{
    if ( da_ ) delete da_;
    if ( fns_ ) delete fns_;
}

void Correct::correct()
{
    if ( demate_ ) cout << "Scrutinising remaing read pairs..." << endl << endl;
    else cout << "Scrutinising read pairs..." << endl << endl;
    double startTime = clock();
    
    string lines[2][4];
    ifstream ifsReads( ifn_ );
    assert( ifsReads.good() );
    if ( df_ ) df_->rewind();
    if ( demate_ && !ir_ ) ir_ = new IndexReader( fns_ );
    if ( deamplify_ && !df_ ) df_ = new DeamplifyFiles( ofns_[1], ofns_[2] );
    time_t my_time = time(NULL);
    printf("%s", ctime(&my_time));
    
    assert( !cap_ || ( cap_ > 50 && cap_ < 500 ) );
    
    uint64_t id = 0, usedCount = 0;
    uint8_t idTable[8] = { 128, 64, 32, 16, 8, 4, 2, 1 };
    int baseLen = readLen_, cutoff = readLen_ * .9;
    while ( readFile( ifsReads, lines[0] ) && readFile( ifsReads, lines[1] ) )
    {
        if ( cap_ ) for ( int i : { 0, 1 } ) for ( int j : { 1, 3 } ) if ( lines[i][j].size() > cap_ )
        {
            lines[i][j].erase( lines[i][j].begin() + cap_, lines[i][j].end() );
        }
        if ( df_ ) df_->overwrite( lines );
        bool used = used_ && ( used_[id/8] & ( idTable[id%8] ) );
        bool adapter = false, connector = false;
        if ( !used && da_ )
        {
            adapter = da_->isOverlap( lines[0][1], lines[1][1], lines[0][3], lines[1][3], !demate_ );
            for ( int d : { 0, 1 } ) if ( da_->isConnected( lines[d][1], lines[d][3] ) && ( connector = true ) ) connectCount_++;
        }
        
        // Correct mated pair
        if ( !used && !adapter && !connector && demate_ )
        {
            bool trimmed[2]{ false, false };
            int lens[2]{0};
            
            for ( int i = 0; i < 2; i++ ) CorrectQuery cq( ir_, lines[i][1], lens[i], trimmed[i], false );
            
            for ( int i = 0; i < 2; i++ )
            {
                if ( trimmed[i] || !lens[i] || lens[i] == lines[i][1].size() ) continue;
                if ( trimmed[!i] && lens[i] > cutoff ) continue;
                CorrectQuery cq( ir_, lines[i][1], lens[i], trimmed[i], false );
            }
            
            for ( int i = 0; i < 2; i++ )
            {
                if ( !trimmed[i] ) lens[i] = baseLen;
                if ( lines[i][1].size() > lens[i] ) lines[i][1] = lines[i][1].substr( 0, lens[i] );
                if ( lines[i][3].size() > lens[i] ) lines[i][3] = lines[i][3].substr( 0, lens[i] );
                assert( lines[i][1].size() == lines[i][3].size() );
                if ( lines[i][1].size() < baseLen ) lines[i][1] += string( baseLen-lines[i][1].size(), 'N' );
                if ( lines[i][3].size() < baseLen ) lines[i][3] += string( baseLen-lines[i][3].size(), char( 35 ) );
            }
            
            if ( trimmed[0] || trimmed[1] ) trimCount_++;
        }
        if ( used ) usedCount++;
        if ( adapter ) ( demate_ || lines[0][1].empty() ? discardCount_ : adpCount_ )++;
        if ( connector ) connectPairCount_++;
        
        if ( !demate_ && !adapter && overlapCorrect( lines[0][1], lines[1][1], lines[0][3], lines[1][3] ) ) overlapCount_++;
        
        if ( !used && ( !demate_ || !adapter ) && !lines[0][1].empty() )
        {
            // Silence verbose headers
            lines[0][0] = fileType_ == 4 ? "@" : ">";
            lines[1][0] = fileType_ == 4 ? "@" : ">";
            lines[0][2] = "+";
            lines[1][2] = "+";
            
            // Silence reverse pair if adapter found
            if ( adapter ) lines[1][1] = string( lines[0][1].size(), 'N' );
            if ( adapter ) lines[1][3] = string( lines[0][1].size(), 35 );
            for ( int i = 0; i < 2; i++ )
            {
                if ( fileType_ > 1 ) for ( int j = 0; j < fileType_; j++ ) ofs_ << lines[i][j] << "\n";
                else if ( fileType_ == 1 ) ofs_ << lines[i][1] << "\n";
            }
        }
        
        id++;
        if ( !( id % 5000000 ) )
        {
            my_time = time(NULL);
            printf("%s", ctime(&my_time));
            cout << to_string( id ) << " read pairs corrected so far in " << getDuration( startTime ) << endl;
        }
    }
    pairCount_ = id;
    
    cout << "Finished parsing " << id << " read pairs in " << getDuration( startTime ) << endl;
    if ( usedCount ) cout << "Note: " << usedCount << " read pairs were already handled during de-amplification." << endl;
    cout << endl;
    
    if ( ir_ ) delete ir_;
    if ( used_ ) delete used_;
    ir_ = NULL;
    used_ = NULL;
}

void Correct::deamplify()
{
    deamplifyConvert();
    deamplifyIdentify();
    deamplifyScrutinise();
}

void Correct::deamplifyConvert()
{
    cout << "Reading in sequences from file... " << endl;
    double startTime = clock();
    
    ifstream ifsReads( ifn_ );
    FILE* ofpSeq = fopen( ofns_[1].c_str(), "wb" ),* ofpStatus = fopen( ofns_[2].c_str(), "wb" );
    int readLen = 0, p = 0, b = 0, bufSize = 65536;
    uint64_t kmers[2];
    uint32_t pairCount = 0, discardCount = 0;
    uint8_t bufStatus[bufSize]{0}, bufSeq[2][1024], status;
    fwrite( &pairCount, 4, 1, ofpSeq );
    fwrite( &readLen_, 1, 1, ofpSeq );
    
    if ( !ir_ ) ir_ = new IndexReader( fns_ );
    
    string lines[2][4];
    
    while ( readFastq( ifsReads, lines[0] ) )
    {
        if ( !readFastq( ifsReads, lines[1] ) ) failure( "Unequal number of read pairs." );
        pairCount++;
        uint8_t lens[2] = { (uint8_t)lines[0][1].size(), (uint8_t)lines[1][1].size() };
        status = lens[0] < 16 || lens[1] < 16 ? 3 : 0;
        for ( int i = 0; i < 2; i++ )
        {
            readLen = max( readLen, (int)lines[i][1].size() );
            kmers[i] = 0;
            size_t it = lines[i][1].find( 'N' );
            if ( it != string::npos && it < 16 && ir_ )
            {
                bool trimmed = false;
                int len = 0;
                CorrectQuery cq( ir_, lines[i][1], len, trimmed, false );
            }
            for ( int j = 0; !status && j < lens[i]; j++ )
            {
                uint8_t c = charToInt[ lines[i][1][j] ];
                if ( j < 16 && c > 3 ) status = 3;
                bufSeq[i][j] = c > 3 ? 255 : c * 63 + lines[i][3][j] - 33;
                if ( j < 16 ) kmers[i] = ( kmers[i] << 2 ) + c;
                assert( bufSeq[i][j] < 252 || bufSeq[i][j] == 255 );
            }
        }
        if ( p == bufSize )
        {
            fwrite( bufStatus, 1, p, ofpStatus );
            memset( bufStatus, 0, p );
            p = 0;
        }
        bufStatus[p] += status << ( ( 3 - b ) * 2 );
        b = b == 3 ? 0 : b+1;
        if ( b == 0 ) p++;
        if ( status ) discardCount++;
        if ( status ) continue;
        
        fwrite( &lens, 1, 2, ofpSeq );
        int d = kmers[1] < kmers[0];
        kmers[0] = ( kmers[d] << 32 ) + kmers[!d];
        fwrite( &kmers[0], 8, 1, ofpSeq );
        fwrite( bufSeq[d], 1, lens[d], ofpSeq );
        fwrite( bufSeq[!d], 1, lens[!d], ofpSeq );
    }
    if ( readFastq( ifsReads, lines[1] ) ) failure( "Error: unequal number of read pairs." );
    if ( readLen > 255 ) failure( "Error: maximum read length is too high. Max supported is 255." );
    readLen_ = readLen;
    pairCount_= pairCount;
    
    if ( p || b ) fwrite( bufStatus, 1, p + bool( b ), ofpStatus );
    
    delete ir_;
    ir_ = NULL;
    fclose( ofpSeq );
    fclose( ofpStatus );
    ofpSeq = fopen( ofns_[1].c_str(), "rb+" );
    fwrite( &pairCount, 4, 1, ofpSeq );
    fwrite( &readLen_, 1, 1, ofpSeq );
    fclose( ofpSeq );
    
    cout << "Read " << to_string( pairCount ) << " read pairs, discarding " << to_string( discardCount ) << " in: " << getDuration( startTime ) << endl << endl;
}

void Correct::deamplifyIdentify()
{
    cout << "Compiling candidate amplicons... " << endl;
    double startTime = clock();
    
    if ( !df_ ) df_ = new DeamplifyFiles( ofns_[1], ofns_[2] );
    pairCount_ = df_->maxPair;
    
    FILE* ofpDump = fopen( ofns_[4].c_str(), "wb" );
    uint64_t kmerSize = 500 * 1000 * 1000, maxMem = 2 * 1000 * 1000 * 1000, ampliconCount = 0;
    bool full = df_->maxPair;
    
    while ( full )
    {
        // Parse as many k-mers as possible
        df_->cycle( ofns_[3], kmerSize, full );
        
        // Collate and de-amplify all read for multiple k-mers
        df_->collate( ofpDump, ofns_[3], kmerSize, maxMem, ampliconCount );
        
        df_->reset( 1 );
    }
    
    fclose( ofpDump );
    
    cout << "Identified " << to_string( ampliconCount ) << " potential amplicons out of " << df_->maxPair << " total read pairs in " << getDuration( startTime ) << endl << endl;
}

void Correct::deamplifyScrutinise()
{
    cout << "Scrutinising candidate amplicons... " << endl;
    double startTime = clock();
    
    assert( fns_ );
    if ( !df_ ) df_ = new DeamplifyFiles( ofns_[1], ofns_[2] );
    if ( !ir_ ) ir_ = new IndexReader( fns_ );
    if ( !pairCount_ ) pairCount_ = df_->maxPair;
    used_ = new uint8_t[ 1 + ( pairCount_ / 8 ) ]{0};
    
    FILE* ifpDump = fopen( ofns_[4].c_str(), "rb" );
    uint64_t lastKmer;
    
    vector<Amp*> amps[2];
    uint64_t candidateCount = 0;
    while ( !feof( ifpDump ) )
    {
        candidateCount++;
        Amp* tmps[2]{ NULL, NULL };
        uint64_t kmer = Amp::create( ifpDump, tmps );
        if ( kmer != lastKmer && !amps[0].empty() )
        {
            for ( Pile* pile : Pile::compile( amps, ir_ ) ) pile->output( ofs_, used_, da_, ir_ );
        }
        if ( tmps[0] ) amps[0].push_back( tmps[0] );
        if ( tmps[1] ) amps[1].push_back( tmps[1] );
        lastKmer = kmer;
        if ( !( candidateCount % 5000000 ) )
        {
            time_t my_time = time(NULL);
            printf("%s", ctime(&my_time));
            cout << to_string( candidateCount ) << " candidate amplicons parsed so far in " << getDuration( startTime ) << endl;
        }
    }
    
    fclose( ifpDump );
    trimCount_ += Pile::trimCount;
    discardCount_ += Pile::adpCount;
    
    cout << "Finished scrutinising " << candidateCount << " candidate amplicons in " << getDuration( startTime ) << endl;
    cout << Amp::ampCount << " candidates were consolidated into " << Pile::pileCount << " unique pairs." << endl << endl;
    
}

bool Correct::overlapCorrect( string &seq1, string &seq2, string &phred1, string &phred2 )
{
    return false;
    string rev = revCompNew( seq2 );
    for ( int i = seq1.size()/2; i+24 < seq1.size(); i+=24 )
    {
        string q = seq1.substr( i, 24 );
        size_t it = rev.find( q );
        if ( it == string::npos ) continue;
        if ( i < (int)it )
        {
            continue;
        }
        
        int off = (int)it - i, l = i, r = i+24;
        while ( min( l, l+off ) > 0 && seq1[l-1] == rev[l+off-1] ) l--;
        while ( r < seq1.size() && r+off < rev.size() && seq1[r] == rev[r+off] ) r++;
        
        if ( r-l < 48 && ( r+12 < seq1.size() ) && l+off > 12 )
        {
            int x = 0;
            break;
        }
        
        bool corrected = false;
        
        for ( ; r < seq1.size() && r+off < rev.size(); r++ ) if ( seq1[r] != rev[r+off] && ( corrected = true ) )
        {
            seq1[r] = rev[r+off];
            if ( r < phred1.size() ) phred1[r] = 35;
        }
        
        for ( ; l-- > 0 && l+off >= 0; ) if ( seq1[l] != rev[l+off] && ( corrected = true ) )
        {
            assert( seq2.end()[-l-off-1] != getComp( seq1[l] ) );
            seq2.end()[-l-off-1] = getComp( seq1[l] );
            if ( l+off < phred2.size() ) phred2.end()[-l-off-1] = 35;
        }
        
        if ( corrected ) reoverlapCount_++;
        
        return true;
        
    }
    
    return false;
}

bool Correct::readFastq( ifstream &ifs, string (&lines)[4] )
{
    if ( !ifs.good() || !getline( ifs, lines[0] ) || !getline( ifs, lines[1] ) || !getline( ifs, lines[2] ) || !getline( ifs, lines[3] ) ) return false;
    return true;
}

bool Correct::readFile( ifstream &ifs, string (&lines)[4] )
{
    if ( !ifs.good() || !fileType_ ) return false;
    if ( fileType_ > 1 && !getline( ifs, lines[0] ) ) return false;
    if ( !getline( ifs, lines[1] ) ) return false;
    if ( fileType_ == 4 && ( !getline( ifs, lines[2] ) || !getline( ifs, lines[3] ) ) ) return false;
    if ( fileType_ != 4 ) lines[3] = lines[1].empty() ? "" : string( lines[1].size(), 63 );
    return true;
}

//void Correct::censorFiles( string fn, ofstream &ofsOut, bool paired, bool chim, bool pyro )
//{
//    assert( paired + pyro + chim < 2 );
//    ifstream ifsSeqs( "/home/glen/Genomes/Lv/" + fn + ".fastq" );
//    MarksFile mf( "/media/glen/ssd/Lv/" + fn + "_marks.bin", "rb" );
//    assert( ofsOut.good() );
//    assert( ifsSeqs.good() );
//    
//    auto reinterval = []( vector< pair<int,int> > &ils )
//    {
//        for ( int i = 1; i < ils.size(); i++ )
//        {
//            int ol = ils[i-1].second - ils[i].first;
//            if ( ol < 0 ) continue;
//            ils[i-1].second -= ol+1;
//            ils[i].first += ol+1;
//        }
//    };
//    
//    string lines[8];
//    vector< pair<int,int> > ils[2];
//    int lineCount = 0, outCount = 0;
//    while ( getline( ifsSeqs, lines[0] ) )
//    {
//        for ( int i = 0; i < ( paired ? 7 : 3 ); i++ ) getline( ifsSeqs, lines[i+1] );
//        ils[0] = mf.getIntervals( lines[1].size() );
//        ils[1] = mf.getIntervals( lines[5].size() );
//        if ( pyro || chim ) reinterval( ils[0] );
//        
//        if ( pyro )
//        {
//            for ( auto it = ils[0].begin(); it != ils[0].end(); it++ )
//            {
//                int len = it->second - it->first - 95;
//                int readCount = len > -21;
//                while ( len > 63 ){ len -= 64; readCount++; };
//                int buffer = len > 0 ? ( len / ( readCount+1 ) ) : 0;
//                for ( int i = 0; i < readCount; i++ )
//                {
//                    int base = it->first + ( i * 64 ) + ( (i+1) * buffer );
//                    ofsOut << lines[1].substr( base, min( 95, it->second - base ) ) + "\n";
//                    outCount++;
//                }
//            }
//        }
//        else if ( chim )
//        {
//            for ( auto it = ils[0].begin(); it != ils[0].end(); it++ )
//            {
//                int len = it->second - it->first;
//                if ( len >= 50 )
//                {
//                    ofsOut << lines[1].substr( it->first, len ) + "\n";
//                    outCount++;
//                }
//            }
//        }
//        else
//        {
//            for ( int i = 0; i < 1+paired; i++ )
//            {
//                int j = 1 + ( i * 4 );
//                int len = lines[j+2].size();
//                while ( len > 0 && lines[j+2][len-1] < 39 ) len--;
//                if ( !ils[i].empty() ) len = max( len, ils[i].back().second );
//                ofsOut << lines[j].substr( 0, max( 1, len ) ) + "\n";
//                outCount++;
//            }
//        }
//        lineCount += 1 + paired;
//    }
//    
//    cout << "Output " << outCount << " sequences from " << lineCount << " reads." << endl;
//}
//
//void Correct::correct( Querier &bwt, string fn, bool paired, bool chim, bool pyro )
//{
//    MarksFile mf( fn + "_marks.bin", "rb" );
//    ifstream ifsSeqs( fn + ".fastq" );
//    ofstream ofsSeqs( fn + "_corrected.fastq" );
//    ofstream ofsDump;
//    if ( chim ) ofsDump.open( fn + "_dump.fastq" );
//    assert( ifsSeqs.good() );
//    double startTime = clock();
//    CharId readCount = 0;
//    Interval::bridgeCount = Interval::capCount = 0;
//    string lines[8];
//    vector<Interval*> ils[4];
//    int lineCount = paired ? 8 : 4;
//    while ( readSequence( ifsSeqs, mf, &lines[0], ils[0], chim, pyro ) )
//    {
//        if ( paired ) assert( readSequence( ifsSeqs, mf, &lines[4], ils[1], chim, pyro ) );
//        
////        if ( readCount < 3517687 )
////        {
////            for ( Interval* il : ils[0] ) delete il;
////            ils[0].clear();
////            readCount++;
////            continue;
////        }
//        Interval::correct( bwt, ils, paired && !chim  && !pyro, chim );
//        
//        Interval::output( ils[0], ofsSeqs, &lines[0], pyro );
//        if ( paired ) Interval::output( ils[1], ofsSeqs, &lines[4], pyro );
//        if ( chim && !ils[2].empty() ) Interval::output( ils[2], ofsDump, &lines[0], false );
//        if ( chim && !ils[3].empty() ) Interval::output( ils[3], ofsDump, &lines[4], false );
//        
//        readCount += 1 + paired;
//        if ( readCount % 100000 ) continue;
//        cout << getDuration( startTime ) << endl;
//        cout << readCount << endl;
//    }
//    
//    ifsSeqs.close();
//    ofsSeqs.close();
//    
//    cout << "Corrected " << Interval::bridgeCount << " bridges and " << Interval::capCount << " caps from " << readCount << " reads in " << getDuration( startTime ) << endl;
//}
//
////void Correct::correct( Querier &bwt, string (&lines)[8], vector<Interval*> (&ils)[3], bool paired, bool chim, bool pyro )
////{
////    Interval::correct( bwt, ils, paired && !chim  && !pyro );
////    
////    if ( chim )
////    {
////        vector<Interval*> split;
////        Interval* splits[2] = { Interval::getSplit( ils[0] ), Interval::getSplit( ils[1] ) };
////        if ( splits[0] && !splits[1] ) splits[0]->split( bwt, ils[0], split );
////        if ( splits[1] && !splits[0] ) splits[1]->split( bwt, ils[1], split );
////        for ( Interval* il : split ) delete il;
////    }
////    
////    Interval::output( ils[0], lines[1], lines[3] );
////    if ( paired ) Interval::output( ils[1], lines[5], lines[7] );
////}
//
//bool Correct::readSequence( ifstream &ifs, MarksFile &mf, string* lines, vector<Interval*> &ils, bool chim, bool pyro )
//{
//    if ( !getline( ifs, lines[0] ) ) return false;
//    if ( lines[0][0] != '@' ) return false;
//    for ( int i = 0; i < 3; i++ ) getline( ifs, lines[i+1] );
//    vector<bool> marks( max( 0, (int)lines[1].size()-31 ), false );
//    for ( int i = 0; i < marks.size(); i++ ) marks[i] = mf.read();
//    ils = Interval::create( lines[1], lines[3], marks, chim || pyro, pyro );
//    return true;
//}
//
////bool Correct::splitJoin( string &seq1, string &seq2, string &phred1, string &phred2 )
////{
////    if ( !seq2.empty() ) return false;
////    
////    static int i = 0;
////    i++;
////    LocalAlignment glocal( seq1, joiner_, true, true );
////    int coords[2];
////    int score = glocal.setCoords( coords );
//////    string align1, align2;
//////    glocal.setAlign( align1, align2 );
////    seq2 = seq1.substr( coords[1] );
////    seq1 = seq1.substr( 0, coords[0] );
////    phred2 = phred1.substr( coords[1] );
////    phred1 = phred1.substr( 0, coords[0] );
//////    if ( score < 41 ) cout << align1 << endl << align2 << endl << endl;
////    while ( seq1.back() == 'N' )
////    {
////        seq1.pop_back();
////        phred1.pop_back();
////    }
////    
////    return true;
////}
//
////void Correct::trimEnd( string &seq, string &phred )
////{
////    int it = seq.find_last_not_of( 'N' );
////    if ( it == seq.npos ) return;
////    it = max( 0, it-70 );
////    string tmp = seq.substr( it );
////    LocalAlignment glocal( tmp, ender_, true, true );
////    int coords[2];
////    int score = glocal.setCoords( coords );
//////    string al = seq.substr( coords[0] + it );
////    seq = seq.substr( 0, coords[0] + it );
////    phred = phred.substr( 0, coords[0] + it );
////    while ( seq.back() == 'N' )
////    {
////        seq.pop_back();
////        phred.pop_back();
////    }
////}