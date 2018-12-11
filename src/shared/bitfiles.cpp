/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "bitfiles.h"
#include "error.h"
#include <cassert>

BitFile::BitFile( string fn, int bits )
: fn_( fn ), fp_( NULL ), bits_( bits ), per_( 8 / max( 1, bits ) ), edited_( false )
{
    if ( bits_ != 1 && bits_ != 2 && bits_ != 4 ) failure( "Invalid bit file type." );
    
    uint8_t limit = bits == 1 ? 1 : ( bits == 2 ? 3 : 15 );
    for ( int i = 0; i < per_; i++ )
    {
        int shift = ( per_ - 1 - i ) * bits_;
        flags_[i] = ~( limit << shift );
        for ( int j = 0; j < 256; j++ ) readTable_[i][j] = ( j >> shift ) & limit;
        for ( int j = 0; j <= limit; j++ ) writeTable_[i][j] = j << shift;
    }
}

BitFile::~BitFile()
{
    close();
}

bool BitFile::open( bool doRead, bool doWrite, bool overwrite )
{
    edited_ = writing_ = false;
    i_ = doRead ? per_ - 1 : 0;
    if ( fp_ ) failure( "Tried to open file before closing." );
    
    if ( doRead && doWrite )
    {
        fp_ = fopen( fn_.c_str(), "rb+" );
        if ( !fp_ ) failure( "Failed to open file for editing." );
    }
    else if ( doRead )
    {
        fp_ = fopen( fn_.c_str(), "rb" );
        if ( !fp_ ) failure( "Failed to open file for reading." );
    }
    else if ( doWrite )
    {
        if ( !overwrite )
        {
            fp_ = fopen( fn_.c_str(), "rb" );
            if ( fp_ )
            {
                fclose( fp_ );
                fp_ = NULL;
                return false;
            }
        }
        fp_ = fopen( fn_.c_str(), "wb" );
        if ( !fp_ ) failure( "Failed to open file for writing: \"" + fn_ + "\"" );
        writing_ = true;
    }
    
    return true;
}

void BitFile::close( bool ignore )
{
    if ( ignore && !fp_ ) return;
    if ( !fp_ ) failure( "Tried to close a file that was not open" );
    
    if ( writing_ && i_ ) fwrite( &byte_, 1, 1, fp_ );
    if ( edited_ )
    {
        fseek( fp_, -1, SEEK_CUR );
        fwrite( &byte_, 1, 1, fp_ );
    }
    
    fclose( fp_ );
    fp_ = NULL;
    edited_ = writing_ = false;
}

void BitFile::edit( uint8_t code )
{
    byte_ = ( byte_ & flags_[i_] ) + writeTable_[i_][code];
    edited_ = true;
}

uint8_t BitFile::read()
{
    
    if ( ++i_ == per_ )
    {
        if ( edited_ )
        {
            fseek( fp_, -1, SEEK_CUR );
            fwrite( &byte_, 1, 1, fp_ );
        }
        i_ = 0;
        fread( &byte_, 1, 1, fp_ );
        edited_ = false;
    }
    return readTable_[i_][byte_];
}

void BitFile::setPos( uint64_t pos, int i )
{
    fseek( fp_, pos - ( i < per_ ), SEEK_SET );
    if ( i < per_ ) fread( &byte_, 1, 1, fp_ );
    i_ = i;
}

void BitFile::tellPos( uint64_t &pos, int &i )
{
    pos = ftell( fp_ );
    i = i_;
}

void BitFile::write( uint8_t code )
{
    if ( code ) byte_ = ( byte_ & flags_[i_] ) + writeTable_[i_][code];
    if ( ++i_ == per_ )
    {
        fwrite( &byte_, 1, 1, fp_ );
        i_ = 0;
        byte_ = 0;
    }
}

