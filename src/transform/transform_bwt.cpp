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

#include "transform_bwt.h"
#include "filenames.h"
#include <cassert>
#include <string.h>
#include <iostream>

BwtCycler::BwtCycler( PreprocessFiles* filenames )
: fns( filenames )
{
    // Create buffers
    inBwtBuff = new uint8_t[BWT_BUFFER];
    outBwtBuff = new uint8_t[BWT_BUFFER];
    inSapBuff = new ReadId[SAP_BUFFER];
    inPosBuff = new CharId[IDS_BUFFER];
    inEndBuff = new ReadId[IDS_BUFFER];
    outEndBuff = new ReadId[IDS_BUFFER];
    for ( int i ( 0 ); i < 4; i++ )
    {
        inIdsBuff[i] = new ReadId[IDS_BUFFER];
        outSapBuff[i] = new ReadId[SAP_BUFFER];
        outPosBuff[i] = new CharId[IDS_BUFFER];
        for ( int j ( 0 ); j < 5; j++ )
        {
            outIdsBuff[i][j] = new ReadId[IDS_BUFFER];
        }
    }
    inIdsBuff[4] = new ReadId[IDS_BUFFER];
    
    inStacksPerSap = outStacksPerSap = 4;
    inSapBytes = outSapBytes = 16;
    
    samePosFlag = (CharId)1 << 63;
    samePosMask = ~samePosFlag;
    sameByteFlag = (uint8_t)1 << 7;
    sameByteMask = ~sameByteFlag;
    
    isFinal = isPenultimate = false;
    
    FILE* bin = fns->getBinary( true, false );
    readEndBwt = writeEndBwt = readEndIds = writeEndIds = false;
    fseek( bin, 1, SEEK_SET );
    fread( &id, 8, 1, bin );
    fclose( bin );
    
    memset( &isRunArray, false, 256 );
    memset( &readArray[0], 0, 64 );
    memset( &readArray[64], 1, 64 );
    memset( &readArray[128], 2, 64 );
    memset( &readArray[192], 3, 64 );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        for ( int j ( 0 ); j < 64; j++ )
        {
            runLenArray[ i * 64 + j] = j + 1;
        }
        isRunArray[ i * 64 + 63 ] = true;
        writeBaseBit[i] = 64 * i;
        writeMaxBase[i] = 63;
        writeFullByte[i] = 64 * i + 63;
    }
    
    for ( int i = 0; i < 8; i++ )
    {
        endBitArray[i] = 1 << ( 7 - i );
    }
}

BwtCycler::~BwtCycler()
{
    if ( inBwtBuff ) delete[] inBwtBuff;
    if ( outBwtBuff ) delete[] outBwtBuff;
    if ( inSapBuff ) delete[] inSapBuff;
    if ( inPosBuff ) delete[] inPosBuff;
    if ( inEndBuff ) delete[] inEndBuff;
    if ( outEndBuff ) delete[] outEndBuff;
    for ( int i ( 0 ); i < 4; i++ )
    {
        if ( inIdsBuff[i] ) delete[] inIdsBuff[i];
        if ( outSapBuff[i] ) delete[] outSapBuff[i];
        if ( outPosBuff[i] ) delete[] outPosBuff[i];
        for ( int j ( 0 ); j < 5; j++ )
        {
            if ( outIdsBuff[i][j] ) delete[] outIdsBuff[i][j];
        }
    }
    if ( inIdsBuff[4] ) delete[] inIdsBuff[4];
}

void BwtCycler::finish( uint8_t cycle )
{
    fns->setCyclerFinal( inBwt, outBwt, inEnd, outEnd, cycle );
    prepIn();
    prepOutFinal();
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fns->setCyclerFinalIter( inPos, inSap, inIds[4], cycle, i );
        prepIter();
        finishIter( i );
    }
    
    // Flush buffers and close write files
    writeLast();
    fclose( inBwt );
    fwrite( outBwtBuff, 1, pOutBwt, outBwt );
    fclose( outBwt );
    fwrite( outEndBuff, 4, pOutEnd, outEnd );
    fclose( outEnd );
    outBwt = fns->getReadPointer( fns->bwt, true );
    uint8_t bwtBegin = 57;
    fwrite( &bwtBegin, 1, 1, outBwt );
    fwrite( &id, 8, 1, outBwt );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &charCounts[4], 8, 1, outBwt );
    fwrite( &charCounts, 8, 4, outBwt );
    fclose( outBwt );
}

void BwtCycler::finishIter( uint8_t i )
{
    ++posLeft;
    readNextPos();
    
    while( posLeft )
    {
        if ( currPos == nextPos )
        {
            writeEnd();
            readNextPos();
        }
        else
        {
            writeBwt();
        }
    }
    
    while ( currPos < charSizes[i] )
    {
        writeBwt();
    }
    assert( currPos >= charSizes[i] );
    currPos -= charSizes[i];
    
    fclose( inPos );
    fclose( inSap );
    fclose( inIds[4] );
}

void BwtCycler::flush( uint8_t cycle )
{
    // Flush buffers and close write files
    writeLast();
    fwrite( outBwtBuff, 1, pOutBwt, outBwt );
    fclose( outBwt );
    fwrite( outEndBuff, 4, pOutEnd, outEnd );
    fclose( outEnd );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fwrite( outPosBuff[i], 8, pOutPos[i], outPos[i] );
        fclose( outPos[i] );
        fwrite( outSapBuff[i], 4, pOutSap[i], outSap[i] );
        fclose( outSap[i] );
        for ( int j ( 0 ); j < 5; j++ )
        {
            fwrite( outIdsBuff[i][j], 4, pOutIds[i][j], outIds[i][j] );
            fclose( outIds[i][j] );
        }
    }
    
    // Edit in counts
    fns->setCyclerUpdate( outBwt, outEnd, outPos, outSap, outIds, cycle );
    fseek( outBwt, 9, SEEK_SET );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &charCounts, 8, 4, outBwt );
    fclose( outBwt );
    fwrite ( &endCount, 4, 1, outEnd );
    fclose( outEnd );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fwrite( &posCounts[i], 4, 1, outPos[i] );
        fclose( outPos[i] );
        fseek( outSap[i], 1, SEEK_SET );
        fwrite( &sapCounts[i], 4, 1, outSap[i] );
        fclose( outSap[i] );
        for ( int j ( 0 ); j < 5; j++ )
        {
            fwrite( &idsCounts[i][j], 4, 1, outIds[i][j] );
            fclose( outIds[i][j] );
        }
    }
}

void BwtCycler::prepIn()
{
    // Read sizes for this cycle
    bool doReadBwtEnds;
    CharId thisId;
    fread( &thisId, 8, 1, inBwt );
    if ( thisId != id )
    {
        cerr << "Error: temporary files do not appear to be from the same session." << endl;
        exit( EXIT_FAILURE );
    }
    fread( &doReadBwtEnds, 1, 1, inBwt );
    fread( &bwtLeft, 8, 1, inBwt );
    fread( &charSizes, 8, 5, inBwt );
    fread( &basePos, 4, 4, inBwt );
    fread( &endLeft, 4, 1, inEnd );
    if ( doReadBwtEnds )
    {
        if ( !readEndBwt ) setReadEnds();
        if ( !writeEndBwt ) setWriteEnds();
    }
    
    // Reset pointers and counts for cycle
    bwtFirst = true;
    currSplit = false;
    lastChar = -1;
    bwtCount = currPos = 0;
    pInBwt = BWT_BUFFER;
    pInEnd = IDS_BUFFER;
    pOutBwt = pOutEnd = 0;
    for ( int i ( 0 ); i < 4; i++ )
    {
        charCounts[i] = basePos[i];
    }
    charCounts[4] = 0;
    memset( inSapCount, 0, 20 );
    memset( outSapCount, 0, 20 );
}

void BwtCycler::prepIter()
{
    fread( &posLeft, 4, 1, inPos );
    fread( &inStacksPerSap, 1, 1, inSap );
    fread( &sapLeft, 4, 1, inSap );
    pInPos = IDS_BUFFER - 1;
    pInSap = SAP_BUFFER;
    inSapBytes = inStacksPerSap * 4;
    sapLeft *= inStacksPerSap;
    readEndIds = inStacksPerSap != 4;
    
    if ( !isFinal )
    {
        for ( int j ( 0 ); j < 4; j++ )
        {
            fread( &idsLeft[j], 4, 1, inIds[j] );
            pInIds[j] = IDS_BUFFER;
        }
    }
    
    if ( readEndIds )
    {
        fread( &idsLeft[4], 4, 1, inIds[4] );
        pInIds[4] = IDS_BUFFER;
    }

}

void BwtCycler::prepOut()
{
    anyEnds = false;
    if ( !chars )
    {
        outStacksPerSap = 1;
        isPenultimate = true;
        nextChar = 4;
    }
    else if ( ends )
    {
        anyEnds = true;
        outStacksPerSap = 5;
        writeEndIds = true;
        if ( !writeEndBwt ) setWriteEnds();
    }
    else
    {
        outStacksPerSap = 4;
    }
    outSapBytes = outStacksPerSap * 4;
    endCount = 0;
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        posCounts[i] = sapCounts[i] = 0;
        pOutPos[i] = pOutSap[i] = 0;
        fwrite( &posCounts[i], 4, 1, outPos[i] );
        fwrite( &outStacksPerSap, 1, 1, outSap[i] );
        fwrite( &sapCounts[i], 4, 1, outSap[i] );
        
        for ( int j ( 0 ); j < 5; j++ )
        {
            idsCounts[i][j] = 0;
            pOutIds[i][j] = 0;
            fwrite( &idsCounts[i][j], 4, 1, outIds[i][j] );
        }
    }
    
    // Set place holders for this cycle size
    fwrite( &id, 8, 1, outBwt );
    fwrite( &writeEndBwt, 1, 1, outBwt );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &charCounts, 8, 5, outBwt );
    fwrite( &basePos, 4, 4, outBwt );
    fwrite( &endCount, 4, 1, outEnd );
}

void BwtCycler::prepOutFinal()
{
    isFinal = true;
    if ( !writeEndBwt ) setWriteEnds();
    
    uint8_t bwtBegin = 57, idsBegin = 9;
    CharId endCount = basePos[0] + basePos[1] + basePos[2] + basePos[3];
    fwrite( &bwtBegin, 1, 1, outBwt );
    fwrite( &id, 8, 1, outBwt );
    fwrite( &bwtCount, 8, 1, outBwt );
    fwrite( &endCount, 8, 1, outBwt );
    fwrite( &charCounts, 8, 4, outBwt );
    
    fwrite( &idsBegin, 1, 1, outEnd );
    fwrite( &id, 8, 1, outEnd );
    
    memset( &charCounts, 0, 40 );
    for ( int i ( 0 ); i < 4; i++ )
    {
        if ( basePos[i] )
        {
            writeRun( i, basePos[i] );
        }
    }
    currPos = 0;
}

void BwtCycler::readBwtIn()
{
    fread( inBwtBuff, 1, min( BWT_BUFFER, bwtLeft), inBwt );
    pInBwt = 0;
}

void BwtCycler::readNextId()
{
    // Refresh IDs buffer if necessary
    if ( pInIds[thisChar] == IDS_BUFFER )
    {
        fread( inIdsBuff[thisChar], 4, min( idsLeft[thisChar], IDS_BUFFER ), inIds[thisChar] );
        pInIds[thisChar] = 0;
    }

    // Set id and next character
    nextId = inIdsBuff[thisChar][ pInIds[thisChar]++ ];
    if ( anyEnds && ( ends[ (nextId / 2) / 8 ] & endBitArray[ (nextId / 2) % 8 ] ) )
    {
        nextChar = 4;
    }
    else if ( !isPenultimate )
    {
        nextChar = byteToInt[ nextId & 0x3 ][ chars[ nextId / 4 ] ];
    }
    --idsLeft[thisChar];
}

void BwtCycler::readNextPos()
{
    // Refresh positions buffer if necessary
    if ( --posLeft && ++pInPos == IDS_BUFFER )
    {
        fread( inPosBuff, 8, min( posLeft, IDS_BUFFER ), inPos );
        pInPos = 0;
    }
    
    if ( !posLeft )
    {
        nextPos = -1;
    }
    
    // Next position include multiple characters with the same suffix
    else if ( inPosBuff[pInPos] & samePosFlag )
    {
        nextSame = true;
        nextPos = inPosBuff[pInPos] & samePosMask;
        
        // Refresh same buffer if necessary
        if ( pInSap == SAP_BUFFER )
        {
            fread( inSapBuff, 4, min( sapLeft, SAP_BUFFER ), inSap );
            pInSap = 0;
        }
        
        // Load same counts
        if ( isFinal )
        {
            inSapCount[4] = inSapBuff[ pInSap++ ];
        }
        else
        {
            memcpy( &inSapCount, &inSapBuff[pInSap], inSapBytes );
            pInSap += inStacksPerSap;
        }
    }
    
    // Next position does not share its suffix with any other character
    else
    {
        // Set which stack and position
        thisChar = inPosBuff[pInPos] & 0x7;
        nextPos = inPosBuff[pInPos] >> 3;
        nextSame = false;
        
        // Set id and next character
        if ( !isFinal ) readNextId();
    }
    
    if ( currSplit && currPos < nextPos )
    {
        if ( currPos + splitRun > nextPos )
        {
            ReadId thisRun = nextPos - currPos;
            splitRun -= thisRun;
            writeRun( splitChar, thisRun );
        }
        else
        {
            writeRun( splitChar, splitRun );
            currSplit = false;
        }
    }
}

void BwtCycler::run( uint8_t* inChars, uint8_t* inEnds, uint8_t cycle )
{
    fns->setCycler( inBwt, outBwt, inEnd, outEnd, outPos, outSap, outIds, cycle );
    chars = inChars;
    ends = inEnds;
    prepIn();
    prepOut();
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        fns->setCyclerIter( inPos, inSap, inIds, cycle, i );
        prepIter();
        runIter( i );
    }
    
    flush( cycle );
}

void BwtCycler::runIter( uint8_t i )
{
    ++posLeft;
    readNextPos();
    
    while ( posLeft )
    {
        if ( currPos == nextPos )
        {
            if ( nextSame )
            {
                writeSame();
            }
            else
            {
                writeNext();
            }
            
            readNextPos();
        }
        else
        {
            writeBwt();
        }
        
    }
    
    while ( currPos < charSizes[i] )
    {
        writeBwt();
    }
    assert( currPos >= charSizes[i] );
    currPos -= charSizes[i];
    
    fclose( inPos );
    fclose( inSap );
    for ( int j ( 0 ); j < 4; j++ )
    {
        fclose( inIds[j] );
    }
}

void BwtCycler::setReadEnds()
{
    memset( &isRunArray, false, 256 );
    memset( &readArray[0], 0, 63 );
    memset( &readArray[63], 1, 63 );
    memset( &readArray[126], 2, 63 );
    memset( &readArray[189], 3, 63 );
    memset( &readArray[252], 4, 4 );
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        runLenArray[ 252 + i ] = i + 1;
        for ( int j ( 0 ); j < 63; j++ )
        {
            runLenArray[ i * 63 + j] = j + 1;
        }
        isRunArray[ i * 63 + 62 ] = true;
    }
    isRunArray[255] = true;
    
    readEndBwt = true;
}

void BwtCycler::setWriteEnds()
{
    for ( int i ( 0 ); i < 4; i++ )
    {
        writeBaseBit[i] = 63 * i;
        writeMaxBase[i] = 62;
        writeFullByte[i] = 63 * i + 62;
    }
    writeBaseBit[4] = 252;
    writeMaxBase[4] = 3;
    writeFullByte[4] = 255;
    
    writeEndBwt = true;
}

void BwtCycler::writeIdsToFile( uint8_t i, uint8_t j )
{
    fwrite( outIdsBuff[i][j], 4, pOutIds[i][j], outIds[i][j] );
    pOutIds[i][j] = 0;
}

void BwtCycler::writeBwt()
{
    if ( pInBwt == BWT_BUFFER ) readBwtIn();
    uint8_t currChar = inBwtBuff[pInBwt];
    uint8_t c = readArray[ currChar ];
    ReadId runLen = runLenArray[ currChar ];
    --bwtLeft;
    ++pInBwt;
    
    // Count run length if greater than 63
    if ( isRunArray[currChar] )
    {
        if ( pInBwt == BWT_BUFFER ) readBwtIn();
        ReadId thisRun = 0;
        uint8_t shiftCount = 0;
        do {
            if ( pInBwt == BWT_BUFFER ) readBwtIn();
            thisRun ^= ( inBwtBuff[pInBwt] & sameByteMask ) << ( shiftCount++ * 7 );
        } while ( inBwtBuff[pInBwt++] & sameByteFlag );
        
        runLen += thisRun;
        bwtLeft -= shiftCount;
    }
    
    if ( currPos + runLen > nextPos )
    {
        splitChar = c;
        splitRun = currPos + runLen - nextPos;
        runLen -= splitRun;
        currSplit = true;
    }
    
    writeRun( c, runLen );
    
    if ( c == 4 )
    {
        while ( runLen-- )
        {
            if ( pInEnd == IDS_BUFFER )
            {
                fread( inEndBuff, 4, min( endLeft, IDS_BUFFER ), inEnd );
                pInEnd = 0;
            }
            if ( pOutEnd == IDS_BUFFER )
            {
                fwrite( outEndBuff, 4, IDS_BUFFER, outEnd );
                pOutEnd = 0;
            }
            outEndBuff[ pOutEnd++ ] = inEndBuff[ pInEnd++ ];
            --endLeft;
        }
    }
}

void BwtCycler::writeBwtByte( uint8_t c )
{
    if ( pOutBwt == BWT_BUFFER )
    {
        fwrite( outBwtBuff, 1, BWT_BUFFER, outBwt );
        pOutBwt = 0;
    }
    
    outBwtBuff[ pOutBwt++ ] = c;
    ++bwtCount;
}

void BwtCycler::writeEnd()
{
    if ( !nextSame )
    {
        inSapCount[4] = 1;
    }
    
    writeRun( 4, inSapCount[4] );
    
    while ( inSapCount[4]-- )
    {
        if ( pInIds[4] == IDS_BUFFER )
        {
            fread( inIdsBuff[4], 4, min( IDS_BUFFER, idsLeft[4] ), inIds[4] );
            pInIds[4] = 0;
        }
        if ( pOutEnd == IDS_BUFFER )
        {
            fwrite( outEndBuff, 4, IDS_BUFFER, outEnd );
            pOutEnd = 0;
        }
        
        outEndBuff[ pOutEnd++ ] = inIdsBuff[4][ pInIds[4]++ ];
    }
}

void BwtCycler::writeLast()
{
    if ( bwtFirst )
    {
        bwtFirst = false;
        return;
    }
    
    uint8_t maxBase = writeMaxBase[lastChar];
    if ( lastRun < maxBase )
    {
        writeBwtByte( writeBaseBit[lastChar] + lastRun );
    }
    // Encode a new byte for every power of 128
    else
    {
        writeBwtByte( writeFullByte[lastChar] );
        lastRun -= maxBase;
        while ( lastRun >= 128 )
        {
            writeBwtByte( sameByteFlag ^ ( lastRun & sameByteMask ) );
            lastRun >>= 7;
        }
        writeBwtByte( lastRun );
    }
}

void BwtCycler::writeNext()
{
    if ( thisChar != 4 )
    {
        // Write positions buffer to file if full
        if ( pOutPos[thisChar] == IDS_BUFFER )
        {
            fwrite( outPosBuff[thisChar], 8, IDS_BUFFER, outPos[thisChar] );
            pOutPos[thisChar] = 0;
        }

        // Write position to buffer
        outPosBuff[thisChar][ pOutPos[thisChar]++ ] = ( charCounts[thisChar] << 3 ) + nextChar;
        ++posCounts[thisChar];
    }
    
    // Write id to IDs stack
    writeNextId();
    writeRun( thisChar, 1 );
}

void BwtCycler::writeNextId()
{
    if ( thisChar == 4 )
    {
        if ( pOutEnd == IDS_BUFFER )
        {
            fwrite( outEndBuff, 4, IDS_BUFFER, outEnd );
            pOutEnd = 0;
        }
        
        outEndBuff[ pOutEnd++ ] = nextId;
        ++endCount;
    }
    else
    {
        // Write IDs buffer to file if full
        if ( pOutIds[thisChar][nextChar] == IDS_BUFFER )
        {
            fwrite( outIdsBuff[thisChar][nextChar], 4, IDS_BUFFER, outIds[thisChar][nextChar] );
            pOutIds[thisChar][nextChar] = 0;
        }

        // Write id to buffer
        outIdsBuff[thisChar][nextChar][ pOutIds[thisChar][nextChar]++ ] = nextId;
        idsCounts[thisChar][nextChar]++;
    }
}

void BwtCycler::writeRun( uint8_t c, ReadId runLen )
{
    if ( c == lastChar )
    {
        lastRun += runLen;
    }
    else
    {
        writeLast();
        lastChar = c;
        lastRun = runLen - 1;
    }
    
    charCounts[c] += runLen;
    currPos += runLen;
}

void BwtCycler::writeSame()
{
    if ( inSapCount[4] )
    {
        thisChar = 4;
        writeRun( thisChar, inSapCount[4] );
        while ( inSapCount[4]-- )
        {
            readNextId();
            writeNextId();
        }
    }
    
    for ( int i ( 0 ); i < 4; i++ )
    {
        if ( inSapCount[i] )
        {
            thisChar = i;
            
            // Write positions buffer to file if full
            if ( pOutPos[i] == IDS_BUFFER )
            {
                fwrite( outPosBuff[i], 8, IDS_BUFFER, outPos[i] );
                pOutPos[i] = 0;
            }
            
            // Write same run buffer to file if full
            if ( pOutSap[i] == SAP_BUFFER )
            {
                fwrite( outSapBuff[i], 4, SAP_BUFFER, outSap[i] );
                pOutSap[i] = 0;
            }
            
            // Next char is part of a run
            if ( inSapCount[i] > 1 )
            {
                // Write this run and position
                outPosBuff[i][ pOutPos[i]++ ] = samePosFlag ^ charCounts[i];
                posCounts[i]++;
                writeRun( thisChar, inSapCount[i] );
                
                // Write IDs and count same runs
                memset( &outSapCount, (uint8_t)0, 20 );
                while ( inSapCount[i]-- )
                {
                    readNextId();
                    writeNextId();
                    outSapCount[nextChar]++;
                }
                
                // Write same run to buffer
                if ( isPenultimate )
                {
                    outSapBuff[i][ pOutSap[i]++ ] = outSapCount[4];
                }
                else
                {
                    memcpy( &outSapBuff[i][ pOutSap[i] ], &outSapCount, outSapBytes );
                    pOutSap[i] += outStacksPerSap;
                }
                
                ++sapCounts[i];
            }
            
            // Next char is not part of a run
            else
            {
                readNextId();
                writeNextId();
                
                outPosBuff[i][ pOutPos[i]++ ] = ( charCounts[i] << 3 ) ^ nextChar;
                posCounts[i]++;
                writeRun( thisChar, 1 );
            }
        }
    }
    
    sapLeft -= inStacksPerSap;
}

