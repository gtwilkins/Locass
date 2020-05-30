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

#ifndef PRUNE_BUBBLE_H
#define PRUNE_BUBBLE_H

#include "node.h"

//struct BubbleBranch
//{
//    BubbleBranch( Node* fork, vector<Node*> branch, int32_t coord, int32_t reads );
//    bool pop( Node* tar, Nodes& fwd, bool drxn );
//    void set( Node* node, vector<Node*>& path, int32_t dist, Node* tar, Nodes& fwd, bool drxn );
//    Node* fork;
//    vector<Node*> branch;
//    vector<BubbleBranch> alts;
//    int32_t coord, reads;
//};

class Bubble
{
    Bubble( Node* fork, Edge& e, Nodes& fwd, bool drxn );
    Bubble( vector<Node*>& path, int32_t dist, bool drxn );
    bool confirm();
    bool pop( Nodes& prunes, bool drxn  );
    static bool resolve( Node* node, Nodes& prunes, bool drxn );
    void set( Node* node, vector<Node*>& path, int32_t dist, Nodes& fwd, bool drxn );
    Node* forks[2];
    vector<Node*> branch;
    vector<Bubble> alts;
    int32_t coord, reads;
    bool branched;
    
public:
    static bool prune( NodeRoll& nodes );
};

#endif /* PRUNE_BUBBLE_H */

